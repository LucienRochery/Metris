// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_quafun_connection
#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "codegen_lag2bez.hxx"
#include "ho_constants.hxx"
#include "linalg/symidx.hxx"
#include "quality/quafun.hxx"
#include "quality/low_metqua.hxx"
#include "quality/low_metqua_d.hxx"
#include "quality/objective_quadrature_sample.hxx"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace Metris{
namespace{

// Push-forward of a constant SPD metric by psi(xi)=(1.5xi-.5xi^2,eta,zeta).
// The off-diagonal entry ensures this is not just a diagonal-norm test.
template<int dim>
void manufactured_metric(const AnaMetCtx*, const double *x, double scale,
                          int idiff, double *m, double *dm){
  constexpr int nm = dim*(dim+1)/2;
  const double d = 2.25-2*x[0];
  std::fill(m,m+nm,0.);
  for(int k = 0; k < dim; k++) m[sym2idx(k,k)] = scale;
  m[0] = scale/d;
  m[1] = 0.5*scale/std::sqrt(d);
  if(idiff){
    std::fill(dm,dm+dim*nm,0.);
    dm[0] = 2*scale/(d*d);
    dm[1] = 0.5*scale/std::pow(d,1.5);
  }
}

template<int dim, int degree, class MFT>
void initialize(Mesh<MFT> &msh, MetrisParameters &parameters){
  msh.idim = dim;
  msh.curdeg = degree;
  msh.strdeg = degree;
  msh.param = &parameters;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.set_npoin(getnnode(dim,degree));
  msh.set_nentt(dim,1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);
  for(int i = 0; i < getnnode(dim,degree); i++) msh.ent2poi(dim)(0,i) = i;
}

template<int dim, int degree>
void check_manufactured(){
  constexpr int nn = getnnode(dim,degree);
  auto ordering = ORDELT(dim);
  // Both triangles on each diagonal, followed by a non-grid affine simplex.
  const double triangles[5][3][2] = {
      {{.1,.1},{.7,.1},{.1,.7}}, {{.7,.7},{.1,.7},{.7,.1}},
      {{.1,.1},{.7,.1},{.7,.7}}, {{.1,.1},{.7,.7},{.1,.7}},
      {{.12,.19},{.69,.27},{.25,.83}}};
  for(const auto &triangle : triangles){
    MetrisParameters parameters;
    parameters.iverb = 0;
    Mesh<MetricFieldAnalytical> msh;
    initialize<dim,degree>(msh,parameters);
    if constexpr(dim == 2) parameters.setAnalyticalMetric(9);
    else parameters.setAnalyticalMetric(AnaMetFun(manufactured_metric<dim>));
    msh.met.setAnalyticalMetric(parameters);
    double vertices[dim+1][dim] = {};
    for(int i = 0; i < 3; i++){
      for(int k = 0; k < 2; k++) vertices[i][k] = triangle[i][k];
    }
    if constexpr(dim == 3){
      vertices[3][0] = .23;
      vertices[3][1] = .34;
      vertices[3][2] = .8;
    }
    for(int i = 0; i < nn; i++){
      for(int k = 0; k < dim; k++){
        double x = 0.;
        for(int j = 0; j <= dim; j++) x += ordering[degree][i][j]*vertices[j][k]/degree;
        msh.coord(i,k) = k == 0 ? 1.5*x-.5*x*x : x;
      }
    }
    // Barycenter, every vertex, and an edge point. Test the actual point,
    // including vertices (analytical derivatives must not be skipped there).
    for(int basis = 0; basis < 2; basis++){
      if(basis){
        std::array<double,nn*dim> bez_data;
        dblAr2 bez(nn,dim,bez_data.data());
        if constexpr(dim == 2) lag2bez2<degree,dim>(msh.ent2poi(dim)[0],msh.coord,bez);
        else lag2bez3<degree,dim>(msh.ent2poi(dim)[0],msh.coord,bez);
        for(int i = 0; i < nn; i++){
          for(int k = 0; k < dim; k++) msh.coord(i,k) = bez(i,k);
        }
        msh.forceBasisFlag(FEBasis::Bezier);
      }
      for(int sample = -2; sample <= dim; sample++){
        double bary[dim+1] = {};
        if(sample == -2) std::fill(bary,bary+dim+1,1./(dim+1));
        else if(sample == -1){bary[0] = .3; bary[1] = .7;}
        else bary[sample] = 1.;
        const double value = quafun_connection<MetricFieldAnalytical,dim>(
            msh,msh.ent2poi(dim)[0],bary);
        if constexpr(degree > 1){
          BOOST_CHECK_SMALL(value,1.e-23);
          for(int ivar = 0; ivar < nn; ivar++){
            double gradient[dim];
            const double val = d_quafun_connection<MetricFieldAnalytical,dim>(
                msh,msh.ent2poi(dim)[0],bary,ivar,msh.getBasis(),gradient);
            BOOST_CHECK_EQUAL(val,value);
            for(double g : gradient) BOOST_CHECK_SMALL(g,2.e-10);
          }
        }else{
          // Affine geometry has zero connection; only Gamma^0_00=1/d
          // of the target remains. Its target norm is 16/(9d).
          double x = 0.;
          for(int i = 0; i <= dim; i++) x += bary[i]*msh.coord(i,0);
          BOOST_CHECK_CLOSE(value,16./(9*(2.25-2*x)),1.e-9);
        }
      }
    }
    if constexpr(degree > 1){
      const double bary[dim+1] = {.2,.3,.5};
      msh.coord(dim+1,0) += .01;
      BOOST_CHECK_GT((quafun_connection<MetricFieldAnalytical,dim>(
          msh,msh.ent2poi(dim)[0],bary)),1.e-8);
    }
  }
}

template<int dim, int degree>
void check_fe_derivatives(){
  // log M = alpha*x I is exactly represented, including on a Q2 psi map.
  // Here Gamma_K^0_00 = gamma, all other components zero. The norm is
  // exp(-alpha*x)*(gamma^2-alpha*gamma+(3*dim-2)*alpha^2/4).
  constexpr double alpha = .7;
  constexpr int nm = dim*(dim+1)/2;
  constexpr int nn = getnnode(dim,degree);
  MetrisParameters parameters;
  parameters.iverb = 0;
  Mesh<MetricFieldFE> msh;
  initialize<dim,degree>(msh,parameters);
  auto ordering = ORDELT(dim);
  for(int i = 0; i < nn; i++){
    for(int k = 0; k < dim; k++){
      const double xi = double(ordering[degree][i][k+1])/degree;
      msh.coord(i,k) = degree > 1 && k == 0 ? 1.5*xi-.5*xi*xi : xi;
    }
  }
  double bary[dim+1];
  std::fill(bary,bary+dim+1,1./(dim+1));
  for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
    msh.met.forceSpaceFlag(space);
    for(int i = 0; i < nn; i++){
      for(int k = 0; k < nm; k++) msh.met(i,k) = 0.;
      const double logm = alpha*msh.coord(i,0);
      for(int k = 0; k < dim; k++){
        msh.met(i,sym2idx(k,k)) = space == MetSpace::Log ? logm : std::exp(logm);
      }
    }
    const double x = degree > 1 ? 1.5*bary[1]-.5*bary[1]*bary[1] : bary[1];
    const double gamma = degree > 1 ? 1./(2.25-2*x) : 0.;
    const double expected = std::exp(-alpha*x)
        *(gamma*gamma-alpha*gamma+(3*dim-2)*alpha*alpha/4);
    std::array<double,nn*dim> original_coordinates, bez_data;
    for(int i = 0; i < nn; i++){
      for(int k = 0; k < dim; k++) original_coordinates[i*dim+k] = msh.coord(i,k);
    }
    for(int basis = 0; basis < 2; basis++){
      if(basis){
        dblAr2 bez(nn,dim,bez_data.data());
        if constexpr(dim == 2) lag2bez2<degree,dim>(msh.ent2poi(dim)[0],msh.coord,bez);
        else lag2bez3<degree,dim>(msh.ent2poi(dim)[0],msh.coord,bez);
        for(int i = 0; i < nn; i++){
          for(int k = 0; k < dim; k++) msh.coord(i,k) = bez(i,k);
        }
        msh.forceBasisFlag(FEBasis::Bezier);
      }
      BOOST_CHECK_CLOSE((quafun_connection<MetricFieldFE,dim>(
          msh,msh.ent2poi(dim)[0],bary)),expected,1.e-10);
      BOOST_CHECK(msh.met.getSpace() == space);
    }
    for(int i = 0; i < nn; i++){
      for(int k = 0; k < dim; k++) msh.coord(i,k) = original_coordinates[i*dim+k];
    }
    msh.forceBasisFlag(FEBasis::Lagrange);
  }
}

template<int dim>
void check_coordinate_invariance(){
  constexpr int nm = dim*(dim+1)/2;
  using Matrix = Eigen::Matrix<double,dim,dim>;
  Matrix jac = Matrix::Identity(), B = Matrix::Identity(), m = Matrix::Identity();
  jac(0,1) = .4; jac(1,0) = -.2; jac(0,0) = 1.3;
  B(0,1) = .7; B(1,0) = -.3; B(0,0) = 1.7;
  m(0,0) = 2.; m(0,1) = m(1,0) = .3;
  std::array<double,dim*dim> j, jnew;
  std::array<double,nm> met, metnew;
  std::array<double,dim*nm> h, hnew, dm, dmnew;
  for(int i = 0; i < dim*nm; i++){
    h[i] = .03*(i+1);
    dm[i] = .02*(i-2);
  }
  const Matrix Binv = B.inverse();
  const Matrix newjac = B*jac, newmet = Binv.transpose()*m*Binv;
  for(int i = 0; i < dim; i++){
    for(int k = 0; k < dim; k++){
      j[i*dim+k] = jac(k,i); jnew[i*dim+k] = newjac(k,i);
      met[sym2idx(i,k)] = m(i,k); metnew[sym2idx(i,k)] = newmet(i,k);
    }
  }
  for(int ab = 0; ab < nm; ab++){
    for(int k = 0; k < dim; k++){
      hnew[ab*dim+k] = 0.;
      for(int l = 0; l < dim; l++) hnew[ab*dim+k] += B(k,l)*h[ab*dim+l];
    }
  }
  for(int r = 0; r < dim; r++){
    Matrix deriv = Matrix::Zero();
    for(int s = 0; s < dim; s++){
      for(int i = 0; i < dim; i++){
        for(int k = 0; k < dim; k++) deriv(i,k) += Binv(s,r)*dm[s*nm+sym2idx(i,k)];
      }
    }
    const Matrix transformed = Binv.transpose()*deriv*Binv;
    for(int i = 0; i < dim; i++){
      for(int k = 0; k <= i; k++) dmnew[r*nm+sym2idx(i,k)] = transformed(i,k);
    }
  }
  const double value = connection_difference_norm2<dim>(j.data(),h.data(),met.data(),dm.data());
  BOOST_CHECK_GT(value,0.);
  BOOST_CHECK_CLOSE(connection_difference_norm2<dim>(jnew.data(),hnew.data(),metnew.data(),dmnew.data()),value,1.e-10);
  // Constant scaling preserves connections but rescales this norm by 1/c.
  for(double &v : met) v *= 7.;
  for(double &v : dm) v *= 7.;
  BOOST_CHECK_CLOSE(connection_difference_norm2<dim>(j.data(),h.data(),met.data(),dm.data()),value/7,1.e-10);
}

// Smooth non-diagonal SPD field with mixed derivatives in every direction.
// M = I + w(x) w(x)^T, w_i = .2*(i+1) + sum_r A_ir x_r.
template<int dim>
double affine_weight(int i, int r){return .1*(i+1)*(r+1)+(i == r ? .2 : 0.);}

template<int dim>
void polynomial_metric(const AnaMetCtx*, const double *x, double scale,
                        int idiff, double *m, double *dm){
  constexpr int nm = dim*(dim+1)/2;
  double w[dim];
  for(int i = 0; i < dim; i++){
    w[i] = .2*(i+1);
    for(int r = 0; r < dim; r++) w[i] += affine_weight<dim>(i,r)*x[r];
  }
  for(int i = 0; i < dim; i++){
    for(int j = 0; j <= i; j++){
      m[sym2idx(i,j)] = scale*((i == j ? 1. : 0.) + w[i]*w[j]);
      if(idiff){
        for(int r = 0; r < dim; r++){
          dm[r*nm+sym2idx(i,j)] = scale*(affine_weight<dim>(i,r)*w[j]
                                                        +w[i]*affine_weight<dim>(j,r));
        }
      }
    }
  }
}

// Express a unit perturbation in dofbas in the actual geometry storage basis.
template<int dim, int degree>
std::array<double,getnnode(dim,degree)> perturbation_coefficients(
    FEBasis storage, FEBasis dofbas, int ivar){
  constexpr int nn = getnnode(dim,degree);
  std::array<double,nn> coefficients{}, result{};
  std::array<int,nn> nodes;
  for(int i = 0; i < nn; i++) nodes[i] = i;
  coefficients[ivar] = 1.;
  if(storage == dofbas) return coefficients;
  dblAr2 input(nn,1,coefficients.data()), output(nn,1,result.data());
  if(storage == FEBasis::Bezier){
    if constexpr(dim == 2) lag2bez2<degree,1>(nodes.data(),input,output);
    else lag2bez3<degree,1>(nodes.data(),input,output);
  }else{
    auto ordering = ORDELT(dim);
    for(int i = 0; i < nn; i++){
      double bary[dim+1];
      for(int k = 0; k <= dim; k++) bary[k] = double(ordering[degree][i][k])/degree;
      if constexpr(dim == 2){
        eval2<1,degree>(input,nodes.data(),dofbas,DifVar::None,DifVar::None,
                        bary,&result[i],nullptr,nullptr);
      }else{
        eval3<1,degree>(input,nodes.data(),dofbas,DifVar::None,DifVar::None,
                        bary,&result[i],nullptr,nullptr);
      }
    }
  }
  return result;
}

template<int dim, int degree, class MFT>
void evaluate_test_geometry(Mesh<MFT> &msh, const double *bary,
                            double *point, double *jac, double *hess){
  constexpr int nn = getnnode(dim,degree);
  std::array<double,nn*dim> data;
  dblAr2 coefficients(nn,dim,data.data());
  const int *nodes = msh.ent2poi(dim)[0];
  if(msh.getBasis() == FEBasis::Lagrange){
    if constexpr(dim == 2) lag2bez2<degree,dim>(nodes,msh.coord,coefficients);
    else lag2bez3<degree,dim>(nodes,msh.coord,coefficients);
  }else{
    for(int i = 0; i < nn; i++) for(int k = 0; k < dim; k++)
      coefficients(i,k) = msh.coord(i,k);
  }
  if constexpr(dim == 2)
    eval2<dim,degree>(coefficients,nodes,FEBasis::Bezier,DifVar::Bary,
                     DifVar::Bary,bary,point,jac,hess);
  else
    eval3<dim,degree>(coefficients,nodes,FEBasis::Bezier,DifVar::Bary,
                     DifVar::Bary,bary,point,jac,hess);
}

template<int dim, int degree, class MFT>
void check_frozen_control_gradients(){
  constexpr int nn = getnnode(dim,degree), nm = dim*(dim+1)/2;
  MetrisParameters parameters;
  parameters.iverb = 0;
  Mesh<MFT> msh;
  initialize<dim,degree>(msh,parameters);
  auto ordering = ORDELT(dim);
  for(int i = 0; i < nn; i++){
    double xi[dim];
    for(int k = 0; k < dim; k++) xi[k] = double(ordering[degree][i][k+1])/degree;
    for(int k = 0; k < dim; k++){
      msh.coord(i,k) = .1 + .6*xi[k] + .05*xi[(k+1)%dim]
                            + .08*xi[0]*xi[1] + .03*xi[k]*xi[k];
    }
  }
  if constexpr(std::is_same_v<MFT,MetricFieldAnalytical>){
    parameters.setAnalyticalMetric(AnaMetFun(polynomial_metric<dim>));
    msh.met.setAnalyticalMetric(parameters);
  }else{
    for(int i = 0; i < nn; i++){
      polynomial_metric<dim>(nullptr,msh.coord[i],1.,0,msh.met[i],nullptr);
    }
  }
  std::array<double,nn*nm> lag_metric_data{};
  if constexpr(std::is_same_v<MFT,MetricFieldFE>){
    for(int i = 0; i < nn; i++){
      for(int m = 0; m < nm; m++) lag_metric_data[i*nm+m] = msh.met(i,m);
      getlogmet_inp<dim,double>(&lag_metric_data[i*nm]);
    }
  }
  for(FEBasis storage : {FEBasis::Lagrange,FEBasis::Bezier}){
    if(storage == FEBasis::Bezier){
      std::array<double,nn*dim> bez_data;
      dblAr2 bez(nn,dim,bez_data.data());
      if constexpr(dim == 2) lag2bez2<degree,dim>(msh.ent2poi(dim)[0],msh.coord,bez);
      else lag2bez3<degree,dim>(msh.ent2poi(dim)[0],msh.coord,bez);
      for(int i = 0; i < nn; i++){
        for(int k = 0; k < dim; k++) msh.coord(i,k) = bez(i,k);
      }
      msh.forceBasisFlag(storage);
    }
    // Check both FE metric storage spaces and bases; conversions preserve
    // the polynomial interpolant. Analytical fields have no FE coefficients.
    for(int representation = 0; representation < (std::is_same_v<MFT,MetricFieldFE> ? 4 : 1); representation++){
      if constexpr(std::is_same_v<MFT,MetricFieldFE>){
        // Construct local coefficients directly: this isolated-element test
        // does not have the tags/neighbors required by whole-mesh conversion.
        auto data = lag_metric_data;
        if(representation >= 2){
          const dblAr2 lag(nn,nm,lag_metric_data.data());
          dblAr2 bez(nn,nm,data.data());
          if constexpr(dim == 2) lag2bez2<degree,nm>(msh.ent2poi(dim)[0],lag,bez);
          else lag2bez3<degree,nm>(msh.ent2poi(dim)[0],lag,bez);
        }
        msh.met.forceBasisFlag(representation < 2 ? FEBasis::Lagrange : FEBasis::Bezier);
        msh.met.forceSpaceFlag(representation%2 ? MetSpace::Log : MetSpace::Exp);
        for(int i = 0; i < nn; i++){
          if(representation%2 == 0) getexpmet_inp<dim,double>(&data[i*nm]);
          for(int m = 0; m < nm; m++) msh.met(i,m) = data[i*nm+m];
        }
      }
      std::array<double,nn*nm> original_metric{};
      if constexpr(std::is_same_v<MFT,MetricFieldFE>){
        for(int i = 0; i < nn; i++){
          for(int m = 0; m < nm; m++) original_metric[i*nm+m] = msh.met(i,m);
        }
      }
      for(int sample = 0; sample < 3; sample++){
        double bary[dim+1] = {};
        if(sample == 0) std::fill(bary,bary+dim+1,1./(dim+1));
        if(sample == 1){bary[0] = .3; bary[1] = .7;}
        if(sample == 2) bary[0] = 1.;
        // Capture target data once, before perturbing the geometry. FE target
        // derivatives use an independent finite-difference value oracle.
        double point[dim], jac[dim*dim], hess[nm*dim], met[nm], dmet[dim*nm]{};
        evaluate_test_geometry<dim,degree>(msh,bary,point,jac,hess);
        if constexpr(std::is_same_v<MFT,MetricFieldAnalytical>){
          msh.met.getMetPhys(DifVar::Phys,MetSpace::Exp,point,met,dmet);
        }else{
          const int *nodes = msh.ent2poi(dim)[0];
          msh.met.getMetBary(AsDeg::Pk,DifVar::None,MetSpace::Exp,
                            nodes,dim,bary,met,nullptr);
          Eigen::Matrix<double,dim,dim> J;
          for(int k = 0; k < dim; k++) for(int a = 0; a < dim; a++)
            J(k,a) = jac[a*dim+k];
          const Eigen::Matrix<double,dim,dim> B = J.inverse();
          constexpr double eps = 2.e-5;
          for(int a = 0; a < dim; a++){
            const auto target_value = [&](double step, double *out){
              double sample_bary[dim+1];
              std::copy(bary,bary+dim+1,sample_bary);
              sample_bary[0] -= step;
              sample_bary[a+1] += step;
              msh.met.getMetBary(AsDeg::Pk,DifVar::None,MetSpace::Exp,
                                nodes,dim,sample_bary,out,nullptr);
            };
            double plus[nm], minus[nm], plus_half[nm], minus_half[nm];
            target_value(eps,plus); target_value(-eps,minus);
            target_value(eps/2,plus_half); target_value(-eps/2,minus_half);
            for(int m = 0; m < nm; m++){
              const double coarse = (plus[m]-minus[m])/(2*eps);
              const double fine = (plus_half[m]-minus_half[m])/eps;
              for(int r = 0; r < dim; r++)
                dmet[r*nm+m] += B(a,r)*(4*fine-coarse)/3;
            }
          }
        }
        const double frozen_value = connection_difference_norm2<dim>(jac,hess,met,dmet);
        BOOST_CHECK_SMALL((frozen_value-quafun_connection<MFT,dim>(
            msh,msh.ent2poi(dim)[0],bary)),1.e-9*(1.+std::abs(frozen_value)));
        for(FEBasis dofbas : {FEBasis::Lagrange,FEBasis::Bezier}){
          double translation_gradient[dim] = {};
          for(int ivar = 0; ivar < nn; ivar++){
            BOOST_TEST_CONTEXT("dim=" << dim << " degree=" << degree
                               << " storage=" << int(storage) << " dofbas=" << int(dofbas)
                               << " sample=" << sample << " node=" << ivar){
              double grad[dim];
              const double value = d_quafun_connection<MFT,dim>(
                  msh,msh.ent2poi(dim)[0],bary,ivar,dofbas,grad);
              BOOST_CHECK_EQUAL(value,(quafun_connection<MFT,dim>(msh,msh.ent2poi(dim)[0],bary)));
              const auto direction = perturbation_coefficients<dim,degree>(storage,dofbas,ivar);
              for(int v = 0; v < dim; v++){
                translation_gradient[v] += grad[v];
                std::array<double,nn> coords;
                for(int i = 0; i < nn; i++) coords[i] = msh.coord(i,v);
                const auto displaced_value = [&](double step){
                  for(int i = 0; i < nn; i++) msh.coord(i,v) = coords[i] + step*direction[i];
                  evaluate_test_geometry<dim,degree>(msh,bary,point,jac,hess);
                  return connection_difference_norm2<dim>(jac,hess,met,dmet);
                };
                // Richardson-extrapolated central differences, used ONLY
                // for verification, never by the production evaluator.
                constexpr double eps = 2.e-5;
                const double coarse = (displaced_value(eps)-displaced_value(-eps))/(2*eps);
                const double fine = (displaced_value(eps/2)-displaced_value(-eps/2))/eps;
                for(int i = 0; i < nn; i++) msh.coord(i,v) = coords[i];
                const double numerical = (4*fine-coarse)/3;
                BOOST_CHECK_SMALL(grad[v]-numerical,2.e-8*(1.+std::abs(numerical)));
              }
            }
          }
          // With target data frozen, a common translation leaves the
          // objective unchanged for both analytical and FE targets.
          for(double v : translation_gradient) BOOST_CHECK_SMALL(v,2.e-10);
        }
      }
      if constexpr(std::is_same_v<MFT,MetricFieldFE>){
        for(int i = 0; i < nn; i++){
          for(int m = 0; m < nm; m++) BOOST_CHECK_EQUAL(msh.met(i,m),original_metric[i*nm+m]);
        }
      }
    }
  }
}

} // namespace

template<int dim>
void check_frozen_connection_gradient(){
  constexpr int nn = getnnode(dim,2), nm = dim*(dim+1)/2;
  Mesh<MetricFieldAnalytical> msh;
  MetrisParameters parameters;
  initialize<dim,2>(msh,parameters);
  parameters.setAnalyticalMetric(AnaMetFun(polynomial_metric<dim>));
  msh.met.setAnalyticalMetric(parameters); // No Hessian callback: frozen needs none.
  auto ordering = ORDELT(dim);
  for(int i = 0; i < nn; i++){
    for(int k = 0; k < dim; k++)
      msh.coord(i,k) = .1+.6*ordering[2][i][k+1]/2.+.007*(i+1)*(k+1);
  }
  msh.forceBasisFlag(FEBasis::Bezier);
  const int *nodes = msh.ent2poi(dim)[0];
  double bary[dim+1];
  std::fill(bary,bary+dim+1,1./(dim+1));
  double point[dim], jac[dim*dim], hess[nm*dim], met[nm], dmet[dim*nm];
  const auto geometry = [&](){
    if constexpr(dim == 2)
      eval2<dim,2>(msh.coord,nodes,FEBasis::Bezier,DifVar::Bary,DifVar::Bary,bary,point,jac,hess);
    else
      eval3<dim,2>(msh.coord,nodes,FEBasis::Bezier,DifVar::Bary,DifVar::Bary,bary,point,jac,hess);
  };
  geometry();
  msh.met.getMetPhys(DifVar::Phys,MetSpace::Exp,point,met,dmet);
  for(int ivar = 0; ivar < nn; ivar++){
    double grad[dim];
    const double value = d_quafun_connection<MetricFieldAnalytical,dim>(
        msh,nodes,bary,ivar,FEBasis::Bezier,grad);
    BOOST_CHECK_EQUAL(value,(quafun_connection<MetricFieldAnalytical,dim>(msh,nodes,bary)));
    for(int k = 0; k < dim; k++){
      const double initial = msh.coord(ivar,k), eps = 1.e-6;
      msh.coord(ivar,k) = initial+eps;
      geometry();
      const double plus = connection_difference_norm2<dim>(jac,hess,met,dmet);
      msh.coord(ivar,k) = initial-eps;
      geometry();
      const double minus = connection_difference_norm2<dim>(jac,hess,met,dmet);
      msh.coord(ivar,k) = initial;
      BOOST_CHECK_SMALL(grad[k]-(plus-minus)/(2*eps),2.e-7*(1.+std::abs(grad[k])));
    }
  }
}

BOOST_AUTO_TEST_CASE(frozen_connection_gradient){
  check_frozen_connection_gradient<2>();
  check_frozen_connection_gradient<3>();
}

template<int dim>
void check_frozen_hessian(){
  constexpr int nm=dim*(dim+1)/2;
  double J[dim*dim], H[nm*dim], M[nm]{}, dM[dim*nm];
  double dphi[dim], hphi[nm];
  for(int a=0;a<dim;a++){
    dphi[a]=.2*(a+1);
    M[sym2idx(a,a)]=1.+a;
    for(int k=0;k<dim;k++) J[a*dim+k]=(a==k ? .7 : .04*(a+k+1));
  }
  M[1]=.17;
  for(int a=0;a<nm;a++){
    hphi[a]=.08*(a+1);
    for(int k=0;k<dim;k++) H[a*dim+k]=.03*(a+1)*(k+1);
  }
  for(int i=0;i<dim*nm;i++) dM[i]=.015*(i+1);
  double g[dim], h[nm], dj[dim*dim*dim]{}, dh[dim*nm*dim]{};
  const double val=d2_connection_difference_norm2_frozen<dim>(J,H,M,dM,dphi,hphi,g,h);
  for(int v=0;v<dim;v++){
    for(int a=0;a<dim;a++) dj[v*dim*dim+a*dim+v]=dphi[a];
    for(int a=0;a<nm;a++) dh[v*nm*dim+a*dim+v]=hphi[a];
  }
  double adgrad[dim];
  const double adval=d_connection_difference_norm2<dim>(J,H,M,dM,dj,dh,adgrad);
  BOOST_CHECK_CLOSE_FRACTION(val,adval,1.e-13);
  for(int u=0;u<dim;u++) BOOST_CHECK_SMALL(g[u]-adgrad[u],1.e-11*(1.+std::abs(g[u])));
  constexpr double eps=1.e-6;
  for(int v=0;v<dim;v++){
    double gp[dim],gm[dim];
    const auto displaced=[&](double step,double *out){
      double jt[dim*dim],ht[nm*dim];
      for(int i=0;i<dim*dim;i++) jt[i]=J[i]+step*dj[v*dim*dim+i];
      for(int i=0;i<nm*dim;i++) ht[i]=H[i]+step*dh[v*nm*dim+i];
      d_connection_difference_norm2<dim>(jt,ht,M,dM,dj,dh,out);
    };
    displaced(eps,gp); displaced(-eps,gm);
    for(int u=0;u<dim;u++){
      const double fd=(gp[u]-gm[u])/(2*eps);
      BOOST_CHECK_SMALL(h[sym2idx(u,v)]-fd,1.e-7*(1.+std::abs(fd)));
    }
  }
}

BOOST_AUTO_TEST_CASE(frozen_connection_hessian){
  check_frozen_hessian<2>();
  check_frozen_hessian<3>();
}

BOOST_AUTO_TEST_CASE(integrated_connection_uses_distance_weights_and_frozen_gradient){
  MetrisParameters parameters;
  Mesh<MetricFieldAnalytical> mesh;
  initialize<2,2>(mesh,parameters);
  parameters.setAnalyticalMetric(9);
  parameters.objective_quadrature_order=5;
  mesh.met.setAnalyticalMetric(parameters);
  mesh.forceBasisFlag(FEBasis::Bezier);
  for(int i=0;i<6;i++) for(int k=0;k<2;k++)
    mesh.coord(i,k)=.1+.6*ordfac.s[2][i][k+1]/2.+.004*(i+1)*(k+1);
  const auto rule=get_objective_quadrature<2>(5);
  const auto vf=get_quafun<MetricFieldAnalytical,2,2>(QuaFun::Connection);
  const auto df=get_d_quafun<MetricFieldAnalytical,2,2>(QuaFun::Connection);
  for(int ivar=0;ivar<6;ivar++) for(FEBasis dofbas : {FEBasis::Lagrange,FEBasis::Bezier}){
    double integral=0., g[2]{}, gradient[2], hessian[3], gradientOnly[2];
    for(int q=0;q<rule.size();q++){
      const auto sample=prepare_objective_quadrature_sample<MetricFieldAnalytical,2,2,2>(
          mesh,AsDeg::Pk,mesh.fac2poi[0],rule[q],standard_objective_quadrature_theta_mode());
      double pointGradient[2];
      const double value=d_quafun_connection<MetricFieldAnalytical,2>(
          mesh,mesh.fac2poi[0],rule[q].bary,ivar,dofbas,pointGradient);
      const double w=rule[q].weight*sample.theta;
      integral+=w*value;
      for(int k=0;k<2;k++) g[k]+=w*pointGradient[k];
    }
    BOOST_CHECK_CLOSE_FRACTION(vf(mesh,AsDeg::Pk,AsDeg::Pk,0,1.),integral,1.e-13);
    BOOST_CHECK_CLOSE_FRACTION(df(mesh,AsDeg::Pk,AsDeg::Pk,0,ivar,dofbas,DifVar::None,gradient,hessian,1.),integral,1.e-13);
    df(mesh,AsDeg::Pk,AsDeg::Pk,0,ivar,dofbas,DifVar::None,gradientOnly,nullptr,1.);
    for(int k=0;k<2;k++){
      BOOST_CHECK_SMALL(gradient[k]-g[k],1.e-10*(1.+std::abs(g[k])));
      BOOST_CHECK_SMALL(gradientOnly[k]-g[k],1.e-10*(1.+std::abs(g[k])));
    }
  }
}

BOOST_AUTO_TEST_CASE(manufactured_triangles){
  check_manufactured<2,1>();
  check_manufactured<2,2>();
}
BOOST_AUTO_TEST_CASE(manufactured_tetrahedra){
  check_manufactured<3,1>();
  check_manufactured<3,2>();
}
BOOST_AUTO_TEST_CASE(fe_log_euclidean_derivatives){
  check_fe_derivatives<2,1>();
  check_fe_derivatives<3,1>();
  check_fe_derivatives<2,2>();
  check_fe_derivatives<3,2>();
}
BOOST_AUTO_TEST_CASE(affine_coordinate_invariance_and_target_scaling){
  check_coordinate_invariance<2>();
  check_coordinate_invariance<3>();
}
BOOST_AUTO_TEST_CASE(affine_constant_and_invalid_input){
  double jac[4] = {1.,.2,.4,1.};
  double hess[6] = {}, met[3] = {2.,.3,1.}, dmet[6] = {};
  BOOST_CHECK_EQUAL(connection_difference_norm2<2>(jac,hess,met,dmet),0.);
  met[0] = -1.;
  BOOST_CHECK_THROW(connection_difference_norm2<2>(jac,hess,met,dmet),MetrisExcept);
  met[0] = 2.;
  jac[0] = jac[1] = 0.;
  BOOST_CHECK_THROW(connection_difference_norm2<2>(jac,hess,met,dmet),MetrisExcept);
  jac[0] = 1.;
  dmet[0] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_THROW(connection_difference_norm2<2>(jac,hess,met,dmet),MetrisExcept);
}

BOOST_AUTO_TEST_CASE(analytical_frozen_control_point_gradient){
  check_frozen_control_gradients<2,1,MetricFieldAnalytical>();
  check_frozen_control_gradients<2,2,MetricFieldAnalytical>();
  check_frozen_control_gradients<3,1,MetricFieldAnalytical>();
  check_frozen_control_gradients<3,2,MetricFieldAnalytical>();
}
BOOST_AUTO_TEST_CASE(fe_frozen_control_point_gradient){
  check_frozen_control_gradients<2,1,MetricFieldFE>();
  check_frozen_control_gradients<2,2,MetricFieldFE>();
  check_frozen_control_gradients<3,1,MetricFieldFE>();
  check_frozen_control_gradients<3,2,MetricFieldFE>();
}

BOOST_AUTO_TEST_CASE(frozen_gradient_rejects_invalid_node){
  Mesh<MetricFieldAnalytical> msh;
  MetrisParameters parameters;
  initialize<2,1>(msh,parameters);
  parameters.setAnalyticalMetric(AnaMetFun(polynomial_metric<2>));
  msh.met.setAnalyticalMetric(parameters);
  for(int i = 0; i < 3; i++){
    for(int k = 0; k < 2; k++) msh.coord(i,k) = i == k+1 ? 1. : 0.;
  }
  const double bary[3] = {.2,.3,.5};
  double grad[2];
  for(int ivar : {-1,3}){
    BOOST_CHECK_THROW((d_quafun_connection<MetricFieldAnalytical,2>(
        msh,msh.ent2poi(2)[0],bary,ivar,FEBasis::Lagrange,grad)),MetrisExcept);
  }
  BOOST_CHECK_NO_THROW((d_quafun_connection<MetricFieldAnalytical,2>(
      msh,msh.ent2poi(2)[0],bary,0,FEBasis::Lagrange,grad)));
  for(int metric : {1,9}){
    parameters.setAnalyticalMetric(metric);
    msh.met.setAnalyticalMetric(parameters);
    BOOST_CHECK_NO_THROW((d_quafun_connection<MetricFieldAnalytical,2>(
        msh,msh.ent2poi(2)[0],bary,0,FEBasis::Lagrange,grad)));
    // Q1 geometry has zero connection and zero frozen connection gradient.
    for(double g : grad) BOOST_CHECK_SMALL(g,1.e-14);
  }
}

} // namespace Metris
