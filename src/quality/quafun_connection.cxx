// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#include "quafun_connection.hxx"
#include "objective_quadrature_sample.hxx"

#include "../Mesh/Mesh.hxx"
#include "../aux_exceptions.hxx"
#include "../codegen_lag2bez.hxx"
#include "../linalg/symidx.hxx"
#include "../low_eval.hxx"
#include "../libs/SANS/Surreal/SurrealS.h"
#include "../utils/CT_loop.hxx"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>

namespace Metris{

template<int gdim>
double connection_difference_norm2(const double *jmat, const double *hmat,
                                   const double *met, const double *dmet){
  static_assert(gdim == 2 || gdim == 3);
  constexpr int nmet = gdim*(gdim+1)/2;
  using Matrix = Eigen::Matrix<double,gdim,gdim>;
  Matrix jac, metric;
  for(int i = 0; i < gdim; i++){
    for(int j = 0; j < gdim; j++){
      jac(i,j) = jmat[j*gdim+i];
      metric(i,j) = met[sym2idx(i,j)];
    }
  }
  METRIS_ENFORCE_MSG(jac.allFinite() && metric.allFinite(),
                    "Nonfinite Jacobian or metric in connection objective");
  for(int i = 0; i < gdim*nmet; i++){
    METRIS_ENFORCE_MSG(std::isfinite(hmat[i]) && std::isfinite(dmet[i]),
                      "Nonfinite derivative in connection objective");
  }
  // J need not be symmetric or positive definite. Pivoted LU checks its
  // numerical invertibility and supplies B = J^-1 for the mesh connection.
  const Eigen::FullPivLU<Matrix> jacLU(jac);
  METRIS_ENFORCE_MSG(jacLU.isInvertible(),
                    "Singular Jacobian in connection objective");
  const Matrix jinv = jacLU.inverse();
  // M must be SPD. Cholesky supplies M^-1 for the target connection and
  // M = U^T U for evaluating the target-metric tensor norm below.
  const Eigen::LLT<Matrix> cholM(metric);
  METRIS_ENFORCE_MSG(cholM.info() == Eigen::Success,
                    "Target metric is not SPD in connection objective");
  const Matrix minv = cholM.solve(Matrix::Identity());

  // deltaChristoffelTensor[k][i][j] = (Gamma_K)^k_ij - (Gamma_M)^k_ij.
  double deltaChristoffelTensor[gdim][gdim][gdim] = {};
  for(int k = 0; k < gdim; k++){
    for(int i = 0; i < gdim; i++){
      for(int j = 0; j < gdim; j++){
        // Add (Gamma_K)^k_ij = -sum_{a,b} (d_a d_b F^k) B^a_i B^b_j.
        for(int a = 0; a < gdim; a++){
          for(int b = 0; b < gdim; b++){
            deltaChristoffelTensor[k][i][j] -=
                hmat[sym2idx(a,b)*gdim+k]*jinv(a,i)*jinv(b,j);
          }
        }
        // Subtract (Gamma_M)^k_ij = (1/2) sum_l M^{kl}
        //                         (d_i M_lj + d_j M_li - d_l M_ij).
        for(int l = 0; l < gdim; l++){
          deltaChristoffelTensor[k][i][j] -= 0.5*minv(k,l)
              *(dmet[i*nmet+sym2idx(l,j)]
                + dmet[j*nmet+sym2idx(l,i)]
                - dmet[l*nmet+sym2idx(i,j)]);
        }
      }
    }
  }

  // Compute the squared tensor norm defined by the target metric.

  // Let U = upperM and E = inverseUpperM. Then M = U^T U,
  // M^-1 = E E^T, and the columns of E form an M-orthonormal basis.
  // For each (a,b,c), compute the tensor component in that basis:
  //   R_abc = sum_{k,i,j} U_ak DeltaGamma^k_ij E_ib E_jc.

  // Thus sum_{a,b,c} R_abc^2 is the squared tensor norm.
  // Summing squares avoids cancellation in the direct contraction and
  // gives a nonnegative result without forcing negative values to zero.

  const Matrix upperM = cholM.matrixU();
  const Matrix inverseUpperM = upperM.template triangularView<Eigen::Upper>()
                                    .solve(Matrix::Identity());
  double value = 0.;
  for(int a = 0; a < gdim; a++){
    for(int b = 0; b < gdim; b++){
      for(int c = 0; c < gdim; c++){
        double orthonormalComponent = 0.;
        for(int k = 0; k < gdim; k++){
          for(int i = 0; i < gdim; i++){
            for(int j = 0; j < gdim; j++){
              orthonormalComponent += upperM(a,k)*deltaChristoffelTensor[k][i][j]
                                      *inverseUpperM(i,b)*inverseUpperM(j,c);
            }
          }
        }
        value += orthonormalComponent*orthonormalComponent;
      }
    }
  }
  METRIS_ENFORCE_MSG(std::isfinite(value),
                    "Nonfinite connection objective");
  return value;
}

namespace{ // Implementation helpers private to this translation unit.

// Evaluate the polynomial sum_A N_A c_A and its first two reference derivatives.
// Used for F_K (ncomponent = gdim) and N_A (ncomponent = 1, with one
// coefficient set to one and the rest to zero).
// Derivative layouts: referenceGradient[a*ncomponent+k] and
// referenceHessian[sym2idx(a,b)*ncomponent+k].
template<int ncomponent, int gdim, int ideg>
void evaluate_reference_polynomial(
    const dblAr2 &coefficients, const int *coefficientIndices, FEBasis basis,
    const double *bary, double *value,
    double *referenceGradient, double *referenceHessian){
  static_assert(gdim == 2 || gdim == 3);
  METRIS_ENFORCE(basis == FEBasis::Lagrange || basis == FEBasis::Bezier);
  const auto evaluateBasis = [&](const dblAr2 &basisCoefficients,
                                  const int *indices, FEBasis evaluationBasis){
    if constexpr(gdim == 2){
      eval2<ncomponent,ideg>(basisCoefficients,indices,evaluationBasis,
                            DifVar::Bary,DifVar::Bary,
                            bary,value,referenceGradient,referenceHessian);
    }else{
      eval3<ncomponent,ideg>(basisCoefficients,indices,evaluationBasis,
                            DifVar::Bary,DifVar::Bary,
                            bary,value,referenceGradient,referenceHessian);
    }
  };
  if constexpr(ideg > 1){
    if(basis == FEBasis::Lagrange){
      // Higher-order eval2/eval3 provide second derivatives in Bernstein
      // basis only, so convert a local copy of the Lagrange coefficients.
      constexpr int nnode = getnnode(gdim,ideg);
      std::array<double,nnode*ncomponent> lag_data, bez_data;
      std::array<int,nnode> local_nodes;
      dblAr2 lag(nnode,ncomponent,lag_data.data());
      dblAr2 bez(nnode,ncomponent,bez_data.data());
      for(int i = 0; i < nnode; i++){
        local_nodes[i] = i;
        for(int k = 0; k < ncomponent; k++)
          lag(i,k) = coefficients(coefficientIndices[i],k);
      }
      if constexpr(gdim == 2){
        lag2bez2<ideg,ncomponent>(local_nodes.data(),lag,bez);
      }else{
        lag2bez3<ideg,ncomponent>(local_nodes.data(),lag,bez);
      }
      evaluateBasis(bez,local_nodes.data(),FEBasis::Bezier);
      return;
    }
  }
  // At degree one, Lagrange and Bernstein are the same
  evaluateBasis(coefficients,coefficientIndices,basis);
}

// Evaluate the SPD target metric M and its first physical derivatives at
// point = F(bary). Analytical targets use point; FE targets use the owning
// entity when supplied.
// met[sym2idx(i,j)] = M_ij; dmet[r*nm+sym2idx(i,j)] = partial_r M_ij.
template<class MFT, int gdim>
void evaluate_target_metric_and_derivatives(
    Mesh<MFT> &msh, const int *nodes, const double *bary, const double *point,
    double *met, double *dmet, int ientt){
  constexpr int nm = gdim*(gdim+1)/2;
  for(int i = 0; i < nm; i++) met[i] = std::numeric_limits<double>::quiet_NaN();
  for(int i = 0; i < gdim*nm; i++) dmet[i] = std::numeric_limits<double>::quiet_NaN();
  if constexpr(std::is_same_v<MFT,MetricFieldAnalytical>){
    msh.met.getMetFullinfo(AsDeg::Pk,DifVar::Phys,MetSpace::Exp,
                           nodes,gdim,bary,point,met,dmet);
  }else if(ientt >= 0){
    msh.met.getMetOnEntity(ientt,DifVar::Phys,bary,met,dmet);
  }else{
    // With L = log M, partial_r M = Dexp_L[partial_r L]. Evaluate this
    // Frechet derivative using exponential divided differences, including
    // their continuous limit at repeated eigenvalues.
    double logm[nm], dlogm[gdim*nm];
    msh.met.getMetFullinfo(AsDeg::Pk,DifVar::Phys,MetSpace::Log,
                           nodes,gdim,bary,point,logm,dlogm);
    getexpmet_frechet<gdim>(logm,dlogm,gdim,met,dmet);
  }
}

} // namespace

template<class MFT, int gdim>
static double quafun_connection_impl(Mesh<MFT> &msh, const int *ent2poi,
                                      const double *bary, int ientt){
  static_assert(gdim == 2 || gdim == 3);
  METRIS_ENFORCE_MSG(msh.idim == gdim,
                    "Connection objective requires full-dimensional elements");
  METRIS_ENFORCE_MSG(msh.curdeg >= 1 && msh.curdeg <= METRIS_MAX_DEG,
                    "Unsupported geometry degree in connection objective");
  METRIS_ENFORCE(msh.getBasis() == FEBasis::Lagrange
                 || msh.getBasis() == FEBasis::Bezier);
  constexpr int nmet = gdim*(gdim+1)/2;
  double point[gdim], jac[gdim*gdim], hess[nmet*gdim];
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){
    if(msh.curdeg == ideg){
      evaluate_reference_polynomial<gdim,gdim,ideg>(msh.coord,ent2poi,msh.getBasis(),
                                      bary,point,jac,hess);
    }
  }CT_FOR1(ideg);
  double metric[nmet], dmetric[gdim*nmet];
  evaluate_target_metric_and_derivatives<MFT,gdim>(msh,ent2poi,bary,point,metric,dmetric,ientt);
  return connection_difference_norm2<gdim>(jac,hess,metric,dmetric);
}

template<class MFT, int gdim>
double quafun_connection(Mesh<MFT> &msh, const int *ent2poi, const double *bary){
  return quafun_connection_impl<MFT,gdim>(msh,ent2poi,bary,-1);
}

template<class MFT, int gdim>
double quafun_connection(Mesh<MFT> &msh, int ientt, const double *bary){
  METRIS_ENFORCE(ientt >= 0 && ientt < msh.nentt(gdim));
  METRIS_ENFORCE(msh.ent2poi(gdim)(ientt,0) >= 0);
  return quafun_connection_impl<MFT,gdim>(msh,msh.ent2poi(gdim)[ientt],bary,ientt);
}

template<int gdim>
double d_connection_difference_norm2(
    const double *jmat, const double *hmat,
    const double *met, const double *dmet,
    const double *d_jmat, const double *d_hmat,
    double *gradient){
  static_assert(gdim == 2 || gdim == 3);
  constexpr int nm = gdim*(gdim+1)/2;
  using AD = SANS::SurrealS<gdim,double>;
  using Matrix = Eigen::Matrix<double,gdim,gdim>;
  METRIS_ENFORCE(gradient != nullptr);
  // Use precisely the same primal evaluation and validity checks as value-only.
  const double value = connection_difference_norm2<gdim>(jmat,hmat,met,dmet);
  const auto seed = [](double primal, const double *derivatives,
                        int stride, int index){
    AD result = primal;
    for(int v = 0; v < gdim; v++){
      const double derivative = derivatives[v*stride+index];
      METRIS_ENFORCE_MSG(std::isfinite(derivative),
                        "Nonfinite input variation in connection gradient");
      result.deriv(v) = derivative;
    }
    return result;
  };
  Matrix jac, metric;
  AD Ji[gdim][gdim], H[nm*gdim];
  AD dJ[gdim][gdim];
  for(int i = 0; i < gdim; i++){
    for(int j = 0; j < gdim; j++){
      jac(i,j) = jmat[j*gdim+i];
      metric(i,j) = met[sym2idx(i,j)];
      dJ[i][j] = seed(jac(i,j),d_jmat,gdim*gdim,j*gdim+i);
    }
  }
  for(int i = 0; i < gdim*nm; i++){
    H[i] = seed(hmat[i],d_hmat,gdim*nm,i);
  }
  // Only geometry varies. Differentiate J^-1 analytically so the pivot/rank
  // decisions are not differentiated; all target-metric factors stay fixed.
  const Eigen::FullPivLU<Matrix> jacLU(jac);
  const Matrix ji = jacLU.inverse();
  const Eigen::LLT<Matrix> cholM(metric);
  const Matrix mi = cholM.solve(Matrix::Identity());
  const Matrix upperM = cholM.matrixU();
  const Matrix inverseUpperM = upperM.template triangularView<Eigen::Upper>()
                                    .solve(Matrix::Identity());
  for(int i = 0; i < gdim; i++){
    for(int j = 0; j < gdim; j++){
      Ji[i][j] = ji(i,j);
      for(int v = 0; v < gdim; v++){
        double dj = 0.;
        for(int a = 0; a < gdim; a++){
          for(int b = 0; b < gdim; b++){
            dj -= ji(i,a)*dJ[a][b].deriv(v)*ji(b,j);
          }
        }
        Ji[i][j].deriv(v) = dj;
      }
    }
  }
  AD T[gdim][gdim][gdim];
  for(int k = 0; k < gdim; k++){
    for(int i = 0; i < gdim; i++){
      for(int j = 0; j < gdim; j++){
        T[k][i][j] = 0.;
        for(int a = 0; a < gdim; a++){
          for(int b = 0; b < gdim; b++){
            T[k][i][j] -= H[sym2idx(a,b)*gdim+k]*Ji[a][i]*Ji[b][j];
          }
        }
        for(int l = 0; l < gdim; l++){
          T[k][i][j] -= 0.5*mi(k,l)
              *(dmet[i*nm+sym2idx(l,j)] + dmet[j*nm+sym2idx(l,i)]
                - dmet[l*nm+sym2idx(i,j)]);
        }
      }
    }
  }
  // Transform the varying connection difference using the fixed metric
  // factors, then differentiate its sum of squared components.
  AD norm2 = 0.;
  for(int a = 0; a < gdim; a++){
    for(int b = 0; b < gdim; b++){
      for(int c = 0; c < gdim; c++){
        AD component = 0.;
        for(int k = 0; k < gdim; k++){
          for(int i = 0; i < gdim; i++){
            for(int j = 0; j < gdim; j++){
              component += upperM(a,k)*T[k][i][j]*inverseUpperM(i,b)*inverseUpperM(j,c);
            }
          }
        }
        norm2 += component*component;
      }
    }
  }
  for(int v = 0; v < gdim; v++){
    METRIS_ENFORCE_MSG(std::isfinite(norm2.deriv(v)),
                      "Nonfinite connection gradient");
    gradient[v] = norm2.deriv(v);
  }
  return value;
}

template<class MFT, int gdim>
static double d_quafun_connection_impl(Mesh<MFT> &msh, const int *ent2poi,
                           const double *bary, int ivar, FEBasis dofbas,
                           double *gradient, int ientt){
  static_assert(gdim == 2 || gdim == 3);
  static_assert(std::is_same_v<MFT,MetricFieldFE>
                 || std::is_same_v<MFT,MetricFieldAnalytical>);
  METRIS_ENFORCE(msh.idim == gdim);
  METRIS_ENFORCE(msh.curdeg >= 1 && msh.curdeg <= METRIS_MAX_DEG);
  METRIS_ENFORCE(ivar >= 0 && ivar < getnnode(gdim,msh.curdeg));
  METRIS_ENFORCE(msh.getBasis() == FEBasis::Lagrange
                 || msh.getBasis() == FEBasis::Bezier);
  METRIS_ENFORCE(dofbas == FEBasis::Lagrange || dofbas == FEBasis::Bezier);
  METRIS_ENFORCE(gradient != nullptr);
  constexpr int nm = gdim*(gdim+1)/2;
  double point[gdim], jac[gdim*gdim], hess[nm*gdim];
  double phi, dphi[gdim], hphi[nm];
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){
    if(msh.curdeg == ideg){
      evaluate_reference_polynomial<gdim,gdim,ideg>(msh.coord,ent2poi,msh.getBasis(),
                                      bary,point,jac,hess);
      constexpr int nn = getnnode(gdim,ideg);
      std::array<double,nn> coefficients{};
      std::array<int,nn> nodes;
      for(int i = 0; i < nn; i++) nodes[i] = i;
      coefficients[ivar] = 1.;
      const dblAr2 shapeCoefficients(nn,1,coefficients.data());
      evaluate_reference_polynomial<1,gdim,ideg>(shapeCoefficients,nodes.data(),dofbas,
                                   bary,&phi,dphi,hphi);
    }
  }CT_FOR1(ideg);
  double met[nm], dm[gdim*nm];
  evaluate_target_metric_and_derivatives<MFT,gdim>(msh,ent2poi,bary,point,met,dm,ientt);
  double d_jac[gdim*gdim*gdim] = {}, d_hess[gdim*nm*gdim] = {};
  for(int v = 0; v < gdim; v++){
    for(int a = 0; a < gdim; a++) d_jac[v*gdim*gdim+a*gdim+v] = dphi[a];
    for(int ab = 0; ab < nm; ab++) d_hess[v*nm*gdim+ab*gdim+v] = hphi[ab];
  }
  return d_connection_difference_norm2<gdim>(jac,hess,met,dm,
                                             d_jac,d_hess,gradient);
}

template<class MFT, int gdim>
double d_quafun_connection(Mesh<MFT> &msh, const int *ent2poi,
                           const double *bary, int ivar, FEBasis dofbas,
                           double *gradient){
  return d_quafun_connection_impl<MFT,gdim>(msh,ent2poi,bary,ivar,dofbas,gradient,-1);
}

template<class MFT, int gdim>
double d_quafun_connection(Mesh<MFT> &msh, int ientt,
                           const double *bary, int ivar, FEBasis dofbas,
                           double *gradient){
  METRIS_ENFORCE(ientt >= 0 && ientt < msh.nentt(gdim));
  METRIS_ENFORCE(msh.ent2poi(gdim)(ientt,0) >= 0);
  return d_quafun_connection_impl<MFT,gdim>(msh,msh.ent2poi(gdim)[ientt],
                                           bary,ivar,dofbas,gradient,ientt);
}

template<int dim>
double d2_connection_difference_norm2_frozen(
    const double *jac, const double *hess, const double *met, const double *dm,
    const double *dphi, const double *hphi, double *gradient, double *hessian){
  constexpr int nm = dim*(dim+1)/2;
  using Matrix = Eigen::Matrix<double,dim,dim>;
  METRIS_ENFORCE(gradient && hessian);
  const double value = connection_difference_norm2<dim>(jac,hess,met,dm);
  Matrix J, M;
  for(int i=0;i<dim;i++) for(int j=0;j<dim;j++){
    J(i,j)=jac[j*dim+i]; M(i,j)=met[sym2idx(i,j)];
  }
  const Matrix P=J.inverse();
  const Eigen::LLT<Matrix> cholM(M);
  const Matrix minv=cholM.solve(Matrix::Identity());
  const Matrix upperM=cholM.matrixU();
  const Matrix inverseUpperM=upperM.template triangularView<Eigen::Upper>().solve(Matrix::Identity());
  Matrix dJ[dim], dP[dim], ddP[dim][dim];
  for(int v=0;v<dim;v++){
    dJ[v].setZero();
    for(int a=0;a<dim;a++) dJ[v](v,a)=dphi[a];
    dP[v]=-P*dJ[v]*P;
  }
  for(int u=0;u<dim;u++) for(int v=0;v<dim;v++)
    ddP[u][v]=P*dJ[u]*P*dJ[v]*P+P*dJ[v]*P*dJ[u]*P;
  constexpr int nt=dim*dim*dim;
  std::array<double,nt> T{};
  std::array<double,nt> dT[dim]{}, ddT[nm]{};
  for(int k=0;k<dim;k++) for(int i=0;i<dim;i++) for(int j=0;j<dim;j++){
    const int index=(k*dim+i)*dim+j;
    for(int l=0;l<dim;l++)
      T[index]-=.5*minv(k,l)*(dm[i*nm+sym2idx(l,j)]+dm[j*nm+sym2idx(l,i)]-dm[l*nm+sym2idx(i,j)]);
    for(int a=0;a<dim;a++) for(int b=0;b<dim;b++){
      const double H=hess[sym2idx(a,b)*dim+k];
      T[index]-=H*P(a,i)*P(b,j);
      for(int u=0;u<dim;u++){
        const double Hu=k==u ? hphi[sym2idx(a,b)] : 0.;
        dT[u][index]-=Hu*P(a,i)*P(b,j)+H*(dP[u](a,i)*P(b,j)+P(a,i)*dP[u](b,j));
        for(int v=0;v<=u;v++){
          const double Hv=k==v ? hphi[sym2idx(a,b)] : 0.;
          ddT[sym2idx(u,v)][index]-=
              Hu*(dP[v](a,i)*P(b,j)+P(a,i)*dP[v](b,j))
             +Hv*(dP[u](a,i)*P(b,j)+P(a,i)*dP[u](b,j))
             +H*(ddP[u][v](a,i)*P(b,j)+dP[u](a,i)*dP[v](b,j)
                 +dP[v](a,i)*dP[u](b,j)+P(a,i)*ddP[u][v](b,j));
        }
      }
    }
  }
  const auto whiten=[&](const std::array<double,nt>& tensor){
    std::array<double,nt> result{};
    for(int a=0;a<dim;a++) for(int b=0;b<dim;b++) for(int c=0;c<dim;c++)
      for(int k=0;k<dim;k++) for(int i=0;i<dim;i++) for(int j=0;j<dim;j++)
        result[(a*dim+b)*dim+c]+=upperM(a,k)*tensor[(k*dim+i)*dim+j]*inverseUpperM(i,b)*inverseUpperM(j,c);
    return result;
  };
  const auto R=whiten(T);
  std::array<double,nt> dR[dim];
  for(int u=0;u<dim;u++){
    dR[u]=whiten(dT[u]); gradient[u]=0.;
    for(int k=0;k<nt;k++) gradient[u]+=2.*R[k]*dR[u][k];
    METRIS_ENFORCE(std::isfinite(gradient[u]));
  }
  for(int u=0;u<dim;u++) for(int v=0;v<=u;v++){
    const int uv=sym2idx(u,v);
    const auto ddR=whiten(ddT[uv]);
    hessian[uv]=0.;
    for(int k=0;k<nt;k++) hessian[uv]+=2.*(dR[u][k]*dR[v][k]+R[k]*ddR[k]);
    METRIS_ENFORCE(std::isfinite(hessian[uv]));
  }
  return value;
}

template<class MFT, int dim>
double d_metqua_connection(Mesh<MFT>& msh, AsDeg asdmsh, AsDeg asdmet,
                           int ientt, int ivar, FEBasis dofbas, DifVar idifmet,
                           double *gradient, double *hessian, double){
  METRIS_ENFORCE(asdmsh==AsDeg::Pk && asdmet==AsDeg::Pk && idifmet==DifVar::None);
  METRIS_ENFORCE(msh.idim==dim && msh.get_tdim()==dim);
  METRIS_ENFORCE(msh.curdeg>=1 && msh.curdeg<=METRIS_MAX_DEG);
  METRIS_ENFORCE(!hessian || gradient);
  constexpr int nm=dim*(dim+1)/2;
  if(gradient){
    METRIS_ENFORCE(ivar>=0 && ivar<getnnode(dim,msh.curdeg));
    METRIS_ENFORCE(dofbas==FEBasis::Lagrange || dofbas==FEBasis::Bezier);
    std::fill(gradient,gradient+dim,0.);
  }
  if(hessian) std::fill(hessian,hessian+nm,0.);
  METRIS_ENFORCE(ientt >= 0 && ientt < msh.nentt(dim));
  METRIS_ENFORCE(msh.ent2poi(dim)(ientt,0) >= 0);
  const int *nodes=msh.ent2poi(dim)[ientt];
  const auto rule=get_objective_quadrature<dim>(msh.param->objective_quadrature_order);
  double integral=0.;
  CT_FOR0_INC(1,METRIS_MAX_DEG,degree){
    if(msh.curdeg==degree){
      for(int q=0;q<rule.size();q++){
        const auto sample=prepare_objective_quadrature_sample<MFT,dim,dim,degree>(
            msh,AsDeg::Pk,ientt,rule[q],standard_objective_quadrature_theta_mode());
        METRIS_ENFORCE(sample.theta_is_valid);
        double value, g[dim], h[nm];
        if(hessian){
          double point[dim], J[dim*dim], H[nm*dim], phi, dphi[dim], hphi[nm];
          evaluate_reference_polynomial<dim,dim,degree>(msh.coord,nodes,msh.getBasis(),rule[q].bary,point,J,H);
          constexpr int nn=getnnode(dim,degree);
          std::array<double,nn> coeff{};
          std::array<int,nn> local;
          for(int i=0;i<nn;i++) local[i]=i;
          coeff[ivar]=1.;
          const dblAr2 shapeCoefficients(nn,1,coeff.data());
          evaluate_reference_polynomial<1,dim,degree>(shapeCoefficients,local.data(),dofbas,rule[q].bary,&phi,dphi,hphi);
          double M[nm], dM[dim*nm];
          evaluate_target_metric_and_derivatives<MFT,dim>(msh,nodes,rule[q].bary,point,M,dM,ientt);
          value=d2_connection_difference_norm2_frozen<dim>(J,H,M,dM,dphi,hphi,g,h);
        }else if(gradient){
          value=d_quafun_connection<MFT,dim>(msh,ientt,rule[q].bary,ivar,dofbas,g);
        }else value=quafun_connection<MFT,dim>(msh,ientt,rule[q].bary);
        const double weight=rule[q].weight*sample.theta;
        integral+=weight*value;
        if(gradient) for(int i=0;i<dim;i++) gradient[i]+=weight*g[i];
        if(hessian) for(int i=0;i<nm;i++) hessian[i]+=weight*h[i];
      }
    }
  }CT_FOR1(degree);
  METRIS_ENFORCE(std::isfinite(integral));
  return integral;
}

template<class MFT, int dim>
double metqua_connection(Mesh<MFT>& msh, AsDeg asdmsh, AsDeg asdmet,
                         int ientt, double difto){
  return d_metqua_connection<MFT,dim>(msh,asdmsh,asdmet,ientt,-1,
      msh.getBasis(),DifVar::None,nullptr,nullptr,difto);
}

#define INSTANTIATE_CONNECTION_INTEGRAL(MFT,dim) \
template double metqua_connection<MFT,dim>(Mesh<MFT>&,AsDeg,AsDeg,int,double); \
template double d_metqua_connection<MFT,dim>(Mesh<MFT>&,AsDeg,AsDeg,int,int,FEBasis,DifVar,double*,double*,double);
INSTANTIATE_CONNECTION_INTEGRAL(MetricFieldAnalytical,2)
INSTANTIATE_CONNECTION_INTEGRAL(MetricFieldAnalytical,3)
INSTANTIATE_CONNECTION_INTEGRAL(MetricFieldFE,2)
INSTANTIATE_CONNECTION_INTEGRAL(MetricFieldFE,3)
#undef INSTANTIATE_CONNECTION_INTEGRAL

#define INSTANTIATE_CONNECTION_ENTITY(MFT,dim) \
template double quafun_connection<MFT,dim>(Mesh<MFT>&,int,const double*); \
template double d_quafun_connection<MFT,dim>(Mesh<MFT>&,int,const double*,int,FEBasis,double*);
INSTANTIATE_CONNECTION_ENTITY(MetricFieldFE,2)
INSTANTIATE_CONNECTION_ENTITY(MetricFieldFE,3)
INSTANTIATE_CONNECTION_ENTITY(MetricFieldAnalytical,2)
INSTANTIATE_CONNECTION_ENTITY(MetricFieldAnalytical,3)
#undef INSTANTIATE_CONNECTION_ENTITY
template double d2_connection_difference_norm2_frozen<2>(const double*,const double*,const double*,const double*,const double*,const double*,double*,double*);
template double d2_connection_difference_norm2_frozen<3>(const double*,const double*,const double*,const double*,const double*,const double*,double*,double*);

#define INSTANTIATE_CONNECTION_GRADIENT(dim) \
template double d_connection_difference_norm2<dim>( \
    const double*,const double*,const double*,const double*, \
    const double*,const double*,double*); \
template double d_quafun_connection<MetricFieldFE,dim>( \
    Mesh<MetricFieldFE>&,const int*,const double*,int,FEBasis,double*); \
template double d_quafun_connection<MetricFieldAnalytical,dim>( \
    Mesh<MetricFieldAnalytical>&,const int*,const double*,int,FEBasis,double*);
INSTANTIATE_CONNECTION_GRADIENT(2)
INSTANTIATE_CONNECTION_GRADIENT(3)
#undef INSTANTIATE_CONNECTION_GRADIENT

template double connection_difference_norm2<2>(const double*,const double*,
                                               const double*,const double*);
template double connection_difference_norm2<3>(const double*,const double*,
                                               const double*,const double*);
template double quafun_connection<MetricFieldFE,2>(
    Mesh<MetricFieldFE>&,const int*,const double*);
template double quafun_connection<MetricFieldFE,3>(
    Mesh<MetricFieldFE>&,const int*,const double*);
template double quafun_connection<MetricFieldAnalytical,2>(
    Mesh<MetricFieldAnalytical>&,const int*,const double*);
template double quafun_connection<MetricFieldAnalytical,3>(
    Mesh<MetricFieldAnalytical>&,const int*,const double*);

} // namespace Metris
