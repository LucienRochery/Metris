// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_metric_discontinuous_evaluation
#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "codegen_lag2bez.hxx"
#include "quality/objective_quadrature_sample.hxx"
#include "quality/quafun_connection.hxx"

#include <cmath>
#include <limits>
#include <vector>

namespace Metris{
namespace{

constexpr double a = 0.3, b = 0.2, c = -0.15;
constexpr double points[9][2] = {
  {0.,0.},{1.,0.},{0.,1.},{1.,1.},
  {0.5,0.5},{0.,0.5},{0.5,0.},{0.5,1.},{1.,0.5}
};
constexpr double offsets[2][3] = {{0.2,0.1,-0.3},{-0.4,0.3,0.6}};
constexpr double dx[2][3] = {{0.4,-0.2,0.7},{0.3,0.5,-0.1}};
constexpr double dy[2][3] = {{-0.5,0.3,0.2},{0.8,-0.4,0.9}};

struct CurvedTriangles{
  MetrisParameters parameters;
  Mesh<MetricFieldFE> mesh;

  CurvedTriangles(){
    parameters.iverb = 0;
    mesh.param = &parameters;
    mesh.idim = 2;
    mesh.curdeg = mesh.strdeg = 2;
    mesh.forceBasisFlag(FEBasis::Lagrange);
    mesh.set_npoin(9);
    mesh.set_nentt(2,2);
    mesh.met.forceBasisFlag(FEBasis::Lagrange);
    mesh.met.forceSpaceFlag(MetSpace::Log);
    const int nodes[2][6] = {{0,1,2,4,5,6},{1,3,2,7,4,8}};
    for(int ientt = 0; ientt < 2; ientt++)
      for(int node = 0; node < 6; node++) mesh.fac2poi(ientt,node) = nodes[ientt][node];
    for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 6; node++){
      const int point = mesh.fac2poi(ientt,node);
      if(node < 3) mesh.set_poi2ent(Vertex{point},2,ientt);
      else mesh.set_poi2ent(CtrlPt{point},2,ientt);
    }
    for(int point = 0; point < 9; point++){
      const double xi = points[point][0], eta = points[point][1];
      mesh.coord(point,0) = xi + a*xi*xi + b*eta;
      mesh.coord(point,1) = eta + c*xi;
      for(int k = 0; k < 3; k++) mesh.met(point,k) = 0.1*(k+1)*(1.+xi*xi-eta);
    }
  }

  void initializeDG(){
    mesh.met.initializeDiscontinuous();
    for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++){
      const int point = mesh.fac2poi(ientt,node);
      for(int k = 0; k < 3; k++)
        mesh.met.elementCoefficient(ientt,node)[k]
            = offsets[ientt][k] + dx[ientt][k]*points[point][0] + dy[ientt][k]*points[point][1];
    }
  }

  void poisonInactivePointBuffer(){
    for(int point = 0; point < mesh.npoin; point++) for(int k = 0; k < 3; k++)
      mesh.met.rfld(point,k) = std::numeric_limits<double>::quiet_NaN();
  }

  void useBezierGeometry(){
    dblAr2 converted(9,2);
    for(int ientt = 0; ientt < 2; ientt++)
      lag2bez2<2,2>(mesh.fac2poi[ientt],mesh.coord,converted);
    converted.copyTo(mesh.coord);
    mesh.forceBasisFlag(FEBasis::Bezier);
  }

  void computationalPoint(int ientt, const double *bary, double &xi, double &eta){
    xi = eta = 0.;
    for(int node = 0; node < 3; node++){
      const int point = mesh.fac2poi(ientt,node);
      xi += bary[node]*points[point][0];
      eta += bary[node]*points[point][1];
    }
  }
};

void inverseMap(int ientt, double x, double y, double *bary){
  const double d = 1.-b*c, z = x-b*y;
  const double xi = 2.*z/(d+std::sqrt(d*d+4.*a*z));
  const double eta = y-c*xi;
  if(ientt == 0){
    bary[0] = 1.-xi-eta; bary[1] = xi; bary[2] = eta;
  }else{
    bary[0] = 1.-eta; bary[1] = xi+eta-1.; bary[2] = 1.-xi;
  }
}

// Closed-form 2x2 symmetric exponential, independent of the spectral evaluator.
void exactMetric(int ientt, double xi, double eta, double *met){
  double L[3];
  for(int k = 0; k < 3; k++) L[k] = offsets[ientt][k]+dx[ientt][k]*xi+dy[ientt][k]*eta;
  const double mean = (L[0]+L[2])/2., halfDifference = (L[0]-L[2])/2.;
  const double radius = std::hypot(halfDifference,L[1]);
  const double scale = std::exp(mean), sinhc = radius == 0. ? 1. : std::sinh(radius)/radius;
  met[0] = scale*(std::cosh(radius)+sinhc*halfDifference);
  met[1] = scale*sinhc*L[1];
  met[2] = scale*(std::cosh(radius)-sinhc*halfDifference);
}

} // namespace

BOOST_FIXTURE_TEST_CASE(dg_log_metric_uses_owning_entity_and_curved_jacobian,CurvedTriangles){
  initializeDG();
  const double samples[][3] = {{0.2,0.3,0.5},{0.,0.4,0.6},{1.,0.,0.},{0.,1.,0.}};
  for(int basis = 0; basis < 2; basis++){
    if(basis == 1) useBezierGeometry();
    for(int ientt = 0; ientt < 2; ientt++) for(const auto &bary : samples){
      double xi, eta, logm[3], deriv[6], reference[6], valueOnly[3];
      computationalPoint(ientt,bary,xi,eta);
      const double determinant = 1.+2.*a*xi-b*c;
      mesh.met.getLogMetOnEntity(ientt,DifVar::Phys,bary,logm,deriv);
      mesh.met.getLogMetOnEntity(ientt,DifVar::None,bary,valueOnly,nullptr);
      mesh.met.getLogMetOnEntity(ientt,DifVar::Bary,bary,valueOnly,reference);
      for(int k = 0; k < 3; k++){
        BOOST_CHECK_SMALL(logm[k]-offsets[ientt][k]-dx[ientt][k]*xi-dy[ientt][k]*eta,2.e-14);
        BOOST_CHECK_SMALL(valueOnly[k]-logm[k],2.e-14);
        BOOST_CHECK_SMALL(deriv[k]-(dx[ientt][k]-c*dy[ientt][k])/determinant,2.e-13);
        BOOST_CHECK_SMALL(deriv[3+k]-((1.+2.*a*xi)*dy[ientt][k]-b*dx[ientt][k])/determinant,2.e-13);
        BOOST_CHECK_SMALL(reference[k]-(ientt == 0 ? dx[ientt][k] : dy[ientt][k]),2.e-14);
        BOOST_CHECK_SMALL(reference[3+k]-(ientt == 0 ? dy[ientt][k] : dy[ientt][k]-dx[ientt][k]),2.e-14);
      }
    }
  }
}

BOOST_FIXTURE_TEST_CASE(dg_physical_derivatives_match_independent_inverse_map_differences,CurvedTriangles){
  initializeDG();
  const double bary[3] = {0.2,0.3,0.5};
  for(int ientt = 0; ientt < 2; ientt++){
    double xi, eta, logm[3], deriv[6];
    computationalPoint(ientt,bary,xi,eta);
    const double point[2] = {xi+a*xi*xi+b*eta,eta+c*xi};
    mesh.met.getLogMetOnEntity(ientt,DifVar::Phys,bary,logm,deriv);
    for(int r = 0; r < 2; r++){
      constexpr double h = 1.e-6;
      double plusBary[3], minusBary[3], plus[3], minus[3];
      inverseMap(ientt,point[0]+(r==0?h:0.),point[1]+(r==1?h:0.),plusBary);
      inverseMap(ientt,point[0]-(r==0?h:0.),point[1]-(r==1?h:0.),minusBary);
      mesh.met.getLogMetOnEntity(ientt,DifVar::None,plusBary,plus,nullptr);
      mesh.met.getLogMetOnEntity(ientt,DifVar::None,minusBary,minus,nullptr);
      for(int k = 0; k < 3; k++)
        BOOST_CHECK_SMALL(deriv[3*r+k]-(plus[k]-minus[k])/(2.*h),3.e-9);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(dg_constant_local_fields_keep_jumps_and_zero_derivatives,CurvedTriangles){
  initializeDG();
  for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
    for(int k = 0; k < 3; k++) mesh.met.elementCoefficient(ientt,node)[k] = offsets[ientt][k];
  // These two reference points are the same physical point on the shared edge.
  const double bary[2][3] = {{0.,0.4,0.6},{0.4,0.,0.6}};
  for(int ientt = 0; ientt < 2; ientt++){
    double logm[3], deriv[6];
    mesh.met.getLogMetOnEntity(ientt,DifVar::Phys,bary[ientt],logm,deriv);
    for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(logm[k]-offsets[ientt][k],2.e-14);
    for(double value : deriv) BOOST_CHECK_SMALL(value,2.e-14);
  }
}

BOOST_FIXTURE_TEST_CASE(dg_edge_motion_changes_physical_derivatives_only_through_the_geometry,CurvedTriangles){
  initializeDG(); poisonInactivePointBuffer();
  const double bary[3] = {.2,.3,.5}, delta[2] = {.04,-.025};
  double before[2][3], beforeDeriv[2][6];
  for(int ientt = 0; ientt < 2; ientt++)
    mesh.met.getLogMetOnEntity(ientt,DifVar::Phys,bary,before[ientt],beforeDeriv[ientt]);
  mesh.coord(4,0) += delta[0]; mesh.coord(4,1) += delta[1];
  for(int ientt = 0; ientt < 2; ientt++){
    double xi, eta; computationalPoint(ientt,bary,xi,eta);
    const double dshapeXi = ientt == 0 ? 4.*eta : -4.*(1.-eta);
    const double dshapeEta = ientt == 0 ? 4.*xi : -4.*(1.-xi);
    const double j00 = 1.+2.*a*xi+delta[0]*dshapeXi;
    const double j01 = b+delta[0]*dshapeEta;
    const double j10 = c+delta[1]*dshapeXi;
    const double j11 = 1.+delta[1]*dshapeEta;
    const double determinant = j00*j11-j01*j10;
    double met[3], deriv[6], reference[6];
    mesh.met.getLogMetOnEntity(ientt,DifVar::Phys,bary,met,deriv);
    mesh.met.getLogMetOnEntity(ientt,DifVar::Bary,bary,met,reference);
    double change = 0.;
    for(int k = 0; k < 3; k++){
      BOOST_CHECK_EQUAL(met[k],before[ientt][k]);
      BOOST_CHECK_SMALL(deriv[k]-(j11*dx[ientt][k]-j10*dy[ientt][k])/determinant,2e-13);
      BOOST_CHECK_SMALL(deriv[3+k]-(-j01*dx[ientt][k]+j00*dy[ientt][k])/determinant,2e-13);
      change += std::abs(deriv[k]-beforeDeriv[ientt][k])+std::abs(deriv[3+k]-beforeDeriv[ientt][3+k]);
    }
    BOOST_CHECK_GT(change,1e-3);
    for(int node = 0; node < 3; node++){
      const int ipoin = mesh.fac2poi(ientt,node);
      for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],
          offsets[ientt][k]+dx[ientt][k]*points[ipoin][0]+dy[ientt][k]*points[ipoin][1]);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(entity_evaluation_matches_continuous_metric_path,CurvedTriangles){
  const double bary[3] = {0.2,0.3,0.5};
  for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
    mesh.met.setSpace(space);
    for(int ientt = 0; ientt < 2; ientt++){
      double logm[3], deriv[6], expected[3], expectedDeriv[6];
      mesh.met.getLogMetOnEntity(ientt,DifVar::Phys,bary,logm,deriv);
      mesh.met.getMetBary(AsDeg::Pk,DifVar::Phys,MetSpace::Log,
                         mesh.fac2poi[ientt],2,bary,expected,expectedDeriv);
      for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(logm[k]-expected[k],2.e-13);
      for(int k = 0; k < 6; k++) BOOST_CHECK_SMALL(deriv[k]-expectedDeriv[k],2.e-13);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(entity_evaluation_checks_owner_output_and_jacobian,CurvedTriangles){
  initializeDG();
  const double bary[3] = {0.2,0.3,0.5};
  double logm[3], deriv[6];
  BOOST_CHECK_THROW(mesh.met.getLogMetOnEntity(-1,DifVar::None,bary,logm,nullptr),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met.getLogMetOnEntity(2,DifVar::None,bary,logm,nullptr),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met.getLogMetOnEntity(0,DifVar::Phys,bary,logm,nullptr),MetrisExcept);
  mesh.fac2poi(1,0) = -1;
  BOOST_CHECK_THROW(mesh.met.getLogMetOnEntity(1,DifVar::None,bary,logm,nullptr),MetrisExcept);
  for(int point = 0; point < mesh.npoin; point++) mesh.coord(point,1) = 0.;
  BOOST_CHECK_THROW(mesh.met.getLogMetOnEntity(0,DifVar::Phys,bary,logm,deriv),MetrisExcept);
}

BOOST_FIXTURE_TEST_CASE(dg_spd_evaluation_matches_closed_form_and_physical_differences,CurvedTriangles){
  initializeDG();
  const double bary[3] = {0.2,0.3,0.5};
  for(int geometryBasis = 0; geometryBasis < 2; geometryBasis++){
    if(geometryBasis == 1) useBezierGeometry();
    for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
      mesh.met.setBasis(basis);
      for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
        mesh.met.setSpace(space);
        for(int ientt = 0; ientt < 2; ientt++){
          double xi, eta, met[3], deriv[6], expected[3], valueOnly[3];
          computationalPoint(ientt,bary,xi,eta);
          const double point[2] = {xi+a*xi*xi+b*eta,eta+c*xi};
          exactMetric(ientt,xi,eta,expected);
          mesh.met.getMetOnEntity(ientt,DifVar::Phys,bary,met,deriv);
          mesh.met.getMetOnEntity(ientt,DifVar::None,bary,valueOnly,nullptr);
          for(int k = 0; k < 3; k++){
            BOOST_CHECK_SMALL(met[k]-expected[k],3.e-13);
            BOOST_CHECK_SMALL(valueOnly[k]-met[k],3.e-13);
          }
          BOOST_CHECK_GT(met[0]*met[2]-met[1]*met[1],0.);
          for(int r = 0; r < 2; r++){
            constexpr double h = 1.e-6;
            double shifted[3], plus[3], minus[3], sx, sy;
            inverseMap(ientt,point[0]+(r==0?h:0.),point[1]+(r==1?h:0.),shifted);
            computationalPoint(ientt,shifted,sx,sy);
            exactMetric(ientt,sx,sy,plus);
            inverseMap(ientt,point[0]-(r==0?h:0.),point[1]-(r==1?h:0.),shifted);
            computationalPoint(ientt,shifted,sx,sy);
            exactMetric(ientt,sx,sy,minus);
            for(int k = 0; k < 3; k++)
              BOOST_CHECK_SMALL(deriv[3*r+k]-(plus[k]-minus[k])/(2.*h),5.e-9);
          }
          mesh.met.getMetOnEntity(ientt,DifVar::Bary,bary,met,deriv);
          for(int r = 0; r < 2; r++){
            constexpr double h = 1.e-6;
            double shifted[3] = {bary[0]-h,bary[1],bary[2]}, plus[3], minus[3], sx, sy;
            shifted[r+1] += h;
            computationalPoint(ientt,shifted,sx,sy);
            exactMetric(ientt,sx,sy,plus);
            shifted[0] += 2.*h; shifted[r+1] -= 2.*h;
            computationalPoint(ientt,shifted,sx,sy);
            exactMetric(ientt,sx,sy,minus);
            for(int k = 0; k < 3; k++)
              BOOST_CHECK_SMALL(deriv[3*r+k]-(plus[k]-minus[k])/(2.*h),5.e-9);
          }
        }
      }
    }
  }
}

BOOST_FIXTURE_TEST_CASE(dg_frechet_derivative_handles_repeated_eigenvalues,CurvedTriangles){
  initializeDG();
  const double bary[3] = {1.,0.,0.};
  for(double gap : {0.,1.e-12}){
    for(int node = 0; node < 3; node++){
      double *coefficient = mesh.met.elementCoefficient(0,node);
      coefficient[0] = 0.4+gap+(node==1?0.3:0.);
      coefficient[1] = node==1 ? 0.2 : 0.;
      coefficient[2] = 0.4+(node==2?-0.2:0.);
    }
    double logm[3], dlogm[6], met[3], dmet[6];
    mesh.met.getLogMetOnEntity(0,DifVar::Phys,bary,logm,dlogm);
    mesh.met.getMetOnEntity(0,DifVar::Phys,bary,met,dmet);
    BOOST_CHECK_SMALL(met[0]-std::exp(0.4+gap),2.e-14);
    BOOST_CHECK_SMALL(met[1],2.e-14);
    BOOST_CHECK_SMALL(met[2]-std::exp(0.4),2.e-14);
    for(int k = 0; k < 6; k++)
      BOOST_CHECK_SMALL(dmet[k]-std::exp(0.4)*dlogm[k],2.e-12);
  }
}

BOOST_FIXTURE_TEST_CASE(dg_transformations_and_copies_preserve_each_trace,CurvedTriangles){
  initializeDG();
  const double bary[2][3] = {{0.,0.4,0.6},{0.4,0.,0.6}};
  for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
    for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
      MetricFieldFE copy(mesh.met);
      BOOST_CHECK(copy.elementCoefficient(0,0) != mesh.met.elementCoefficient(0,0));
      copy.setSpace(space);
      copy.setBasis(basis);
      BOOST_CHECK(copy.getSpace() == space);
      BOOST_CHECK(copy.getBasis() == basis);
      CurvedTriangles destination;
      destination.mesh.met = copy;
      for(int ientt = 0; ientt < 2; ientt++){
        double original[3], originalD[6], actual[3], actualD[6];
        mesh.met.getMetOnEntity(ientt,DifVar::Phys,bary[ientt],original,originalD);
        destination.mesh.met.getMetOnEntity(ientt,DifVar::Phys,bary[ientt],actual,actualD);
        for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(actual[k]-original[k],3.e-13);
        for(int k = 0; k < 6; k++) BOOST_CHECK_SMALL(actualD[k]-originalD[k],3.e-13);
      }
      copy.normalize(2.);
      for(int ientt = 0; ientt < 2; ientt++){
        double original[3], originalD[6], actual[3], actualD[6], unchanged[3];
        mesh.met.getMetOnEntity(ientt,DifVar::Phys,bary[ientt],original,originalD);
        copy.getMetOnEntity(ientt,DifVar::Phys,bary[ientt],actual,actualD);
        destination.mesh.met.getMetOnEntity(ientt,DifVar::None,bary[ientt],unchanged,nullptr);
        for(int k = 0; k < 3; k++){
          BOOST_CHECK_SMALL(actual[k]-original[k]/4.,3.e-13);
          BOOST_CHECK_SMALL(unchanged[k]-original[k],3.e-13);
        }
        for(int k = 0; k < 6; k++) BOOST_CHECK_SMALL(actualD[k]-originalD[k]/4.,3.e-13);
      }
      destination.mesh.met.setSpace(MetSpace::Log);
      destination.mesh.met.setBasis(FEBasis::Lagrange);
      for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
        for(int k = 0; k < 3; k++)
          BOOST_CHECK_SMALL(destination.mesh.met.elementCoefficient(ientt,node)[k]
                            -mesh.met.elementCoefficient(ientt,node)[k],3.e-13);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(dg_failed_space_conversion_preserves_all_coefficients,CurvedTriangles){
  initializeDG();
  mesh.met.elementCoefficient(1,2)[0] = 1000.;
  MetricFieldFE saved(mesh.met);
  BOOST_CHECK_THROW(mesh.met.setSpace(MetSpace::Exp),MetrisExcept);
  BOOST_CHECK(mesh.met.getSpace() == MetSpace::Log);
  for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
    for(int k = 0; k < 3; k++)
      BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],saved.elementCoefficient(ientt,node)[k]);
  BOOST_CHECK_THROW(mesh.met.normalize(std::numeric_limits<double>::quiet_NaN()),MetrisExcept);

  mesh.met.elementCoefficient(1,2)[0] = 0.2;
  mesh.met.setSpace(MetSpace::Exp);
  double *invalid = mesh.met.elementCoefficient(1,2);
  invalid[0] = invalid[2] = -1.; invalid[1] = 0.;
  saved = mesh.met;
  BOOST_CHECK_THROW(mesh.met.setSpace(MetSpace::Log),MetrisExcept);
  BOOST_CHECK(mesh.met.getSpace() == MetSpace::Exp);
  for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
    for(int k = 0; k < 3; k++)
      BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],saved.elementCoefficient(ientt,node)[k]);
}

BOOST_FIXTURE_TEST_CASE(dg_connection_samples_use_local_traces_in_spd_space,CurvedTriangles){
  initializeDG();
  poisonInactivePointBuffer();
  const double samples[2][3][3] = {
    {{0.,1.,0.},{0.,0.4,0.6},{0.2,0.3,0.5}},
    {{1.,0.,0.},{0.4,0.,0.6},{0.2,0.3,0.5}}
  };
  for(int basis = 0; basis < 2; basis++){
    if(basis == 1) useBezierGeometry();
    for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
      mesh.met.setSpace(space);
      for(int ientt = 0; ientt < 2; ientt++) for(const auto &bary : samples[ientt]){
        double xi, eta, expected[3];
        computationalPoint(ientt,bary,xi,eta);
        exactMetric(ientt,xi,eta,expected);
        const auto sample = prepare_objective_quadrature_sample<MetricFieldFE,2,2,2>(
            mesh,AsDeg::Pk,ientt,{bary,1.},ObjectiveQuadratureTheta::PhysicalMetricMeasure);
        BOOST_CHECK(sample.theta_is_valid);
        for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(sample.metric[k]-expected[k],5.e-13);
        const double measure = 0.5*(1.+2.*a*xi-b*c);
        BOOST_CHECK_SMALL(sample.theta-measure*std::sqrt(expected[0]*expected[2]-expected[1]*expected[1]),5.e-13);

        // Independent geometry derivatives and physical finite differences of M.
        double J[4] = {1.+2.*a*xi,c,b,1.};
        double H[6] = {2.*a,0.,0.,0.,0.,0.}, dM[6];
        if(ientt == 1){
          J[0] = b; J[1] = 1.; J[2] = b-1.-2.*a*xi; J[3] = 1.-c;
          H[0] = 0.; H[4] = 2.*a;
        }
        const double x = xi+a*xi*xi+b*eta, y = eta+c*xi;
        for(int r = 0; r < 2; r++){
          constexpr double h = 1.e-6;
          double shifted[3], plus[3], minus[3], sx, sy;
          inverseMap(ientt,x+(r==0?h:0.),y+(r==1?h:0.),shifted);
          computationalPoint(ientt,shifted,sx,sy);
          exactMetric(ientt,sx,sy,plus);
          inverseMap(ientt,x-(r==0?h:0.),y-(r==1?h:0.),shifted);
          computationalPoint(ientt,shifted,sx,sy);
          exactMetric(ientt,sx,sy,minus);
          for(int k = 0; k < 3; k++) dM[3*r+k] = (plus[k]-minus[k])/(2.*h);
        }
        const double expectedValue = connection_difference_norm2<2>(J,H,expected,dM);
        const double value = quafun_connection<MetricFieldFE,2>(mesh,ientt,bary);
        BOOST_CHECK_SMALL(value-expectedValue,2.e-8*(1.+expectedValue));
        double gradient[2];
        const double differentiated = d_quafun_connection<MetricFieldFE,2>(
            mesh,ientt,bary,4,mesh.getBasis(),gradient);
        BOOST_CHECK_SMALL(differentiated-value,2.e-13*(1.+value));
      }
    }
  }
}

BOOST_FIXTURE_TEST_CASE(dg_integrated_connection_matches_frozen_target_differences,CurvedTriangles){
  initializeDG();
  poisonInactivePointBuffer();
  useBezierGeometry();
  parameters.objective_quadrature_order = 5;
  const auto rule = get_objective_quadrature<2>(5);
  for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
    mesh.met.setSpace(space);
    for(int ientt = 0; ientt < 2; ientt++) for(int node : {0,4}){
      struct FrozenSample { double metric[3], derivative[6], weight; };
      std::vector<FrozenSample> frozen(rule.size());
      for(int q = 0; q < rule.size(); q++){
        auto &sample = frozen[q];
        mesh.met.getMetOnEntity(ientt,DifVar::Phys,rule[q].bary,sample.metric,sample.derivative);
        double xi, eta, expected[3];
        computationalPoint(ientt,rule[q].bary,xi,eta);
        exactMetric(ientt,xi,eta,expected);
        sample.weight = rule[q].weight*0.5*(1.+2.*a*xi-b*c);
        if(standard_objective_quadrature_theta_mode() == ObjectiveQuadratureTheta::PhysicalMetricMeasure)
          sample.weight *= std::sqrt(expected[0]*expected[2]-expected[1]*expected[1]);
      }
      const int point = mesh.fac2poi(ientt,node);
      const double original[2] = {mesh.coord(point,0),mesh.coord(point,1)};
      const auto frozenValue = [&](double shiftX, double shiftY){
        mesh.coord(point,0) = original[0]+shiftX;
        mesh.coord(point,1) = original[1]+shiftY;
        double result = 0.;
        for(int q = 0; q < rule.size(); q++){
          double position[2], J[4], H[6];
          eval2<2,2>(mesh.coord,mesh.fac2poi[ientt],FEBasis::Bezier,
                     DifVar::Bary,DifVar::Bary,rule[q].bary,position,J,H);
          result += frozen[q].weight*connection_difference_norm2<2>(
              J,H,frozen[q].metric,frozen[q].derivative);
        }
        mesh.coord(point,0) = original[0];
        mesh.coord(point,1) = original[1];
        return result;
      };
      double gradient[2], hessian[3], gradientOnly[2];
      const double value = metqua_connection<MetricFieldFE,2>(mesh,AsDeg::Pk,AsDeg::Pk,ientt);
      const double differentiated = d_metqua_connection<MetricFieldFE,2>(
          mesh,AsDeg::Pk,AsDeg::Pk,ientt,node,FEBasis::Bezier,DifVar::None,gradient,hessian);
      const double firstOrder = d_metqua_connection<MetricFieldFE,2>(
          mesh,AsDeg::Pk,AsDeg::Pk,ientt,node,FEBasis::Bezier,DifVar::None,gradientOnly,nullptr);
      BOOST_CHECK_SMALL(value-frozenValue(0.,0.),2.e-12*(1.+value));
      BOOST_CHECK_SMALL(differentiated-value,2.e-13*(1.+value));
      BOOST_CHECK_SMALL(firstOrder-value,2.e-13*(1.+value));
      for(int r = 0; r < 2; r++){
        BOOST_CHECK_SMALL(gradient[r]-gradientOnly[r],2.e-11*(1.+std::abs(gradient[r])));
        constexpr double h = 1.e-5;
        const double plus = frozenValue(r==0?h:0.,r==1?h:0.);
        const double minus = frozenValue(r==0?-h:0.,r==1?-h:0.);
        BOOST_CHECK_SMALL(gradient[r]-(plus-minus)/(2.*h),2.e-7*(1.+std::abs(gradient[r])));
        BOOST_CHECK_SMALL(hessian[2*r]-(plus-2.*value+minus)/(h*h),5.e-5*(1.+std::abs(hessian[2*r])));
      }
      constexpr double h = 1.e-5;
      const double mixed = (frozenValue(h,h)-frozenValue(h,-h)
                          -frozenValue(-h,h)+frozenValue(-h,-h))/(4.*h*h);
      BOOST_CHECK_SMALL(hessian[1]-mixed,5.e-5*(1.+std::abs(hessian[1])));
    }
  }
}

BOOST_FIXTURE_TEST_CASE(dg_connection_rejects_owner_free_and_invalid_calls,CurvedTriangles){
  initializeDG();
  const double samples[2][3] = {{1.,0.,0.},{0.2,0.3,0.5}};
  double gradient[2];
  for(const auto &bary : samples){
    BOOST_CHECK_THROW((quafun_connection<MetricFieldFE,2>(mesh,mesh.fac2poi[0],bary)),MetrisExcept);
    BOOST_CHECK_THROW((d_quafun_connection<MetricFieldFE,2>(
        mesh,mesh.fac2poi[0],bary,0,FEBasis::Lagrange,gradient)),MetrisExcept);
    BOOST_CHECK_THROW((prepare_objective_quadrature_sample<MetricFieldFE,2,2,2>(
        mesh,AsDeg::Pk,mesh.fac2poi[0],{bary,1.},ObjectiveQuadratureTheta::PhysicalMetricMeasure)),MetrisExcept);
    for(int ientt : {-1,2}){
      BOOST_CHECK_THROW((quafun_connection<MetricFieldFE,2>(mesh,ientt,bary)),MetrisExcept);
      BOOST_CHECK_THROW((d_quafun_connection<MetricFieldFE,2>(
          mesh,ientt,bary,0,FEBasis::Lagrange,gradient)),MetrisExcept);
      BOOST_CHECK_THROW((metqua_connection<MetricFieldFE,2>(mesh,AsDeg::Pk,AsDeg::Pk,ientt)),MetrisExcept);
      BOOST_CHECK_THROW((prepare_objective_quadrature_sample<MetricFieldFE,2,2,2>(
          mesh,AsDeg::Pk,ientt,{bary,1.},ObjectiveQuadratureTheta::PhysicalMetricMeasure)),MetrisExcept);
    }
  }
}

} // namespace Metris
