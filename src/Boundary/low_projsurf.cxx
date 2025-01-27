//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_projsurf.hxx"

#include "../Mesh/MeshBase.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../low_geo.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_EGADSprinterr.hxx"
#include "../opt_generic.hxx"
#include "../low_eval.hxx"

namespace Metris{

/* -----------------------------------------------------------------
   FUNCTION1: P1 edge proj, no CAD
   -----------------------------------------------------------------*/


template <int gdim, int ideg>
int projptedg(MeshBase &msh, const double*__restrict__ coop, 
              int iedge, 
              double*__restrict__ bary,
              double*__restrict__ coopr){
  int ierro = 0;
  const int iverb = msh.param->iverb;

  int ipoi1 = msh.edg2poi(iedge,0);
  int ipoi2 = msh.edg2poi(iedge,1);

  if constexpr(ideg == 1){

    double tan[gdim];
    for(int ii = 0; ii < gdim; ii++) tan[ii] = msh.coord(ipoi2,ii)
                                             - msh.coord(ipoi1,ii);

    double x0 = getprdl2<gdim>(msh.coord[ipoi1], tan);
    double x1 = getprdl2<gdim>(msh.coord[ipoi2], tan);
    double tp = (getprdl2<gdim>(coop, tan) - x0) / (x1 - x0);

    // Raw, potentially negative barys
    bary[0] = 1.0 - tp;
    bary[1] =       tp;

    ierro = bary[0] < -Constants::baryTol || bary[0] > 1 + Constants::baryTol
          ||bary[1] < -Constants::baryTol || bary[1] > 1 + Constants::baryTol;

    // Truncate for projection
    tp = MAX(0.0,MIN(1.0,tp));
    for(int ii = 0; ii < gdim; ii++) coopr[ii] = (1.0 - tp) * msh.coord(ipoi1,ii)
                                               +        tp  * msh.coord(ipoi2,ii);

  }else{ // if ideg > 1

    double tol0 = geterrl2<gdim>(msh.coord[ipoi1],msh.coord[ipoi2]);
    tol0 = sqrt(tol0);
    newton_drivertype_args<1> args;
    args.stpmin = 1.0e-12;
    args.iprt = msh.param->iverb - 1;
    int iflag = 0, ihess;
    double xcur, fcur, gcur, hcur;

    double d1F[gdim],d2F[gdim];

    xcur = 0.5;
    while(true){

      ierro = optim_newton_drivertype<1>(args, &xcur, &fcur, &gcur, &hcur, &iflag, &ihess);
      
      if(iverb >= 3) printf("     - newton ret ierro %d iflag %d xcur %f \n",ierro,iflag,xcur);

      if(ierro > 0){
        if(iverb >= 1) printf("  ## optim_newton_drivertype error %d\n",ierro);
        return ierro;
      }
      if(iflag <= 0) {
        if(iverb >= 3) printf("   - iflag = 0 termination\n");
        break;
      }

      bary[0] = xcur;
      bary[1] = 1 - xcur;


      DifVar d2flag = ihess > 0 ? DifVar::Bary : DifVar::None;
      eval1<gdim, ideg>(msh.coord, msh.edg2poi[iedge], msh.getBasis(), DifVar::Bary, d2flag,
                        bary, coopr, d1F, d2F);
      fcur = geterrl2<gdim>(coopr,coop) / 2;
      gcur = getprdl2<gdim>(d1F, coop)  
           - getprdl2<gdim>(d1F, coopr);
      if(abs(gcur) < tol0 * Constants::projedgTol){
        if(iverb >= 3) printf("    - grad = %15.7e < %15.7e = tol\n",
                              abs(gcur), tol0*Constants::projedgTol);
        return 0;
      }
      if(ihess > 0) hcur = getprdl2<gdim>(d2F, coopr)
                         - getprdl2<gdim>(d2F, coop)
                         + getnrml2<gdim>(d1F);
      if(iverb >= 3) printf("     - projptedg bary %f dist = %15.7e gcur %f \n",xcur,fcur,gcur);

    }

    bary[0] = args.xopt[0];
    bary[1] = 1 - args.xopt[0];

    // Projected is outside
    if(bary[0] <  - Constants::baryTol){
      for(int ii = 0; ii < gdim; ii++) coopr[ii] = msh.coord(ipoi2,ii);
      return 0;
    }
    if(bary[0] > 1 + Constants::baryTol){
      for(int ii = 0; ii < gdim; ii++) coopr[ii] = msh.coord(ipoi1,ii);
      return 0;
    }

    eval1<gdim, ideg>(msh.coord, msh.edg2poi[iedge], msh.getBasis(), DifVar::None, DifVar::None,
                      bary, coopr, NULL, NULL);
  }

  return ierro;
}
#define BOOST_PP_LOCAL_MACRO(n)\
template \
int projptedg<2,n>(MeshBase &msh, const double*__restrict__ coop, \
                   int iedge, \
                   double*__restrict__ bary,\
                   double*__restrict__ coopr);\
template \
int projptedg<3,n>(MeshBase &msh, const double*__restrict__ coop, \
                   int iedge, \
                   double*__restrict__ bary,\
                   double*__restrict__ coopr);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()






/* -----------------------------------------------------------------
FUNCTION2: P1 edge proj, with CAD
   -----------------------------------------------------------------*/


// Project point coop on edge iedge, return barycentric of proj and CAD param
// Requires msh to be linked to CAD object. 
template <int gdim>
int projptedgCAD(MeshBase &msh, const double*__restrict__ coop, double tol, 
                 int iedge, 
                 double*__restrict__ param, double*__restrict__ bary,
                 double*__restrict__ coopr){
  int iref = msh.edg2ref[iedge];
  // Parameters from above
  const int iverb = msh.param->iverb;
  // Specific parameters
  const int miter = 100;

  int ipoi1 = msh.edg2poi(iedge,0);
  int ipoi2 = msh.edg2poi(iedge,1);

  int ibpo1 = msh.poi2ebp(ipoi1,1,iedge,-1);
  int ibpo2 = msh.poi2ebp(ipoi2,1,iedge,-1);

  double u1 = msh.bpo2rbi(ibpo1,0);
  double u2 = msh.bpo2rbi(ibpo2,0);

  bool iflip = false;
  if(u1 > u2){
    iflip = true;
    double tmp = u1;
    u1 = u2;
    u2 = tmp;
  }

  ego obj = msh.CAD.cad2edg[iref];
  METRIS_ASSERT(obj != NULL); 

  double nrm0 = geterrl2<gdim>(msh.coord[ipoi1], msh.coord[ipoi2]);

  double tan[gdim];
  for(int ii = 0; ii < gdim; ii++) tan[ii] = msh.coord(ipoi2,ii)
                                           - msh.coord(ipoi1,ii);
  //double dp[gdim];
  //for(int ii = 0; ii < gdim; ii++) tan[ii] = msh.coord(ipoi2,ii)
  //                                         - coop[ii];


  double x0 = getprdl2<gdim>(msh.coord[ipoi1], tan);
  double x1 = getprdl2<gdim>(msh.coord[ipoi2], tan);
  double tp = (getprdl2<gdim>(coop, tan) - x0) / (x1 - x0);
  if(iverb >= 3) printf("   - projptedg ipoi1 %d ipoi2 %d x0 %f x1 %f tp %f\n",
                        ipoi1,ipoi2,x0,x1,tp);
  // If outside, return index of nearest +1:
  // t = 1 iff coop = coord(ipoi1,:)
  bary[0] = 1.0 - tp;
  bary[1] = tp;
  if(tp >= 1+tol) return 1;
  if(tp <   -tol) return 2;



  // This is double the number of iterations to get to machine precision while 
  // dividing 1 by 2 each iter. 
  for(int niter = 0; niter < miter; niter++){

    // t was initialized outside 
    double result[18];
    //*param =        tp *u1 
    //       + (1.0 - tp)*u2;

    *param = (u1 + u2) / 2;

    int ierro = EG_evaluate(obj, param, result);
    if(ierro > 0){
      print_EGADS_error("EG_evaluate", ierro);
      METRIS_THROW_MSG(GeomExcept(), "EGADS error");
    }

    for(int ii = 0; ii < gdim; ii++) coopr[ii] = result[ii];
    //double dist = geterrl2<gdim>(result, coop);

    double t = (getprdl2<gdim>(result, tan) - x0) / (x1 - x0);
    if(abs(t-tp) < tol){
      if(iverb >= 3) 
        printf("  -- projptedg success |t - tp| = %15.7e < tol %15.7e * %15.7e\n",
               abs(t-tp),tol,nrm0);
      //bary[0] = tp;
      //bary[1] = 1.0 - tp;
      return 0;
    }
    if(iverb >= 3) printf("   - projptedg %d / %d param = %f t = %f tp %f u1 %f u2 %f\n",
                          niter,miter,*param,t,tp,u1,u2);
    if(t >= 1+tol) return 1;
    if(t <   -tol) return 2;

    // If EG evaluated point is closer to ipoi1n update u2
    printf("Debug t %f tp %f t > tp %d \n ",t,tp,t>tp);
    if(t > tp && iflip){
      u2 = *param;
    }else{
      u1 = *param;
    }
    //tp = t;

  }


  // Add corner for visualization 
  // Since this is suicide, just add a point (forbidden on MeshBase otherwise)
  int ipoin = msh.npoin;
  msh.set_npoin(msh.npoin + 1);
  msh.template newbpotopo<0>(ipoin,ipoin);
  for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipoin,ii) = coop[ii];
  writeMesh("dbg_projptedg.meshb",msh);

  printf("Failed to proj on edge %d %d \n",msh.edg2poi(iedge,0),msh.edg2poi(iedge,1));

  METRIS_THROW_MSG(AlgoExcept(), "Too many iterations projptedg")
}
template 
int projptedgCAD<2>(MeshBase &msh, const double*__restrict__ coop, double tol, 
                    int iedge, 
                    double*__restrict__ param, double*__restrict__ bary,
                    double*__restrict__ coopr);
template 
int projptedgCAD<3>(MeshBase &msh, const double*__restrict__ coop, double tol, 
                    int iedge, 
                    double*__restrict__ param, double*__restrict__ bary,
                    double*__restrict__ coopr);


}
