//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "low_localization.hxx"

//#include <nlopt.h>

#include "../low_geo.hxx"
#include "../opt_generic.hxx"
#include "../ho_constants.hxx"
#include "../low_eval.hxx"

#include "../linalg/eigen.hxx"
#include "../linalg/matprods.hxx"
#include "../Mesh/MeshBase.hxx"

#include "../aux_pp_inc.hxx"


#include <random>

namespace Metris{
  
// Tess covering tessellations of ielem for whether the point is in there. 
// Return 1 if point inside tess, 0 otherwise. 
// Put bary guess in bary
// Negative tolerance is advised (ensure outside is outside)
// NOTE: this needs revamping. The principle itself is probably wrong. Also useless in 1D.
template<int gdim, int ideg>
int inveval_tessreject(MeshBase &msh, int ientt, const double* coor0, double*__restrict__ bary, double tol, int iprt){
  METRIS_THROW_MSG(TODOExcept(), "Rework inveval_tessreject");
  return 0;
  #if 0

  static_assert(ideg > 1);
//  static_assert(gdim >= 1 && gdim <= 3);
  static_assert(gdim == 3);

  double tolP1 = 1.0e-2;
  if(msh.poi2rwk.size() < 3*msh.npoin) METRIS_THROW_MSG(DMemExcept(), "poi2rwk insufficient size = "<<msh.poi2rwk.size())
  dblAr2 rfld2(msh.npoin,3,&msh.poi2rwk[0]);

  intAr2 &ent2poi = gdim == 1 ? msh.edg2poi : 
                    gdim == 2 ? msh.fac2poi : msh.tet2poi;

  // We're not going to be clever. Just run it on Lagrange then Bézier. 
  // Point must be outside both tessellations. 
  // Note: this could perhaps still fail in degree >= 3 where we'd have preferred
  // () - Lag - Béz - ()
  // This is not very important as the point of this is to accelerate rejection. 
  // Standard inveval3 can still handle rejecting a point if it comes to that. 
  // In the future, devise selection strategy depending on local "concavity". 
  if(msh.getBasis() == FEBasis::Lagrange){
    if(gdim == 1){
      lag2bez1<ideg,gdim>(ent2poi[ientt],msh.coord,rfld2);
    }else if(gdim == 2){
      lag2bez2<ideg,gdim>(ent2poi[ientt],msh.coord,rfld2);
    }else{
      lag2bez3<ideg,gdim>(ent2poi[ientt],msh.coord,rfld2);
    }
  }else{
    if(gdim == 1){
      bez2lag1<ideg,gdim>(ent2poi[ientt],msh.coord,rfld2);
    }else if(gdim == 2){
      bez2lag2<ideg,gdim>(ent2poi[ientt],msh.coord,rfld2);
    }else{
      bez2lag3<ideg,gdim>(ent2poi[ientt],msh.coord,rfld2);
    }
  }


  double meas0 = gdim*getmeasentP1<gdim>(ent2poi[ientt],msh.coord);
  if(gdim > 1) meas0 *= gdim - 1;
  meas0 = abs(meas0);

  // The lower the minimum, the least confident we are in saying the point is outside. 
  double ccoef[tetnpps[gdim*(ideg-1)]];
  getccoef<gdim,gdim,ideg>(msh, ientt, NULL, ccoef);
  double jmin = 1.0e30;
  double sumj = 0.0;
  constexpr int npps = gdim == 1 ? edgnpps[gdim*(ideg-1)] : 
                       gdim == 2 ? facnpps[gdim*(ideg-1)] : tetnpps[gdim*(ideg-1)];

  for(int ii = 0; ii < npps; ii++){
    jmin = MIN(jmin,ccoef[ii]);
    sumj += ccoef[ii];
  }
  jmin /= sumj;
  if(iprt > 0) printf("Using jmin = %15.7e\n",jmin);


  double baryl[gdim+1];
  double bary1[gdim+1];
  double coopr[gdim];
  bool outP1 = true; 
  for(int ii = 0; ii < gdim + 1; ii++) bary[ii] = 1.0 / (gdim + 1);
  double nrmi = 1.0e30;
  // Do Lagrange and Bézier. 
  for(int ifld = 0; ifld <= 1 ; ifld ++){
    dblAr2 &crd = ifld == 0 ? msh.coord : rfld2;

    int nsub = subref.ntet[ideg];

    for(int isub = 0; isub < nsub; isub++){
      int irn[gdim+1];
      int irn1 = subref.tet[ideg][isub][0];
      int irn2 = subref.tet[ideg][isub][1];
      int irn3 = subref.tet[ideg][isub][2];
      int irn4 = subref.tet[ideg][isub][3];

      if(iprt > 0) printf("  isub %d irn1, irn2, irn3, irn4 %d %d %d %d \n",isub,irn1, irn2, irn3, irn4);

      assert(irn1 >= 0 && irn1 < tetnpps[ideg]);
      assert(irn2 >= 0 && irn2 < tetnpps[ideg]);
      assert(irn3 >= 0 && irn3 < tetnpps[ideg]);
      assert(irn4 >= 0 && irn4 < tetnpps[ideg]);

      inventP1<gdim>(ent2poi[ientt],crd,coor0,baryl);
      //invtetP1(crd[ent2poi(ientt,irn1)],crd[ent2poi(ientt,irn2)],
      //         crd[ent2poi(ientt,irn3)],crd[ent2poi(ientt,irn4)],
      //         coor0, baryl);
      if(iprt > 0) printf(" bary in sub (%d) %f %f %f %f  \n",ifld,baryl[0],baryl[1],baryl[2],baryl[3]);
      for(int ii = 0; ii < 4; ii++){
        bary1[ii] = ( baryl[0]*ordtet.s[ideg][irn1][ii] + baryl[1]*ordtet.s[ideg][irn2][ii] 
                    + baryl[2]*ordtet.s[ideg][irn3][ii] + baryl[3]*ordtet.s[ideg][irn4][ii])/ideg;
      }
      if(iprt > 0) printf("Corresponding to master bary %f %f %f %f \n",bary1[0],bary1[1],bary1[2],bary1[3]);
      bool iinside = true; 
      for(int ii = 0; ii < 4; ii++){
        if(bary1[ii] < -1.0e-16 || 1 - bary1[ii] < -1.0e-16) iinside = false;
      }
      if(iinside){
        eval3<3,ideg>(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,DifVar::None,
                     bary,coopr,NULL,NULL);
        double nrm = geterrl2<3>(coopr,coor0);
        if(iprt > 0) printf("Distance of evaluated %15.7e \n",sqrt(nrm));
        if(nrm < tol*tol){
          if(iprt > 0) printf("Found within tolerance\n");
          bary[0] = bary1[0];
          bary[1] = bary1[1];
          bary[2] = bary1[2];
          bary[3] = bary1[3];
          return 1;
        }else if(nrm < nrmi){
          nrmi = nrm;
          bary[0] = bary1[0];
          bary[1] = bary1[1];
          bary[2] = bary1[2];
          bary[3] = bary1[3];
        }
      }
      // Even if the guess has not convinced us, we need to see if the bary within the P1 sub elt is decent. 
      // With negative tolerance to cut some slack. 
      double jmat[3][3];
      for(int i = 0; i < 3; i++){
        jmat[0][i] = crd[ent2poi(ientt,irn2)][i] - crd[ent2poi(ientt,irn1)][i];
        jmat[1][i] = crd[ent2poi(ientt,irn3)][i] - crd[ent2poi(ientt,irn1)][i];
        jmat[2][i] = crd[ent2poi(ientt,irn4)][i] - crd[ent2poi(ientt,irn1)][i];
      }
      double meas1 = detmat<gdim>(jmat[0]);
      meas1 = abs(meas1);
      
      // Small volume sub-elements should be given a higher tolerance. 
      // These are the ones we're the least sure of. 
      double fac = (meas0/meas1)*(meas0/meas1)/jmin/jmin;
      const double tol0 = -0.00001;
      if(iprt > 0)printf("Effective tolerance for inP1 test %15.7e\n",tol0*fac);
      bool in = true; 
      for(int ii = 0; ii < 4; ii++){
        if(baryl[ii] < tol0*fac || 1 - baryl[ii] < tol0*fac) in = false;
      }
      if(in) outP1 = false;
    }
  }

  // None produced a < tol guess, but at least some were valid. 
  // We may or may not have improved bary; we have not if no guess translated to viable barys. 
  // In that case, the init 0.25 stands. 
  if(!outP1){
    if(iprt > 0) printf("No guess good enough but cannot reject.\n");
    return 1; 
  }

  if(iprt > 1){
    int nelem = subref.ntet[ideg];
    intAr2 tet2pol(nelem,4,&msh.poi2iwk[0]);
    for(int isub=0;isub<nelem;isub++){
      tet2pol[isub][0] = ent2poi[ientt][subref.tet[ideg][isub][0]];
      tet2pol[isub][1] = ent2poi[ientt][subref.tet[ideg][isub][1]];
      tet2pol[isub][2] = ent2poi[ientt][subref.tet[ideg][isub][2]];
      tet2pol[isub][3] = ent2poi[ientt][subref.tet[ideg][isub][3]];
    }
    writeMesh("tess1.meshb",1,1,msh.npoin,msh.coord,
              nelem,tet2pol,msh.tet2ref,
              0,msh.fac2poi,msh.fac2ref,
              0,msh.edg2poi,msh.edg2ref);
    writeMesh("tess2.meshb",1,1,msh.npoin,rfld2,
              nelem,tet2pol,msh.tet2ref,
              0,msh.fac2poi,msh.fac2ref,
              0,msh.edg2poi,msh.edg2ref);
  }
  return 0;
  #endif
}

#define BOOST_PP_LOCAL_MACRO(n)\
template int inveval_tessreject<3, n >(MeshBase &msh, int ielem, const double* coor0, double*__restrict__ bary,double,int);
#define BOOST_PP_LOCAL_LIMITS     (2, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


/* Second-order method; Newton's on min_xi ||F_K(xi) - X_0||^2  */

/*
Return values:
 -  0 found point in element
 - -1 (if ilazy): point not found but outside with correct bary signs
*/
template<int ideg> 
int inveval3_q2_Newton(MeshBase &msh, int ielem, const double* coor0, 
	            double* __restrict__ coopr, double* __restrict__ bary,
	            double tol, int *ncall, bool ilazy, int iprt){

  double bary00[4] = {0.25,0.25,0.25,0.25};

  *ncall = 0;

  double xtol = -1,stpmin = 1.0e-3, 
         wlfc1 = 0.1 , wlfc2 = 10.0,ratnew = 0.5;
  int niter = 0, maxit = 30, iflgp, iflag = 0, ihess, ierro;
  int nrwrk = 34, niwrk = 3;
  double rwrkN[nrwrk];
  int    iwrkN[niwrk];
  double xopt[3], fopt = -1;
  double xcur[4], fcur, rhs[3];
  double jmat[3][3], bar0[4], hmat[6][3];
  double hess[6];
  double rwork[10];
  double eigval[3], eigvec[9];

  bool icstr[4] = {false,false,false,false};

  int irestart = 0;
  const int mrestart = 50;


  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(1);

  restart: 
  iflag = 0;

  bary[0] = bary00[0];
	bary[1] = bary00[1];
	bary[2] = bary00[2];
	bary[3] = bary00[3];

  bar0[0] = bary[0];
  bar0[1] = bary[1];
  bar0[2] = bary[2];
  bar0[3] = bary[3];

  xcur[0] = 0;
  xcur[1] = 0;
  xcur[2] = 0;
  xcur[3] = 0;


	if(iprt > 0) printf("--- \nNewton 2 start ideg = %d ilag = %d bary = %f %f %f %f\n" ,
                          ideg,(int)msh.getBasis(),bary[0],bary[1],bary[2],bary[3]);

  fcur = 1;
  while(true){

    iflgp = iflag;
    optim_newton_drivertype(3    ,
                            &xcur[1] ,&fcur  ,rhs   ,hess ,
                            xtol ,stpmin, 1 , 
                            wlfc1,wlfc2 ,ratnew ,
                            &niter,maxit ,iprt > 1 ? 2 : 0 ,
                            &iflag,&ihess , 
                            nrwrk,rwrkN ,
                            niwrk,iwrkN ,
                            xopt ,&fopt ,&ierro);
    if(ierro > 0){
      printf("Newton ierro %d\n",ierro);
      return ierro;
    }
    if(iflag <= 0) {
      if(iprt > 0) printf("iflag = 0 termination\n");
      break;
    }

    xcur[0] = - xcur[1] - xcur[2] - xcur[3];


    if(iflgp == 2 && iflag == 1){
      bar0[0] = bary[0];
      bar0[1] = bary[1];
      bar0[2] = bary[2];
      bar0[3] = bary[3];
      if(iprt > 0) printf("Catch update: set bar0 = %f %f %f %f \n",bar0[0],bar0[1],bar0[2],bar0[3]);
                   xcur[0] = 0;
      rwrkN[3+0] = xcur[1] = 0;
      rwrkN[3+1] = xcur[2] = 0;
      rwrkN[3+2] = xcur[3] = 0;
      icstr[0] = false;
      icstr[1] = false;
      icstr[2] = false;
      icstr[3] = false;
      rwrkN[1] = ratnew; // Reset to 1, it will be redivided by ratnew in Newton. 
    }

    double dnrm = abs(xcur[0]) + abs(xcur[1]) + abs(xcur[2]) + abs(xcur[3]);
    if(iprt >= 2) printf("  - Newton descent norm = %23.15e \n",dnrm);

    if(dnrm >= 1.0e-15){

      // Since we recentre at each update, xcur is only ever the descent direction (multiplied by step).
      // Only at the first step can it lead us outside the domain. 
      if(ideg > 1 && iflag == 2 && iflgp == 1){
        ierro = constrain_bary_desc(3,icstr,bary,bar0,xcur, iprt);
        if(ierro > 0){
          if(iprt > 0) printf("All constrained iterate, restart.\n");
          break; 
        }
        rwrkN[2*3+3] = xcur[1];
        rwrkN[2*3+4] = xcur[2];
        rwrkN[2*3+5] = xcur[3];
      }


      bary[0] = bar0[0] + xcur[0];
      bary[1] = bar0[1] + xcur[1];
      bary[2] = bar0[2] + xcur[2];
      bary[3] = bar0[3] + xcur[3];

      METRIS_ASSERT_MSG(bary[0] >= -MAX(tol/10,1.0e-14) && bary[0]-1 <= MAX(tol/10,1.0e-14), "bary[0] within bounds " << bary[0]);
      METRIS_ASSERT_MSG(bary[1] >= -MAX(tol/10,1.0e-14) && bary[1]-1 <= MAX(tol/10,1.0e-14), "bary[1] within bounds " << bary[1]);
      METRIS_ASSERT_MSG(bary[2] >= -MAX(tol/10,1.0e-14) && bary[2]-1 <= MAX(tol/10,1.0e-14), "bary[2] within bounds " << bary[2]);
      METRIS_ASSERT_MSG(bary[3] >= -MAX(tol/10,1.0e-14) && bary[3]-1 <= MAX(tol/10,1.0e-14), "bary[3] within bounds " << bary[3]);

    }

    if(ihess > 0){
  		eval3<3,ideg>(msh.coord,msh.tet2poi[ielem],msh.getBasis(),DifVar::Bary,DifVar::Bary,bary,coopr,jmat[0],hmat[0]);
    }else{
  		eval3<3,ideg>(msh.coord,msh.tet2poi[ielem],msh.getBasis(),DifVar::Bary,DifVar::None,bary,coopr,jmat[0],NULL);
    }
    (*ncall)++;

    double nrm = geterrl2<3>(coopr,coor0);
    fcur = nrm;

    if(iprt > 0) printf("  - Newton %d dist = %15.7e bary %f %f %f %f \n",
                                niter,sqrt(nrm),bary[0],bary[1],bary[2],bary[3]);

    if(nrm < tol*tol){
      if(iprt > 0) printf("Newton terminate with low dist\n");
      return 0;
    } 


		matvdft(3,jmat[0],coor0,coopr,rhs);

    if(ihess > 0){
      // d_11 
      hess[0] = jmat[0][0]*jmat[0][0] 
              + jmat[0][1]*jmat[0][1] 
              + jmat[0][2]*jmat[0][2] 
              + hmat[0][0]*(coopr[0] - coor0[0]) 
              + hmat[0][1]*(coopr[1] - coor0[1])
              + hmat[0][2]*(coopr[2] - coor0[2]);
      // d_12 
      hess[1] = jmat[0][0]*jmat[1][0] 
              + jmat[0][1]*jmat[1][1] 
              + jmat[0][2]*jmat[1][2] 
              + hmat[1][0]*(coopr[0] - coor0[0]) 
              + hmat[1][1]*(coopr[1] - coor0[1])
              + hmat[1][2]*(coopr[2] - coor0[2]);
      // d_22 
      hess[2] = jmat[1][0]*jmat[1][0] 
              + jmat[1][1]*jmat[1][1] 
              + jmat[1][2]*jmat[1][2] 
              + hmat[2][0]*(coopr[0] - coor0[0]) 
              + hmat[2][1]*(coopr[1] - coor0[1])
              + hmat[2][2]*(coopr[2] - coor0[2]);
      // d_13 
      hess[3] = jmat[0][0]*jmat[2][0] 
              + jmat[0][1]*jmat[2][1] 
              + jmat[0][2]*jmat[2][2] 
              + hmat[3][0]*(coopr[0] - coor0[0]) 
              + hmat[3][1]*(coopr[1] - coor0[1])
              + hmat[3][2]*(coopr[2] - coor0[2]);
      // d_23 
      hess[4] = jmat[1][0]*jmat[2][0] 
              + jmat[1][1]*jmat[2][1] 
              + jmat[1][2]*jmat[2][2] 
              + hmat[4][0]*(coopr[0] - coor0[0]) 
              + hmat[4][1]*(coopr[1] - coor0[1])
              + hmat[4][2]*(coopr[2] - coor0[2]);
      // d_33 
      hess[5] = jmat[2][0]*jmat[2][0] 
              + jmat[2][1]*jmat[2][1] 
              + jmat[2][2]*jmat[2][2] 
              + hmat[5][0]*(coopr[0] - coor0[0]) 
              + hmat[5][1]*(coopr[1] - coor0[1])
              + hmat[5][2]*(coopr[2] - coor0[2]);

      // Some notes:
      // The mapping hessians need not be positive definite. 
      // Even in 1D, an under-linear mapping is perfectly possible and not convex (in fact, concave)
      // While the 2-norm itself is convex, there is no guarantee the composition with F_K
      // remains so. 

      geteigsym<3>(hess,10,rwork,eigval,eigvec);
      bool ineg = false;
      for(int ii = 0; ii < 3; ii++){
        if(eigval[ii] < 1.0e-6){
          ineg = true;
          eigval[ii] = MAX(1.0e-6, -eigval[ii]);
        }
      }
      if(ineg) eig2met<3>(eigval,eigvec,hess);
    }

  }


  if constexpr (ideg > 1){
    if(irestart == 0){
      irestart ++;

      // Could possibly do cheaper. 
      int imin = -1;
      double dmin = 1.0e30;
      for(int ii = 0; ii < tetnpps[ideg]; ii++){
        bary[0] = ordtet.s[ideg][ii][0]/(double)ideg;
        bary[1] = ordtet.s[ideg][ii][1]/(double)ideg;
        bary[2] = ordtet.s[ideg][ii][2]/(double)ideg;
        bary[3] = ordtet.s[ideg][ii][3]/(double)ideg;
  		  eval3<3,ideg>(msh.coord,msh.tet2poi[ielem],msh.getBasis(),DifVar::None,DifVar::None,bary,coopr,NULL,NULL);
        double nrm = geterrl2<3>(coor0,coopr);
        if(nrm < dmin){
          dmin = nrm;
          imin = ii;
        }
        if(nrm < tol*tol){
          if(iprt > 0) printf("Found optimum while sampling.\n");
          return 0;
        }
      }
      bary00[0] = ordtet.s[ideg][imin][0]/(double)ideg;
      bary00[1] = ordtet.s[ideg][imin][1]/(double)ideg;
      bary00[2] = ordtet.s[ideg][imin][2]/(double)ideg;
      bary00[3] = ordtet.s[ideg][imin][3]/(double)ideg;
      if(iprt > 0) printf("Try again after sampling got mindist = %15.7e at %f %f %f %f \n",
            dmin, bary00[0], bary00[1], bary00[2], bary00[3]);
      goto restart;
    }else if (irestart < mrestart){
      irestart++;
      double sum = 0;
      for(int jj = 0; jj < 4; jj++){
        bary00[jj] = unif(rng);
        sum += bary00[jj];
      }
      for(int jj = 0; jj < 4; jj++){
        bary00[jj] /= sum;
      }
      if(iprt > 0)printf("Restart %d try random bary %f %f %f %f \n",irestart,
                              bary00[0], bary00[1], bary00[2], bary00[3]);
      goto restart;
    }
  }


  // If a distance had been low enough before, we would have exited already. 
  return 1;
}
#define BOOST_PP_LOCAL_MACRO(n)\
template int inveval3_q2_Newton< n >(MeshBase&, int, const double*, double* __restrict__, double* __restrict__, double, int*, bool, int);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



#if 0
// Not necessary: we have inventP1; also these don't evaluate coop. 
/* First-order methods; Newton's F_K = X_0 */

// Specializations for degree 1 should be very simple 
template<>
int inveval_linesearch<1,1>(MeshBase &msh, int ientt, const double* coor0, 
              double* __restrict__ coopr, double* __restrict__ bary,
              double tol, int *ncall, bool ilazy, int iprt){
  METRIS_THROW_MSG(TODOExcept(),"inveval_linesearch not implemented for gdim = 1");
}

template<>
int inveval_linesearch<2,1>(MeshBase &msh, int ientt, const double* coor0, 
              double* __restrict__ coopr, double* __restrict__ bary,
              double tol, int *ncall, bool ilazy, int iprt){

  const int ip1 = msh.fac2poi[ientt][0];
  const int ip2 = msh.fac2poi[ientt][1];
  const int ip3 = msh.fac2poi[ientt][2];

  double meas, meas0 = abs(det2_vdif(msh.coord[ip2],msh.coord[ip1],
                                     msh.coord[ip3],msh.coord[ip1]));

  bary[0] = det2_vdif(msh.coord[ip2],coor0,
                      msh.coord[ip3],coor0) / meas0;

  bary[1] = det2_vdif(coor0         ,msh.coord[ip1],
                      msh.coord[ip3],msh.coord[ip1]) / meas0;

  bary[2] = det2_vdif(msh.coord[ip2],msh.coord[ip1],
                      coor0         ,msh.coord[ip1]) / meas0;

  //Unstable: often times, one of the barys is almost 0, can't get accurate 
  // 1.0e-15 order value from 1 - ~1/2 - ~1/2
  //bary[2] = 1.0 - bary[0] - bary[1];

  if(bary[0] < -tol || bary[1] < -tol || bary[2] < -tol) return -1;
  return 0;
}

template<>
int inveval_linesearch<3,1>(MeshBase &msh, int ientt, const double* coor0, 
              double* __restrict__ coopr, double* __restrict__ bary,
              double tol, int *ncall, bool ilazy, int iprt){

  const int ip1 = msh.tet2poi[ientt][0];
  const int ip2 = msh.tet2poi[ientt][1];
  const int ip3 = msh.tet2poi[ientt][2];
  const int ip4 = msh.tet2poi[ientt][3];

  double meas, meas0 = abs(det3_vdif(msh.coord[ip2],msh.coord[ip1],
                                     msh.coord[ip3],msh.coord[ip1],
                                     msh.coord[ip4],msh.coord[ip1]));

  bary[0] = det3_vdif(msh.coord[ip2],coor0,
                      msh.coord[ip3],coor0,
                      msh.coord[ip4],coor0) / meas0;

  bary[1] = det3_vdif(coor0         ,msh.coord[ip1],
                      msh.coord[ip3],msh.coord[ip1],
                      msh.coord[ip4],msh.coord[ip1]) / meas0;

  bary[2] = det3_vdif(msh.coord[ip2],msh.coord[ip1],
                      coor0         ,msh.coord[ip1],
                      msh.coord[ip4],msh.coord[ip1]) / meas0;

  bary[3] = det3_vdif(msh.coord[ip2],msh.coord[ip1],
                      msh.coord[ip3],msh.coord[ip1],
                      coor0         ,msh.coord[ip1]) / meas0;

  //bary[3] = 1.0 - bary[0] - bary[1] - bary[2];


  if(bary[0] < -tol || bary[1] < -tol || bary[2] < -tol || bary[3] < -tol) return -1;
  return 0;
}
#endif

/*
(originally)
Tetrahedron inversion using Newton's (secant) method, 
i.e. Jacobian linearizations
Return values:
 -  0 found point in element
 - -1 (if ilazy): point not found but outside with correct bary signs
 Update: 
 element of dimension gdim. Geometric == topo dim (invertibility)
*/
template<int gdim, int ideg> 
int inveval_linesearch(MeshBase &msh, int ientt, const double* coor0, 
	            double* __restrict__ coopr, double* __restrict__ bary,
	            double tol, int *ncall, bool ilazy, int iprt){
  static_assert(gdim == 1 || gdim == 2 || gdim == 3);
  constexpr int tdim = gdim;
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  constexpr auto evalf = tdim == 1 ? eval1<gdim,ideg> : 
                         tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;


  if constexpr(ideg == 1){
    inventP1<gdim>(ent2poi[ientt], msh.coord, coor0, bary);
    if(iprt >= 3) printf("    - inveval_linesearch inventP1 got bary = ");
    if(iprt >= 3) dblAr1(gdim + 1, bary).print();
    
    //invtetP1(msh.coord[ent2poi(ientt,0)],msh.coord[ent2poi(ientt,1)],
    //         msh.coord[ent2poi(ientt,2)],msh.coord[ent2poi(ientt,3)],
    //         coor0, bary);
    bool iout = false;
    for(int ii = 0 ;ii < gdim+1; ii++){
      if(bary[ii] >   - Constants::baryTol 
      && bary[ii] < 1 + Constants::baryTol) continue;
      iout = true;
    }
    if(gdim == 1){
      eval1<gdim,1>(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                    DifVar::None,bary,coopr,NULL,NULL);
    }else if(gdim == 2){
      eval2<gdim,1>(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                    DifVar::None,bary,coopr,NULL,NULL);
    }else{
      eval3<gdim,1>(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                    DifVar::None,bary,coopr,NULL,NULL);
    }
    return iout;
  }

  if(iprt >= 3){
    printf("  -- Start inveval_linesearch ientt %d \n",ientt);

    printf("   - Vertices : ");
    intAr1(msh.nnode(gdim),ent2poi[ientt]).print();
    for(int ii = 0; ii < msh.nnode(gdim); ii++){
      printf("     %d : ",ent2poi[ientt][ii]);
      dblAr1(gdim,msh.coord[ent2poi[ientt][ii]]).print();
    }
  }

//	if(iprt >= 3) printf("Start 0 bary = %f %f %f %f\n" ,bary[0],bary[1],bary[2],bary[3]);
  //if constexpr(ideg == 1){
  //  return inveval_linesearch<3,1>(msh,ientt,coor0,coopr,bary,tol,ilazy,iprt);
  //}
  *ncall = 0;

  double bary00[gdim+1];
  //for(int ii = 0; ii < gdim + 1; ii++) bary00[ii] = 1.0 / (gdim + 1);
  for(int ii = 0; ii < gdim + 1; ii++) bary00[ii] = bary[ii];
  if constexpr(ideg > 1 && gdim == 3){
    int iout = inveval_tessreject<gdim,ideg>(msh, ientt, coor0, bary00, tol, iprt);
    if(iout == 0) return 1;
  }

  bool icstr[gdim+1]; 
  double bar0[gdim+1];
  //for(int ii = 0; ii < gdim + 1; ii++) bar0[ii] = 1.0 / (gdim + 1);
  for(int ii = 0; ii < gdim + 1; ii++) bar0[ii] = bary[ii];

  double xcur[gdim+1];
  for(int ii = 0; ii < gdim + 1; ii++) xcur[ii] = bary00[ii] - bar0[ii];
  constrain_bary_desc(gdim,icstr,bary00,bar0,xcur,iprt);
  for(int ii = 0 ;ii < gdim + 1; ii++) bary00[ii] = bar0[ii] + xcur[ii];

  if(iprt >= 3){
    printf("   - Initial bary guess = ");
    dblAr1(gdim+1,bary00).print();
  }


  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(1);


  double xtol = -1,stpmin = 1.0e-3, 
         wlfc1 = 0.1 , wlfc2 = 10.0,ratnew = 0.5;
  int niter = 0, maxit = 50, iflgp, iflag = 0, ijaco, ierro;
  int nrwrk = 34, niwrk = 3;
  double rwrkN[nrwrk];
  int    iwrkN[niwrk];
  double xopt[gdim], fopt = -1;
  double fcur, rhs[gdim];
  double jmat[gdim][gdim];

  int irestart = 0;
  const int mrestart = 50;



  restart: 
  iflag = 0;

  for(int ii = 0; ii < gdim + 1; ii++){
    bary[ii] = bary00[ii];
    bar0[ii] = bary[ii];
    xcur[ii] = 0;
  }

	if(iprt >= 3){
    printf("   - Start ideg = %d ilag = %d bary = \n" ,
                          ideg,(int)msh.getBasis());
    dblAr1(gdim+1,bary).print();
  }

  fcur = 1;
  while(true){

    iflgp = iflag;
    optim_newton_drivertype(gdim     ,
                            &xcur[1] ,&fcur  ,rhs   ,jmat[0] ,
                            xtol ,stpmin, -1 , 
                            wlfc1,wlfc2 ,ratnew ,
                            &niter,maxit ,iprt > 1 ? 2 : 0 ,
                            &iflag,&ijaco , 
                            nrwrk,rwrkN ,
                            niwrk,iwrkN ,
                            xopt ,&fopt ,&ierro);
    if(ierro > 0){
      printf("   - optim_newton_drivertype error %d\n",ierro);
      return ierro;
    }
    if(iflag <= 0) {
      if(iprt >= 3) printf("   - iflag = 0 termination\n");
      break;
    }


    xcur[0] = - xcur[1];
    for(int ii = 2; ii < gdim + 1; ii++) xcur[0] -= xcur[ii];
    //xcur[0] = - xcur[1] - xcur[2] - xcur[3];

    if(iflgp == 2 && iflag == 1){
      for(int ii = 0; ii < gdim + 1; ii++){
        bar0[ii] = bary[ii];
        // Reset iterate to 0, we prefer to see it as a descent direction.
        xcur[ii] = 0;
        icstr[ii] = false;
      }
                   
      // Previous iterate: only used to compute next iterate. We can set it to zero, 
      // as we're manually updating the current iterate. 
      for(int ii = 0; ii < gdim; ii++){
        rwrkN[gdim+ii] = 0;
      }

      if(iprt >= 3){
        printf("   - update: set bar0 = \n");
        dblAr1(gdim+1,bar0).print();
      }

      rwrkN[1] = ratnew; // Reset to 1, it will be redivided by ratnew in Newton. 
    }




    // Since we recentre at each update, xcur is only ever the descent direction (multiplied by step).
    // Only at the first step can it lead us outside the domain. 
    if(ideg > 1 && iflag == 2 && iflgp == 1){

      ierro = constrain_bary_desc(gdim, icstr,bary,bar0,xcur, iprt);
      if(ierro > 0){
        if(iprt >= 3) printf("   - All constrained iterate, restart.\n");
        break; 
      }
      for(int ii = 0; ii < gdim; ii++){
        rwrkN[2*gdim+3+ii] = xcur[1+ii];
      }
      //rwrkN[2*3+3] = xcur[1];
      //rwrkN[2*3+4] = xcur[2];
      //rwrkN[2*3+5] = xcur[3];
    }
    for(int ii = 0; ii < gdim + 1; ii++){
      bary[ii] = bar0[ii] + xcur[ii];
      METRIS_ASSERT_MSG(bary[ii] >= -MAX(tol/10,1.0e-14) && bary[ii]-1 <= MAX(tol/10,1.0e-14), "bary[ii] within bounds " << bary[ii]);
    }


    if(ijaco > 0){
  		evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::Bary,DifVar::None,bary,coopr,jmat[0],NULL);
      if(iprt >= 3){
        printf("   - Computed at bary = %f %f %f jmat = \n",bary[0],bary[1],bary[2]);
        dblAr2(gdim,gdim,jmat[0]).print();
        printf("   - Computed coopr = ");
        dblAr1(gdim,coopr).print();
        printf("   - Using coor0 = ");
        dblAr1(gdim,coor0).print();
      }
    }else{
  		evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,DifVar::None,bary,coopr,NULL,NULL);
    }
    (*ncall)++;
    double nrm = geterrl2<gdim>(coopr,coor0);
    fcur = sqrt(nrm);

    if(iprt >= 3){
      printf("   - Newton %d dist = %15.7e (^2 = %15.7e) bary = \n",
                                niter,sqrt(nrm),nrm);
      dblAr1(gdim+1,bary).print();
    }

    if(nrm < tol*tol){
      if(iprt >= 3) printf("   - Newton terminate with low dist\n");
      return 0;
    } 

    //    fcur /= 2; // Trick Newton into accepting every step. 
    for(int ii = 0; ii < gdim; ii++) rhs[ii] = coopr[ii] - coor0[ii];
  }


  if constexpr (ideg > 1){
    if(irestart == 0){
      irestart ++;

      // Could possibly do cheaper. 
      constexpr int nn = gdim == 1 ? edgnpps[ideg] 
                       : gdim == 2 ? facnpps[ideg] : tetnpps[ideg];
      auto ordent = ORDELT(gdim);

      int imin = -1;
      double dmin = 1.0e30;
      for(int ii = 0; ii < nn; ii++){
        for(int jj = 0; jj < gdim + 1; jj++) bary[jj] = ordent[ideg][ii][jj]/(double)ideg;
  		  evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,DifVar::None,bary,coopr,NULL,NULL);
        double nrm = geterrl2<gdim>(coor0,coopr);
        if(nrm < dmin){
          dmin = nrm;
          imin = ii;
        }
        if(nrm < tol*tol){
          if(iprt >= 3) printf("   - Found optimum while sampling.\n");
          return 0;
        }
      }
      for(int jj = 0; jj < gdim + 1; jj++) bary00[jj] = ordent[ideg][imin][jj]/(double)ideg;
      if(iprt >= 3){
        printf("   - Try again after sampling got mindist = %15.7e at \n",dmin);
        dblAr1(gdim+1,bary00).print();
      }
      goto restart;
    }else if (irestart < mrestart){
      irestart++;
      double sum = 0;
      for(int jj = 0; jj < gdim+1; jj++){
        bary00[jj] = unif(rng);
        sum += bary00[jj];
      }
      for(int jj = 0; jj < gdim+1; jj++){
        bary00[jj] /= sum;
      }
      if(iprt >= 3){
        printf("   - Restart %d try random bary \n",irestart);
        dblAr1(gdim+1,bary).print();
      }
      goto restart;
    }
  }
  // If a distance had been low enough before, we would have exited already. 
  return 1;
}


// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int inveval_linesearch<1, n>(MeshBase&, int, const double*, double* __restrict__, double* __restrict__, double, int*, bool, int);\
template int inveval_linesearch<2, n>(MeshBase&, int, const double*, double* __restrict__, double* __restrict__, double, int*, bool, int);\
template int inveval_linesearch<3, n>(MeshBase&, int, const double*, double* __restrict__, double* __restrict__, double, int*, bool, int);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




int constrain_bary_desc(int gdim, bool *icstr, double* __restrict__ bary,
                        double* __restrict__ bar0, 
                        double* __restrict__ desc, int iprt){
  for(int ii = 0; ii < gdim + 1; ii++) icstr[ii] = false;
  if(iprt > 0){
    printf("     - test constraints dir = ");
    dblAr1(gdim+1,desc).print();
    printf("     - with bar0 = ");
    dblAr1(gdim+1,bar0).print();
    printf("     - constraints = ");
    MeshArray1D<bool>(gdim+1,icstr).print();

  }

  int ncstr = 0, ncst0;
  do{
    ncst0 = ncstr;
    for(int ii = 0; ii < gdim+1; ii++) bary[ii] = bar0[ii] + desc[ii];

    #ifndef NDEBUG
      double mm1 = abs(desc[0]);
      for(int ii = 1; ii < gdim + 1; ii++) mm1 = MAX(mm1,abs(desc[ii]));
      mm1 = MAX(mm1,1);
    #endif

    // Accrue overflow by hard restraining variables (can be positive or negative)
    double overflow = 0;
    for(int ii = 0; ii < gdim+1; ii++){
      if(icstr[ii]) continue;

      // Either small, or small resulting coefficient. 
      if(bar0[ii] < 1.0e-15*MAX(desc[ii],1) && desc[ii] < 1.0e-15){
        ncstr++;
        icstr[ii] = true;
        overflow += desc[ii];
        desc[ii] = 0;
        if(iprt > 0) printf("     - New constraint (-) %d overflow %23.15e \n",ii,overflow);
      }else if(1 - bar0[ii] < 1.0e-15*MAX(desc[ii],1) && desc[ii] > 1.0e-15){
        ncstr++;
        icstr[ii] = true;
        overflow += desc[ii];
        desc[ii] = 0;
        if(iprt > 0) printf("     - New constraint (+) %d overflow %23.15e \n",ii,overflow);
      } 
    }

    if(ncstr == ncst0) break;

    // The sum of the xcur must remain 0. For this we distribute the overflow among 
    // remaining unconstrained. 
    if(ncstr == gdim+1){
      // This is actually a legitimate situation. See localize_in_this.jpg... 
      // Starting from the top vertex, looking for a point on the (concave) face,
      // the first descent direction has you stuck in all directions. 
      return 1;
    }
    overflow /= (gdim + 1 - ncstr);
    for(int ii = 0; ii < gdim + 1 ; ii++){
      if(icstr[ii]) continue;
      desc[ii] += overflow;
    }
    if(iprt > 0){
      printf("     - desc reajusted to ");
      dblAr1(gdim+1,desc).print();
    }
    #ifndef NDEBUG 
    double sum = 0;
    for(int ii = 0; ii < gdim + 1; ii++) sum+=desc[ii];
    #endif
    METRIS_ASSERT_MSG((abs(sum) < 4.0e-12*mm1), "Sum desc = 0 got = " << abs(sum));
    if(iprt > 0) printf("     - New constraints %d -> %d restart.\n",ncst0,ncstr);
  }while(ncstr > ncst0);
        

  for(int ii = 0; ii < gdim + 1; ii++) bary[ii] = bar0[ii] + desc[ii];

  bool icc = false;


  // Check if new bary (i.e. + desc) exits domain, then apply mult
  bool check = bary[0] > 0 ;
  for(int ii = 1; ii < gdim + 1; ii++) check = check && (bary[ii] > 0);
  if(!check){
    double ccmin = 1.0e30;
    for(int ii = 0; ii < gdim+1; ii++){
      if(icstr[ii]) continue;

      if(iprt > 0) printf("     - Testing bary[%d] = %15.7e for scaling \n",ii,bary[ii]);

      if(bary[ii] < -1.0e-16){
        // If the original iterate is on the boundary, and the descent direction 
        // component is negative, simply constrain this variable. 
        icc = true;
        if(abs(desc[ii]) < 1.0e-16) return 1;
        double cc = - bar0[ii] / desc[ii];
        if(iprt > 0) printf("     - (1) bary %d %23.15e -> %23.15e unconst cc = %15.7e \n",
                                                        ii,bar0[ii],bary[ii],cc);
        ccmin = MIN(ccmin,cc);
      }
      if(1 - bary[ii] < -1.0e-16){
        // If the original iterate is on the boundary, and the descent direction 
        // component is negative, simply constrain this variable. 
        icc = true;
        if(abs(desc[ii]) < 1.0e-16) return 1;
        double cc = (1 - bar0[ii]) / desc[ii];
        if(iprt > 0) printf("     - (2) bary %d %23.15e -> %23.15e unconst cc = %15.7e \n",
                                                        ii,bar0[ii],bary[ii],cc);
        ccmin = MIN(ccmin,cc);
      }
    }
    //if(iprt > 0) printf("  Found cstr: %d %d %d %d coeff = %15.7e apply ? %d\n",
    //                                icstr[0],icstr[1],icstr[2],icstr[3],ccmin,icc);


    // Skip if correction coefficient very close to 1. 
    if(icc && abs(ccmin-1.0) > 1.0e-16){    
      if(ccmin < -1.0e-16) METRIS_THROW_MSG(AlgoExcept(), 
      "Negative or nearly vanishing correction coeff "<<ccmin);
      //if(ccmin - 1.0 > 1.0e-16) METRIS_THROW_MSG(AlgoExcept(), 
      //"> 1 or almost correction coeff "<<ccmin);

      for(int ii = 0; ii < gdim + 1; ii++) desc[ii] *= ccmin;
      if(iprt > 0){
        printf("     - Rescale descent direction by %f -> ",ccmin);
        dblAr1(gdim+1,desc).print();
      }

      // Rescale descent direction 
      //if(iprt > 0) printf("Desc becomes %f %f %f %f new step %15.7e\n", - xcur[0] - xcur[1] - xcur[2],xcur[0],xcur[1],xcur[2],rwrkN[1]);
    }
  }

  return 0;
}



template<int gdim>
int locMeshQuick(MeshBase &msh, const double *coor0){
	for(int ii = 0; ii < gdim; ++ii){
		if(coor0[ii] > msh.bb[ii][1] || coor0[ii] < msh.bb[ii][0]){
			return 1;
		}
	}
	return 0;
}

template int locMeshQuick<1>(MeshBase &msh, const double *coor0);
template int locMeshQuick<2>(MeshBase &msh, const double *coor0);
template int locMeshQuick<3>(MeshBase &msh, const double *coor0);


} // End namespace

