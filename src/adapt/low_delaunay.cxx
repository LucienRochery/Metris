//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "low_delaunay.hxx"

#include "../Mesh/MeshBase.hxx"
#include "../linalg/matprods.hxx"
#include "../linalg/invmat.hxx"
#include "../linalg/det.hxx"
#include "../aux_exceptions.hxx"
#include "../low_normal.hxx"
#include "../low_geo.hxx"
#include "../utils/mprintf.hxx"


namespace Metris{

// If gdim > tdim, we first project into tdim using the triangle normal. 
template <int gdim, int tdim>
bool indelsphere(const MeshBase &msh, const double *coop, const double *metl, 
                 const int *ent2pol){
  static_assert(gdim == tdim || gdim == 3 && tdim == 2);
  //if(gdim > tdim){
  //  METRIS_ASSERT(nrmal != NULL);
  //}

  GETVDEPTH(msh.param);

  const double orthTol = 1.0e-15; // VecNrmTol = 1.0e-16 doesn't cut it

  double mat[tdim][tdim];
  double rhs[tdim];
  double centr[tdim]; // Sized tdim as expressed in Frénet frame if gdim > tdim
  double nrm1, nrm2, r, r1;

  // Strange as it seems, this is always sized tdim as, if gdim > tdim, we project. 
  double buf[tdim];

  // Only for case gdim > tdim
  double tau1[3], tau2[3], met2[3], buf1[3];

  CPRINTF1("-- Start indelsphere ent2pol = ");
  if(DOPRINTS1()){
    intAr1(tdim+1,ent2pol).print();
  }


  // Matrix is (P2-P1)^TM
  //           (P3-P1)^TM
  //   (if 3D) (P4-P1)^TM
  if constexpr(tdim == gdim){
    for(int ii = 0; ii < tdim; ii++){
      for(int jj = 0; jj < gdim; jj++) 
        buf[jj] = msh.coord(ent2pol[ii+1],jj) - msh.coord(ent2pol[0],jj);
      symXvec<gdim>(metl,buf,mat[ii]);
      //CPRINTF1(" - buf %f %f \n",buf[0],buf[1]);
    }
  }else{
    static_assert(gdim == 3);
    // In this case, we're going to compute the Frénet frame associated to the 
    // normal. Then write the metric in that frame, as well as project the 
    // above buf vector. 

    // 1. Build the rest of the frame: tau1, tau2
    double nrmal[3];
    getnorfacP1(ent2pol, msh.coord, nrmal);

    int imax = -1;
    double rmax = -1;
    for(int ii = 0; ii < gdim; ii++){
      if(abs(nrmal[ii]) > rmax){
        rmax = abs(nrmal[ii]);
        imax = ii;
      }
    }
    // Use the max value as pivot
    tau1[(imax+1)%3] =  nrmal[imax];
    tau1[ imax     ] = -nrmal[(imax+1)%3];
    tau1[(imax+2)%3] =  0;

    double nrm = sqrt(getnrml2<gdim>(tau1));
    tau1[0] /= nrm;
    tau1[1] /= nrm;
    tau1[2] /= nrm;
    #ifndef NDEBUG
    double nrmn = sqrt(getnrml2<gdim>(nrmal));
    METRIS_ASSERT(nrmn > orthTol);
    if(getprdl2<gdim>(tau1,nrmal) >= orthTol*nrmn){
      printf("## Normal norm %15.7e prod w tau1 %15.7e  \n",nrmn,getprdl2<gdim>(tau1,nrmal));
      printf(" normal: ");
      dblAr1(3,nrmal).print();
      printf(" tau1: ");
      dblAr1(3,tau1).print();
      METRIS_THROW(GeomExcept());
    }
    #endif

    vecprod(nrmal,tau1,tau2);
    double nrm2 = sqrt(getnrml2<gdim>(tau2));
    tau2[0] /= nrm2;
    tau2[1] /= nrm2;
    tau2[2] /= nrm2;
    METRIS_ASSERT(getprdl2<gdim>(tau1,tau2) < orthTol);
    #ifndef NDEBUG
    if(getprdl2<gdim>(tau2,nrmal) >= orthTol*nrmn){
      printf("## Normal norm %15.7e prod w tau2 %15.7e rat %15.7e \n",
        nrmn,getprdl2<gdim>(tau2,nrmal),getprdl2<gdim>(tau2,nrmal)/nrmn);
      printf(" normal: ");
      dblAr1(3,nrmal).print();
      printf(" tau2: ");
      dblAr1(3,tau2).print();

      double dtprd = getprdl2<gdim>(tau2,nrmal);
      for(int ii = 0; ii < gdim ;ii++) tau2[ii] -= dtprd*nrmal[ii]/nrmn;
      printf("## -> prod w tau2 %15.7e  \n",getprdl2<gdim>(tau2,nrmal));

      METRIS_THROW(GeomExcept());
    }
    #endif

    CPRINTF1(" - indelsphere nrmal %f %f %f tau1 = %f %f %f tau2 %f %f %f\n",
      nrmal[0],nrmal[1],nrmal[2],tau1[0],tau1[1],tau1[2],tau2[0],tau2[1],tau2[2]);


    // Modify the metric as well. This is such that (a b) Q (a b)^T 
    // = (a tau_1^T + b tau_2^T) M (a tau_1 + b tau_2) 
    // for all a,b
    met2[0] = tvecXsymXvec<gdim>(tau1,tau1,metl);
    met2[1] = tvecXsymXvec<gdim>(tau1,tau2,metl);
    met2[2] = tvecXsymXvec<gdim>(tau2,tau2,metl);


    CPRINTF1(" - indelsphere metl %f %f %f %f %f %f det = %15.7e met2 %f %f %f det = %15.7e \n",
      metl[0],metl[1],metl[2],metl[3],metl[4],metl[5],detsym<3>(metl),met2[0],met2[1],met2[2],detsym<2>(met2));

    // Compute t1^T P as the components instead of straight P
    for(int ii = 0; ii < tdim; ii++){
      for(int jj = 0; jj < gdim; jj++) 
        buf1[jj] = msh.coord(ent2pol[ii+1],jj) - msh.coord(ent2pol[0],jj);
      buf[0] = getprdl2<gdim>(tau1,buf1);
      buf[1] = getprdl2<gdim>(tau2,buf1);
      symXvec<tdim>(met2,buf,mat[ii]);
      CPRINTF1(" - buf1 = %f %f buf %f %f \n",buf1[0],buf1[1],buf[0],buf[1]);
    }

  }

  if(DOPRINTS1()){
    if constexpr (tdim == 2){
      CPRINTF1(" - indelsphere mat %f %f %f %f det = %15.7e\n",
        mat[0][0],mat[0][1],mat[1][0],mat[1][1],detmat<2>(mat[0]));
    }else{
      CPRINTF1(" - indelsphere mat %f %f %f ; %f %f %f ; %f %f %f det = %15.7e\n",
        mat[0][0],mat[0][1],mat[0][2],
        mat[1][0],mat[1][1],mat[1][2],
        mat[2][0],mat[2][1],mat[2][2],detmat<3>(mat[0]));
    }
  }


  //rhs is (||P2||^2 - ||P1||^2)/2
  //        ...
  if constexpr(tdim == gdim){
    nrm1 = tvecXsymXvec<tdim>(msh.coord[ent2pol[0]],msh.coord[ent2pol[0]],metl);
    for(int ii = 0; ii < tdim; ii++){
      nrm2 = tvecXsymXvec<tdim>(msh.coord[ent2pol[ii+1]],msh.coord[ent2pol[ii+1]],metl);
      rhs[ii] = (nrm2 - nrm1)/2;
    }
  }else{
    // Note store P1 proj into buf1
    buf1[0] = getprdl2<gdim>(tau1,msh.coord[ent2pol[0]]);
    buf1[1] = getprdl2<gdim>(tau2,msh.coord[ent2pol[0]]);
    nrm1 = tvecXsymXvec<tdim>(buf1,buf1,met2);
    for(int ii = 0; ii < tdim; ii++){
      buf[0] = getprdl2<gdim>(tau1,msh.coord[ent2pol[ii+1]]);
      buf[1] = getprdl2<gdim>(tau2,msh.coord[ent2pol[ii+1]]);
      nrm2 = tvecXsymXvec<tdim>(buf,buf,met2);
      rhs[ii] = (nrm2 - nrm1)/2;
    }
  }

  // C is solution of mat*C = rhs
  //invmat(gdim, mat[0]);
  if(invmat<tdim>(mat[0])){
    #ifndef NDEBUG
    METRIS_THROW_MSG(GeomExcept(), "Invmat failed Delaunay")
    #endif
    return false;
  }
  matXvec<tdim>(mat[0], rhs, centr);

  if constexpr(tdim == gdim){
    for(int jj = 0; jj < gdim; jj++) buf[jj] = msh.coord(ent2pol[0],jj) - centr[jj];
    r = tvecXsymXvec<tdim>(buf,buf,metl);

    for(int jj = 0; jj < gdim; jj++) buf[jj] = coop[jj] - centr[jj];
    r1 = tvecXsymXvec<tdim>(buf,buf,metl);
  }else{
    // We stored P1 proj into buf1 previously. centr is already in frame tau1 tau2
    for(int jj = 0; jj < tdim; jj++) buf[jj] = buf1[jj] - centr[jj];
    r = tvecXsymXvec<tdim>(buf,buf,metl);

    // Project coop as well
    buf1[0] = getprdl2<gdim>(tau1,coop);
    buf1[1] = getprdl2<gdim>(tau2,coop);
    for(int jj = 0; jj < tdim; jj++) buf[jj] = buf1[jj] - centr[jj];
    r1 = tvecXsymXvec<tdim>(buf,buf,metl);
  }

  CPRINTF1("-- END indelsphere %d %d  r = %f r1 = %f centr ",gdim,tdim,r,r1);
  if(DOPRINTS1()){
    dblAr1(tdim,centr).print();
  }

  return r1 < r;
}


template bool indelsphere<2,2>(const MeshBase &msh, const double *coop, const double *metl, 
                               const int *ent2pol);
template bool indelsphere<3,2>(const MeshBase &msh, const double *coop, const double *metl, 
                               const int *ent2pol);
template bool indelsphere<3,3>(const MeshBase &msh, const double *coop, const double *metl, 
                               const int *ent2pol);

}//end namespace