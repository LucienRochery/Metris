//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../adapt/low_delaunay.hxx"

#include "../linalg/matprods.hxx"
#include "../linalg/invmat.hxx"
#include "../linalg/det.hxx"
#include "../aux_exceptions.hxx"
#include "../low_geo.hxx"


namespace Metris{

// If gdim > tdim, we first project into tdim using the normal. 
template <int gdim, int tdim>
bool indelsphere(const double *coop, const double *metl, 
                 const dblAr2 &coord, const int *ent2pol, const double *nrmal){
  static_assert(gdim == tdim || gdim == 3 && tdim == 2);
  if(gdim > tdim){
    METRIS_ASSERT(nrmal != NULL);
  }
  if(tdim != 2) METRIS_THROW_MSG(TODOExcept(), "Implement Delaunay tdim 3");

  double mat[tdim][tdim];
  double rhs[tdim];
  double centr[tdim]; // Sized tdim as expressed in Frénet frame if gdim > tdim
  double nrm1, nrm2, r, r1;

  // Strange as it seems, this is always sized tdim as, if gdim > tdim, we project. 
  double buf[tdim];

  // Only for case gdim > tdim
  double tau1[3], tau2[3], met2[3], buf1[3];


  // Matrix is (P2-P1)^TM
  //           (P3-P1)^TM
  //   (if 3D) (P4-P1)^TM
  if constexpr(tdim == gdim){
    for(int ii = 0; ii < tdim; ii++){
      for(int jj = 0; jj < gdim; jj++) buf[jj] = coord[ent2pol[ii+1]][jj] - coord[ent2pol[0]][jj];
      symXvec<gdim>(metl,buf,mat[ii]);
    }
  }else{
    // In this case, we're going to compute the Frénet frame associated to the 
    // normal. Then write the metric in that frame, as well as project the 
    // above buf vector. 

    // 1. Build the rest of the frame: tau1, tau2
    int imax = -1;
    double rmax = -1;
    for(int ii = 0; ii < gdim; ii++){
      if(abs(nrmal[ii]) > rmax){
        rmax = abs(nrmal[ii]);
        imax = ii;
      }
    }
    // Use the max value as pivot
    tau1[(imax+1)%1] =  nrmal[imax];
    tau1[ imax     ] = -nrmal[(imax+1)%1];
    tau1[(imax+2)%2] =  nrmal[(imax+2)%2];

    double nrm = sqrt(getnrml2<gdim>(nrmal));
    METRIS_ASSERT(nrm > Constants::vecNrmTol);
    tau1[0] /= nrm;
    tau1[1] /= nrm;
    tau1[2] /= nrm;
    METRIS_ASSERT(getprdl2<gdim>(tau1,nrmal) < Constants::vecNrmTol*nrm);

    vecprod(tau1,nrmal,tau2);
    double nrm2 = sqrt(getnrml2<gdim>(tau2));
    tau2[0] /= nrm2;
    tau2[1] /= nrm2;
    tau2[2] /= nrm2;
    METRIS_ASSERT(getprdl2<gdim>(tau1,tau2) < Constants::vecNrmTol);
    METRIS_ASSERT(getprdl2<gdim>(tau2,nrmal) < Constants::vecNrmTol*nrm);

    // Modify the metric as well. This is such that (a b) Q (a b)^T 
    // = (a tau_1^T + b tau_2^T) M (a tau_1 + b tau_2) 
    // for all a,b
    met2[0] = tvecXsymXvec<gdim>(tau1,tau1,metl);
    met2[1] = tvecXsymXvec<gdim>(tau1,tau2,metl);
    met2[2] = tvecXsymXvec<gdim>(tau2,tau2,metl);

    // Compute t1^T P as the components instead of straight P
    for(int ii = 0; ii < tdim; ii++){
      for(int jj = 0; jj < tdim; jj++) 
        buf1[jj] = coord[ent2pol[ii+1]][jj] - coord[ent2pol[0]][jj];
      buf[0] = getprdl2<gdim>(tau1,buf1);
      buf[1] = getprdl2<gdim>(tau2,buf1);
      symXvec<gdim>(met2,buf,mat[ii]);
    }

  }


  //rhs is (||P2||^2 - ||P1||^2)/2
  //        ...
  if constexpr(tdim == gdim){
    nrm1 = tvecXsymXvec<tdim>(coord[ent2pol[0]],coord[ent2pol[0]],metl);
    for(int ii = 0; ii < tdim; ii++){
      nrm2 = tvecXsymXvec<tdim>(coord[ent2pol[ii+1]],coord[ent2pol[ii+1]],metl);
      rhs[ii] = (nrm2 - nrm1)/2;
    }
  }else{
    // Note store P1 proj into buf1
    buf1[0] = getprdl2<gdim>(tau1,coord[ent2pol[0]]);
    buf1[1] = getprdl2<gdim>(tau2,coord[ent2pol[0]]);
    nrm1 = tvecXsymXvec<tdim>(buf1,buf1,met2);
    for(int ii = 0; ii < tdim; ii++){
      buf[0] = getprdl2<gdim>(tau1,coord[ent2pol[ii+1]]);
      buf[1] = getprdl2<gdim>(tau2,coord[ent2pol[ii+1]]);
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
    for(int jj = 0; jj < gdim; jj++) buf[jj] = coord[ent2pol[0]][jj] - centr[jj];
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

  return r1 < r;
}


template bool indelsphere<2,2>(const double *coop, const double *metl, 
                 const dblAr2 &coord, const int *ent2pol, const double* nrmal);
template bool indelsphere<3,2>(const double *coop, const double *metl, 
                 const dblAr2 &coord, const int *ent2pol, const double* nrmal);
template bool indelsphere<3,3>(const double *coop, const double *metl, 
                 const dblAr2 &coord, const int *ent2pol, const double* nrmal);

}//end namespace