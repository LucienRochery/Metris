//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../adapt/low_delaunay.hxx"

#include "../linalg/matprods.hxx"
#include "../linalg/invmat.hxx"
#include "../aux_exceptions.hxx"
#include "../low_geo.hxx"


namespace Metris{

template <int gdim>
bool indelsphere(const double *coop, const double *metl, 
                 const dblAr2 &coord, const int *ent2pol){
  if(gdim != 2) METRIS_THROW_MSG(TODOExcept(), "Implement dim 3");


  double mat[gdim][gdim];
  double rhs[gdim];
  double centr[gdim];
  double nrm1, nrm2, r, r1;

  double buf[gdim];

  #ifndef NDEBUG
  double buf0[gdim];
  try{
  #endif
  // Matrix is (P2-P1)^TM
  //           (P3-P1)^TM
  //   (if 3D) (P4-P1)^TM
  for(int ii = 0; ii < gdim; ii++){
    for(int jj = 0; jj < gdim; jj++) buf[jj] = coord[ent2pol[ii+1]][jj] - coord[ent2pol[0]][jj];
    symXvec<gdim>(metl,buf,mat[ii]);
  }


  #ifndef NDEBUG
    for(int ii = 0; ii < gdim; ii++) buf0[ii] = buf[ii];
  #endif

  //rhs is (||P2||^2 - ||P1||^2)/2
  //        ...
  symXvec<gdim>(metl,coord[ent2pol[0]],buf);
  nrm1 = getprdl2<gdim>(buf,coord[ent2pol[0]]);
  for(int ii = 0; ii < gdim; ii++){
    symXvec<gdim>(metl,coord[ent2pol[ii+1]],buf);
    nrm2 = getprdl2<gdim>(buf,coord[ent2pol[ii+1]]);
    rhs[ii] = (nrm2 - nrm1)/2;
  }

  // C is solution of mat*C = rhs
  //invmat(gdim, mat[0]);
  if(invmat<gdim>(mat[0])){
    #ifndef NDEBUG
    METRIS_THROW_MSG(GeomExcept(), "Invmat failed Delaunay")
    #endif
    return false;
  }
  matXvec<gdim>(mat[0], rhs, centr);

  for(int jj = 0; jj < gdim; jj++) buf[jj] = coord[ent2pol[0]][jj] - centr[jj];
  symXvec<gdim>(metl,buf,mat[0]);
  r = getprdl2<gdim>(buf,mat[0]);


  for(int jj = 0; jj < gdim; jj++) buf[jj] = coop[jj] - centr[jj];
  symXvec<gdim>(metl,buf,mat[0]);
  r1 = getprdl2<gdim>(buf,mat[0]);

  //#ifndef NDEBUG
  //  for(int ii = 1; ii <= gdim; ii++){
  //    for(int jj = 0; jj < gdim; jj++) buf[jj] = coord[ent2pol[ii]][jj] - centr[jj];
  //    symXvec<gdim>(metl,buf,mat[0]);
  //    double r2 = getprdl2<gdim>(buf,mat[0]);
  //    METRIS_ASSERT(abs(r2 - r) <= 1.0e-12);
  //  }
  //#endif


  #ifndef NDEBUG
  }catch(const MetrisExcept &e){
    constexpr int nnmet = (gdim * (gdim + 1)) / 2;
    printf("indelsphere exception raised\n");
    printf("Inp coop = ");
    dblAr1(gdim, coop).print();

    printf("Inp metl = ");
    dblAr1(nnmet, metl).print();

    printf("buf0 = ");
    dblAr1(gdim, buf0).print();

    printf("nrm1 = %15.7e rhs = ",nrm1);
    dblAr1(gdim, rhs).print();

    throw(e);
  }
  #endif

  return r1 < r;
}


template bool indelsphere<2>(const double *coop, const double *metl, 
                 const dblAr2 &coord, const int *ent2pol);
template bool indelsphere<3>(const double *coop, const double *metl, 
                 const dblAr2 &coord, const int *ent2pol);

}//end namespace