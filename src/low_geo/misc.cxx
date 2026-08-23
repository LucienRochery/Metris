//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "misc.hxx"
#include "lenedg.hxx"
#include "normal.hxx"

#include "../linalg/matprods.hxx"
#include "../linalg/invmat.hxx"
#include "../linalg/det.hxx"
#include "../aux_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../utils/mprintf.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"


namespace Metris{



void vdiff_perp(const double* a, const double* b, double *res){
  res[0] = a[1] - b[1];
  res[1] = b[0] - a[0];
}
void vdiff_perp(const double* a, const double* b, int up, int lo, double *res){
  res[0] = up*(a[1] - b[1])/lo; 
  res[1] = up*(b[0] - a[0])/lo;
}
void vdiff_perp_sum(const double* a, const double* b, int up, int lo, double *res){
  res[0] += up*(a[1] - b[1])/lo;
  res[1] += up*(b[0] - a[0])/lo;
}

//double vdiff_perp_x(const double* a, const double* b){
//  return a[1]-b[1]; 
//}
//
//double vdiff_perp_y(const double* a, const double* b){
//  return b[0]-a[0]; 
//}


bool isintetP1(const double *p1, const double *p2,
               const double *p3, const double *p4,
               const double *pp, double tol){
  double vol0 = detvdif3(p2,p1
                         ,p3,p1
                         ,p4,p1);
  vol0 = abs(vol0);
  double vol;

  vol = detvdif3(p2,pp
                 ,p3,pp
                 ,p4,pp);
  if(vol < -tol*vol0) return false;
  vol = detvdif3(pp,p1
                 ,p3,p1
                 ,p4,p1);
  if(vol < -tol*vol0) return false;
  vol = detvdif3(p2,p1
                 ,pp,p1
                 ,p4,p1);
  if(vol < -tol*vol0) return false;
  vol = detvdif3(p2,p1
                 ,p3,p1
                 ,pp,p1);
  if(vol < -tol*vol0) return false;

  return true;
}

bool isinfacP1(const double *p1, const double *p2,
               const double *p3, const double *pp, double tol){
  double are0 = detvdif2(p2,p1,p3,p1);
  are0 = abs(are0);
  double are;

  are = detvdif2(p2,pp,p3,pp);
  if(are < -tol*are0) return false;

  are = detvdif2(pp,p1,p3,p1);
  if(are < -tol*are0) return false;

  are = detvdif2(p2,p1,pp,p1);
  if(are < -tol*are0) return false;

  return true;
}


// gdim == tdim
template<int gdim>
void inventP1(const int*__restrict__ ent2pol, const dblAr2 &coord, 
              const double*__restrict__ pp, 
              double*__restrict__ bary){
  double jmat[gdim][gdim]; 

  for(int ii = 0; ii < gdim; ii++){
    for(int jj = 0; jj < gdim; jj++){
      jmat[ii][jj] = coord(ent2pol[1+ii],jj) - coord(ent2pol[0],jj); 
    }
  }
  //invmat(gdim,jmat[0]);
  METRIS_ENFORCE(!invmat<gdim>(jmat[0]));

  matvdft(gdim,jmat[0],pp,coord[ent2pol[0]],&bary[1]);
//  bary[0] = 1 - bary[1] - bary[2] - bary[3]; 
  bary[0] = 1;
  for(int ii = 1; ii <= gdim; ii++) bary[0] -= bary[ii];
}
template void inventP1<1>(const int*__restrict__ ent2pol, const dblAr2 &coord,const double*__restrict__ pp,
                          double*__restrict__ bary);
template void inventP1<2>(const int*__restrict__ ent2pol, const dblAr2 &coord,const double*__restrict__ pp,
                          double*__restrict__ bary);
template void inventP1<3>(const int*__restrict__ ent2pol, const dblAr2 &coord,const double*__restrict__ pp,
                          double*__restrict__ bary);






// Intrinsic metric
// Metric stored in order 
// 1 2 4 
// * 3 5
// * * 6
// This makes it easier to store 2D metrics if needed
template<int gdim, int tdim, int ideg>
int getintmetxi(const dblAr2 &coord, const int* __restrict__ ent2pol, FEBasis ibasis, 
                const double* bary, double* __restrict__ met){
  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim);

  int ierro = 0;

  double eval[gdim], jmat[gdim*tdim];

  constexpr auto evalf = tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;

  //  Jacobian J verifies: JJ_0^{-1} = M^{-1/2} R 
  //  Multiplying right by the transpose,  J J_0^{-1}J_0^{-T} J^T = M^{-1}
  //  Thus M = J^{-T} J_0^T J_0 J^{-1} 
  //  Now, J_0^T J_0 is none other than (1 1/2 1/2 
  //                                   *  1  1/2) etc 
  //  So its ij-th term is (d_i F, d_j F) = 1 - (i!=j)... 
  
  //  Let's compute instead M^{-1} which is cheaper to invert than J_0, which is not symmetric
  //  For this, compute J J_0^{-1}J_0^{-T} J^T = M^{-1}
  //  
  //  J_0^{-1} J_0^{-T} = 3/2 -1/2 -1/2 ... in 3D, 
  //  4/3 , -2/3 in 2D (-> in 5D, 1/0 +1/0 ? :))
  //  Hence (J_0^TJ_0)^{-1} = hardcoded in the product  

  // 

  evalf(coord,ent2pol,ibasis,DifVar::Bary,DifVar::None,bary,eval,jmat,NULL);
  // ATTENTION Jmat is stored transposed !


  if constexpr(gdim == 2 && tdim == 2){
    met[0] = 4*(jmat[2*0+0]*jmat[2*0+0] + jmat[2*1+0]*jmat[2*1+0])/3
           - 4* jmat[2*0+0]*jmat[2*1+0]/3;
    
    met[1] = 4*(jmat[2*0+1]*jmat[2*0+0] + jmat[2*1+1]*jmat[2*1+0])/3
           - 2*(jmat[2*0+1]*jmat[2*1+0] + jmat[2*0+0]*jmat[2*1+1])/3;

    met[2] = 4*(jmat[2*0+1]*jmat[2*0+1] + jmat[2*1+1]*jmat[2*1+1])/3
           - 4* jmat[2*0+1]*jmat[2*1+1]/3;
    ierro = invspd<gdim>(met);
  
  // This constexpr would seem spurious but I get clang warnings about the size of jmat
  // if I don't add it explicitely:
  }else if constexpr(gdim == 3 && tdim == 3){
    met[0] = 3*(jmat[3*0+0]*jmat[3*0+0] + jmat[3*1+0]*jmat[3*1+0] + jmat[3*2+0]*jmat[3*2+0])/2
           -    jmat[3*0+0]*jmat[3*1+0]
           -    jmat[3*0+0]*jmat[3*2+0]
           -    jmat[3*1+0]*jmat[3*2+0];
    
    met[1] = 3*(jmat[3*0+0]*jmat[3*0+1] + jmat[3*1+0]*jmat[3*1+1] + jmat[3*2+0]*jmat[3*2+1])/2
           -   (jmat[3*0+0]*jmat[3*1+1] + jmat[3*0+1]*jmat[3*1+0])/2
           -   (jmat[3*0+0]*jmat[3*2+1] + jmat[3*0+1]*jmat[3*2+0])/2
           -   (jmat[3*1+0]*jmat[3*2+1] + jmat[3*1+1]*jmat[3*2+0])/2;

    met[2] = 3*(jmat[3*0+1]*jmat[3*0+1] + jmat[3*1+1]*jmat[3*1+1] + jmat[3*2+1]*jmat[3*2+1])/2
           -    jmat[3*0+1]*jmat[3*1+1]
           -    jmat[3*0+1]*jmat[3*2+1]
           -    jmat[3*1+1]*jmat[3*2+1];

    met[3] = 3*(jmat[3*0+0]*jmat[3*0+2] + jmat[3*1+0]*jmat[3*1+2] + jmat[3*2+0]*jmat[3*2+2])/2
           -   (jmat[3*0+0]*jmat[3*1+2] + jmat[3*0+2]*jmat[3*1+0])/2
           -   (jmat[3*0+0]*jmat[3*2+2] + jmat[3*0+2]*jmat[3*2+0])/2
           -   (jmat[3*1+0]*jmat[3*2+2] + jmat[3*1+2]*jmat[3*2+0])/2;

    met[4] = 3*(jmat[3*0+1]*jmat[3*0+2] + jmat[3*1+1]*jmat[3*1+2] + jmat[3*2+1]*jmat[3*2+2])/2
           -   (jmat[3*0+1]*jmat[3*1+2] + jmat[3*0+2]*jmat[3*1+1])/2
           -   (jmat[3*0+1]*jmat[3*2+2] + jmat[3*0+2]*jmat[3*2+1])/2
           -   (jmat[3*1+1]*jmat[3*2+2] + jmat[3*1+2]*jmat[3*2+1])/2;

    met[5] = 3*(jmat[3*0+2]*jmat[3*0+2] + jmat[3*1+2]*jmat[3*1+2] + jmat[3*2+2]*jmat[3*2+2])/2
           -    jmat[3*0+2]*jmat[3*1+2]
           -    jmat[3*0+2]*jmat[3*2+2]
           -    jmat[3*1+2]*jmat[3*2+2];
    ierro = invspd<gdim>(met);
  }else if(gdim == 3 && tdim == 2){
    // This case is a mess! There must be a more elegant way but this seems to work. 
    // The columns of J are T1, T2 (stored transposed here, don't forget)
    // Let t1 = T1 / ||T1||, t2 <- T2 - (T2,t1)t1, t2 <- t2 / ||t2||
    // nu = t1 x t2 (not involved in computations)
    // Seek M = (t1 t2 nu) D (t1 t2 nu)'
    // with D = (d11 d12 0 
    //           d12 d22 0 
    //           0   0   a)
    // st J^T M J = J_0^{-1}J_0 = (4/3 -2/3
    //                             -2/3 4/3)
    // With T1 aligned with t1 hence T1 orthogonal to t2, then 
    // T1^T M T1 = 4/3 = d11 T1^T t1 + 0 -> d11 = 4/(3||T1||^2)
    // T1^TMT2 = d11 (T1^T t1)(T2^T t1)
    //         + d12 (T1^T t1)(T2^T t2) + 0 | <- (T1^T t2 = 0)
    //         + d22 * 0                    |  
    //         = -2/3
    // and T1^T t1 = ||T1||
    // => d12 = ( -2/3 - d11 ||T1|| (T2^T t1) ) / ( ||T1|| (T2^T t2) )

    // Finally T2^T M T2 = d11 (T2^T t1)^2 
    //                   + 2d12 (T2^T t1)(T2^T t2)
    //                   + d22 (T2^T t2)^2  = 4/3
    // thus d22 =  [4/3 - d11 (T2^T t1)^2 - 2d12 (T2^T t1)(T2^T t2)] / 
    //           (T2^T t2)^2 
    double *T1 = &jmat[3*0+0];
    double *T2 = &jmat[3*1+0];

    double sqnrm1 = getnrml2<3>(T1); // don't forget, squared
    double d11 = 4.0/(3.0*sqnrm1);

    double nrm1 = sqrt(sqnrm1); // Now ||T1||

    double t1[3], t2[3]; 
    for(int ii = 0; ii < 3; ii++) t1[ii] = T1[ii]/nrm1;
    double dtprd_t1_T2 = getprdl2<3>(t1,T2);
    for(int ii = 0; ii < 3; ii++) t2[ii] = T2[ii] - dtprd_t1_T2*t1[ii];
    double sqnrm2 = getnrml2<3>(t2);
    double nrm2 = sqrt(sqnrm2);
    for(int ii = 0 ; ii < 3; ii++) t2[ii] /= nrm2;
    // t1, t2 now form an orthonormal basis

    // Compute d12
    double dtprd_t2_T2 = getprdl2<3>(t2,T2);
    double d12 = -2.0/3.0 - d11 * nrm1 * dtprd_t1_T2; 
    d12 /= nrm1*dtprd_t2_T2;

    // Compute d22 
    double d22 = 4.0/3.0 -   d11 * dtprd_t1_T2*dtprd_t1_T2 
                         - 2*d12 * dtprd_t1_T2*dtprd_t2_T2;
    d22 /= dtprd_t2_T2*dtprd_t2_T2;

    // And finally compute the damn matrix 
    // We'll just encode D as a 3x3 sym matrix and use A M A^T

    double nrsz = 10000; // 1 /  sqrt(nrsz) size, do something more clever some day

    double Dm[6] = {0};
    Dm[sym2idx(0,0)] = d11;
    Dm[sym2idx(0,1)] = d12;
    Dm[sym2idx(1,0)] = d12;
    Dm[sym2idx(1,1)] = d22;
    Dm[sym2idx(2,2)] = nrsz; // Arbitrarily, size 1 along the normal

    double nu[3]; 
    vecprod(t1,t2,nu);
    double P[3][3];
    for(int ii = 0; ii < 3; ii++){
      P[ii][0] = t1[ii];
      P[ii][1] = t2[ii];
      P[ii][2] = nu[ii];
    }

    matXsymXtmat<3,3,double,double,double>(Dm,P[0],met);


  }

  return ierro;
}


// This guy is not used and we probably should be using getMetBary / getMetPhys now 
//// This one computes physical derivatives of the metric as well (of \circ F_K^{-1}). 
//template<int ndim, int ideg>
//void getintmetxi(const dblAr2 &coord, const int* __restrict__ tet2pol, FEBasis ibasis, 
//                 const double* bary, Metris::SurrealS<3,double>* __restrict__ metS){
//  static_assert(ndim == 3);
//  double eval[3], jmat[9], hmat[18], djmat[27], invjmat[9];
//
//  //  eval3_bezier<3,ideg,1,0>(coord,tet2pol,bary,eval, jmat,NULL);
//  eval3<3,ideg>(coord,tet2pol,ibasis,DifVar::Bary,DifVar::Bary,bary,eval,jmat,hmat);
//
//
//  // Physical derivative of jmat = dF \circ F_K^{-1} = dF_K^{-1} x_1 d^2 F_K
//  for(int i = 0; i < 9; i++) invjmat[i] = jmat[i];
//  invmat(3,invjmat);
//  mat3X1tens3sym3(hmat,invjmat,djmat);
//
//  Metris::SurrealS<3,double> jmatS[9];
//  for(int ii = 0; ii < 9; ii++){
//    jmatS[ii].value() = jmat[ii];
//  }
//  for(int jj = 0; jj < 3; jj++){
//    for(int ii = 0; ii < 9; ii++){
//      jmatS[ii].deriv(jj) = djmat[9*jj+ii];
//    }
//  }
//  
//  Metris::SurrealS<3,double> invmetS[6];
//
//  invmetS[0] = 3*(jmatS[3*0+0]*jmatS[3*0+0] + jmatS[3*1+0]*jmatS[3*1+0] + jmatS[3*2+0]*jmatS[3*2+0])/2
//             -    jmatS[3*0+0]*jmatS[3*1+0]
//             -    jmatS[3*0+0]*jmatS[3*2+0]
//             -    jmatS[3*1+0]*jmatS[3*2+0];
//  
//  invmetS[1] = 3*(jmatS[3*0+0]*jmatS[3*0+1] + jmatS[3*1+0]*jmatS[3*1+1] + jmatS[3*2+0]*jmatS[3*2+1])/2
//             -   (jmatS[3*0+0]*jmatS[3*1+1] + jmatS[3*0+1]*jmatS[3*1+0])/2
//             -   (jmatS[3*0+0]*jmatS[3*2+1] + jmatS[3*0+1]*jmatS[3*2+0])/2
//             -   (jmatS[3*1+0]*jmatS[3*2+1] + jmatS[3*1+1]*jmatS[3*2+0])/2;
//
//  invmetS[2] = 3*(jmatS[3*0+1]*jmatS[3*0+1] + jmatS[3*1+1]*jmatS[3*1+1] + jmatS[3*2+1]*jmatS[3*2+1])/2
//             -    jmatS[3*0+1]*jmatS[3*1+1]
//             -    jmatS[3*0+1]*jmatS[3*2+1]
//             -    jmatS[3*1+1]*jmatS[3*2+1];
//
//  invmetS[3] = 3*(jmatS[3*0+0]*jmatS[3*0+2] + jmatS[3*1+0]*jmatS[3*1+2] + jmatS[3*2+0]*jmatS[3*2+2])/2
//             -   (jmatS[3*0+0]*jmatS[3*1+2] + jmatS[3*0+2]*jmatS[3*1+0])/2
//             -   (jmatS[3*0+0]*jmatS[3*2+2] + jmatS[3*0+2]*jmatS[3*2+0])/2
//             -   (jmatS[3*1+0]*jmatS[3*2+2] + jmatS[3*1+2]*jmatS[3*2+0])/2;
//
//  invmetS[4] = 3*(jmatS[3*0+1]*jmatS[3*0+2] + jmatS[3*1+1]*jmatS[3*1+2] + jmatS[3*2+1]*jmatS[3*2+2])/2
//             -   (jmatS[3*0+1]*jmatS[3*1+2] + jmatS[3*0+2]*jmatS[3*1+1])/2
//             -   (jmatS[3*0+1]*jmatS[3*2+2] + jmatS[3*0+2]*jmatS[3*2+1])/2
//             -   (jmatS[3*1+1]*jmatS[3*2+2] + jmatS[3*1+2]*jmatS[3*2+1])/2;
//
//  invmetS[5] = 3*(jmatS[3*0+2]*jmatS[3*0+2] + jmatS[3*1+2]*jmatS[3*1+2] + jmatS[3*2+2]*jmatS[3*2+2])/2
//             -    jmatS[3*0+2]*jmatS[3*1+2]
//             -    jmatS[3*0+2]*jmatS[3*2+2]
//             -    jmatS[3*1+2]*jmatS[3*2+2];
//
//  inv3sym<Metris::SurrealS<3,double>>(invmetS, metS);
//
//}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int getintmetxi<2, 2, n>(const dblAr2 &coord, const int* __restrict__ tet2pol,FEBasis ibasis,  \
                                const double* bary,double* __restrict__ met);\
template int getintmetxi<3, 2, n>(const dblAr2 &coord, const int* __restrict__ tet2pol,FEBasis ibasis,  \
                                const double* bary,double* __restrict__ met);\
template int getintmetxi<3, 3, n>(const dblAr2 &coord, const int* __restrict__ tet2pol,FEBasis ibasis,  \
                                const double* bary,double* __restrict__ met);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




// Characteristic element length for tolerance scaling
// Minimum edge length for now
template<int gdim>
double getepsent(MeshBase &msh, int tdim, int ientt){
  // Replace with Frobenius of jmat. Controlling this times eps in dx norm
  // controls dx by eps. 
  double eps;

  //if(tdim == gdim){ // We can do this even in that case but gotta change matprods
  //  // And for now this is for localization
  //  double jmat[gdim*tdim], dum[gdim];
  //  if(tdim == 3){
  //    eval3<gdim,1>(msh.coord,msh.tet2poi[ientt],msh.getBasis(),DifVar::Bary,
  //                  DifVar::None,dum,jmat,NULL);

  //  }else if(tdim == 2){
  //    eval2<gdim,1>(msh.coord,msh.fac2poi[ientt],msh.getBasis(),DifVar::Bary,
  //                  DifVar::None,dum,jmat,NULL);
  //  }else{
  //    eval1<gdim,1>(msh.coord,msh.edg2poi[ientt],msh.getBasis(),DifVar::Bary,
  //                  DifVar::None,dum,jmat,NULL);
  //  }
  //  double jmat2[tdim][tdim];
  //  matXtmat<gdim>(jmat,jmat,jmat2[0]);
  //  eps = 0;
  //  for(int ii = 0; ii < tdim; ii++) eps += jmat2[ii][ii];

  //}else{
    double x2;
    if(tdim == 3){
      eps = geterrl2<gdim>(msh.coord[msh.tet2poi(ientt,0)],msh.coord[msh.tet2poi(ientt,1)]);
      x2 = geterrl2<gdim>(msh.coord[msh.tet2poi(ientt,0)],msh.coord[msh.tet2poi(ientt,2)]);
      eps = eps < x2 ? eps : x2;
      x2 = geterrl2<gdim>(msh.coord[msh.tet2poi(ientt,0)],msh.coord[msh.tet2poi(ientt,3)]);
      eps = eps < x2 ? eps : x2;
      x2 = geterrl2<gdim>(msh.coord[msh.tet2poi(ientt,1)],msh.coord[msh.tet2poi(ientt,2)]);
      eps = eps < x2 ? eps : x2;
      x2 = geterrl2<gdim>(msh.coord[msh.tet2poi(ientt,1)],msh.coord[msh.tet2poi(ientt,3)]);
      eps = eps < x2 ? eps : x2;
      x2 = geterrl2<gdim>(msh.coord[msh.tet2poi(ientt,2)],msh.coord[msh.tet2poi(ientt,3)]);
      eps = eps < x2 ? eps : x2;
    }else if(tdim == 2){
      eps = geterrl2<gdim>(msh.coord[msh.fac2poi(ientt,0)],msh.coord[msh.fac2poi(ientt,1)]);
      x2 = geterrl2<gdim>(msh.coord[msh.fac2poi(ientt,0)],msh.coord[msh.fac2poi(ientt,2)]);
      eps = eps < x2 ? eps : x2;
      x2 = geterrl2<gdim>(msh.coord[msh.fac2poi(ientt,1)],msh.coord[msh.fac2poi(ientt,2)]);
      eps = eps < x2 ? eps : x2;
    }else{
      eps = geterrl2<gdim>(msh.coord[msh.edg2poi(ientt,0)],msh.coord[msh.edg2poi(ientt,1)]);
    }
  //}

  return sqrt(eps); // in both cases
}

template double getepsent<1>(MeshBase &msh, int tdim, int ientt);
template double getepsent<2>(MeshBase &msh, int tdim, int ientt);
template double getepsent<3>(MeshBase &msh, int tdim, int ientt);


} // End namespace




