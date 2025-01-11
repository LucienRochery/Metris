//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_EVAL__
#define __LOW_EVAL__

#include "types.hxx"
#include "ho_constants.hxx"
#include "metris_constants.hxx"

#include "codegen_lagrange.hxx"


namespace Metris{

  
template<int n> inline double idpow(double x);

/*
~~ low_eval.hxx: polynomial evaluation routines ~~ 

~~ Main functions:
  - eval1, eval2, eval3: evaluate topo dim 1,2,3 polynomial
  - eval_lagrangefunc, eval_bezierfunc: evaluate single Lagrange or Bernstein poly

~~ Common arguments (or template params) are:
 - szfld: field size, e.g. R (scalar field), R^3 (point field e.g. control points), or R^6 (metric field)
 - ideg: degree 
 - idif1: differentiate once (Jacobian matrix) 
 - idif2: differentiate twice (D2 order 3 tensor). Not available in Lagrange yet (TODO for consistency but not need identified yet). 
 - rfld: arbitrary size unordered coefficients field (e.g. control points)
 - lfld: points onto rfld, sized depending on topo dim and degree
 - bary: tdim+1 barycentric coordinates
 - eval: array of size szfld, the evaluation at bary
 - jmat: array of size tdim*szfld. 1D but (C convention)  jmat[i][j] = d_i eval[j]
 - hmat: array of size [(tdim+1)*(tdim+2)/2]*szfld. 1D but (C convention)  hmat[i][j] = s_i eval[j] 
 where s_i designates some d_jd_k ordered the same as other symmetric mats (e.g. metrics):
    1: d_11
    2: d_12
    3: d_22
    4: d_13
    5: d_23
    6: d_33


~~ Note on derivatives: 
A tdim poly P is of tdim+1 barycentric variables. Only tdim physical derivatives are given:
(phys) d_i  | d_{i+1} - d_1 (bary)
This corresponds to the reference element:
P_1 = 0
P_i = e_{i+1}
with e_i the canonical basis vectors.

~~ Performance:
First implementation involved repeated calls to eval_lagrangefunc and eval_bezierfunc. 
These are runtime evaluated. 
Similarly, runtime recursive Casteljau was implemented, it is at least 10 times slower 
as a recursive template implementation. See docs/Performance.md for more details. 
*/

template <int szfld, int ideg>
void eval3(const dblAr2& __restrict__ rfld, 
           const int*__restrict__ lfld,
           FEBasis ibasis, DifVar idif1, DifVar idif2, 
           const double*__restrict__ bary, double*__restrict__ eval, 
           double*__restrict__ jmat, double*__restrict__ hmat);
template <int szfld, int ideg>
void eval2(const dblAr2& __restrict__ rfld, 
           const int*__restrict__ lfld,
           FEBasis ibasis, DifVar idif1, DifVar idif2, 
           const double*__restrict__ bary, double*__restrict__ eval, 
           double*__restrict__ jmat, double*__restrict__ hmat);
template <int szfld, int ideg>
void eval1(const dblAr2& __restrict__ rfld, 
           const int*__restrict__ lfld,
           FEBasis ibasis, DifVar idif1, DifVar idif2, 
           const double*__restrict__ bary, double*__restrict__ eval, 
           double*__restrict__ jmat, double*__restrict__ hmat);


// These shouldn't be called directly. For evaluation of a single Lag or Béz function. 
template <int ideg,int tdim>
double eval_lagrangefunc(const int*__restrict__ idx, const double*__restrict__ bary,int ider,double*__restrict__ dlag);
template <int ideg,int tdim>
double eval_bezierfunc(const int*__restrict__ idx, const double*__restrict__ bary,int ider,double*__restrict__ dlag);



/*
With Hessian computed from Jmats
The Hessian is compressed as column first as usual:
1 2 4
  3 5
    6
However, note that each entry is, itself, comprised of szfld values. 
*/
template <int szfld, int ideg, int di=0,int dj=0, int dk=0, int dl=0>
struct eval3_bezier{
  eval3_bezier(const dblAr2 & __restrict__  rfld,  
               const int * __restrict__  lfld,
               DifVar idif1, DifVar idif2,
               const double  * __restrict__  bary,  
               double * __restrict__  eval,
               double * __restrict__ jmat,
               double * __restrict__ hmat){
    double eva1[szfld], eva2[szfld], eva3[szfld], eva4[szfld];
    double jmat1[3*szfld];
    double jmat2[3*szfld];
    double jmat3[3*szfld];
    double jmat4[3*szfld];


    eval3_bezier<szfld,ideg-1,di+1,dj+0,dk+0,dl+0>(rfld,lfld,idif2,DifVar::None,bary,eva1,jmat1,NULL);
    eval3_bezier<szfld,ideg-1,di+0,dj+1,dk+0,dl+0>(rfld,lfld,idif2,DifVar::None,bary,eva2,jmat2,NULL);
    eval3_bezier<szfld,ideg-1,di+0,dj+0,dk+1,dl+0>(rfld,lfld,idif2,DifVar::None,bary,eva3,jmat3,NULL);
    eval3_bezier<szfld,ideg-1,di+0,dj+0,dk+0,dl+1>(rfld,lfld,idif2,DifVar::None,bary,eva4,jmat4,NULL);


    for(int i = 0; i < szfld; i++){
      eval[i] = bary[0] * eva1[i]
              + bary[1] * eva2[i]
              + bary[2] * eva3[i]
              + bary[3] * eva4[i];
    }

    if(idif1 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        double d1 = eva1[i];
        jmat[0*szfld+i] = ideg * (eva2[i] - d1);
        jmat[1*szfld+i] = ideg * (eva3[i] - d1);
        jmat[2*szfld+i] = ideg * (eva4[i] - d1);
      }
    }

    if(idif2 == DifVar::Bary){
      for(int i=0; i<szfld; i++){
                  // d11 (i.e. (d2-d1)(d2-d1))
        hmat[0*szfld + i] = ideg * (jmat2[0*szfld+i] - jmat1[0*szfld+i]);
                  // d21 
        hmat[1*szfld + i] = ideg * (jmat3[0*szfld+i] - jmat1[0*szfld+i]);
                  // d31 
        hmat[3*szfld + i] = ideg * (jmat4[0*szfld+i] - jmat1[0*szfld+i]);
  
                  // d22 
        hmat[2*szfld + i] = ideg * (jmat3[1*szfld+i] - jmat1[1*szfld+i]);
                  // d32 
        hmat[4*szfld + i] = ideg * (jmat4[1*szfld+i] - jmat1[1*szfld+i]);
  
                  // d33 
        hmat[5*szfld + i] = ideg * (jmat4[2*szfld+i] - jmat1[2*szfld+i]);
      }
    }
  }
};

template <int szfld, int di,int dj, int dk, int dl>
struct eval3_bezier<szfld,1,di,dj,dk,dl>{
  eval3_bezier(const dblAr2 & __restrict__  rfld,  
               const int * __restrict__  lfld,
               DifVar idif1, DifVar idif2, 
               const double  * __restrict__  bary,  
               double * __restrict__  eval,
               double * __restrict__ jmat,
               double * __restrict__ hmat){

   // constexpr int ideg = 1 + di + dj + dk + dl; // The true degree. 

    constexpr int i1000 = mul2nod(1+di,0+dj,0+dk,0+dl);
    constexpr int i0100 = mul2nod(0+di,1+dj,0+dk,0+dl);
    constexpr int i0010 = mul2nod(0+di,0+dj,1+dk,0+dl);
    constexpr int i0001 = mul2nod(0+di,0+dj,0+dk,1+dl);

//          printf("  - ideg = %d with (di,dj,dk,dl) = (%d,%d,%d,%d)\n",ideg,di,dj,dk,dl);
//          printf("  - Vertex indices %d %d %d %d \n",i1000,i0100,i0010,i0001);

    for(int i = 0; i < szfld; i++){
      eval[i] = bary[0]*rfld[lfld[i1000]][i]
              + bary[1]*rfld[lfld[i0100]][i]
              + bary[2]*rfld[lfld[i0010]][i]
              + bary[3]*rfld[lfld[i0001]][i];
    }

    //if(eval[0] > 1.0e10){
    //  printf("## 1 DEBUG EVAL3_BEZIER VERY LARGE VALUE %20.16e \n",eval[0]);
    //  printf("rfld 1000 0100 0010 0001 %20.16e %20.16e %20.16e %20.16e \n",rfld[lfld[i1000]][0],
    //    rfld[lfld[i0100]][0],rfld[lfld[i0010]][0],rfld[lfld[i0001]][0]);
    //  wait();
    //}

    if(idif1 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        jmat[0*szfld+i] = rfld[lfld[i0100]][i] - rfld[lfld[i1000]][i];
        jmat[1*szfld+i] = rfld[lfld[i0010]][i] - rfld[lfld[i1000]][i];
        jmat[2*szfld+i] = rfld[lfld[i0001]][i] - rfld[lfld[i1000]][i];
      }
    }

    if(idif2 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        hmat[0*szfld+i] = 0.0;
        hmat[1*szfld+i] = 0.0;
        hmat[2*szfld+i] = 0.0;
        hmat[3*szfld+i] = 0.0;
        hmat[4*szfld+i] = 0.0;
        hmat[5*szfld+i] = 0.0;
      }
    }

  }
};


template <int szfld,int di,int dj, int dk, int dl>
  struct eval3_bezier<szfld,2,di,dj,dk,dl>{
    eval3_bezier(const dblAr2 & __restrict__  rfld,  
                 const int * __restrict__ lfld,
                 DifVar idif1,DifVar idif2,
                 const double  * __restrict__ bary,  
                 double * __restrict__  eval,
                 double * __restrict__ jmat,
                 double * __restrict__ hmat){
    //constexpr int ideg = 2 + di + dj + dk + dl; // The true degree. 
    
    if(idif1 == DifVar::Bary){


        // The Jacobian is quicker to compute recursively as it only uses the previous degree
        // evals. 

      double eva1[szfld], eva2[szfld], eva3[szfld], eva4[szfld];
      double jmat1[3*szfld];
      double jmat2[3*szfld];
      double jmat3[3*szfld];
      double jmat4[3*szfld];
    
      eval3_bezier<szfld,1,di+1,dj+0,dk+0,dl+0>(rfld,lfld,idif2,DifVar::None,bary,eva1,jmat1,NULL);
      eval3_bezier<szfld,1,di+0,dj+1,dk+0,dl+0>(rfld,lfld,idif2,DifVar::None,bary,eva2,jmat2,NULL);
      eval3_bezier<szfld,1,di+0,dj+0,dk+1,dl+0>(rfld,lfld,idif2,DifVar::None,bary,eva3,jmat3,NULL);
      eval3_bezier<szfld,1,di+0,dj+0,dk+0,dl+1>(rfld,lfld,idif2,DifVar::None,bary,eva4,jmat4,NULL);

      for(int i = 0; i < szfld; i++){
        eval[i] = bary[0] * eva1[i]
                + bary[1] * eva2[i]
                + bary[2] * eva3[i]
                + bary[3] * eva4[i];
      }


      for(int i = 0; i < szfld;i++){
        double d1 = eva1[i];
        jmat[0*szfld+i] = 2 * (eva2[i] - d1);
        jmat[1*szfld+i] = 2 * (eva3[i] - d1);
        jmat[2*szfld+i] = 2 * (eva4[i] - d1);
      }


      if(idif2 == DifVar::Bary){
        for(int i=0; i<szfld; i++){
                    // d11 (i.e. (d2-d1)(d2-d1))
          hmat[0*szfld + i] = 2 * (jmat2[0*szfld+i] - jmat1[0*szfld+i]);
                    // d21 
          hmat[1*szfld + i] = 2 * (jmat3[0*szfld+i] - jmat1[0*szfld+i]);
                    // d31 
          hmat[3*szfld + i] = 2 * (jmat4[0*szfld+i] - jmat1[0*szfld+i]);
    
                    // d22 
          hmat[2*szfld + i] = 2 * (jmat3[1*szfld+i] - jmat1[1*szfld+i]);
                    // d32 
          hmat[4*szfld + i] = 2 * (jmat4[1*szfld+i] - jmat1[1*szfld+i]);
    
                    // d33 
          hmat[5*szfld + i] = 2 * (jmat4[2*szfld+i] - jmat1[2*szfld+i]);
        }
      }


    }else{

      constexpr int i2000 = mul2nod(2+di,0+dj,0+dk,0+dl);
      constexpr int i0200 = mul2nod(0+di,2+dj,0+dk,0+dl);
      constexpr int i0020 = mul2nod(0+di,0+dj,2+dk,0+dl);
      constexpr int i0002 = mul2nod(0+di,0+dj,0+dk,2+dl);
      constexpr int i1100 = mul2nod(1+di,1+dj,0+dk,0+dl);
      constexpr int i0110 = mul2nod(0+di,1+dj,1+dk,0+dl);
      constexpr int i1010 = mul2nod(1+di,0+dj,1+dk,0+dl);
      constexpr int i1001 = mul2nod(1+di,0+dj,0+dk,1+dl);
      constexpr int i0101 = mul2nod(0+di,1+dj,0+dk,1+dl);
      constexpr int i0011 = mul2nod(0+di,0+dj,1+dk,1+dl);

      for(int i = 0; i < szfld; i++){
        eval[i] =   bary[0]*bary[0]*rfld[lfld[i2000]][i] 
                + 2*bary[0]*bary[1]*rfld[lfld[i1100]][i] 
                + 2*bary[0]*bary[2]*rfld[lfld[i1010]][i] 
                + 2*bary[0]*bary[3]*rfld[lfld[i1001]][i]
                +   bary[1]*bary[1]*rfld[lfld[i0200]][i] 
                + 2*bary[1]*bary[2]*rfld[lfld[i0110]][i] 
                + 2*bary[1]*bary[3]*rfld[lfld[i0101]][i]
                +   bary[2]*bary[2]*rfld[lfld[i0020]][i] 
                + 2*bary[2]*bary[3]*rfld[lfld[i0011]][i]
                +   bary[3]*bary[3]*rfld[lfld[i0002]][i] ;
      }


    }
  }
};
//------------------------------------------------------------------------------------------------
//------------------------------         TRIANGLE         ----------------------------------------
//------------------------------------------------------------------------------------------------


//// ideg is used externally to signify the degree
// internally, it is decreased as a stop condition equivalent to ideg == di + dj + dk + dl 
// which cannot be directly enforced by template specialization.
template <int szfld, int ideg, int di=0,int dj=0, int dk=0>
struct eval2_bezier{
  eval2_bezier(const dblAr2 & __restrict__  rfld,  
               const int * __restrict__  lfld,
               DifVar idif1, DifVar idif2,
               const double  * __restrict__  bary,  
               double * __restrict__  eval,
               double * __restrict__ jmat,
               double * __restrict__ hmat){
    double eva1[szfld], eva2[szfld], eva3[szfld];
    double jmat1[2*szfld], jmat2[2*szfld], jmat3[2*szfld];

    eval2_bezier<szfld,ideg-1,di+1,dj+0,dk+0>(rfld,lfld,idif2,DifVar::None,bary,eva1,jmat1,NULL);
    eval2_bezier<szfld,ideg-1,di+0,dj+1,dk+0>(rfld,lfld,idif2,DifVar::None,bary,eva2,jmat2,NULL);
    eval2_bezier<szfld,ideg-1,di+0,dj+0,dk+1>(rfld,lfld,idif2,DifVar::None,bary,eva3,jmat3,NULL);

    for(int i = 0; i < szfld; i++){
      eval[i] = bary[0] * eva1[i]
              + bary[1] * eva2[i]
              + bary[2] * eva3[i];
    }
    if(idif1 == DifVar::Bary){
      for(int i = 0; i < szfld;i++){
        double d1 = eva1[i];
        jmat[0*szfld+i] = ideg * (eva2[i] - d1);
        jmat[1*szfld+i] = ideg * (eva3[i] - d1);
      }
    }

    if(idif2 == DifVar::Bary){
      for(int i=0; i<szfld; i++){
                  // d11 (i.e. (d2-d1)(d2-d1))
        hmat[0*szfld + i] = ideg * (jmat2[0*szfld+i] - jmat1[0*szfld+i]);
                  // d21 
        hmat[1*szfld + i] = ideg * (jmat3[0*szfld+i] - jmat1[0*szfld+i]);
                  // d22 
        hmat[2*szfld + i] = ideg * (jmat3[1*szfld+i] - jmat1[1*szfld+i]);
      }
    }

  }
};

template <int szfld,int di,int dj, int dk>
struct eval2_bezier<szfld, 1,di,dj,dk>{
  eval2_bezier(const dblAr2 & __restrict__  rfld,  
               const int * __restrict__  lfld,
               DifVar idif1, DifVar idif2,
               const double  * __restrict__  bary,  
               double * __restrict__  eval,
               double * __restrict__ jmat,
               double * __restrict__ hmat){
    //constexpr int ideg = 1 + di + dj + dk; // The true degree. 

    constexpr int i100 = mul2nod(1+di,0+dj,0+dk);
    constexpr int i010 = mul2nod(0+di,1+dj,0+dk);
    constexpr int i001 = mul2nod(0+di,0+dj,1+dk);

    for(int i = 0; i < szfld; i++){
      eval[i] = bary[0]*rfld[lfld[i100]][i]
              + bary[1]*rfld[lfld[i010]][i]
              + bary[2]*rfld[lfld[i001]][i];
    }


    if(idif1 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        jmat[0*szfld+i] = rfld[lfld[i010]][i] - rfld[lfld[i100]][i];
        jmat[1*szfld+i] = rfld[lfld[i001]][i] - rfld[lfld[i100]][i];
      }
    }

    if(idif2 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        hmat[0*szfld+i] = 0.0;
        hmat[1*szfld+i] = 0.0;
        hmat[2*szfld+i] = 0.0;
      }
    }

  }
};


template <int szfld, int di,int dj, int dk>
struct eval2_bezier<szfld,2,di,dj,dk>{
  eval2_bezier(const dblAr2 & __restrict__  rfld,  
               const int * __restrict__  lfld,
               DifVar idif1, DifVar idif2,
               const double  * __restrict__  bary,  
               double * __restrict__  eval,
               double * __restrict__ jmat,
               double * __restrict__ hmat){
    //constexpr int ideg = 2 + di + dj + dk ; // The true degree. 
    
    if(idif1 == DifVar::Bary){

        // The Jacobian is quicker to compute recursively as it only uses the previous degree
        // evals. 

      double eva1[szfld], eva2[szfld], eva3[szfld];
      double jmat1[2*szfld], jmat2[2*szfld], jmat3[2*szfld];

      eval2_bezier<szfld,1,di+1,dj+0,dk+0>(rfld,lfld,idif2,DifVar::None,bary,eva1,jmat1,NULL);
      eval2_bezier<szfld,1,di+0,dj+1,dk+0>(rfld,lfld,idif2,DifVar::None,bary,eva2,jmat2,NULL);
      eval2_bezier<szfld,1,di+0,dj+0,dk+1>(rfld,lfld,idif2,DifVar::None,bary,eva3,jmat3,NULL);

      for(int i = 0; i < szfld; i++){
        eval[i] = bary[0] * eva1[i]
                + bary[1] * eva2[i]
                + bary[2] * eva3[i];
      }

      for(int i = 0; i < szfld;i++){
        double d1 = eva1[i];
        jmat[0*szfld+i] = 2 * (eva2[i] - d1);
        jmat[1*szfld+i] = 2 * (eva3[i] - d1);
      }

      if(idif2 == DifVar::Bary){
        for(int i=0; i<szfld; i++){
                    // d11 (i.e. (d2-d1)(d2-d1))
          hmat[0*szfld + i] = 2 * (jmat2[0*szfld+i] - jmat1[0*szfld+i]);
                    // d21 
          hmat[1*szfld + i] = 2 * (jmat3[0*szfld+i] - jmat1[0*szfld+i]);
                    // d22 
          hmat[2*szfld + i] = 2 * (jmat3[1*szfld+i] - jmat1[1*szfld+i]);
        }
      }
      
    }else{

      constexpr int i200 = mul2nod(2+di,0+dj,0+dk);
      constexpr int i020 = mul2nod(0+di,2+dj,0+dk);
      constexpr int i002 = mul2nod(0+di,0+dj,2+dk);
      constexpr int i110 = mul2nod(1+di,1+dj,0+dk);
      constexpr int i011 = mul2nod(0+di,1+dj,1+dk);
      constexpr int i101 = mul2nod(1+di,0+dj,1+dk);


      for(int i = 0; i < szfld; i++){
        eval[i] =   bary[0]*bary[0]*rfld[lfld[i200]][i] 
                + 2*bary[0]*bary[1]*rfld[lfld[i110]][i] 
                + 2*bary[0]*bary[2]*rfld[lfld[i101]][i] 
                +   bary[1]*bary[1]*rfld[lfld[i020]][i] 
                + 2*bary[1]*bary[2]*rfld[lfld[i011]][i] 
                +   bary[2]*bary[2]*rfld[lfld[i002]][i] ;
      }
    }
  }
};






//------------------------------------------------------------------------------------------------
//------------------------------            EDGE          ----------------------------------------
//------------------------------------------------------------------------------------------------

//// ideg is used externally to signify the degree
// internally, it is decreased as a stop condition equivalent to ideg == di + dj + dk + dl 
// which cannot be directly enforced by template specialization.
template <int szfld, int ideg, int di=0,int dj=0>
struct eval1_bezier{
  eval1_bezier(const dblAr2 & __restrict__  rfld,  
               const int * __restrict__  lfld,
               DifVar idif1, DifVar idif2,
               const double  * __restrict__  bary,  
               double * __restrict__  eval,
               double * __restrict__ jmat,
               double * __restrict__ hmat){
    double eva1[szfld], eva2[szfld];
    double jmat1[szfld];
    double jmat2[szfld];

    eval1_bezier<szfld,ideg-1,di+1,dj+0>(rfld,lfld,idif2,DifVar::None,bary,eva1,jmat1,NULL);
    eval1_bezier<szfld,ideg-1,di+0,dj+1>(rfld,lfld,idif2,DifVar::None,bary,eva2,jmat2,NULL);

    for(int i = 0; i < szfld; i++){
      eval[i] = bary[0]*eva1[i] + bary[1]*eva2[i];
    }

    if(idif1 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        jmat[i] = ideg * (eva2[i] - eva1[i]);
      }
    }
    if(idif2 == DifVar::Bary){
      for(int i=0; i<szfld; i++){
                  // d11 (i.e. (d2-d1)(d2-d1))
        hmat[0*szfld + i] = ideg * (jmat2[0*szfld+i] - jmat1[0*szfld+i]);
      }
    }

  }
};

template <int szfld, int di,int dj>
struct eval1_bezier<szfld,1,di,dj>{
  eval1_bezier(const dblAr2 & __restrict__  rfld,  
               const int * __restrict__  lfld,
               DifVar idif1, DifVar idif2,
               const double  * __restrict__  bary,  
               double * __restrict__  eval,
               double * __restrict__ jmat,
               double * __restrict__ hmat){

    constexpr int i10 = mul2nod(1+di,0+dj);
    constexpr int i01 = mul2nod(0+di,1+dj);

    for(int i = 0; i < szfld; i++){
      eval[i] = bary[0]*rfld[lfld[i10]][i] 
              + bary[1]*rfld[lfld[i01]][i];
    }

    if(idif1 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        jmat[i] = rfld[lfld[i01]][i] - rfld[lfld[i10]][i];
      }
    }

    if(idif2 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        hmat[0*szfld+i] = 0.0;
      }
    }

  }
};


template <int szfld, int di,int dj>
struct eval1_bezier<szfld,2,di,dj>{
  eval1_bezier(const dblAr2 & __restrict__  rfld,  
               const int * __restrict__  lfld,
               DifVar idif1, DifVar idif2,
               const double  * __restrict__  bary,  
               double * __restrict__  eval,
               double * __restrict__ jmat,
               double * __restrict__ hmat){
    constexpr int i20 = mul2nod(2+di,0+dj);
    constexpr int i02 = mul2nod(0+di,2+dj);
    constexpr int i11 = mul2nod(1+di,1+dj);
    
    if(idif1 == DifVar::Bary){
        // The Jacobian is quicker to compute recursively as it only uses the previous degree
        // evals. 
      double eva1[szfld], eva2[szfld];
      double jmat1[szfld],jmat2[szfld];

      eval1_bezier<szfld,1,di+1,dj+0>(rfld,lfld,idif2,DifVar::None,bary,eva1,jmat1,NULL);
      eval1_bezier<szfld,1,di+0,dj+1>(rfld,lfld,idif2,DifVar::None,bary,eva2,jmat2,NULL);

      for(int i = 0; i < szfld; i++){
        eval[i] = bary[0]*eva1[i] + bary[1]*eva2[i];
        jmat[i] = 2 * (eva2[i] - eva1[i]);
      }

      if(idif2 == DifVar::Bary){
        for(int i=0; i<szfld; i++){
                    // d11 (i.e. (d2-d1)(d2-d1))
          hmat[0*szfld + i] = 2 * (jmat2[0*szfld+i] - jmat1[0*szfld+i]);
        }
      }

    }else{

      for(int i = 0; i < szfld; i++){
        eval[i] =   bary[0]*bary[0]*rfld[lfld[i20]][i] 
                + 2*bary[0]*bary[1]*rfld[lfld[i11]][i] 
                +   bary[1]*bary[1]*rfld[lfld[i02]][i];
      }

    }
  }
};






} // End namespace


#endif