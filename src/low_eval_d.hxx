//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __LOW_EVAL_SURREALS__
#define __LOW_EVAL_SURREALS__

#include "types.hxx"
#include "ho_constants.hxx"
#include "aux_utils.hxx"
#include "codegen_lagrange.hxx"
//#include "low_eval_d_SurrealS.hxx"

#include "../SANS/Surreal/SurrealS.h"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
#include <boost/hana.hpp> 
namespace hana = boost::hana;
using namespace hana::literals;



namespace Metris{


constexpr int tdim3 = 3;
constexpr int sdim3 = 6;
/*
~~ low_eval_d.hxx: polynomial evaluation routines diff'd wrt poly coeffs ~~ 

~~ Main functions:
 - eval3_d: idem as eval3 with added ivar template param. Diff eval3 wrt poly coeff of 
 rank ivar. 
 - eval3_d_SurrealS0: Not to be called directly. Expects a tuple containing double* and Surreal<double,...>* 
 types determining the variable. eval3_d wraps around this for ease of use. 

~~ Inputs:
 - inherited from eval3: rfld,lfld,bary
 - ivar (template param): rank of the DoF in the element, comprised in [[0,tetnpps[ideg]-1]].
 - dfld (optional argument): derivatives of rfld[ivar][:] (default Id mat i.e. rfld[ivar] itself)
    dfld(i,j) = d_i rfld(idof,j) 

~~ Outputs
Derivatives w.r.t. coeff are D_i, w.r.t. bary are d_i
 - deval: szfld*szfld array. deval(i,j) = D_i eval[j]
 - djmat: szfld*szfld*tdim3 array. djmat(i,j)[k] = D_i d_j eval[k]
 - dhmat: szfld*[(szfld+1)*(szfld+2)/2]*tdim3 array. dhmat(i,j)[k] = D_i s_j eval[k]
 where s_j is a D^2 (bary) op ordered 1,1 1,2 etc. see low_eval.hxx. 

See low_eval.hxx header for other parameters and outputs. 
*/


/*
eval3_d_SurrealS0 is not designed to be called directly: call instead eval3_d. 

T designed to be a heterogeneous Boost::hana container such as a hana::tuple
Each type in the tuple can be a szfld allocated array. 
This tuple may hold SANS::SurrealS<szfld>* and no other n. 
Ideally, it should also hold exactly 1 SANS::SurrealS<szfld>*, doing otherwise
would be wasteful. 

This container can only be accessed by compile-time constants, which is why
every loop is replaced by a hana::while_ loop. 

Global indexing via lfld is no longer supported. This is because rfld is accessed
(in eval3) at lfld[i], but this cannot be made compile-time known (it depends on the mesh...)
As such, rfld is assumed ordered as given by ordedg, ordfac, ordtet. 
*/



/*
This is the wrapper that should be called. 
It sets up the tuple, gathers the derivatives.
- deval: derivatives of eval w.r.t. DoFs
  deval[N*i + j] with 0 <= i < szfld, N = szfld*szfld, 0 <= j < N
                 is the derivative w.r.t. to the j-th global (sized szfld) DoF
                 of the i-th evaluated coordinat
- djmat: similarly, deval[N*szfld*i + N*j + k] = d_k jmat(i,j)
- dhmat: idem, hmat (4D)

The plan is to constitute a pack of types using a hana::tuple
For instance hana::tuple<double, double, SANS::SurrealS, double>
Now, we'd want to populate it with values from rfld, or from our allocated SANS::SurrealS types. 
But we can't, because a hana::tuple is a purely compile-time construct. 
So we'll have to unpack its type template parameter into an std::tuple which can be 
modified at runtime! 
This is what we'll pass to eval3_d_SurrealS0. 

- hana offers hana::replicate to build a "constexpr dynamically sized" tuple
- we can "modify" types of this hana::tuple with the replace_at_c helper 
-> all to (double*) then ivar to (SANS::SurrealS<N,double>*)
- Problem: hana::tuple is fully constexpr so the double* etc values cannot be set
-> convert to std::tuple and pass that to the low-level routine
- Problem: std::tuple does not offer [] so we'll have to use std::get... which is not 
 overloadable to a basic dblAr2... 
-> wrap std::tuple in basic struct overloading []
- semi-Problem: [] argument must be compile-time known type; we can't deduce
 non-type template parameters... so we'll pass a hana::integral_constant
 and parameter deduction will catch its "value" (non-type parameter)
-> overload [] in MeshArrays to handle integral constants (simple cast to int)
*/

template <int szfld, int ideg, int ivar, int nvar = szfld>
void eval3_d(const dblAr2 & __restrict__ rfld,  
             const int * __restrict__ lfld,
             FEBasis ibasis, DifVar idif1, DifVar idif2, 
             const double * __restrict__ bary,
             double * __restrict__  eval,
             double * __restrict__  jmat,
             double * __restrict__  hmat,
             double * __restrict__  deval,
             double * __restrict__  djmat,
             double * __restrict__  dhmat,  
             const double * __restrict__ dfld = NULL);

template <int szfld, int ideg, int ivar, int nvar = szfld>
void eval2_d(const dblAr2 & __restrict__ rfld,  
             const int * __restrict__ lfld,
             FEBasis ibasis, DifVar idif1, DifVar idif2, 
             const double * __restrict__ bary,
             double * __restrict__  eval,
             double * __restrict__  jmat,
             double * __restrict__  hmat,
             double * __restrict__  deval,
             double * __restrict__  djmat,
             double * __restrict__  dhmat,  
             const double * __restrict__ dfld = NULL);






/* Secondary routines
  1. 3D */
//template <int szfld, int ideg, int ivar, int nvar = szfld>
//void eval3_d_direct(const dblAr2 & __restrict__ rfld,  
//                    const int * __restrict__ lfld,
//                    FEBasis ibasis, DifVar idif1, DifVar idif2, 
//                    const double * __restrict__ bary,  
//                    double * __restrict__  eval,
//                    double * __restrict__  jmat,
//                    double * __restrict__  hmat,
//                    double * __restrict__  deval,
//                    double * __restrict__  djmat,
//                    double * __restrict__  dhmat,
//                    const double * __restrict__  dfld = NULL);
//template <int szfld, int ideg, int ivar, int nvar = szfld>
//void eval3_d_SurrealS(const dblAr2 & __restrict__ rfld,  
//                      const int * __restrict__ lfld,
//                      FEBasis ibasis, DifVar idif1, DifVar idif2, 
//                      const double * __restrict__ bary,  
//                      double * __restrict__  eval,
//                      double * __restrict__  jmat,
//                      double * __restrict__  hmat,
//                      double * __restrict__  deval,
//                      double * __restrict__  djmat,
//                      double * __restrict__  dhmat,
//                      const double * __restrict__  dfld = NULL);
//// Attention: pointers to arrays; this is to handle NULLs (idif1/2 == 0).
//template <typename T, int szfld, int ideg, int nvar = szfld>
//void eval3_d_SurrealS0(const       T& __restrict__  rfld,  
//                       FEBasis ibasis, DifVar idif1, DifVar idif2, 
//                       const double * __restrict__  bary, 
//                       SANS::DLA::VectorS<  szfld,SANS::SurrealS<nvar,double>> *eval, 
//                       SANS::DLA::MatrixS<3,szfld,SANS::SurrealS<nvar,double>> *jmat, 
//                       SANS::DLA::MatrixS<6,szfld,SANS::SurrealS<nvar,double>> *hmat);


template <int szfld, int tdim, int ideg, int ivar, int nvar = szfld>
void eval_d_direct(const dblAr2 & __restrict__ rfld,  
                   const int * __restrict__ lfld,
                   FEBasis ibasis, DifVar idif1, DifVar idif2, 
                   const double * __restrict__ bary,  
                   double * __restrict__  eval,
                   double * __restrict__  jmat,
                   double * __restrict__  hmat,
                   double * __restrict__  deval,
                   double * __restrict__  djmat,
                   double * __restrict__  dhmat,
                   const double * __restrict__  dfld = NULL);
//// Attention: pointers to arrays; this is to handle NULLs (idif1/2 == 0).
//template <typename T, int szfld, int tdim, int ideg, int nvar = szfld>
//void eval_d_SurrealS0(const       T& __restrict__  rfld,  
//                      FEBasis ibasis, DifVar idif1, DifVar idif2, 
//                      const double * __restrict__  bary, 
//                      SANS::DLA::VectorS<                  szfld,SANS::SurrealS<nvar,double>> *eval, 
//                      SANS::DLA::MatrixS< tdim,            szfld,SANS::SurrealS<nvar,double>> *jmat, 
//                      SANS::DLA::MatrixS<(tdim*(tdim+1))/2,szfld,SANS::SurrealS<nvar,double>> *hmat);







/* ------------

Bezier eval 3D

------------ */

/*
With Hessian computed from Jmats
The Hessian is compressed as column first as usual:
1 2 4
  3 5
    6
However, note that each entry is, itself, comprised of szfld values. 

Eval, jmat and hmat are pointers in order to pass NULL when necessary. 
*/
template <typename T, int szfld, int ideg, int nvar, int di=0,int dj=0, int dk=0, int dl=0>
struct eval3_d_bezier{
  eval3_d_bezier(const T& __restrict__  rfld,  
                 DifVar idif1, DifVar idif2, 
                 const double  * __restrict__  bary,  
                 SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> *eval, 
                 SANS::DLA::MatrixS<tdim3,szfld,SANS::SurrealS<nvar,double>> *jmat, 
                 SANS::DLA::MatrixS<sdim3,szfld,SANS::SurrealS<nvar,double>> *hmat){


    SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> eva1, eva2, eva3, eva4;
    SANS::DLA::MatrixS<tdim3,szfld,SANS::SurrealS<nvar,double>> jmat1,jmat2,jmat3,jmat4;

    if(idif2 == DifVar::None){
      eval3_d_bezier<T,szfld,ideg-1,nvar,di+1,dj+0,dk+0,dl+0>(rfld,DifVar::None,DifVar::None,bary,&eva1,NULL,NULL);
      eval3_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+1,dk+0,dl+0>(rfld,DifVar::None,DifVar::None,bary,&eva2,NULL,NULL);
      eval3_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+0,dk+1,dl+0>(rfld,DifVar::None,DifVar::None,bary,&eva3,NULL,NULL);
      eval3_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+0,dk+0,dl+1>(rfld,DifVar::None,DifVar::None,bary,&eva4,NULL,NULL);
    }else if(idif2 == DifVar::Bary){
      eval3_d_bezier<T,szfld,ideg-1,nvar,di+1,dj+0,dk+0,dl+0>(rfld,DifVar::Bary,DifVar::None,bary,&eva1,&jmat1,NULL);
      eval3_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+1,dk+0,dl+0>(rfld,DifVar::Bary,DifVar::None,bary,&eva2,&jmat2,NULL);
      eval3_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+0,dk+1,dl+0>(rfld,DifVar::Bary,DifVar::None,bary,&eva3,&jmat3,NULL);
      eval3_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+0,dk+0,dl+1>(rfld,DifVar::Bary,DifVar::None,bary,&eva4,&jmat4,NULL);
    }

    for(int i = 0; i < szfld; i++){
      (*eval)[i] = bary[0] * eva1[i]
                 + bary[1] * eva2[i]
                 + bary[2] * eva3[i]
                 + bary[3] * eva4[i];
    }

    if(idif1 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        auto d1 = eva1[i];
        (*jmat)(0,i) = ideg * (eva2[i] - d1);
        (*jmat)(1,i) = ideg * (eva3[i] - d1);
        (*jmat)(2,i) = ideg * (eva4[i] - d1);
      }
    }

    if(idif2 == DifVar::Bary){
      for(int i=0; i<szfld; i++){
        // d11 (i.e. (d2-d1)(d2-d1))
        (*hmat)(0,i) = ideg*(jmat2(0,i) - jmat1(0,i));
        // d21 
        (*hmat)(1,i) = ideg*(jmat3(0,i) - jmat1(0,i));
        // d31 
        (*hmat)(3,i) = ideg*(jmat4(0,i) - jmat1(0,i));
  
        // d22 
        (*hmat)(2,i) = ideg*(jmat3(1,i) - jmat1(1,i));
        // d32 
        (*hmat)(4,i) = ideg*(jmat4(1,i) - jmat1(1,i));
  
        // d33 
        (*hmat)(5,i) = ideg*(jmat4(2,i) - jmat1(2,i));
      }
    }

  }
};

template <typename T, int szfld, int nvar, int di,int dj, int dk, int dl>
struct eval3_d_bezier<T,szfld,1,nvar,di,dj,dk,dl>{
  eval3_d_bezier(const T& __restrict__  rfld,  
                 DifVar idif1, DifVar idif2, 
                 const double  * __restrict__  bary,  
                 SANS::DLA::VectorS<  szfld,SANS::SurrealS<nvar,double>> *eval, 
                 SANS::DLA::MatrixS<3,szfld,SANS::SurrealS<nvar,double>> *jmat, 
                 SANS::DLA::MatrixS<6,szfld,SANS::SurrealS<nvar,double>> *hmat){

    //constexpr int ideg = 1 + di + dj + dk + dl; // The true degree. 

    constexpr int i1000 = mul2nod(1+di,0+dj,0+dk,0+dl);
    constexpr int i0100 = mul2nod(0+di,1+dj,0+dk,0+dl);
    constexpr int i0010 = mul2nod(0+di,0+dj,1+dk,0+dl);
    constexpr int i0001 = mul2nod(0+di,0+dj,0+dk,1+dl);

    for(int i = 0; i < szfld; i++){
      (*eval)[i] = bary[0]*rfld[hana::int_c<i1000>][i]
                 + bary[1]*rfld[hana::int_c<i0100>][i]
                 + bary[2]*rfld[hana::int_c<i0010>][i]
                 + bary[3]*rfld[hana::int_c<i0001>][i];
    }

    if(idif1 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        (*jmat)(0,i) = rfld[hana::int_c< i0100 >][i] 
                     - rfld[hana::int_c< i1000 >][i];

        (*jmat)(1,i) = rfld[hana::int_c< i0010 >][i] 
                     - rfld[hana::int_c< i1000 >][i];

        (*jmat)(2,i) = rfld[hana::int_c< i0001 >][i] 
                     - rfld[hana::int_c< i1000 >][i];
      }
    }

    if(idif2 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        (*hmat)(0,i) = 0;
        (*hmat)(1,i) = 0;
        (*hmat)(2,i) = 0;
        (*hmat)(3,i) = 0;
        (*hmat)(4,i) = 0;
        (*hmat)(5,i) = 0;
      }
    }

  }
};


template <typename T, int szfld,int nvar, int di,int dj, int dk, int dl>
  struct eval3_d_bezier<T,szfld,2,nvar,di,dj,dk,dl>{
    eval3_d_bezier(const T& __restrict__  rfld,  
                   DifVar idif1, DifVar idif2, 
                   const double  * __restrict__  bary,  
                   SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> *eval, 
                   SANS::DLA::MatrixS<tdim3,szfld,SANS::SurrealS<nvar,double>> *jmat, 
                   SANS::DLA::MatrixS<sdim3,szfld,SANS::SurrealS<nvar,double>> *hmat){
    //constexpr int ideg = 2 + di + dj + dk + dl; // The true degree. 

    //printf("Debug in this spec print rfld bary = %f %f %f %f\n",bary[0]
    //,bary[1],bary[2],bary[3]);
    //std::cout<<"0: "<<rfld[0_c][0]<<" "<<rfld[0_c][1]<<" "
    //                <<rfld[0_c][2]<<std::endl;
    //std::cout<<"1: "<<rfld[1_c][0]<<" "<<rfld[1_c][1]<<" "
    //                <<rfld[1_c][2]<<std::endl;
    //std::cout<<"2: "<<rfld[2_c][0]<<" "<<rfld[2_c][1]<<" "
    //                <<rfld[2_c][2]<<std::endl;
    //std::cout<<"3: "<<rfld[3_c][0]<<" "<<rfld[3_c][1]<<" "
    //                <<rfld[3_c][2]<<std::endl;


    if (idif1 == DifVar::Bary){
      // The Jacobian is quicker to compute recursively as it only uses the previous degree
      // evals. 

      SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> eva1, eva2, eva3, eva4;
      SANS::DLA::MatrixS<tdim3,szfld,SANS::SurrealS<nvar,double>> jmat1,jmat2,jmat3,jmat4;

      if(idif2 == DifVar::None){
        eval3_d_bezier<T,szfld,1,nvar,di+1,dj+0,dk+0,dl+0>(rfld,DifVar::None,DifVar::None,bary,&eva1,NULL,NULL);
        eval3_d_bezier<T,szfld,1,nvar,di+0,dj+1,dk+0,dl+0>(rfld,DifVar::None,DifVar::None,bary,&eva2,NULL,NULL);
        eval3_d_bezier<T,szfld,1,nvar,di+0,dj+0,dk+1,dl+0>(rfld,DifVar::None,DifVar::None,bary,&eva3,NULL,NULL);
        eval3_d_bezier<T,szfld,1,nvar,di+0,dj+0,dk+0,dl+1>(rfld,DifVar::None,DifVar::None,bary,&eva4,NULL,NULL);
      }else if(idif2 == DifVar::Bary){
        eval3_d_bezier<T,szfld,1,nvar,di+1,dj+0,dk+0,dl+0>(rfld,DifVar::Bary,DifVar::None,bary,&eva1,&jmat1,NULL);
        eval3_d_bezier<T,szfld,1,nvar,di+0,dj+1,dk+0,dl+0>(rfld,DifVar::Bary,DifVar::None,bary,&eva2,&jmat2,NULL);
        eval3_d_bezier<T,szfld,1,nvar,di+0,dj+0,dk+1,dl+0>(rfld,DifVar::Bary,DifVar::None,bary,&eva3,&jmat3,NULL);
        eval3_d_bezier<T,szfld,1,nvar,di+0,dj+0,dk+0,dl+1>(rfld,DifVar::Bary,DifVar::None,bary,&eva4,&jmat4,NULL);
      }


      for(int i = 0; i < szfld; i++){
        (*eval)[i] = bary[0]*eva1[i]
                   + bary[1]*eva2[i]
                   + bary[2]*eva3[i]
                   + bary[3]*eva4[i];
      }

      for(int i = 0; i < szfld;i++){
        auto d1 = eva1[i];
        (*jmat)(0,i) = 2*(eva2[i] - d1);
        (*jmat)(1,i) = 2*(eva3[i] - d1);
        (*jmat)(2,i) = 2*(eva4[i] - d1);
      }


      if (idif2 == DifVar::Bary){
        for(int i=0; i<szfld; i++){
          // d11 (i.e. (d2-d1)(d2-d1))
          (*hmat)(0,i) = 2*(jmat2(0,i) - jmat1(0,i));
          // d21 
          (*hmat)(1,i) = 2*(jmat3(0,i) - jmat1(0,i));
          // d31 
          (*hmat)(3,i) = 2*(jmat4(0,i) - jmat1(0,i));
    
          // d22 
          (*hmat)(2,i) = 2*(jmat3(1,i) - jmat1(1,i));
          // d32 
          (*hmat)(4,i) = 2*(jmat4(1,i) - jmat1(1,i));
    
          // d33 
          (*hmat)(5,i) = 2*(jmat4(2,i) - jmat1(2,i));
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
        (*eval)[i] =   bary[0]*bary[0]*rfld[hana::int_c<i2000>][i] 
                   + 2*bary[0]*bary[1]*rfld[hana::int_c<i1100>][i] 
                   + 2*bary[0]*bary[2]*rfld[hana::int_c<i1010>][i] 
                   + 2*bary[0]*bary[3]*rfld[hana::int_c<i1001>][i]
                   +   bary[1]*bary[1]*rfld[hana::int_c<i0200>][i] 
                   + 2*bary[1]*bary[2]*rfld[hana::int_c<i0110>][i] 
                   + 2*bary[1]*bary[3]*rfld[hana::int_c<i0101>][i]
                   +   bary[2]*bary[2]*rfld[hana::int_c<i0020>][i] 
                   + 2*bary[2]*bary[3]*rfld[hana::int_c<i0011>][i]
                   +   bary[3]*bary[3]*rfld[hana::int_c<i0002>][i] ;
      }
    }
  }
};




/* ------------

Bezier eval 2D

------------ */

template <typename T, int szfld, int ideg, int nvar, int di=0,int dj=0, int dk=0>
struct eval2_d_bezier{
  eval2_d_bezier(const T& __restrict__  rfld,  
                 DifVar idif1, DifVar idif2, 
                 const double  * __restrict__  bary,  
                 SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> *eval, 
                 SANS::DLA::MatrixS<2,szfld,SANS::SurrealS<nvar,double>> *jmat, 
                 SANS::DLA::MatrixS<3,szfld,SANS::SurrealS<nvar,double>> *hmat){


    SANS::DLA::VectorS<  szfld,SANS::SurrealS<nvar,double>> eva1, eva2, eva3;
    SANS::DLA::MatrixS<2,szfld,SANS::SurrealS<nvar,double>> jmat1,jmat2,jmat3;

    if(idif2 == DifVar::None){
      eval2_d_bezier<T,szfld,ideg-1,nvar,di+1,dj+0,dk+0>(rfld,DifVar::None,DifVar::None,bary,&eva1,NULL,NULL);
      eval2_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+1,dk+0>(rfld,DifVar::None,DifVar::None,bary,&eva2,NULL,NULL);
      eval2_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+0,dk+1>(rfld,DifVar::None,DifVar::None,bary,&eva3,NULL,NULL);
    }else if(idif2 == DifVar::Bary){
      eval2_d_bezier<T,szfld,ideg-1,nvar,di+1,dj+0,dk+0>(rfld,DifVar::Bary,DifVar::None,bary,&eva1,&jmat1,NULL);
      eval2_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+1,dk+0>(rfld,DifVar::Bary,DifVar::None,bary,&eva2,&jmat2,NULL);
      eval2_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+0,dk+1>(rfld,DifVar::Bary,DifVar::None,bary,&eva3,&jmat3,NULL);
    }

    for(int i = 0; i < szfld; i++){
      (*eval)[i] = bary[0] * eva1[i]
                 + bary[1] * eva2[i]
                 + bary[2] * eva3[i];
    }

    if(idif1 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        auto d1 = eva1[i];
        (*jmat)(0,i) = ideg * (eva2[i] - d1);
        (*jmat)(1,i) = ideg * (eva3[i] - d1);
      }
    }

    if(idif2 == DifVar::Bary){
      for(int i=0; i<szfld; i++){
        // d11 (i.e. (d2-d1)(d2-d1))
        (*hmat)(0,i) = ideg*(jmat2(0,i) - jmat1(0,i));
        // d21 
        (*hmat)(1,i) = ideg*(jmat3(0,i) - jmat1(0,i));
  
        // d22 
        (*hmat)(2,i) = ideg*(jmat3(1,i) - jmat1(1,i));
  
      }
    }

  }
};

template <typename T, int szfld, int nvar, int di,int dj, int dk>
struct eval2_d_bezier<T,szfld,1,nvar,di,dj,dk>{
  eval2_d_bezier(const T& __restrict__  rfld,  
                 DifVar idif1, DifVar idif2, 
                 const double  * __restrict__  bary,  
                 SANS::DLA::VectorS<  szfld,SANS::SurrealS<nvar,double>> *eval, 
                 SANS::DLA::MatrixS<2,szfld,SANS::SurrealS<nvar,double>> *jmat, 
                 SANS::DLA::MatrixS<3,szfld,SANS::SurrealS<nvar,double>> *hmat){

    //constexpr int ideg = 1 + di + dj + dk + dl; // The true degree. 

    constexpr int i100 = mul2nod(1+di,0+dj,0+dk);
    constexpr int i010 = mul2nod(0+di,1+dj,0+dk);
    constexpr int i001 = mul2nod(0+di,0+dj,1+dk);

    for(int i = 0; i < szfld; i++){
      (*eval)[i] = bary[0]*rfld[hana::int_c<i100>][i]
                 + bary[1]*rfld[hana::int_c<i010>][i]
                 + bary[2]*rfld[hana::int_c<i001>][i];
    }

    if(idif1 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        (*jmat)(0,i) = rfld[hana::int_c< i010 >][i] 
                     - rfld[hana::int_c< i100 >][i];

        (*jmat)(1,i) = rfld[hana::int_c< i001 >][i] 
                     - rfld[hana::int_c< i100 >][i];
      }
    }

    if(idif2 == DifVar::Bary){
      for(int i = 0; i < szfld; i++){
        (*hmat)(0,i) = 0;
        (*hmat)(1,i) = 0;
        (*hmat)(2,i) = 0;
      }
    }

  }
};

template <typename T, int szfld,int nvar, int di,int dj, int dk>
  struct eval2_d_bezier<T,szfld,2,nvar,di,dj,dk>{
    eval2_d_bezier(const T& __restrict__  rfld,  
                   DifVar idif1, DifVar idif2, 
                   const double  * __restrict__  bary,  
                   SANS::DLA::VectorS<  szfld,SANS::SurrealS<nvar,double>> *eval, 
                   SANS::DLA::MatrixS<2,szfld,SANS::SurrealS<nvar,double>> *jmat, 
                   SANS::DLA::MatrixS<3,szfld,SANS::SurrealS<nvar,double>> *hmat){
    //constexpr int ideg = 2 + di + dj + dk; // The true degree. 

    //printf("Debug in this spec print rfld bary = %f %f %f %f\n",bary[0]
    //,bary[1],bary[2],bary[3]);
    //std::cout<<"0: "<<rfld[0_c][0]<<" "<<rfld[0_c][1]<<" "
    //                <<rfld[0_c][2]<<std::endl;
    //std::cout<<"1: "<<rfld[1_c][0]<<" "<<rfld[1_c][1]<<" "
    //                <<rfld[1_c][2]<<std::endl;
    //std::cout<<"2: "<<rfld[2_c][0]<<" "<<rfld[2_c][1]<<" "
    //                <<rfld[2_c][2]<<std::endl;
    //std::cout<<"3: "<<rfld[3_c][0]<<" "<<rfld[3_c][1]<<" "
    //                <<rfld[3_c][2]<<std::endl;


    if (idif1 == DifVar::Bary){
      // The Jacobian is quicker to compute recursively as it only uses the previous degree
      // evals. 

      SANS::DLA::VectorS<  szfld,SANS::SurrealS<nvar,double>> eva1, eva2, eva3;
      SANS::DLA::MatrixS<2,szfld,SANS::SurrealS<nvar,double>> jmat1,jmat2,jmat3;

      if(idif2 == DifVar::None){
        eval2_d_bezier<T,szfld,1,nvar,di+1,dj+0,dk+0>(rfld,DifVar::None,DifVar::None,bary,&eva1,NULL,NULL);
        eval2_d_bezier<T,szfld,1,nvar,di+0,dj+1,dk+0>(rfld,DifVar::None,DifVar::None,bary,&eva2,NULL,NULL);
        eval2_d_bezier<T,szfld,1,nvar,di+0,dj+0,dk+1>(rfld,DifVar::None,DifVar::None,bary,&eva3,NULL,NULL);
      }else if(idif2 == DifVar::Bary){
        eval2_d_bezier<T,szfld,1,nvar,di+1,dj+0,dk+0>(rfld,DifVar::Bary,DifVar::None,bary,&eva1,&jmat1,NULL);
        eval2_d_bezier<T,szfld,1,nvar,di+0,dj+1,dk+0>(rfld,DifVar::Bary,DifVar::None,bary,&eva2,&jmat2,NULL);
        eval2_d_bezier<T,szfld,1,nvar,di+0,dj+0,dk+1>(rfld,DifVar::Bary,DifVar::None,bary,&eva3,&jmat3,NULL);
      }


      for(int i = 0; i < szfld; i++){
        (*eval)[i] = bary[0]*eva1[i]
                   + bary[1]*eva2[i]
                   + bary[2]*eva3[i];
      }

      for(int i = 0; i < szfld;i++){
        auto d1 = eva1[i];
        (*jmat)(0,i) = 2*(eva2[i] - d1);
        (*jmat)(1,i) = 2*(eva3[i] - d1);
      }


      if (idif2 == DifVar::Bary){
        for(int i=0; i<szfld; i++){
          // d11 (i.e. (d2-d1)(d2-d1))
          (*hmat)(0,i) = 2*(jmat2(0,i) - jmat1(0,i));
          // d21 
          (*hmat)(1,i) = 2*(jmat3(0,i) - jmat1(0,i));
    
          // d22 
          (*hmat)(2,i) = 2*(jmat3(1,i) - jmat1(1,i));
    
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
        (*eval)[i] =   bary[0]*bary[0]*rfld[hana::int_c<i200>][i] 
                   + 2*bary[0]*bary[1]*rfld[hana::int_c<i110>][i] 
                   + 2*bary[0]*bary[2]*rfld[hana::int_c<i101>][i] 
                   +   bary[1]*bary[1]*rfld[hana::int_c<i020>][i] 
                   + 2*bary[1]*bary[2]*rfld[hana::int_c<i011>][i] 
                   +   bary[2]*bary[2]*rfld[hana::int_c<i002>][i]  ;
      }
    }
  }
};


template <typename T, int szfld, int tdim, int ideg,  int nvar>
void eval_d_SurrealS0(const T& __restrict__  rfld,   
                      FEBasis ibasis, DifVar idif1, DifVar idif2, 
                      const double * __restrict__  bary, 
                      SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> *eval, 
                      SANS::DLA::MatrixS<tdim,szfld,SANS::SurrealS<nvar,double>> *jmat, 
                      SANS::DLA::MatrixS<(tdim*(tdim+1))/2,szfld,SANS::SurrealS<nvar,double>> *hmat){
  //printf("Debug surreal0 bary = %f %f %f %f field:\n",bary[0]
  //  ,bary[1],bary[2],bary[3]);
  //std::cout<<"0: "<<rfld[0_c][0]<<" "<<rfld[0_c][1]<<" "
  //                <<rfld[0_c][2]<<std::endl;
  //std::cout<<"1: "<<rfld[1_c][0]<<" "<<rfld[1_c][1]<<" "
  //                <<rfld[1_c][2]<<std::endl;
  //std::cout<<"2: "<<rfld[2_c][0]<<" "<<rfld[2_c][1]<<" "
  //                <<rfld[2_c][2]<<std::endl;
  //std::cout<<"3: "<<rfld[3_c][0]<<" "<<rfld[3_c][1]<<" "
  //                <<rfld[3_c][2]<<std::endl;
  if constexpr (ideg == 1){
    //hana::while_(hana::less.than(szfld_c), 0_c, [&](auto i_c){
    //  constexpr int i = i_c;
    
    for(int icmp = 0; icmp < szfld; icmp++){
      (*eval)[icmp] = bary[0]*rfld[0_c][icmp];
      // Note bounds included
      CT_FOR0_INC(1,tdim,ii){  
        // Macro creates c_ii
        (*eval)[icmp] += bary[ii]*rfld[c_ii][icmp];
      }CT_FOR1(ii);
      //(*eval)[icmp] = bary[0]*rfld[0_c][icmp]
      //              + bary[1]*rfld[1_c][icmp]
      //              + bary[2]*rfld[2_c][icmp]
      //              + bary[3]*rfld[3_c][icmp];
    }
    //return i_c+1_c;});

    if (idif1 == DifVar::Bary){
      
      for(int i = 0; i < szfld; i++){
        // Note upper bound excluded
        CT_FOR0_INC(1,tdim,jj){
          (*jmat)(jj-1,i) = rfld[c_jj][i] - rfld[0_c][i];

        }CT_FOR1(jj);
        //(*jmat)(0,i) = rfld[1_c][i] - rfld[0_c][i];
        //(*jmat)(1,i) = rfld[2_c][i] - rfld[0_c][i];
        //(*jmat)(2,i) = rfld[3_c][i] - rfld[0_c][i];
      }
    }
    if (idif2 == DifVar::Bary){
      constexpr int nnsym = tdim*(tdim+1)/2;
      
      for(int i = 0; i < szfld; i++){
        for(int jj = 0; jj < nnsym; jj++) (*hmat)(jj,i) = 0;
        //(*hmat)(0,i) = 0;
        //(*hmat)(1,i) = 0;
        //(*hmat)(2,i) = 0;
        //(*hmat)(3,i) = 0;
        //(*hmat)(4,i) = 0;
        //(*hmat)(5,i) = 0;
      }
    }
    return;
  } 

  //std::cout<<"Debug standard:\n";
  //constexpr int npps = tdim == 2 ? facnpps[ideg] : tetnpps[ideg];
  //CT_FOR0_INC(0,npps-1,ii){
  //  std::cout<<"rfld["<<ii<<"] = "<<(*rfld[c_ii])<<"\n";
  //}CT_FOR1(ii);

  if(ibasis == FEBasis::Bezier){
    if constexpr(tdim == 2){
      eval2_d_bezier<T,szfld,ideg,nvar>(rfld, idif1, idif2, bary, eval, jmat, hmat);
    }else{
      eval3_d_bezier<T,szfld,ideg,nvar>(rfld, idif1, idif2, bary, eval, jmat, hmat);
    }
  }else{
    METRIS_THROW_MSG(TODOExcept(),"Implement eval3_d_LAGRANGE");
  } 
  return;
}





template <int szfld, int tdim, int ideg,  int nvar>
void eval_d_SurrealS0_simple(MeshArray2D<SANS::SurrealS<nvar,double>> &rfld,   
                             FEBasis ibasis, DifVar idif1, DifVar idif2, 
                             const double * __restrict__  bary, 
                             SANS::DLA::VectorS<     szfld,SANS::SurrealS<nvar,double>> *eval, 
                             SANS::DLA::MatrixS<tdim,szfld,SANS::SurrealS<nvar,double>> *jmat, 
                             SANS::DLA::MatrixS<(tdim*(tdim+1))/2,szfld,SANS::SurrealS<nvar,double>> *hmat){
  //printf("Debug surreal0 bary = %f %f %f %f field:\n",bary[0]
  //  ,bary[1],bary[2],bary[3]);
  //std::cout<<"0: "<<rfld[0_c][0]<<" "<<rfld[0_c][1]<<" "
  //                <<rfld[0_c][2]<<std::endl;
  //std::cout<<"1: "<<rfld[1_c][0]<<" "<<rfld[1_c][1]<<" "
  //                <<rfld[1_c][2]<<std::endl;
  //std::cout<<"2: "<<rfld[2_c][0]<<" "<<rfld[2_c][1]<<" "
  //                <<rfld[2_c][2]<<std::endl;
  //std::cout<<"3: "<<rfld[3_c][0]<<" "<<rfld[3_c][1]<<" "
  //                <<rfld[3_c][2]<<std::endl;
  if constexpr (ideg == 1){
    //hana::while_(hana::less.than(szfld_c), 0_c, [&](auto i_c){
    //  constexpr int i = i_c;
    
    for(int icmp = 0; icmp < szfld; icmp++){
      (*eval)[icmp] = bary[0]*rfld[0_c][icmp];
      // Note bounds included
      for(int ii = 1; ii <= tdim; ii++){
        (*eval)[icmp] += bary[ii]*rfld[ii][icmp];
      }
      //CT_FOR0_INC(1,tdim,ii){  
      //  // Macro creates c_ii
      //  (*eval)[icmp] += bary[ii]*rfld[c_ii][icmp];
      //}CT_FOR1(ii);
      //(*eval)[icmp] = bary[0]*rfld[0_c][icmp]
      //              + bary[1]*rfld[1_c][icmp]
      //              + bary[2]*rfld[2_c][icmp]
      //              + bary[3]*rfld[3_c][icmp];
    }
    //return i_c+1_c;});

    if (idif1 == DifVar::Bary){
      
      for(int i = 0; i < szfld; i++){
        // Note upper bound excluded
        for(int jj = 1;jj <= tdim; jj++){
          (*jmat)(jj-1,i) = rfld[jj][i] - rfld[0][i];
        }
        //CT_FOR0_INC(1,tdim,jj){
        //  (*jmat)(jj-1,i) = rfld[c_jj][i] - rfld[0_c][i];
        //}CT_FOR1(jj);
        //(*jmat)(0,i) = rfld[1_c][i] - rfld[0_c][i];
        //(*jmat)(1,i) = rfld[2_c][i] - rfld[0_c][i];
        //(*jmat)(2,i) = rfld[3_c][i] - rfld[0_c][i];
      }
    }
    if (idif2 == DifVar::Bary){
      constexpr int nnsym = tdim*(tdim+1)/2;
      
      for(int i = 0; i < szfld; i++){
        for(int jj = 0; jj < nnsym; jj++) (*hmat)(jj,i) = 0;
        //(*hmat)(0,i) = 0;
        //(*hmat)(1,i) = 0;
        //(*hmat)(2,i) = 0;
        //(*hmat)(3,i) = 0;
        //(*hmat)(4,i) = 0;
        //(*hmat)(5,i) = 0;
      }
    }
    return;
  } 

  //std::cout<<"Debug simple:\n";
  //constexpr int npps = tdim == 2 ? facnpps[ideg] : tetnpps[ideg];
  //CT_FOR0_INC(0,npps-1,ii){
  //  //std::cout<<"rfld["<<ii<<"] = "<<rfld[c_ii]<<"\n";
  //  std::cout<<"rfld["<<ii<<"] = "<<(*rfld[c_ii])<<"\n";
  //}CT_FOR1(ii);

  if(ibasis == FEBasis::Bezier){
    if constexpr(tdim == 2){
      eval2_d_bezier<MeshArray2D<SANS::SurrealS<nvar,double>>,szfld,ideg,nvar>(rfld, idif1, idif2, bary, eval, jmat, hmat);
    }else{
      eval3_d_bezier<MeshArray2D<SANS::SurrealS<nvar,double>>,szfld,ideg,nvar>(rfld, idif1, idif2, bary, eval, jmat, hmat);
    }
  }else{
    METRIS_THROW_MSG(TODOExcept(),"Implement eval3_d_LAGRANGE");
  } 
  return;
}




// SANS::SurrealS-based version can only handle Bézier for now. 
// Slower in all cases except P2 + Jacobian. 
template <int szfld, int tdim, int ideg,  int ivar, int nvar = szfld>
void eval_d_SurrealS(const dblAr2 & __restrict__ rfld,  
                     const int * __restrict__ lfld,
                     FEBasis ibasis, DifVar idif1, DifVar idif2, 
                     const double * __restrict__ bary,
                     double * __restrict__  eval,
                     double * __restrict__  jmat,
                     double * __restrict__  hmat,
                     double * __restrict__  deval,
                     double * __restrict__  djmat,
                     double * __restrict__  dhmat,  
                     const double * __restrict__ dfld = NULL){

  constexpr auto entnpps = ENTNPPS(tdim);

  // -- For clarity. 
  // Number of (r|l)fld entries
  constexpr int nrfld = entnpps[ideg];
  static_assert(nrfld < 31);
  // -

  static_assert(ivar >= 0 && ivar < entnpps[ideg]);

  auto rfld_0 = hana::replicate<hana::tuple_tag>((const double *) 1, hana::size_c<nrfld>);

  auto rfld_1 = to_std_tuple(replace_at_c<ivar>(rfld_0,(SANS::SurrealS<nvar,double>*) 0));
  // Store non-dof rfld entries
//  double rmem[szfld*(nrfld-1)];
  // Store the dof as SANS::SurrealS
  SANS::SurrealS<nvar,double> smem[szfld];

  if(dfld == NULL){
    assert(nvar == szfld);
    for(int icmp = 0; icmp < szfld; icmp++){
      smem[icmp].value() = rfld[lfld[ivar]][icmp];
      for(int j = 0; j < nvar; j++){
        smem[icmp].deriv(j) = 0;
      }
      smem[icmp].deriv(icmp) = 1;
    }
  }else{
    for(int icmp = 0; icmp < szfld; icmp++){
      smem[icmp].value() = rfld[lfld[ivar]][icmp];
      for(int idof = 0; idof < nvar; idof++){
        smem[icmp].deriv(idof) = dfld[idof*szfld + icmp];
      }
    }
  }

  //int imem = 0;
  //for(int ipoly = 0; ipoly < entnpps[ideg]; ipoly++){
  //  if(ipoly == ivar) continue;
  //  for(int icmp = 0; icmp < szfld; icmp++){
  //    rmem[imem*szfld + icmp] = rfld[lfld[ipoly]][icmp];
  //  }
  //  imem++;
  //}

  tuple_wrapper w_op(rfld_1);

  // Populate w_op tuple with appropriate types. 
  // The idea is to get something like a tuple<double*, SANS::SurrealS*, double*, double*...>
  // Storing, for the double*'s, the non-dof rfld entries and, for the one SANS::SurrealS*, 
  // the DoF SANS::SurrealS array. 
  //imem = 0;
  hana::while_(hana::less.than(hana::int_c<nrfld>), 0_c, [&](auto i_c){
    constexpr int i = i_c;
    if constexpr(i == ivar){
      w_op[i_c] = smem;
    }else{
//      w_op[i_c] = &rmem[szfld*imem];
      w_op[i_c] = rfld[lfld[i]];
//      imem++;
    }
  return i_c+1_c;});  
//  SANS::SurrealS<nvar,double> seval[szfld];
  constexpr int nnsym = (tdim*(tdim+1))/2;

  SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> seval;
  SANS::DLA::MatrixS<tdim ,szfld,SANS::SurrealS<nvar,double>> sjmat;
//  SANS::SurrealS<nvar,double> sjmat[tdim3*szfld];
  SANS::DLA::MatrixS<nnsym,szfld,SANS::SurrealS<nvar,double>> shmat;
//  SANS::SurrealS<nvar,double> shmat[sdim3*szfld];


  //hana::while_(hana::less.than(hana::int_c<nrfld>), 0_c, [&](auto i_c){
  //  constexpr int i = i_c;
  //  printf("Debug rfld[%d] = \n",i);
  //  std::cout<<"0:"<<w_op[i_c][0]<<std::endl;
  //  std::cout<<"1:"<<w_op[i_c][1]<<std::endl;
  //  std::cout<<"2:"<<w_op[i_c][2]<<std::endl;
  //return i_c+1_c;});  


  eval_d_SurrealS0<decltype(w_op),szfld,tdim,ideg>
                     (w_op,ibasis,idif1,idif2,bary,&seval,&sjmat,&shmat);

  for(int icmp=0; icmp<szfld; icmp++){
    eval[icmp] = seval[icmp].value();

    for(int idof = 0; idof < nvar; idof++){
      deval[szfld*idof + icmp] = seval[icmp].deriv(idof);
    }

    if(idif1 == DifVar::Bary){
      for(int itdim = 0; itdim < tdim; itdim++){
        jmat[itdim*szfld + icmp] = sjmat(itdim,icmp).value();
        for(int idof = 0; idof < nvar; idof++){
          djmat[idof*tdim*szfld + itdim*szfld + icmp] = sjmat(itdim,icmp).deriv(idof);
//          djmat[itdim*szfld*3 + icmp*3 + j] = sjmat[itdim*szfld + icmp].deriv(nvar);
        } 
      }
    }

    if(idif2 == DifVar::Bary){
      for(int ii = 0; ii < nnsym; ii++) hmat[ii*szfld + icmp] = shmat(ii,icmp).value();
      //hmat[0*szfld + icmp] = shmat(0,icmp).value();
      //hmat[1*szfld + icmp] = shmat(1,icmp).value();
      //hmat[2*szfld + icmp] = shmat(2,icmp).value();
      //hmat[3*szfld + icmp] = shmat(3,icmp).value();
      //hmat[4*szfld + icmp] = shmat(4,icmp).value();
      //hmat[5*szfld + icmp] = shmat(5,icmp).value();
      for(int idof = 0; idof < tdim; idof++){
        for(int ii = 0; ii < nnsym; ii++) dhmat[idof*nnsym*szfld + 0*szfld + icmp] = shmat(ii,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 0*szfld + icmp] = shmat(0,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 1*szfld + icmp] = shmat(1,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 2*szfld + icmp] = shmat(2,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 3*szfld + icmp] = shmat(3,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 4*szfld + icmp] = shmat(4,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 5*szfld + icmp] = shmat(5,icmp).deriv(idof);
      }
    }

  }
}




// Same as eval_d_SurrealS 
// but instead of computing szfld, compute size 1, then notice whole vector is repeated
template <int szfld, int tdim, int ideg,  int ivar, int nvar>
void eval_d_SurrealS_bcast(const dblAr2 & __restrict__ rfld,  
                     const int * __restrict__ lfld,
                     FEBasis ibasis, DifVar idif1, DifVar idif2, 
                     const double * __restrict__ bary,
                     double * __restrict__  eval,
                     double * __restrict__  jmat,
                     double * __restrict__  hmat,
                     double * __restrict__  deval,
                     double * __restrict__  djmat,
                     double * __restrict__  dhmat,  
                     const double * __restrict__ dfld = NULL){

  constexpr auto entnpps = ENTNPPS(tdim);

  // -- For clarity. 
  // Number of (r|l)fld entries
  constexpr int nrfld = entnpps[ideg];
  static_assert(nrfld < 31);
  // -

  static_assert(ivar >= 0 && ivar < entnpps[ideg]);

  auto rfld_0 = hana::replicate<hana::tuple_tag>((const double *) 1, hana::size_c<nrfld>);

  auto rfld_1 = to_std_tuple(replace_at_c<ivar>(rfld_0,(SANS::SurrealS<nvar,double>*) 0));
  // Store non-dof rfld entries
//  double rmem[szfld*(nrfld-1)];
  // Store the dof as SANS::SurrealS
  SANS::SurrealS<nvar,double> smem[szfld];

  if(dfld == NULL){
    assert(nvar == 1);
    for(int icmp = 0; icmp < szfld; icmp++){
      smem[icmp].value() = rfld[lfld[ivar]][icmp];
      smem[icmp].deriv(0) = 0;
    }
    smem[0].deriv(0) = 1; // They will all be the same
  }else{
    for(int icmp = 0; icmp < szfld; icmp++){
      smem[icmp].value() = rfld[lfld[ivar]][icmp];
      for(int idof = 0; idof < nvar; idof++){
        smem[icmp].deriv(idof) = dfld[idof*szfld + icmp];
      }
    }
  }

  //int imem = 0;
  //for(int ipoly = 0; ipoly < entnpps[ideg]; ipoly++){
  //  if(ipoly == ivar) continue;
  //  for(int icmp = 0; icmp < szfld; icmp++){
  //    rmem[imem*szfld + icmp] = rfld[lfld[ipoly]][icmp];
  //  }
  //  imem++;
  //}

  tuple_wrapper w_op(rfld_1);

  // Populate w_op tuple with appropriate types. 
  // The idea is to get something like a tuple<double*, SANS::SurrealS*, double*, double*...>
  // Storing, for the double*'s, the non-dof rfld entries and, for the one SANS::SurrealS*, 
  // the DoF SANS::SurrealS array. 
  //imem = 0;
  hana::while_(hana::less.than(hana::int_c<nrfld>), 0_c, [&](auto i_c){
    constexpr int i = i_c;
    if constexpr(i == ivar){
      w_op[i_c] = smem;
    }else{
//      w_op[i_c] = &rmem[szfld*imem];
      w_op[i_c] = rfld[lfld[i]];
//      imem++;
    }
  return i_c+1_c;});  
//  SANS::SurrealS<nvar,double> seval[szfld];
  constexpr int nnsym = (tdim*(tdim+1))/2;

  SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> seval;
  SANS::DLA::MatrixS<tdim ,szfld,SANS::SurrealS<nvar,double>> sjmat;
//  SANS::SurrealS<nvar,double> sjmat[tdim3*szfld];
  SANS::DLA::MatrixS<nnsym,szfld,SANS::SurrealS<nvar,double>> shmat;
//  SANS::SurrealS<nvar,double> shmat[sdim3*szfld];


  //hana::while_(hana::less.than(hana::int_c<nrfld>), 0_c, [&](auto i_c){
  //  constexpr int i = i_c;
  //  printf("Debug rfld[%d] = \n",i);
  //  std::cout<<"0:"<<w_op[i_c][0]<<std::endl;
  //  std::cout<<"1:"<<w_op[i_c][1]<<std::endl;
  //  std::cout<<"2:"<<w_op[i_c][2]<<std::endl;
  //return i_c+1_c;});  


  eval_d_SurrealS0<decltype(w_op),szfld,tdim,ideg>
                     (w_op,ibasis,idif1,idif2,bary,&seval,&sjmat,&shmat);

  //std::cout<<"Debug seval = \n"<<seval<<"\n";

  //std::cout<<"Debug sjmat = \n"<<sjmat<<"\n";

  //std::cout<<"debug ideg = "<<ideg<< " bary = "; 
  //dblAr1(tdim+1,bary).print();



  int nvareff = dfld == NULL ? szfld : nvar;

  for(int icmp=0; icmp<szfld; icmp++){
    eval[icmp] = seval[icmp].value();

    if(dfld != NULL){
      for(int idof = 0; idof < nvareff; idof++){
        deval[szfld*idof + icmp] = seval[icmp].deriv(idof);
      }
      if(idif1 == DifVar::Bary){
        for(int itdim = 0; itdim < tdim; itdim++){
          jmat[itdim*szfld + icmp] = sjmat(itdim,icmp).value();
          for(int idof = 0; idof < nvareff; idof++){
            djmat[idof*tdim*szfld + itdim*szfld + icmp] = sjmat(itdim,icmp).deriv(idof);
          } 
        }
      }
    }else{
      double de = seval[0].deriv(0);
      for(int idof = 0; idof < nvareff; idof++){
        deval[szfld*idof + icmp] = de * (idof == icmp);
      }
      if(idif1 == DifVar::Bary){
        double dj = sjmat(0,0).deriv(0);
        for(int itdim = 0; itdim < tdim; itdim++){
          jmat[itdim*szfld + icmp] = sjmat(itdim,icmp).value();
          for(int idof = 0; idof < nvareff; idof++){
            djmat[idof*tdim*szfld + itdim*szfld + icmp] = dj * (idof == icmp);
          } 
        }
      }
    }



    if(idif2 == DifVar::Bary){
      for(int ii = 0; ii < nnsym; ii++) hmat[ii*szfld + icmp] = shmat(ii,icmp).value();
      //hmat[0*szfld + icmp] = shmat(0,icmp).value();
      //hmat[1*szfld + icmp] = shmat(1,icmp).value();
      //hmat[2*szfld + icmp] = shmat(2,icmp).value();
      //hmat[3*szfld + icmp] = shmat(3,icmp).value();
      //hmat[4*szfld + icmp] = shmat(4,icmp).value();
      //hmat[5*szfld + icmp] = shmat(5,icmp).value();
      for(int idof = 0; idof < tdim; idof++){
        for(int ii = 0; ii < nnsym; ii++) dhmat[idof*nnsym*szfld + 0*szfld + icmp] = shmat(ii,icmp).deriv(0);
        //dhmat[idof*sdim3*szfld + 0*szfld + icmp] = shmat(0,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 1*szfld + icmp] = shmat(1,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 2*szfld + icmp] = shmat(2,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 3*szfld + icmp] = shmat(3,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 4*szfld + icmp] = shmat(4,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 5*szfld + icmp] = shmat(5,icmp).deriv(idof);
      }
    }

  }
}






// This version just makes one large SurrealS buffer instead of trying to limit it to 1 variable
template <int szfld, int tdim, int ideg, int nvar>
void eval_d_SurrealS_simple(const dblAr2 & __restrict__ rfld,  
                            const int * __restrict__ lfld,
                            FEBasis ibasis, DifVar idif1, DifVar idif2, 
                            const double * __restrict__ bary,
                            int ivar, 
                            double * __restrict__  eval,
                            double * __restrict__  jmat,
                            double * __restrict__  hmat,
                            double * __restrict__  deval,
                            double * __restrict__  djmat,
                            double * __restrict__  dhmat,  
                            const double * __restrict__ dfld = NULL){

  constexpr auto entnpps = ENTNPPS(tdim);

  // -- For clarity. 
  // Number of (r|l)fld entries
  constexpr int nrfld = entnpps[ideg];
  // -


  SANS::SurrealS<nvar,double> buffer[nrfld][szfld];
  MeshArray2D<SANS::SurrealS<nvar,double>> rfllS(nrfld,szfld,buffer[0]);

  if(dfld == NULL){
    METRIS_ASSERT(nvar == szfld);
  }
  for(int ii = 0; ii < nrfld; ii++){
    for(int jj = 0; jj < szfld; jj++){
      rfllS[ii][jj].value() = rfld[lfld[ii]][jj];
      for(int kk = 0; kk < nvar; kk++){
        rfllS[ii][jj].deriv(kk) = 0;
      }
      if(ii == ivar) rfllS[ii][jj].deriv(jj) = 1;
    }
  }
  if(dfld != NULL){
    for(int icmp = 0; icmp < szfld; icmp++){
      rfllS[ivar][icmp].value() = rfld[lfld[ivar]][icmp];
      for(int idof = 0; idof < nvar; idof++){
        rfllS[ivar][icmp].deriv(idof) = dfld[idof*szfld + icmp];
      }
    }
  }




  constexpr int nnsym = (tdim*(tdim+1))/2;

  SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> seval;
  SANS::DLA::MatrixS<tdim ,szfld,SANS::SurrealS<nvar,double>> sjmat;
  SANS::DLA::MatrixS<nnsym,szfld,SANS::SurrealS<nvar,double>> shmat;

  eval_d_SurrealS0_simple<szfld,tdim,ideg>
                           (rfllS,ibasis,idif1,idif2,bary,&seval,&sjmat,&shmat);

  for(int icmp=0; icmp<szfld; icmp++){
    eval[icmp] = seval[icmp].value();

    for(int idof = 0; idof < nvar; idof++){
      deval[szfld*idof + icmp] = seval[icmp].deriv(idof);
    }

    if(idif1 == DifVar::Bary){
      for(int itdim = 0; itdim < tdim; itdim++){
        jmat[itdim*szfld + icmp] = sjmat(itdim,icmp).value();
        for(int idof = 0; idof < nvar; idof++){
          djmat[idof*tdim*szfld + itdim*szfld + icmp] = sjmat(itdim,icmp).deriv(idof);
        } 
      }
    }

    if(idif2 == DifVar::Bary){
      for(int ii = 0; ii < nnsym; ii++) hmat[ii*szfld + icmp] = shmat(ii,icmp).value();
      for(int idof = 0; idof < tdim; idof++){
        for(int ii = 0; ii < nnsym; ii++) dhmat[idof*nnsym*szfld + 0*szfld + icmp] = shmat(ii,icmp).deriv(idof);
      }
    }

  }
}











/*
Main eval3_d routine. Calls either eval3_d_SurrealS or eval3_d_direct 
depending on arguments. Direct version can handle both Bézier and Lagrange. 
SANS::SurrealS can only handle Bézier (TODO: implement Lagrange evals with SANS::SurrealS)
*/
/* --- dfld
dfld of size (szfld,szfld) specifies the derivatives of the input field. 
This holds the derivatives of the ivar-th coeff in rfld w.r.t. considered
variable. 

Example, when evaluating metric at xi using metrics at control points, 
metrics at control points depend on the variable in some manner 
(determined by back mesh). 

We restrict ourselves to the case where only the ivar-th rfld entry (e.g. metric)
is non-constant, though we could imagine cases where all polynomial coeffs
depend on a common variable. 
*/
template <int szfld, int ideg, int ivar, int nvar>
void eval3_d(const dblAr2 & __restrict__ rfld,  
             const int * __restrict__ lfld,
             FEBasis ibasis, DifVar idif1, DifVar idif2, 
             const double * __restrict__ bary,
             double * __restrict__  eval,
             double * __restrict__  jmat,
             double * __restrict__  hmat,
             double * __restrict__  deval,
             double * __restrict__  djmat,
             double * __restrict__  dhmat,  
             const double * __restrict__ dfld){

  assert(! (dfld == NULL && nvar != szfld) && "If nvar != szfld, Jacobian matrix must be specified");
  if(idif2 == DifVar::Bary && ibasis == FEBasis::Bezier){
    printf("## EITHER IMPLEMENT HESSIAN IN DIRECT OR IMPLEMENT LAGRANGE IN SANS::SurrealS\n");
    exit(1);
  }

  // Only case where SANS::SurrealS benchmarked faster than direct method
  // But also the only one with a Hessian implementation...
  if(idif2 == DifVar::Bary || (ideg == 2 && idif1 == DifVar::Bary && ibasis == FEBasis::Bezier)){
    eval_d_SurrealS<szfld, 3, ideg,  ivar, nvar>
                 (rfld,lfld,ibasis, idif1, idif2,bary,eval,jmat,hmat,deval,djmat,dhmat, dfld);
    //eval3_d_SurrealS<szfld, ideg,  ivar, nvar>
    //              (rfld,lfld,ibasis, idif1, idif2,bary,eval,jmat,hmat,deval,djmat,dhmat, dfld);
  }else{
    eval_d_direct<szfld, 3, ideg,  ivar, nvar>
                 (rfld,lfld,ibasis, idif1, idif2,bary,eval,jmat,hmat,deval,djmat,dhmat, dfld);
    //eval3_d_direct<szfld, ideg,  ivar, nvar>
    //              (rfld,lfld,ibasis, idif1, idif2,bary,eval,jmat,hmat,deval,djmat,dhmat, dfld);
  }
}



template <int szfld, int ideg, int ivar, int nvar>
void eval2_d(const dblAr2 & __restrict__ rfld,  
             const int * __restrict__ lfld,
             FEBasis ibasis, DifVar idif1, DifVar idif2, 
             const double * __restrict__ bary,
             double * __restrict__  eval,
             double * __restrict__  jmat,
             double * __restrict__  hmat,
             double * __restrict__  deval,
             double * __restrict__  djmat,
             double * __restrict__  dhmat,  
             const double * __restrict__ dfld){

  assert(! (dfld == NULL && nvar != szfld) && "If nvar != szfld, Jacobian matrix must be specified");
  if(idif2 == DifVar::Bary && ibasis == FEBasis::Bezier){
    printf("## EITHER IMPLEMENT HESSIAN IN DIRECT OR IMPLEMENT LAGRANGE IN SANS::SurrealS\n");
    exit(1);
  }

  // Only case where SANS::SurrealS benchmarked faster than direct method
  // But also the only one with a Hessian implementation...
  if(idif2 == DifVar::Bary || (ideg == 2 && idif1 == DifVar::Bary && ibasis == FEBasis::Bezier)){
    eval_d_SurrealS<szfld,2,  ideg,  ivar, nvar>
                 (rfld,lfld,ibasis, idif1, idif2,bary,eval,jmat,hmat,deval,djmat,dhmat, dfld);
  }else{
    eval_d_direct<szfld, 2, ideg,  ivar, nvar>
                 (rfld,lfld,ibasis, idif1, idif2,bary,eval,jmat,hmat,deval,djmat,dhmat, dfld);
  }
}


/* --------------------------------------------------------------------
  Secondary functions
   -------------------------------------------------------------------- */


/* Exploit F_K = sum f_\alpha P_\alpha with f either Bernstein or Lagrange
   and fill in derivatives directly. If we denote J_idof the Jacobian of 
   rfld[idof] w.r.t. variables, then the final Jacobian is simply 
   fun * J_idof where fun is either the Bernstein or Lagrange ivar-th 
   basis function. 
   Similar reasoning for jmat and hmat (not implemented). */
template <int szfld, int tdim, int ideg,  int ivar, int nvar>
void eval_d_direct(const dblAr2 & __restrict__ rfld,  
             const int * __restrict__ lfld,
             FEBasis ibasis, DifVar idif1, DifVar idif2, 
             const double * __restrict__ bary,  
             double * __restrict__   eval,
             double * __restrict__   jmat,
             double * __restrict__   hmat,
             double * __restrict__  deval,
             double * __restrict__  djmat,
             double * __restrict__  dhmat,
             const double * __restrict__  dfld){

  if constexpr(tdim == 3){
    eval3<szfld,ideg>(rfld,lfld,ibasis,idif1,idif2,bary,eval,jmat,hmat);
  }else{
    eval2<szfld,ideg>(rfld,lfld,ibasis,idif1,idif2,bary,eval,jmat,hmat);
  }

  auto ordent = ORDELT(tdim);

  double fun, dfun[tdim];
  int ideriv = idif1 == DifVar::None ? 0 : 1;
  if(ibasis == FEBasis::Bezier){
    fun = eval_bezierfunc<ideg,tdim>(ordent[ideg][ivar],bary,ideriv,dfun);
  }else{
    fun = eval_lagrangefunc<ideg,tdim>(ordent[ideg][ivar],bary,ideriv,dfun);
  } 

  //-- Evaluation and derivatives
  if(dfld == NULL){
    /* In this case, the Jacobian matrix (of ivar(:) wrt variables)
       is the identity. The gradient w.r.t. idof is simply 
       the basis function times (1, ..., 1). */
    for(int idof = 0; idof < nvar; idof++){
      for(int icmp = 0; icmp < szfld; icmp++){
        deval[szfld*idof + icmp] = 0;
      }
      deval[szfld*idof + idof] = fun;
    }
  }else{
    /* In this case, we may have nvar != szfld 
       The Jacobian of F_K w.r.t. idof is held in dfld:
       dfld(i,j) = d_i rfld(idof,j)  */
    for(int idof = 0; idof < nvar; idof++){
      for(int icmp = 0; icmp < szfld; icmp++){
        deval[szfld*idof + icmp] = fun*dfld[idof*szfld + icmp];
      }
    }
  }
  //-


  //-- Jacobian matrix and derivatives 
  /* Denote D_\alpha^i the i-th component (a tensor) of the derivative w.r.t. 
     P_\alpha. Then:
     D_\alpha^i J_{jk} = D_\alpha^i (d_k J_K^j) 
                       = (i==j) d_k B_\alpha   */
  if(idif1 == DifVar::Bary){
    if(dfld == NULL){
      for(int ii = 0; ii < nvar*tdim*szfld; ii++){
        djmat[ii] = 0;
      }
      for(int idof = 0; idof < nvar; idof++){
        for(int itdim = 0; itdim < tdim; itdim++){
          djmat[tdim*szfld*idof + itdim*szfld + idof] = dfun[itdim];
        }
      }
    }else{
      for(int idof = 0; idof < nvar; idof++){
        for(int itdim = 0; itdim < tdim; itdim++){
          for(int icmp = 0; icmp < szfld; icmp++){
            djmat[idof*tdim*szfld + itdim*szfld + icmp] = dfun[itdim]*dfld[idof*szfld + icmp];
          }
        }
      }
    }
  }
  //-

  //-- Hessian (bary) and derivatives
  // not implemented
  if(idif2 == DifVar::Bary){
    METRIS_THROW_MSG(TODOExcept(),"Unsupported diff2 in eval_d_direct");
    //for(int i=0; i<szfld; i++){
    //  for(int j = 0; j < 3; j++){
    //    dhmat[0*szfld*3 + i*3 + j] = shmat[0*szfld + i].deriv(j);
    //    dhmat[1*szfld*3 + i*3 + j] = shmat[1*szfld + i].deriv(j);
    //    dhmat[2*szfld*3 + i*3 + j] = shmat[2*szfld + i].deriv(j);
    //    dhmat[3*szfld*3 + i*3 + j] = shmat[3*szfld + i].deriv(j);
    //    dhmat[4*szfld*3 + i*3 + j] = shmat[4*szfld + i].deriv(j);
    //    dhmat[5*szfld*3 + i*3 + j] = shmat[5*szfld + i].deriv(j);
    //  }
    //}
    //}
  }
  //-
}












} // End namespace




#endif