//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_EVAL_D_BEZIER__
#define __LOW_EVAL_D_BEZIER__


#include "SANS/Surreal/SurrealS.h"
#include "low_eval_d_helper.hxx"
#include "ho_constants.hxx"
#include "types.hxx"


namespace Metris{

/* ------------

Bezier eval 2D

------------ */

template <typename T, int szfld, int ideg, int nvar, int di=0,int dj=0, int dk=0>
struct eval2_d_bezier{
  eval2_d_bezier(const T& __restrict__  rfld,  
                 DifVar idif1, DifVar idif2, 
                 const double  * __restrict__  bary,  
                 Metris::DLA::VectorS<      szfld,Metris::SurrealS<nvar,double>> *eval,
                 Metris::DLA::MatrixS<2,szfld,Metris::SurrealS<nvar,double>> *jmat,
                 Metris::DLA::MatrixS<3,szfld,Metris::SurrealS<nvar,double>> *hmat){


    Metris::DLA::VectorS<  szfld,Metris::SurrealS<nvar,double>> eva1, eva2, eva3;
    Metris::DLA::MatrixS<2,szfld,Metris::SurrealS<nvar,double>> jmat1,jmat2,jmat3;

    eval2_d_bezier<T,szfld,ideg-1,nvar,di+1,dj+0,dk+0>(rfld,idif2,DifVar::None,bary,&eva1,&jmat1,NULL);
    eval2_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+1,dk+0>(rfld,idif2,DifVar::None,bary,&eva2,&jmat2,NULL);
    eval2_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+0,dk+1>(rfld,idif2,DifVar::None,bary,&eva3,&jmat3,NULL);

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
                 Metris::DLA::VectorS<  szfld,Metris::SurrealS<nvar,double>> *eval,
                 Metris::DLA::MatrixS<2,szfld,Metris::SurrealS<nvar,double>> *jmat,
                 Metris::DLA::MatrixS<3,szfld,Metris::SurrealS<nvar,double>> *hmat){

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
                   Metris::DLA::VectorS<  szfld,Metris::SurrealS<nvar,double>> *eval,
                   Metris::DLA::MatrixS<2,szfld,Metris::SurrealS<nvar,double>> *jmat,
                   Metris::DLA::MatrixS<3,szfld,Metris::SurrealS<nvar,double>> *hmat){

    if (idif1 == DifVar::Bary){
      // The Jacobian is quicker to compute recursively as it only uses the previous degree
      // evals. 

      Metris::DLA::VectorS<  szfld,Metris::SurrealS<nvar,double>> eva1, eva2, eva3;
      Metris::DLA::MatrixS<2,szfld,Metris::SurrealS<nvar,double>> jmat1,jmat2,jmat3;

      eval2_d_bezier<T,szfld,1,nvar,di+1,dj+0,dk+0>(rfld,idif2,DifVar::None,bary,&eva1,&jmat1,NULL);
      eval2_d_bezier<T,szfld,1,nvar,di+0,dj+1,dk+0>(rfld,idif2,DifVar::None,bary,&eva2,&jmat2,NULL);
      eval2_d_bezier<T,szfld,1,nvar,di+0,dj+0,dk+1>(rfld,idif2,DifVar::None,bary,&eva3,&jmat3,NULL);

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
                 Metris::DLA::VectorS<  szfld,Metris::SurrealS<nvar,double>> *eval,
                 Metris::DLA::MatrixS<3,szfld,Metris::SurrealS<nvar,double>> *jmat,
                 Metris::DLA::MatrixS<6,szfld,Metris::SurrealS<nvar,double>> *hmat){


    Metris::DLA::VectorS<  szfld,Metris::SurrealS<nvar,double>> eva1, eva2, eva3, eva4;
    Metris::DLA::MatrixS<3,szfld,Metris::SurrealS<nvar,double>> jmat1,jmat2,jmat3,jmat4;

    eval3_d_bezier<T,szfld,ideg-1,nvar,di+1,dj+0,dk+0,dl+0>(rfld,idif2,DifVar::None,bary,&eva1,&jmat1,NULL);
    eval3_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+1,dk+0,dl+0>(rfld,idif2,DifVar::None,bary,&eva2,&jmat2,NULL);
    eval3_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+0,dk+1,dl+0>(rfld,idif2,DifVar::None,bary,&eva3,&jmat3,NULL);
    eval3_d_bezier<T,szfld,ideg-1,nvar,di+0,dj+0,dk+0,dl+1>(rfld,idif2,DifVar::None,bary,&eva4,&jmat4,NULL);

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
                 Metris::DLA::VectorS<  szfld,Metris::SurrealS<nvar,double>> *eval,
                 Metris::DLA::MatrixS<3,szfld,Metris::SurrealS<nvar,double>> *jmat,
                 Metris::DLA::MatrixS<6,szfld,Metris::SurrealS<nvar,double>> *hmat){

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
                   Metris::DLA::VectorS<  szfld,Metris::SurrealS<nvar,double>> *eval,
                   Metris::DLA::MatrixS<3,szfld,Metris::SurrealS<nvar,double>> *jmat,
                   Metris::DLA::MatrixS<6,szfld,Metris::SurrealS<nvar,double>> *hmat){

    if (idif1 == DifVar::Bary){
      // The Jacobian is quicker to compute recursively as it only uses the previous degree
      // evals. 

      Metris::DLA::VectorS<  szfld,Metris::SurrealS<nvar,double>> eva1, eva2, eva3, eva4;
      Metris::DLA::MatrixS<3,szfld,Metris::SurrealS<nvar,double>> jmat1,jmat2,jmat3,jmat4;

      eval3_d_bezier<T,szfld,1,nvar,di+1,dj+0,dk+0,dl+0>(rfld,idif2,DifVar::None,bary,&eva1,&jmat1,NULL);
      eval3_d_bezier<T,szfld,1,nvar,di+0,dj+1,dk+0,dl+0>(rfld,idif2,DifVar::None,bary,&eva2,&jmat2,NULL);
      eval3_d_bezier<T,szfld,1,nvar,di+0,dj+0,dk+1,dl+0>(rfld,idif2,DifVar::None,bary,&eva3,&jmat3,NULL);
      eval3_d_bezier<T,szfld,1,nvar,di+0,dj+0,dk+0,dl+1>(rfld,idif2,DifVar::None,bary,&eva4,&jmat4,NULL);


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




}//namespace
#endif