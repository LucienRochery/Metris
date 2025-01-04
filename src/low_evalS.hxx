//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php




#ifdef THIS_IS_DEPRECATED



#ifndef __LOW_EVALS__
#define __LOW_EVALS__

#include "types.hxx"
#include "ho_constants.hxx"
#include "aux_utils.hxx"


#include "codegen_lagrange.hxx"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"

#include <boost/hana.hpp> 
namespace hana = boost::hana;
using namespace hana::literals;

namespace Metris{
// The main eval routine. Can be used for scalar fields (szfld = 1), element geometry (szfld = 3)
// or metric interpolation (szfld = 6). 
// eval3S: tetra
// eval2: triangles
// eval1: edges
template <int szfld, int ideg, int ilag, int idif1, int idif2>
void eval3S(const dblAr2  &rfld, 
            const SANS::DLA::VectorS<tetnpps[ideg],int> &lfld, 
            const SANS::DLA::VectorS<4,double> &bary, 
            SANS::DLA::VectorS<  szfld,double> &eval, 
            SANS::DLA::VectorS<3*szfld,double> &jmat,  
            SANS::DLA::VectorS<9*szfld,double> &hmat);


/*
With Hessian computed from Jmats
The Hessian is compressed as column first as usual:
1 2 4
  3 5
    6
However, note that each entry is, itself, comprised of szfld values. 
*/
template <int szfld, int ideg, int idif1, int idif2, int di=0,int dj=0, int dk=0, int dl=0>
struct eval3S_bezier{
  eval3S_bezier(const dblAr2 & __restrict__  rfld,  
               const SANS::DLA::VectorS<tetnpps[ideg],int> &lfld,
               const SANS::DLA::VectorS<4,double> &bary,  
               SANS::DLA::VectorS<  szfld,double> &eval,
               SANS::DLA::VectorS<3*szfld,double> &jmat,
               SANS::DLA::VectorS<9*szfld,double> &hmat){
    SANS::DLA::VectorS<  szfld,double> eva1, eva2, eva3, eva4;
    SANS::DLA::VectorS<3*szfld,double> jmat1, jmat2, jmat3, jmat4;
    SANS::DLA::VectorS<9*szfld,double> hdum;

    if constexpr(idif2 == 0){
      eval3S_bezier<szfld,ideg-1,0,0,di+1,dj+0,dk+0,dl+0>(rfld,lfld,bary,eva1,jmat1,hdum);
      eval3S_bezier<szfld,ideg-1,0,0,di+0,dj+1,dk+0,dl+0>(rfld,lfld,bary,eva2,jmat2,hdum);
      eval3S_bezier<szfld,ideg-1,0,0,di+0,dj+0,dk+1,dl+0>(rfld,lfld,bary,eva3,jmat3,hdum);
      eval3S_bezier<szfld,ideg-1,0,0,di+0,dj+0,dk+0,dl+1>(rfld,lfld,bary,eva4,jmat4,hdum);
    }else{
      eval3S_bezier<szfld,ideg-1,1,0,di+1,dj+0,dk+0,dl+0>(rfld,lfld,bary,eva1,jmat1,hdum);
      eval3S_bezier<szfld,ideg-1,1,0,di+0,dj+1,dk+0,dl+0>(rfld,lfld,bary,eva2,jmat2,hdum);
      eval3S_bezier<szfld,ideg-1,1,0,di+0,dj+0,dk+1,dl+0>(rfld,lfld,bary,eva3,jmat3,hdum);
      eval3S_bezier<szfld,ideg-1,1,0,di+0,dj+0,dk+0,dl+1>(rfld,lfld,bary,eva4,jmat4,hdum);
    }

    for(int i = 0; i < szfld; i++){
      eval[i] = bary[0] * eva1[i]
              + bary[1] * eva2[i]
              + bary[2] * eva3[i]
              + bary[3] * eva4[i];
    }


    if constexpr(idif1 > 0){
      for(int i = 0; i < szfld; i++){
        double d1 = eva1[i];
        jmat[0*szfld+i] = ideg * (eva2[i] - d1);
        jmat[1*szfld+i] = ideg * (eva3[i] - d1);
        jmat[2*szfld+i] = ideg * (eva4[i] - d1);
      }
    }

    if constexpr (idif2 > 0){
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

template <int szfld, int idif1, int idif2, int di,int dj, int dk, int dl>
struct eval3S_bezier<szfld,1,idif1,idif2,di,dj,dk,dl>{
  eval3S_bezier(const dblAr2 & __restrict__  rfld,  
                const SANS::DLA::VectorS<tetnpps[1 + di + dj + dk + dl],int> &lfld,
                const SANS::DLA::VectorS<4,double> &bary,  
                SANS::DLA::VectorS<  szfld,double> &eval,
                SANS::DLA::VectorS<3*szfld,double> &jmat,
                SANS::DLA::VectorS<9*szfld,double> &hmat){

    constexpr int ideg = 1 + di + dj + dk + dl; // The true degree. 

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

    if constexpr(idif1 == 1){
      for(int i = 0; i < szfld; i++){
        jmat[0*szfld+i] = rfld[lfld[i0100]][i] - rfld[lfld[i1000]][i];
        jmat[1*szfld+i] = rfld[lfld[i0010]][i] - rfld[lfld[i1000]][i];
        jmat[2*szfld+i] = rfld[lfld[i0001]][i] - rfld[lfld[i1000]][i];
      }
    }

    if constexpr(idif2 == 1){
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


template <int szfld, int idif1,int idif2,int di,int dj, int dk, int dl>
  struct eval3S_bezier<szfld,2,idif1,idif2,di,dj,dk,dl>{
    eval3S_bezier(const dblAr2 & __restrict__  rfld,  
                 const SANS::DLA::VectorS<tetnpps[2 + di + dj + dk + dl],int> &lfld,
                 const SANS::DLA::VectorS<4,double> &bary,
                 SANS::DLA::VectorS<  szfld,double> &eval,
                 SANS::DLA::VectorS<3*szfld,double> &jmat,
                 SANS::DLA::VectorS<9*szfld,double> &hmat){
    constexpr int ideg = 2 + di + dj + dk + dl; // The true degree. 

    
    if constexpr (idif1 > 0){


        // The Jacobian is quicker to compute recursively as it only uses the previous degree
        // evals. 

      SANS::DLA::VectorS<  szfld,double> eva1, eva2, eva3, eva4;
      SANS::DLA::VectorS<3*szfld,double> jmat1, jmat2, jmat3, jmat4;
      SANS::DLA::VectorS<9*szfld,double> hdum;
    
      if constexpr(idif2 == 0){
        eval3S_bezier<szfld,ideg-1,0,0,di+1,dj+0,dk+0,dl+0>(rfld,lfld,bary,eva1,jmat1,hdum);
        eval3S_bezier<szfld,ideg-1,0,0,di+0,dj+1,dk+0,dl+0>(rfld,lfld,bary,eva2,jmat2,hdum);
        eval3S_bezier<szfld,ideg-1,0,0,di+0,dj+0,dk+1,dl+0>(rfld,lfld,bary,eva3,jmat3,hdum);
        eval3S_bezier<szfld,ideg-1,0,0,di+0,dj+0,dk+0,dl+1>(rfld,lfld,bary,eva4,jmat4,hdum);
      }else{
        eval3S_bezier<szfld,ideg-1,1,0,di+1,dj+0,dk+0,dl+0>(rfld,lfld,bary,eva1,jmat1,hdum);
        eval3S_bezier<szfld,ideg-1,1,0,di+0,dj+1,dk+0,dl+0>(rfld,lfld,bary,eva2,jmat2,hdum);
        eval3S_bezier<szfld,ideg-1,1,0,di+0,dj+0,dk+1,dl+0>(rfld,lfld,bary,eva3,jmat3,hdum);
        eval3S_bezier<szfld,ideg-1,1,0,di+0,dj+0,dk+0,dl+1>(rfld,lfld,bary,eva4,jmat4,hdum);
      }


      for(int i = 0; i < szfld; i++){
        eval[i] = bary[0] * eva1[i]
                + bary[1] * eva2[i]
                + bary[2] * eva3[i]
                + bary[3] * eva4[i];
      }


      for(int i = 0; i < szfld;i++){
        double d1 = eva1[i];
        jmat[0*szfld+i] = ideg * (eva2[i] - d1);
        jmat[1*szfld+i] = ideg * (eva3[i] - d1);
        jmat[2*szfld+i] = ideg * (eva4[i] - d1);
      }


      if constexpr (idif2 > 0){
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
template <int szfld, int ideg, int ilag, int idif1, int idif2>
void eval3S(const dblAr2 & __restrict__  rfld,  
           const SANS::DLA::VectorS<tetnpps[ideg],int> &lfld,
           const SANS::DLA::VectorS<4,double> &bary, 
           SANS::DLA::VectorS<  szfld,double> &eval, 
           SANS::DLA::VectorS<3*szfld,double> &jmat, 
           SANS::DLA::VectorS<9*szfld,double> &hmat){
  if constexpr (ideg == 1){
    
    for(int i=0;i<szfld;i++){
      eval[i] = bary[0]*rfld[lfld[0]][i]
              + bary[1]*rfld[lfld[1]][i]
              + bary[2]*rfld[lfld[2]][i]
              + bary[3]*rfld[lfld[3]][i];
    }
    if constexpr(idif1 > 0){
      
      for(int i=0;i<szfld;i++){
        jmat[0*szfld+i] = rfld[lfld[1]][i] - rfld[lfld[0]][i];
        jmat[1*szfld+i] = rfld[lfld[2]][i] - rfld[lfld[0]][i];
        jmat[2*szfld+i] = rfld[lfld[3]][i] - rfld[lfld[0]][i];
      }
    }
    if constexpr(idif2 > 0){
      
      for(int i=0;i<szfld;i++){
        hmat[0*szfld+i] = 0.0;
        hmat[1*szfld+i] = 0.0;
        hmat[2*szfld+i] = 0.0;
        hmat[3*szfld+i] = 0.0;
        hmat[4*szfld+i] = 0.0;
        hmat[5*szfld+i] = 0.0;
      }
    }
    return;
  } 

  if constexpr (ilag <= 0){
    eval3S_bezier<szfld,ideg,idif1,idif2>(rfld, lfld, bary, eval, jmat, hmat);
  }else{
    //eval3S_lagrange<szfld,ideg,idif1>(rfld, lfld, bary, eval, jmat);
  } 
  return;
}






} // End namespace



#endif

#endif // deprecated