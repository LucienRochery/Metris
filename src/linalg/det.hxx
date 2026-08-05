//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LINALG_DET__
#define __METRIS_LINALG_DET__

#include "symmat.hxx"
#include "Metris_LAPACK.hxx"

#include <type_traits>
#include <cmath>

namespace Metris{

// -----------------------------------------------------------------------------
template<int ndimn, typename T = double>
inline T detmat(const T mat[]){
	static_assert(ndimn == 2 || ndimn == 3);
	if constexpr(ndimn == 2){
		return mat[2*0+0]*mat[2*1+1] - mat[2*1+0]*mat[2*0+1];
	}else{
		return mat[3*0+0]*(mat[3*1+1]*mat[3*2+2] - mat[3*1+2]*mat[3*2+1])
         + mat[3*1+0]*(mat[3*2+1]*mat[3*0+2] - mat[3*2+2]*mat[3*0+1])
         + mat[3*2+0]*(mat[3*0+1]*mat[3*1+2] - mat[3*0+2]*mat[3*1+1]);
	}
}
// With MatrixS type.
template<int ndimn, typename T = double>
inline T detmat(const Metris::DLA::MatrixS<ndimn,ndimn,T> &mat){
  static_assert(ndimn == 2 || ndimn == 3);
  if constexpr(ndimn == 2){
    return mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1);
  }else{
    return mat(0,0)*(mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1))
         + mat(1,0)*(mat(2,1)*mat(0,2) - mat(2,2)*mat(0,1))
         + mat(2,0)*(mat(0,1)*mat(1,2) - mat(0,2)*mat(1,1));
  }
}

// Determinant of matrix whose columns or lines are given by v1, v2(, v3)
template<typename T = double>
inline T detvec2(const T* __restrict__ v1, 
                 const T* __restrict__ v2){
  return v1[0]*v2[1] - v1[1]*v2[0];
}
template<typename T = double>
inline T detvec3(const T* __restrict__ v1, 
                 const T* __restrict__ v2, 
                 const T* __restrict__ v3){
  return v1[0]*(v2[1]*v3[2]-v3[1]*v2[2])
       + v1[1]*(v2[2]*v3[0]-v3[2]*v2[0])
       + v1[2]*(v2[0]*v3[1]-v3[0]*v2[1]);
}

// determinant of matrix with columns x1-x2, y1-y2
template<typename T=double>
inline T detvdif2(const T* x1,const T* x2, 
            const T* y1,const T* y2){
  return (x1[0] - x2[0])*(y1[1] - y2[1]) 
       - (x1[1] - x2[1])*(y1[0] - y2[0]);
}

// determinant of matrix with columns x1-x2, y1-y2, z1-z2
template <typename T=double>
inline T detvdif3(const T* x1,const T* x2,
                   const T* y1,const T* y2,
                   const T* z1,const T* z2){
  return (x1[0] - x2[0])*( (y1[1] - y2[1])*(z1[2] - z2[2]) - (z1[1] - z2[1])*(y1[2] - y2[2]))
       + (x1[1] - x2[1])*( (y1[2] - y2[2])*(z1[0] - z2[0]) - (z1[2] - z2[2])*(y1[0] - y2[0]))
       + (x1[2] - x2[2])*( (y1[0] - y2[0])*(z1[1] - z2[1]) - (z1[0] - z2[0])*(y1[1] - y2[1]));
}




// Same as detvec2/3 but get sub-determinant obtained by extracting i-th line,
// j-th column. Now the vs are necessarily columns (or permute i,j). 
// Note the sub-determinants are SIGNED : i.e. multiplied by (-1)^{i+j}. 
// In other words, an expansion by a column is simply a sum of these subdets,
// without need to apply the signs. 

// Second note: v1 and v2 are lines ! 
template<typename T = double>
inline T subdetvec2(const T* __restrict__ v1, 
                    const T* __restrict__ v2,
                    int i, int j){
  if      (i == 0 && j == 0){
    // return 2,2
    return v2[1];
  }else if(i == 0 && j == 1){
    // return 2,1
    return -v2[0];
  }else if(i == 1 && j == 0){
    // return 1,2
    return -v1[1];
  }
  // return 1,1
  return v1[0];
}
template<typename T = double>
inline T subdetvec3(const T* __restrict__ v1, 
                    const T* __restrict__ v2, 
                    const T* __restrict__ v3,
                    int i, int j){
  // Note the signs are naturally obtained by applying rotations to the indices
  // exclude v1 for 3: 
  if      (i == 0 && j == 0){
    return v2[1]*v3[2] 
         - v2[2]*v3[1];
  }else if(i == 0 && j == 1){
    return v2[2]*v3[0] 
         - v2[0]*v3[2];
  }else if(i == 0 && j == 2){
    return v2[0]*v3[1] 
         - v2[1]*v3[0];
  // exclude v2 for 3: 
  }else if(i == 1 && j == 0){
    return v3[1]*v1[2] 
         - v3[2]*v1[1]; 
  }else if(i == 1 && j == 1){
    return v3[2]*v1[0] 
         - v3[0]*v1[2]; 
  }else if(i == 1 && j == 2){
    return v3[0]*v1[1] 
         - v3[1]*v1[0];
  // exclude v3 for 3: 
  }else if(i == 2 && j == 0){
    return v1[1]*v2[2] 
         - v1[2]*v2[1]; 
  }else if(i == 2 && j == 1){
    return v1[2]*v2[0] 
         - v1[0]*v2[2]; 
  }
  //else if(i == 2 && j == 2){ // comment to suppress warnings
  return v1[0]*v2[1] 
       - v1[1]*v2[0];
  //}
}


template<int ndim, typename ftype>
ftype detsym_Eigen_LLT(const ftype *met);

#ifdef METRIS_USE_LAPACK
template<int ndim>
double detsym_LAPACK(const double* met);
#endif


// -----------------------------------------------------------------------------
template<int ndimn, typename T = double>
inline T detsym(const T met[]){
	static_assert(ndimn == 2 || ndimn == 3);
	if constexpr(ndimn == 2){
    return met[0]*met[2] - met[1]*met[1];
  }else{
    return met[0]*(met[2]*met[5]-met[4]*met[4])
         + met[1]*(met[4]*met[3]-met[5]*met[1])
         + met[3]*(met[1]*met[4]-met[3]*met[2]);
  }
}

// -----------------------------------------------------------------------------
template<int ndimn, typename T = double>
inline T detsym2(const T met[]){
  static_assert(ndimn == 2 || ndimn == 3);
  if constexpr(ndimn == 2){
    return met[0]*met[2] - met[1]*met[1];
  }else{
    #ifdef METRIS_USE_LAPACK
    if constexpr (std::is_same<T,double>::value){
      return detsym_LAPACK<ndimn>(met);
    }else{
    #endif
      return detsym_Eigen_LLT<ndimn>(met);
    #ifdef METRIS_USE_LAPACK
    }
    #endif
  }
}
// -----------------------------------------------------------------------------
template<int ndimn, typename T = double>
inline T detsym3(const T met[]){
  static_assert(ndimn == 2 || ndimn == 3);
  if constexpr(ndimn == 2){
    T mx = std::abs(met[0]); 
    mx = mx > std::abs(met[1]) ? mx : std::abs(met[1]);
    mx = mx > std::abs(met[2]) ? mx : std::abs(met[2]);
    T met2[3] ;
    for(int ii = 0; ii < 3 ;ii++) met2[ii] = met[ii] / mx;
    double det = met2[0]*met2[2] - met2[1]*met2[1];
    return det*mx*mx;
  }else{
    T mx = std::abs(met[0]); 
    mx = mx > std::abs(met[1]) ? mx : std::abs(met[1]);
    mx = mx > std::abs(met[2]) ? mx : std::abs(met[2]);
    mx = mx > std::abs(met[3]) ? mx : std::abs(met[3]);
    mx = mx > std::abs(met[4]) ? mx : std::abs(met[4]);
    mx = mx > std::abs(met[5]) ? mx : std::abs(met[5]);
    T met2[6] ;
    for(int ii =0; ii < 6 ;ii++) met2[ii] = met[ii] / mx;
    double det = met2[0]*(met2[2]*met2[5]-met2[4]*met2[4])
               + met2[1]*(met2[4]*met2[3]-met2[5]*met2[1])
               + met2[3]*(met2[1]*met2[4]-met2[3]*met2[2]);
    return det*mx*mx*mx;
  }
}

// MatrixS version
template<int ndimn, typename T = double>
inline T detsym2(const Metris::DLA::MatSymS<ndimn,T> &met){
  static_assert(ndimn == 2 || ndimn == 3);
  if constexpr(ndimn == 2){
    return met[0]*met[2] - met[1]*met[1];
  }else{
    T mx = std::abs(met[0]); 
    mx = mx > std::abs(met[1]) ? mx : std::abs(met[1]);
    mx = mx > std::abs(met[2]) ? mx : std::abs(met[2]);
    mx = mx > std::abs(met[3]) ? mx : std::abs(met[3]);
    mx = mx > std::abs(met[4]) ? mx : std::abs(met[4]);
    mx = mx > std::abs(met[5]) ? mx : std::abs(met[5]);
    T met2[6] ;
    for(int ii = 0; ii < 6 ;ii++) met2[ii] = met[ii] / mx;
    T det = met2[0]*(met2[2]*met2[5]-met2[4]*met2[4])
          + met2[1]*(met2[4]*met2[3]-met2[5]*met2[1])
          + met2[3]*(met2[1]*met2[4]-met2[3]*met2[2]);
    return det*mx*mx*mx;
  }
}


//// -----------------------------------------------------------------------------
//template<int ndimn, typename T = double>
//inline T detsym3(const T met[]){
//  static_assert(ndimn == 2 || ndimn == 3);
//  if constexpr(ndimn == 2){
//    return met[0]*met[2] - met[1]*met[1];
//  }else{
//    double avg[3] = {met[0] + met[1] + met[3];
//    return met[0]*(met[2]*met[5]-met[4]*met[4])
//         + met[1]*(met[4]*met[3]-met[5]*met[1])
//         + met[3]*(met[1]*met[4]-met[3]*met[2]);
//  }
//}

inline void vecprod(const double*__restrict__ u, const double*__restrict__ v, 
                    double*__restrict__ out){
  out[0] = u[1]*v[2] - u[2]*v[1];
  out[1] = u[2]*v[0] - u[0]*v[2];
  out[2] = u[0]*v[1] - u[1]*v[0];
}
#ifdef USE_MULTIPRECISION
inline void vecprod(const float4*__restrict__ u, const float4*__restrict__ v, 
                    float4*__restrict__ out){
  out[0] = u[1]*v[2] - u[2]*v[1];
  out[1] = u[2]*v[0] - u[0]*v[2];
  out[2] = u[0]*v[1] - u[1]*v[0];
}
inline void vecprod(const float8*__restrict__ u, const float8*__restrict__ v, 
                    float8*__restrict__ out){
  out[0] = u[1]*v[2] - u[2]*v[1];
  out[1] = u[2]*v[0] - u[0]*v[2];
  out[2] = u[0]*v[1] - u[1]*v[0];
}
#endif

// Compute (u1 - u2) x (v1 - v2)
inline void vecprod_vdif(const double*__restrict__ u1, const double*__restrict__ u2, 
                         const double*__restrict__ v1, const double*__restrict__ v2, 
                         double*__restrict__ out){
  out[0] = (u1[1] - u2[1])*(v1[2] - v2[2]) - (u1[2] - u2[2])*(v1[1] - v2[1]);
  out[1] = (u1[2] - u2[2])*(v1[0] - v2[0]) - (u1[0] - u2[0])*(v1[2] - v2[2]);
  out[2] = (u1[0] - u2[0])*(v1[1] - v2[1]) - (u1[1] - u2[1])*(v1[0] - v2[0]);
}





} // End namespace

#endif
