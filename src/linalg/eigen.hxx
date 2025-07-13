//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php
#ifndef __METRIS_EIGEN__
#define __METRIS_EIGEN__

#include "../utils/aux_misc.hxx"

namespace Metris{


#ifdef METRIS_USE_LAPACK
template<int ndim>
void geteigsym_LAPACK(const double* met,int nwork,double* rwork,double* eigval,double* eigvec);
#endif


// This function can take SANS::SurrealS as input.
// Returns UNSORTED eigenvalues, call sorteig if needed
template<int ndim, typename T>
void geteigsym(const T* __restrict__ met,T* __restrict__ eigval,T* __restrict__ eigvec);

// This function can take SANS::SurrealS as input.
template<int ndim, typename T>
int geteigsym_Eigen(const T* __restrict__ met,T* __restrict__ eigval,T* __restrict__ eigvec);


//inline void geteigsym(const double* __restrict__ met,double* __restrict__ eigval,double* __restrict__ eigvec){
//	geteigsym<3,double>(met,eigval,eigvec);
//}

// Sort ascending
template<int ndim, typename T>
void sorteig(T* __restrict__ eigval,T* __restrict__ eigvec){

  static_assert(ndim == 2 || ndim == 3);

  if constexpr(ndim == 3){
   
    if(eigval[1] > eigval[2]){
      swi(eigval[1],eigval[2]);
      for(int ii = 0; ii < ndim; ii++) swi(eigvec[ndim*1 + ii], eigvec[ndim*2 + ii]);
    }
    
    if(eigval[0] > eigval[2]){
      swi(eigval[0],eigval[2]);
      for(int ii = 0; ii < ndim; ii++) swi(eigvec[ndim*0 + ii], eigvec[ndim*2 + ii]);
    }
    
    if(eigval[0] > eigval[1]){
      swi(eigval[0],eigval[1]);
      for(int ii = 0; ii < ndim; ii++) swi(eigvec[ndim*0 + ii], eigvec[ndim*1 + ii]);
    }

  }else if(ndim == 2){

    if(eigval[0] > eigval[1]){
      swi(eigval[0],eigval[1]);
      for(int ii = 0; ii < ndim; ii++) swi(eigvec[ndim*0 + ii], eigvec[ndim*1 + ii]);
    }

  }
}


// -----------------------------------------------------------------------------
// Compute R^T D R with R eigvec matrix
template<int ndim, typename T>
inline void eig2met(const T* __restrict__ eigval, const T* __restrict__ eigvec, T* __restrict__ met){
	static_assert(ndim == 2 || ndim == 3);
	if(ndim == 2){
		met[0] = eigval[0]*eigvec[2*0+0]*eigvec[2*0+0] 
		       + eigval[1]*eigvec[2*1+0]*eigvec[2*1+0]; 
		met[1] = eigval[0]*eigvec[2*0+0]*eigvec[2*0+1]
		       + eigval[1]*eigvec[2*1+0]*eigvec[2*1+1];
		met[2] = eigval[0]*eigvec[2*0+1]*eigvec[2*0+1] 
		       + eigval[1]*eigvec[2*1+1]*eigvec[2*1+1];
	}else{
		met[0] = eigval[0]*eigvec[0]*eigvec[0] 
		       + eigval[1]*eigvec[3]*eigvec[3] 
		       + eigval[2]*eigvec[6]*eigvec[6];
		met[1] = eigval[0]*eigvec[0]*eigvec[1]
		       + eigval[1]*eigvec[3]*eigvec[4]
		       + eigval[2]*eigvec[6]*eigvec[7];
		met[2] = eigval[0]*eigvec[1]*eigvec[1] 
		       + eigval[1]*eigvec[4]*eigvec[4] 
		       + eigval[2]*eigvec[7]*eigvec[7];
		met[3] = eigval[0]*eigvec[0]*eigvec[2]
		       + eigval[1]*eigvec[3]*eigvec[5]
		       + eigval[2]*eigvec[6]*eigvec[8];
		met[4] = eigval[0]*eigvec[1]*eigvec[2]
		       + eigval[1]*eigvec[4]*eigvec[5]
		       + eigval[2]*eigvec[7]*eigvec[8];
		met[5] = eigval[0]*eigvec[2]*eigvec[2] 
		       + eigval[1]*eigvec[5]*eigvec[5] 
		       + eigval[2]*eigvec[8]*eigvec[8];
	}
}

template<int ndim, typename T>
inline void eig2met_sum(const T* __restrict__ eigval, const T* __restrict__ eigvec, T* __restrict__ met){
	static_assert(ndim == 2 || ndim == 3);
	if(ndim == 2){
		met[0] += eigval[0]*eigvec[2*0+0]*eigvec[2*0+0] 
		        + eigval[1]*eigvec[2*1+0]*eigvec[2*1+0]; 
		met[1] += eigval[0]*eigvec[2*0+0]*eigvec[2*0+1]
		        + eigval[1]*eigvec[2*1+0]*eigvec[2*1+1];
		met[2] += eigval[0]*eigvec[2*0+1]*eigvec[2*0+1] 
		        + eigval[1]*eigvec[2*1+1]*eigvec[2*1+1];
	}else{
		met[0] += eigval[0]*eigvec[0]*eigvec[0] 
		        + eigval[1]*eigvec[3]*eigvec[3] 
		        + eigval[2]*eigvec[6]*eigvec[6];
		met[1] += eigval[0]*eigvec[0]*eigvec[1]
		        + eigval[1]*eigvec[3]*eigvec[4]
		        + eigval[2]*eigvec[6]*eigvec[7];
		met[2] += eigval[0]*eigvec[1]*eigvec[1] 
		        + eigval[1]*eigvec[4]*eigvec[4] 
		        + eigval[2]*eigvec[7]*eigvec[7];
		met[3] += eigval[0]*eigvec[0]*eigvec[2]
		        + eigval[1]*eigvec[3]*eigvec[5]
		        + eigval[2]*eigvec[6]*eigvec[8];
		met[4] += eigval[0]*eigvec[1]*eigvec[2]
		        + eigval[1]*eigvec[4]*eigvec[5]
		        + eigval[2]*eigvec[7]*eigvec[8];
		met[5] += eigval[0]*eigvec[2]*eigvec[2] 
		        + eigval[1]*eigvec[5]*eigvec[5] 
		        + eigval[2]*eigvec[8]*eigvec[8];
	}
}




 




} // End namespace








#endif