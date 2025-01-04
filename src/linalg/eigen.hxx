//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_EIGEN__
#define __METRIS_EIGEN__

namespace Metris{


template<int ndim, typename T>
void geteigsym(const T* met,int nwork,double* rwork,T* eigval,T* eigvec);


template<int ndim, typename T>
void geteigsym(const T* __restrict__ met,T* __restrict__ eigval,T* __restrict__ eigvec);


//inline void geteigsym(const double* __restrict__ met,double* __restrict__ eigval,double* __restrict__ eigvec){
//	geteigsym<3,double>(met,eigval,eigvec);
//}



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