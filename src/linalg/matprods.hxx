//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_MATPRODS__
#define __METRIS_MATPRODS__

#include "../linalg/sym3idx.hxx"


namespace Metris{

// mat[i][j] = mat[3*i+j]
// mat1 * mat2
template<int ndimn, typename T1 = double, typename T2 = double, typename T3 = double>
inline void matXmat(const T1* mat1, const T2* mat2, T3* mat3){
	for(int ii = 0; ii < ndimn; ii++){
		for(int jj  = 0; jj < ndimn; jj++){
			mat3[ndimn*ii + jj] = mat1[ndimn*ii + 0]*mat2[ndimn*0 + jj];
			for(int kk = 1; kk < ndimn; kk++){
				mat3[ndimn*ii + jj] += mat1[ndimn*ii + kk]*mat2[ndimn*kk + jj];
			}
		}
	}
}

// mat1 is size n1 x n2 
// mat2 is size n2 x n3  
template<int n1, int n2, int n3, typename T1 = double, typename T2 = double, typename T3 = double>
inline void matXmat(const T1* mat1, const T2* mat2, T3* mat3){
  for(int ii = 0; ii < n1; ii++){
    for(int jj  = 0; jj < n3; jj++){
      mat3[n3*ii + jj] = mat1[n2*ii + 0]*mat2[n3*0 + jj];
      for(int kk = 1; kk < n2; kk++){
        mat3[n3*ii + jj] += mat1[n2*ii + kk]*mat2[n3*kk + jj];
      }
    }
  }
}


template<int ndimn, typename T1 = double, typename T2 = double, typename T3 = double>
inline void matXtmat(const T1* mat1, const T2* mat2, T3* mat3){
	for(int ii = 0; ii < ndimn; ii++){
		for(int jj  = 0; jj < ndimn; jj++){
			mat3[ndimn*ii + jj] = mat1[ndimn*ii + 0]*mat2[ndimn*jj + 0];
			for(int kk = 1; kk < ndimn; kk++){
				mat3[ndimn*ii + jj] += mat1[ndimn*ii + kk]*mat2[ndimn*jj + kk];
			}
		}
	}
}

// mat1 is size n1 x n2 
// mat2 is size n3 x n2 (then transposed) 
template<int n1, int n2, int n3, typename T1 = double, typename T2 = double, typename T3 = double>
inline void matXtmat(const T1* mat1, const T2* mat2, T3* mat3){
  for(int ii = 0; ii < n1; ii++){
    for(int jj  = 0; jj < n3; jj++){
      mat3[n3*ii + jj] = mat1[n2*ii + 0]*mat2[n2*jj + 0];
      for(int kk = 1; kk < n2; kk++){
        mat3[n3*ii + jj] += mat1[n2*ii + kk]*mat2[n2*jj + kk];
      }
    }
  }
}


// mat1*mat2 both symmetric order 1 2 4
//                                  3 5
//                                    6
template<int ndimn, typename T = double>
inline void symXsymsub_fac(const T* mat1, const T* mat2, 
                              const double fac, T* mat3){
	for(int ii = 0; ii < ndimn; ii++){
		for(int jj  = 0; jj < ndimn; jj++){
			mat3[sym3idx(ii,jj)] = mat1[sym3idx(ii,0)]*mat2[sym3idx(0,jj)] * fac;
			for(int kk = 1; kk < ndimn; kk++){
				mat3[sym3idx(ii,jj)] += mat1[sym3idx(ii,kk)]*mat2[sym3idx(kk,jj)] * fac;
			}
		}
	}
}
template<int ndimn, typename T = double>
inline void symXsymsub(const T* mat1, const T* mat2, T* mat3){
	for(int ii = 0; ii < ndimn; ii++){
		for(int jj  = 0; jj < ndimn; jj++){
			mat3[sym3idx(ii,jj)] = mat1[sym3idx(ii,0)]*mat2[sym3idx(0,jj)];
			for(int kk = 1; kk < ndimn; kk++){
				mat3[sym3idx(ii,jj)] += mat1[sym3idx(ii,kk)]*mat2[sym3idx(kk,jj)];
			}
		}
	}
}

// Same as symXsymsub_fac<3> but add into mat3 rather than write over
template<int ndimn, typename T = double>
inline void symXsymadd_fac(const T* mat1, const T* mat2, 
                           const double fac, T* mat3){
	for(int ii = 0; ii < ndimn; ii++){
		for(int jj  = 0; jj < ndimn; jj++){
			for(int kk = 0; kk < ndimn; kk++){
				mat3[sym3idx(ii,jj)] += mat1[sym3idx(ii,kk)]*mat2[sym3idx(kk,jj)]*fac;
			}
		}
	}
}

inline void ABpBA3symsub(const double mat1[], const double mat2[], 
                         double mat3[]){
  mat3[0] = 2*(mat2[0]*mat1[0] + mat2[1]*mat1[1] + mat2[3]*mat1[3]);
  mat3[1] = mat1[1]*(mat2[0]+mat2[2]) + mat2[1]*(mat1[0]+mat1[2]) + mat2[3]*mat1[4] + mat1[3]*mat2[4];
  mat3[3] = mat1[3]*(mat2[0]+mat2[5]) + mat2[3]*(mat1[0]+mat1[5]) + mat2[1]*mat1[4] + mat1[1]*mat2[4];
  mat3[2] = 2*(mat2[1]*mat1[1] + mat2[2]*mat1[2] + mat2[4]*mat1[4]);
  mat3[4] = mat1[4]*(mat2[2]+mat2[5]) + mat2[4]*(mat1[2]+mat1[5]) + mat2[1]*mat1[3] + mat1[1]*mat2[3];
  mat3[5] = 2*(mat2[3]*mat1[3] + mat2[4]*mat1[4] + mat2[5]*mat1[5]);
}

// mat * inp store in out
// Mat stored line first (C)
inline void mat3vec(const double mat[], const double inp[], double out[]){
	out[0] = mat[0]*inp[0] + mat[1]*inp[1] + mat[2]*inp[2];
	out[1] = mat[3]*inp[0] + mat[4]*inp[1] + mat[5]*inp[2];
	out[2] = mat[6]*inp[0] + mat[7]*inp[1] + mat[8]*inp[2];
}

// mat * (in1 - in2) store in out
// Mat stored line first (C)
inline void mat3vdf(const double mat[], const double in1[], const double in2[], double out[]){
	out[0] = mat[0]*(in1[0]-in2[0]) + mat[1]*(in1[1]-in2[1]) + mat[2]*(in1[2]-in2[2]);
	out[1] = mat[3]*(in1[0]-in2[0]) + mat[4]*(in1[1]-in2[1]) + mat[5]*(in1[2]-in2[2]);
	out[2] = mat[6]*(in1[0]-in2[0]) + mat[7]*(in1[1]-in2[1]) + mat[8]*(in1[2]-in2[2]);
}

// mat^T * inp store in out
// Mat stored line first (C)
inline void mat3vect(const double mat[], const double inp[], double out[]){
	out[0] = mat[0]*inp[0] + mat[3]*inp[1] + mat[6]*inp[2];
	out[1] = mat[1]*inp[0] + mat[4]*inp[1] + mat[7]*inp[2];
	out[2] = mat[2]*inp[0] + mat[5]*inp[1] + mat[8]*inp[2];
}


// mat^T * (in1 - in2) store in out
// Mat stored line first (C)
inline void matvdft(int n, const double mat[], const double in1[], const double in2[], double out[]){
  if(n == 3){
  	out[0] = mat[0]*(in1[0]-in2[0]) + mat[3]*(in1[1]-in2[1]) + mat[6]*(in1[2]-in2[2]);
  	out[1] = mat[1]*(in1[0]-in2[0]) + mat[4]*(in1[1]-in2[1]) + mat[7]*(in1[2]-in2[2]);
    out[2] = mat[2]*(in1[0]-in2[0]) + mat[5]*(in1[1]-in2[1]) + mat[8]*(in1[2]-in2[2]);
  }else if(n == 2){
    out[0] = mat[0]*(in1[0]-in2[0]) + mat[2]*(in1[1]-in2[1]);
    out[1] = mat[1]*(in1[0]-in2[0]) + mat[3]*(in1[1]-in2[1]);
  }else{
    out[0] = mat[0]*(in1[0]-in2[0]); 
  }
}

// M*v with M symmetric such as the others. 
inline void sym3vec(const double mat[], const double inp[], double out[]){
	out[0] = mat[0]*inp[0] + mat[1]*inp[1] + mat[3]*inp[2];
	out[1] = mat[1]*inp[0] + mat[2]*inp[1] + mat[4]*inp[2];
	out[2] = mat[3]*inp[0] + mat[4]*inp[1] + mat[5]*inp[2];
}

inline void sym3tmat(const double met[], const double mat[], double out[]){
	out[0] = met[0]*mat[0] + met[1]*mat[1] + met[3]*mat[2];
	out[1] = met[0]*mat[3] + met[1]*mat[4] + met[3]*mat[5];
	out[2] = met[0]*mat[6] + met[1]*mat[7] + met[3]*mat[8];

	out[3] = met[1]*mat[0] + met[2]*mat[1] + met[4]*mat[2];
	out[4] = met[1]*mat[3] + met[2]*mat[4] + met[4]*mat[5];
	out[5] = met[1]*mat[6] + met[2]*mat[7] + met[4]*mat[8];

	out[6] = met[3]*mat[0] + met[4]*mat[1] + met[5]*mat[2];
	out[7] = met[3]*mat[3] + met[4]*mat[4] + met[5]*mat[5];
	out[8] = met[3]*mat[6] + met[4]*mat[7] + met[5]*mat[8];
}
inline void sym3mat(const double met[], const double mat[], double out[]){
	out[0] = met[0]*mat[0] + met[1]*mat[3] + met[3]*mat[6];
	out[1] = met[0]*mat[1] + met[1]*mat[4] + met[3]*mat[7];
	out[2] = met[0]*mat[2] + met[1]*mat[5] + met[3]*mat[8];

	out[3] = met[1]*mat[0] + met[2]*mat[3] + met[4]*mat[6];
	out[4] = met[1]*mat[1] + met[2]*mat[4] + met[4]*mat[7];
	out[5] = met[1]*mat[2] + met[2]*mat[5] + met[4]*mat[8];

	out[6] = met[3]*mat[0] + met[4]*mat[3] + met[5]*mat[6];
	out[7] = met[3]*mat[1] + met[4]*mat[4] + met[5]*mat[7];
	out[8] = met[3]*mat[2] + met[4]*mat[5] + met[5]*mat[8];
}

template <int ndimn>
inline void symXvec(const double*__restrict__ met,
	     				      const double*__restrict__ ve1,
	           				     	double*__restrict__ ve2){
	static_assert(ndimn == 2 || ndimn == 3);
	for(int ii = 0; ii < ndimn; ii++){
		ve2[ii] = met[sym3idx(ii,0)]*ve1[0];
		for(int jj = 1; jj < ndimn; jj++){
			ve2[ii] += met[sym3idx(ii,jj)]*ve1[jj];
		}
	}
}
template <int ndimn>
inline void matXvec(const double*__restrict__ mat,
                    const double*__restrict__ ve1,
                          double*__restrict__ ve2){
  for(int ii = 0; ii < ndimn; ii++){
    ve2[ii] = mat[ndimn*ii+0]*ve1[0];
    for(int jj = 1; jj < ndimn; jj++){
      ve2[ii] += mat[ndimn*ii+jj]*ve1[jj];
    }
  }
}
template <int ndimn>
inline void tmatXvec(const double*__restrict__ mat,
                     const double*__restrict__ ve1,
                           double*__restrict__ ve2){
  for(int ii = 0; ii < ndimn; ii++){
    ve2[ii] = mat[ndimn*0+ii]*ve1[0];
    for(int jj = 1; jj < ndimn; jj++){
      ve2[ii] += mat[ndimn*jj+ii]*ve1[jj];
    }
  }
}

template<int ndimn, typename T1 = double, typename T2 = double, typename T3 = double>
void matXsym(const T1* __restrict__ mat, 
             const T2* __restrict__ met, 
             T3* __restrict__ out){

  static_assert(ndimn==2 || ndimn==3);

  for(int ii = 0; ii < ndimn; ii++){
    for(int jj = 0; jj < ndimn; jj++){
      out[ii*ndimn + jj] = mat[ii*ndimn + 0] * met[sym3idx(0,jj)];
      for(int kk = 1; kk < ndimn; kk++){
        out[ii*ndimn + jj] += mat[ii*ndimn + kk] * met[sym3idx(kk,jj)];
      }
    }
  }

}


// -----------------------------------------------------------------------------
// Compute A^T M A with M sym
// out is symmetric stored column first, like sym: 1 2 4  0 1 3
//                                                   3 5    2 4 
//                                                     6      5
template <typename T1, typename T2, typename T3>
inline void tAMA3(const T1* sym, const T2* mat, T3* out){
	out[0] =   mat[3*0+0]*mat[3*0+0]*sym[0] 
	       +   mat[3*1+0]*mat[3*1+0]*sym[2]
	       +   mat[3*2+0]*mat[3*2+0]*sym[5]
	       + 2*mat[3*0+0]*mat[3*1+0]*sym[1]
	       + 2*mat[3*0+0]*mat[3*2+0]*sym[3]
	       + 2*mat[3*1+0]*mat[3*2+0]*sym[4];

	out[1] =   mat[3*0+0]*mat[3*0+1]*sym[0] 
	       +   mat[3*1+0]*mat[3*1+1]*sym[2]
	       +   mat[3*2+0]*mat[3*2+1]*sym[5]
	       +   mat[3*0+0]*mat[3*1+1]*sym[1]
	       +   mat[3*1+0]*mat[3*0+1]*sym[1]
	       +   mat[3*0+0]*mat[3*2+1]*sym[3]
	       +   mat[3*2+0]*mat[3*0+1]*sym[3]
	       +   mat[3*1+0]*mat[3*2+1]*sym[4]
	       +   mat[3*2+0]*mat[3*1+1]*sym[4];

	out[2] =   mat[3*0+1]*mat[3*0+1]*sym[0] 
	       +   mat[3*1+1]*mat[3*1+1]*sym[2]
	       +   mat[3*2+1]*mat[3*2+1]*sym[5]
	       + 2*mat[3*0+1]*mat[3*1+1]*sym[1]
	       + 2*mat[3*0+1]*mat[3*2+1]*sym[3]
	       + 2*mat[3*1+1]*mat[3*2+1]*sym[4];

	out[3] =   mat[3*0+0]*mat[3*0+2]*sym[0] 
	       +   mat[3*1+0]*mat[3*1+2]*sym[2]
	       +   mat[3*2+0]*mat[3*2+2]*sym[5]
	       +   mat[3*0+0]*mat[3*1+2]*sym[1]
	       +   mat[3*1+0]*mat[3*0+2]*sym[1]
	       +   mat[3*0+0]*mat[3*2+2]*sym[3]
	       +   mat[3*2+0]*mat[3*0+2]*sym[3]
	       +   mat[3*1+0]*mat[3*2+2]*sym[4]
	       +   mat[3*2+0]*mat[3*1+2]*sym[4];

	out[4] =   mat[3*0+1]*mat[3*0+2]*sym[0] 
	       +   mat[3*1+1]*mat[3*1+2]*sym[2]
	       +   mat[3*2+1]*mat[3*2+2]*sym[5]
	       +   mat[3*0+1]*mat[3*1+2]*sym[1]
	       +   mat[3*1+1]*mat[3*0+2]*sym[1]
	       +   mat[3*0+1]*mat[3*2+2]*sym[3]
	       +   mat[3*2+1]*mat[3*0+2]*sym[3]
	       +   mat[3*1+1]*mat[3*2+2]*sym[4]
	       +   mat[3*2+1]*mat[3*1+2]*sym[4];

	out[5] =   mat[3*0+2]*mat[3*0+2]*sym[0] 
	       +   mat[3*1+2]*mat[3*1+2]*sym[2]
	       +   mat[3*2+2]*mat[3*2+2]*sym[5]
	       + 2*mat[3*0+2]*mat[3*1+2]*sym[1]
	       + 2*mat[3*0+2]*mat[3*2+2]*sym[3]
	       + 2*mat[3*1+2]*mat[3*2+2]*sym[4];
}

// -----------------------------------------------------------------------------
// Compute A M A^T with M sym
// out is symmetric stored column first, like sym: 1 2 4  0 1 3
//                                                   3 5    2 4 
//                                                     6      5
// A is size n1 x n2
// M is size n2 x n2 
// out is size n1 x n1 
template <int n1, int n2, typename T1, typename T2, typename T3>
inline void matXsymXtmat(const T1* sym, const T2* mat, T3* out){

  static_assert(n1 == 2 || n1 == 3);
	static_assert(n2 == 2 || n2 == 3);

	if constexpr(n1 == n2 && n1 == 2){
    constexpr int ndimn = n1;
		out[0] =   mat[ndimn*0+0]*mat[ndimn*0+0]*sym[0] 
		       +   mat[ndimn*0+1]*mat[ndimn*0+1]*sym[2]
		       + 2*mat[ndimn*0+0]*mat[ndimn*0+1]*sym[1];
	
		out[1] =   mat[ndimn*0+0]*mat[ndimn*1+0]*sym[0] 
		       +   mat[ndimn*0+1]*mat[ndimn*1+1]*sym[2]
		       +   (mat[ndimn*0+0]*mat[ndimn*1+1] + mat[ndimn*0+1]*mat[ndimn*1+0])*sym[1];
	
		out[2] =   mat[ndimn*1+0]*mat[ndimn*1+0]*sym[0] 
		       +   mat[ndimn*1+1]*mat[ndimn*1+1]*sym[2]
		       + 2*mat[ndimn*1+0]*mat[ndimn*1+1]*sym[1];

	}else if (n1 == n2 && n1 == 3){
    constexpr int ndimn = n1;
		out[0] =   mat[ndimn*0+0]*mat[ndimn*0+0]*sym[0] 
		       +   mat[ndimn*0+1]*mat[ndimn*0+1]*sym[2]
		       +   mat[ndimn*0+2]*mat[ndimn*0+2]*sym[5]
		       + 2*mat[ndimn*0+0]*mat[ndimn*0+1]*sym[1]
		       + 2*mat[ndimn*0+0]*mat[ndimn*0+2]*sym[3]
		       + 2*mat[ndimn*0+1]*mat[ndimn*0+2]*sym[4];
	
		out[1] =   mat[ndimn*0+0]*mat[ndimn*1+0]*sym[0] 
		       +   mat[ndimn*0+1]*mat[ndimn*1+1]*sym[2]
		       +   mat[ndimn*0+2]*mat[ndimn*1+2]*sym[5]
		       +   (mat[ndimn*0+0]*mat[ndimn*1+1] + mat[ndimn*0+1]*mat[ndimn*1+0])*sym[1]
		       +   (mat[ndimn*0+0]*mat[ndimn*1+2] + mat[ndimn*0+2]*mat[ndimn*1+0])*sym[3]
		       +   (mat[ndimn*0+1]*mat[ndimn*1+2] + mat[ndimn*0+2]*mat[ndimn*1+1])*sym[4];
	
		out[2] =   mat[ndimn*1+0]*mat[ndimn*1+0]*sym[0] 
		       +   mat[ndimn*1+1]*mat[ndimn*1+1]*sym[2]
		       +   mat[ndimn*1+2]*mat[ndimn*1+2]*sym[5]
		       + 2*mat[ndimn*1+0]*mat[ndimn*1+1]*sym[1]
		       + 2*mat[ndimn*1+0]*mat[ndimn*1+2]*sym[3]
		       + 2*mat[ndimn*1+1]*mat[ndimn*1+2]*sym[4];
	
		out[3] =   mat[ndimn*0+0]*mat[ndimn*2+0]*sym[0] 
		       +   mat[ndimn*0+1]*mat[ndimn*2+1]*sym[2]
		       +   mat[ndimn*0+2]*mat[ndimn*2+2]*sym[5]
		       +  (mat[ndimn*0+0]*mat[ndimn*2+1] + mat[ndimn*0+1]*mat[ndimn*2+0])*sym[1]
		       +  (mat[ndimn*0+0]*mat[ndimn*2+2] + mat[ndimn*0+2]*mat[ndimn*2+0])*sym[3]
		       +  (mat[ndimn*0+1]*mat[ndimn*2+2] + mat[ndimn*0+2]*mat[ndimn*2+1])*sym[4];
	
		out[4] =   mat[ndimn*1+0]*mat[ndimn*2+0]*sym[0] 
		       +   mat[ndimn*1+1]*mat[ndimn*2+1]*sym[2]
		       +   mat[ndimn*1+2]*mat[ndimn*2+2]*sym[5]
		       +  (mat[ndimn*1+0]*mat[ndimn*2+1] + mat[ndimn*1+1]*mat[ndimn*2+0])*sym[1]
		       +  (mat[ndimn*1+0]*mat[ndimn*2+2] + mat[ndimn*1+2]*mat[ndimn*2+0])*sym[3]
		       +  (mat[ndimn*1+1]*mat[ndimn*2+2] + mat[ndimn*1+2]*mat[ndimn*2+1])*sym[4];
	
		out[5] =   mat[ndimn*2+0]*mat[ndimn*2+0]*sym[0] 
		       +   mat[ndimn*2+1]*mat[ndimn*2+1]*sym[2]
		       +   mat[ndimn*2+2]*mat[ndimn*2+2]*sym[5]
		       + 2*mat[ndimn*2+0]*mat[ndimn*2+1]*sym[1]
		       + 2*mat[ndimn*2+0]*mat[ndimn*2+2]*sym[3]
		       + 2*mat[ndimn*2+1]*mat[ndimn*2+2]*sym[4];
	}else{
    // This if constexpr should not be necessary, but it is 
    if constexpr(n1 != n2)    static_assert(n1 == 2 && n2 == 3);

    out[0] =   mat[n2*0+0]*mat[n2*0+0]*sym[sym3idx(0,0)] 
           +   mat[n2*0+1]*mat[n2*0+1]*sym[sym3idx(1,1)]
           +   mat[n2*0+2]*mat[n2*0+2]*sym[sym3idx(2,2)]
           + 2*mat[n2*0+0]*mat[n2*0+1]*sym[sym3idx(0,1)]
           + 2*mat[n2*0+0]*mat[n2*0+2]*sym[sym3idx(0,2)]
           + 2*mat[n2*0+1]*mat[n2*0+2]*sym[sym3idx(1,2)];
  
    out[1] =   mat[n2*0+0]*mat[n2*1+0]*sym[sym3idx(0,0)] 
           +   mat[n2*0+1]*mat[n2*1+1]*sym[sym3idx(1,1)]
           +   mat[n2*0+2]*mat[n2*1+2]*sym[sym3idx(2,2)]
           +   (mat[n2*0+0]*mat[n2*1+1] + mat[n2*0+1]*mat[n2*1+0])*sym[sym3idx(0,1)]
           +   (mat[n2*0+0]*mat[n2*1+2] + mat[n2*0+2]*mat[n2*1+0])*sym[sym3idx(0,2)]
           +   (mat[n2*0+1]*mat[n2*1+2] + mat[n2*0+2]*mat[n2*1+1])*sym[sym3idx(1,2)];
  
    out[2] =   mat[n2*1+0]*mat[n2*1+0]*sym[sym3idx(0,0)] 
           +   mat[n2*1+1]*mat[n2*1+1]*sym[sym3idx(1,1)]
           +   mat[n2*1+2]*mat[n2*1+2]*sym[sym3idx(2,2)]
           + 2*mat[n2*1+0]*mat[n2*1+1]*sym[sym3idx(0,1)]
           + 2*mat[n2*1+0]*mat[n2*1+2]*sym[sym3idx(0,2)]
           + 2*mat[n2*1+1]*mat[n2*1+2]*sym[sym3idx(1,2)];

  }
}

// Compute the trace of A Sym A^T 
template <int ndimn, typename T1, typename T2, typename T3>
//inline void matXsymXtmat_diag(const T1* sym, const T2* mat, T3* out){
inline T3 tra_matXsymXtmat(const T1* sym, const T2* mat){
	static_assert(ndimn == 2 || ndimn == 3);
  T3 out; 

  #if 0
  // Previously computed the diagonal:
	if constexpr(ndimn == 2){
		out[0] =   mat[ndimn*0+0]*mat[ndimn*0+0]*sym[sym3idx(0,0)] 
		       +   mat[ndimn*0+1]*mat[ndimn*0+1]*sym[sym3idx(1,1)]
		       + 2*mat[ndimn*0+0]*mat[ndimn*0+1]*sym[sym3idx(0,1)];
	
		out[1] =   mat[ndimn*1+0]*mat[ndimn*1+0]*sym[sym3idx(0,0)] 
		       +   mat[ndimn*1+1]*mat[ndimn*1+1]*sym[sym3idx(1,1)]
		       + 2*mat[ndimn*1+0]*mat[ndimn*1+1]*sym[sym3idx(0,1)];
	}else{
		out[0] =   mat[ndimn*0+0]*mat[ndimn*0+0]*sym[sym3idx(0,0)] 
		       +   mat[ndimn*0+1]*mat[ndimn*0+1]*sym[sym3idx(1,1)]
		       +   mat[ndimn*0+2]*mat[ndimn*0+2]*sym[sym3idx(2,2)]
		       + 2*mat[ndimn*0+0]*mat[ndimn*0+1]*sym[sym3idx(0,1)]
		       + 2*mat[ndimn*0+0]*mat[ndimn*0+2]*sym[sym3idx(0,2)]
		       + 2*mat[ndimn*0+1]*mat[ndimn*0+2]*sym[sym3idx(1,2)];
	
		out[1] =   mat[ndimn*1+0]*mat[ndimn*1+0]*sym[sym3idx(0,0)] 
		       +   mat[ndimn*1+1]*mat[ndimn*1+1]*sym[sym3idx(1,1)]
		       +   mat[ndimn*1+2]*mat[ndimn*1+2]*sym[sym3idx(2,2)]
		       + 2*mat[ndimn*1+0]*mat[ndimn*1+1]*sym[sym3idx(0,1)]
		       + 2*mat[ndimn*1+0]*mat[ndimn*1+2]*sym[sym3idx(0,2)]
		       + 2*mat[ndimn*1+1]*mat[ndimn*1+2]*sym[sym3idx(1,2)];
	
		out[2] =   mat[ndimn*2+0]*mat[ndimn*2+0]*sym[sym3idx(0,0)] 
		       +   mat[ndimn*2+1]*mat[ndimn*2+1]*sym[sym3idx(1,1)]
		       +   mat[ndimn*2+2]*mat[ndimn*2+2]*sym[sym3idx(2,2)]
		       + 2*mat[ndimn*2+0]*mat[ndimn*2+1]*sym[sym3idx(0,1)]
		       + 2*mat[ndimn*2+0]*mat[ndimn*2+2]*sym[sym3idx(0,2)]
		       + 2*mat[ndimn*2+1]*mat[ndimn*2+2]*sym[sym3idx(1,2)];
	}
  #endif

  if constexpr(ndimn == 2){
    out =   (mat[ndimn*0+0]*mat[ndimn*0+0] + mat[ndimn*1+0]*mat[ndimn*1+0])*sym[sym3idx(0,0)] 
        +   (mat[ndimn*0+1]*mat[ndimn*0+1] + mat[ndimn*1+1]*mat[ndimn*1+1])*sym[sym3idx(1,1)]
        + 2*(mat[ndimn*0+0]*mat[ndimn*0+1] + mat[ndimn*1+0]*mat[ndimn*1+1])*sym[sym3idx(0,1)];
  }else{
    out =   (mat[ndimn*0+0]*mat[ndimn*0+0] 
           + mat[ndimn*1+0]*mat[ndimn*1+0] 
           + mat[ndimn*2+0]*mat[ndimn*2+0])*sym[sym3idx(0,0)] 
        +   (mat[ndimn*0+1]*mat[ndimn*0+1] 
           + mat[ndimn*1+1]*mat[ndimn*1+1] 
           + mat[ndimn*2+1]*mat[ndimn*2+1])*sym[sym3idx(1,1)]
        +   (mat[ndimn*0+2]*mat[ndimn*0+2] 
           + mat[ndimn*1+2]*mat[ndimn*1+2] 
           + mat[ndimn*2+2]*mat[ndimn*2+2])*sym[sym3idx(2,2)]
        + 2*(mat[ndimn*0+0]*mat[ndimn*0+1] 
           + mat[ndimn*1+0]*mat[ndimn*1+1] 
           + mat[ndimn*2+0]*mat[ndimn*2+1])*sym[sym3idx(0,1)]
        + 2*(mat[ndimn*0+0]*mat[ndimn*0+2] 
           + mat[ndimn*1+0]*mat[ndimn*1+2] 
           + mat[ndimn*2+0]*mat[ndimn*2+2])*sym[sym3idx(0,2)]
        + 2*(mat[ndimn*0+1]*mat[ndimn*0+2] 
           + mat[ndimn*1+1]*mat[ndimn*1+2] 
           + mat[ndimn*2+1]*mat[ndimn*2+2])*sym[sym3idx(1,2)];
  }
  return out;
}

// sym is size nsym x nsym
// mat is size nlin x nsym
// out is naturally nlin x nlin, and is a symmetric matrix ; 
// note: we only compute the diagonal here !
template <int nsym, int nlin, typename T1, typename T2, typename T3>
inline void matXsymXtmat_diag(const T1* sym, const T2* mat, T3* out){
  static_assert(nsym > 0 && nlin > 0);
  static_assert(nsym <= 3 && nlin <= 3); // Not so much that it wouldnt work, but why though? Probably a bug in the caller. 
  for(int ii = 0 ; ii < nlin; ii++){
    out[ii] = mat[ii*nsym + 0]*mat[ii*nsym + 0]*sym[sym3idx(0,0)]; 
    for(int jj = 1; jj < nsym; jj++){
      out[ii] += mat[ii*nsym + jj]*mat[ii*nsym + jj]*sym[sym3idx(jj,jj)];
      for(int kk = 0; kk < jj; kk++){
        out[ii] += 2*mat[ii*nsym + jj]*mat[ii*nsym + kk]*sym[sym3idx(jj,kk)];
      }
    }
  }
}






/*
tens3sym: for metric derivatives, morally a [i][j][k] with each restriction to i[:][:] symmetric. 
mat3sym: as usual, a symmetric matrix.
Symmetric matrix indices: 0 1 3
                            2 4
                              5
Define the products: 

T X1 M = sum_s T_{sjk} M_{si}     | tens3sym1X1mat3sym
T X2 M = sum_s T_{isk} M_{sj}     | tens3symX2mat3sym
T X3 M = sum_s T_{ijs} M_{sk}     | tens3symX3mat3sym

If T is symmetric outside of i, the resulting tensor (by Xi) is symmetric outside of i. 
We only consider symmetric outisde of 1 tensors. 
*/

/* Against symmetric matrices */

// This returns a symmetric tensor (outside of 1)
template<int ndimn>
void tens3sym1X1mat3sym(const double* __restrict__ tens3sym1, 
                        const double* __restrict__ mat3sym, 
                        double* __restrict__ tens3symo);

// Compute A x_2 T x_1 A^T with T 1-symmetryc
// This is not a symmetric tensor
void mat3X2tens3sym1X1tmat3(const double* __restrict__ tens3sym1, 
                           const double* __restrict__ mat3,
                           double* __restrict__ tens3);

// Input 3-symmetric, output not sym
void mat3X1tens3sym3(const double* __restrict__ tens3sym, 
                     const double* __restrict__ mat3,
                     double* __restrict__ tens3);


// Input 1-symmetric, output also 1-symmetric
void mat3X1tens3sym1(const double* __restrict__ tens3sym, 
                     const double* __restrict__ mat3,
                     double* __restrict__ tens3symo);


} // End namespace

#endif