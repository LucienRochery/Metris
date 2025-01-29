//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../linalg/matprods.hxx"
#include <cstdio>

namespace Metris{


//template<int ndimn, typename T1 = double, typename T2 = double, typename T3 = double>
//void matXsym(const double* __restrict__ mat, 
//             const double* __restrict__ met, 
//             double* __restrict__ out){
//
//  static_assert(ndimn==2 || ndimn==3);
//
//  for(int ii = 0; ii < ndimn; ii++){
//    for(int jj = 0; jj < ndimn; jj++){
//      out[ii*ndimn + jj] = mat[ii*ndimn + 0] * met[sym2idx(0,jj)];
//      for(int kk = 1; kk < ndimn; kk++){
//        out[ii*ndimn + jj] += mat[ii*ndimn + kk] * met[sym2idx(kk,jj)];
//      }
//    }
//  }
//
//}
//template void matXsym<2>(const double* __restrict__ mat, const double* __restrict__ met, double* __restrict__ out);
//template void matXsym<3>(const double* __restrict__ mat, const double* __restrict__ met, double* __restrict__ out);



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
void tens3sym1X1mat3sym(const double* __restrict__ tenssym1, 
                        const double* __restrict__ matsym, 
                              double* __restrict__ tenssym2){
  static_assert(ndimn == 2 || ndimn == 3);
  constexpr int nnmet = (ndimn*(ndimn+1))/2;

  for(int ii = 0; ii < ndimn; ii++){
    for(int jj = 0 ; jj < ndimn; jj++){
      for(int kk = jj; kk < ndimn ; kk++){
        tenssym2[nnmet*ii + sym2idx(jj,kk)] = tenssym1[nnmet*0 + sym2idx(jj,kk)]*matsym[sym2idx(0,ii)];
        for(int s = 1; s < ndimn; s++){
          tenssym2[nnmet*ii + sym2idx(jj,kk)] += tenssym1[nnmet*s + sym2idx(jj,kk)]*matsym[sym2idx(s,ii)];
        }
      }
    }
  }
}

template void tens3sym1X1mat3sym<2>(const double* __restrict__ tens3sym1, 
                                    const double* __restrict__ mat3sym, 
                                          double* __restrict__ tens3symo);
template void tens3sym1X1mat3sym<3>(const double* __restrict__ tens3sym1, 
                                    const double* __restrict__ mat3sym, 
                                          double* __restrict__ tens3symo);



// Compute A x_2 T x_1 A^T with T 1-symmetric (but A not)
// This is not a symmetric tensor 
void mat3X2tens3sym1X1tmat3(const double* __restrict__ tens3sym1, 
                           const double* __restrict__ mat3,
                                 double* __restrict__ tens3){
  printf("dbg %f \n",tens3sym1[0]);
  printf("dbg %f \n",tens3sym1[1]);
  printf("dbg %f \n",tens3sym1[2]);
  printf("dbg %f \n",tens3sym1[3]);
  printf("dbg %f \n",tens3sym1[4]);
  printf("dbg %f \n",tens3sym1[5]);
  printf("dbg %f \n",tens3sym1[6]);
  printf("dbg %f \n",tens3sym1[7]);
  printf("dbg %f \n",tens3sym1[8]);
  printf("dbg %f \n",tens3sym1[9]);
  printf("dbg %f \n",tens3sym1[10]);
  printf("dbg %f \n",tens3sym1[11]);
  printf("dbg %f \n",tens3sym1[12]);
  printf("dbg %f \n",tens3sym1[13]);
  printf("dbg %f \n",tens3sym1[14]);
  printf("dbg %f \n",tens3sym1[15]);
  printf("dbg %f \n",tens3sym1[16]);
  printf("dbg %f \n",tens3sym1[17]);

  printf("dbg %f \n",mat3[0]);
  printf("dbg %f \n",mat3[1]);
  printf("dbg %f \n",mat3[2]);
  printf("dbg %f \n",mat3[3]);
  printf("dbg %f \n",mat3[4]);
  printf("dbg %f \n",mat3[5]);
  printf("dbg %f \n",mat3[6]);
  printf("dbg %f \n",mat3[7]);
  printf("dbg %f \n",mat3[8]);

  
  for(int ii = 0; ii < 3; ii++){
    for(int jj = 0; jj < 3; jj++){
      for(int kk = 0; kk < 3; kk++){
        tens3[9*ii + 3*jj + kk] = 0;
        for(int ss = 0; ss < 3; ss++){
          for(int tt = 0; tt < 3; tt++){  
            tens3[9*ii + 3*jj + kk] += 
                                 tens3sym1[6*ss + sym2idx(tt,kk)]*mat3[3*jj + ss]*mat3[3*ii + tt];
          }
        }
      }
    }
  }
}

// T is 3-sym (w.r.t. first 2 indices)
// Output is nothing at all 
void mat3X1tens3sym3(const double* __restrict__ tens3sym, 
                     const double* __restrict__ mat3 ,
                           double* __restrict__ tens3){
  for(int ii = 0; ii < 3; ii++){
    for(int jj = 0 ; jj < 3; jj++){
      for(int kk = 0; kk < 3 ; kk++){
        tens3[9*ii + 3*jj + kk] = tens3sym[3*sym2idx(0,jj) + kk] * mat3[3*ii + 0];
        for(int s = 1; s < 3; s++){
          tens3[9*ii + 3*jj + kk] += tens3sym[3*sym2idx(s,jj) + kk] * mat3[3*ii + s];
        }
      }
    }
  }
}


// T is 3-sym (w.r.t. first 2 indices)
// Output is nothing at all 
void mat3X1tens3sym1(const double* __restrict__ tens3sym, 
                     const double* __restrict__ mat3 ,
                           double* __restrict__ tens3symo){
  for(int ii = 0; ii < 3; ii++){
    for(int jj = 0 ; jj < 3; jj++){
      for(int kk = jj; kk < 3 ; kk++){
        tens3symo[6*ii + sym2idx(jj,kk)] = tens3sym[6*0 + sym2idx(jj,kk)] * mat3[3*ii + 0];
        for(int s = 1; s < 3; s++){
          tens3symo[6*ii + sym2idx(jj,kk)] += tens3sym[6*s + sym2idx(jj,kk)] * mat3[3*ii + s];
        }
      }
    }
  }
}


} // End namespace
