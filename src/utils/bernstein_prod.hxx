//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_AUX_BERNSTEIN__
#define __METRIS_AUX_BERNSTEIN__

#include "../types.hxx"
#include "../types.hxx"
#include "aux_misc.hxx"
#include "../fact_mulidx.hxx"
#include "../utils/binom.hxx"
#include "../linalg/symidx.hxx"
#include "../utils/CT_loop.hxx"


//#define USE_BINOM
//#define CT_LOOPS

namespace Metris{

// Compute the square of a polynomial written in Bernstein basis 
// i.e. get its coeffs in the deg 2ideg basis
// Derivatives (with gdim entries) can be provided, up to second, returning as many. 
template<int gdim, int tdim, int ideg>
void square_bernstein(const dblAr2 &coef_inp, dblAr2 &coef_out, 
                      std::initializer_list<const dblAr2*> dcoef_inp = {}, 
                      std::initializer_list<      dblAr2*> dcoef_out = {});


template<int gdim, int tdim, int ideg1, int ideg2>
void prod_bernstein(const dblAr2 &coef_inp1, const dblAr2 &coef_inp2, dblAr2 &coef_out,
                    std::initializer_list<const dblAr2*> dcoef_inp1 = {}, 
                    std::initializer_list<const dblAr2*> dcoef_inp2 = {}, 
                    std::initializer_list<      dblAr2*> dcoef_out  = {}){
  METRIS_ASSERT(  dcoef_inp1.size() == dcoef_out.size() 
               && dcoef_inp2.size() == dcoef_out.size() );

  //static bool iinit = false;
  //if(!iinit){
  //  printf("instantiate prod %d %d %d %d \n",gdim,tdim,ideg1,ideg2);
  //  iinit = true;
  //}

  int ndiff = dcoef_inp1.size();

  constexpr int nhess = (gdim*(gdim+1))/2;
  constexpr int nnod1 = getnnode(tdim,ideg1);
  constexpr int nnod2 = getnnode(tdim,ideg2);
  constexpr int nnods = getnnode(tdim,ideg1+ideg2);

  coef_out.set_n(nnods);
  coef_out.fill(0);


  // Note symmetry w replacing ideg1 by ideg2
  constexpr int binom1 = ct_binom.get[ideg1+ideg2][ideg1];
  METRIS_ENFORCE(binom1 == binom(MIN(ideg1,ideg2),ideg1+ideg2))
  //int binom1 = binom(MIN(ideg1,ideg2),ideg1+ideg2); 

  constexpr auto ordent = ORDELT(tdim);


  const dblAr2 *d1coef_inp1,*d2coef_inp1;
  const dblAr2 *d1coef_inp2,*d2coef_inp2;
        dblAr2 *d1coef_out,*d2coef_out;
  if(ndiff > 0){
    auto iter_inp1 = dcoef_inp1.begin();
    auto iter_inp2 = dcoef_inp2.begin();
    auto iter_out = dcoef_out.begin();
    d1coef_inp1 = *iter_inp1;
    d1coef_inp2 = *iter_inp2;
    d1coef_out = *iter_out;
    d1coef_out->fill(0);

    if(ndiff > 1){
      iter_inp1++;
      iter_inp2++;
      iter_out++;
      d2coef_inp1 = *iter_inp1;
      d2coef_inp2 = *iter_inp2;
      d2coef_out = *iter_out;
      d2coef_out->fill(0);
    }
  }

  
  #ifdef CT_LOOPS
  CT_FOR0_EXC(0,nnod1,inod1){
  #else
  for(int inod1 = 0; inod1 < nnod1; inod1++){
  #endif
    #ifndef USE_BINOM
    double fact_idx1 = (double) fact_mulidx::get(tdim,ideg1,inod1);
    #endif

    #ifdef CT_LOOPS
    CT_FOR0_EXC(0,nnod2,inod2){
    #else
    for(int inod2 = 0; inod2 < nnod2; inod2++){
    #endif
      #ifndef USE_BINOM
      double fact_idx2 = (double) fact_mulidx::get(tdim,ideg2,inod2);
      #endif

      int inods = tdim == 1 ? mul2nod(ordent[ideg1][inod1][0] + ordent[ideg2][inod2][0],
                                      ordent[ideg1][inod1][1] + ordent[ideg2][inod2][1])
                : tdim == 2 ? mul2nod(ordent[ideg1][inod1][0] + ordent[ideg2][inod2][0],
                                      ordent[ideg1][inod1][1] + ordent[ideg2][inod2][1],
                                      ordent[ideg1][inod1][2] + ordent[ideg2][inod2][2])
                            : mul2nod(ordent[ideg1][inod1][0] + ordent[ideg2][inod2][0],
                                      ordent[ideg1][inod1][1] + ordent[ideg2][inod2][1],
                                      ordent[ideg1][inod1][2] + ordent[ideg2][inod2][2],
                                      ordent[ideg1][inod1][3] + ordent[ideg2][inod2][3]); 

      #ifndef USE_BINOM
        double fact_idxs = (double) fact_mulidx::get(tdim,ideg1+ideg2,inods);
        double mulfac = fact_idxs/(fact_idx1*fact_idx2*binom1);
      #else
        #ifndef CT_LOOPS
        int binom2 = binom(ordent[ideg1][inod1][0], ordent[ideg2][inod2][0] + ordent[ideg1][inod1][0])
                   * binom(ordent[ideg1][inod1][1], ordent[ideg2][inod2][1] + ordent[ideg1][inod1][1])
                   * binom(ordent[ideg1][inod1][2], ordent[ideg2][inod2][2] + ordent[ideg1][inod1][2]);
        #else
        constexpr int binom2 = ct_binom.get[ordent[ideg1][inod1][0]+ordent[ideg2][inod2][0]][ordent[ideg1][inod1][0]]
                             * ct_binom.get[ordent[ideg1][inod1][1]+ordent[ideg2][inod2][1]][ordent[ideg1][inod1][1]]
                             * ct_binom.get[ordent[ideg1][inod1][2]+ordent[ideg2][inod2][2]][ordent[ideg1][inod1][2]];
        METRIS_ASSERT(binom2 == binom(ordent[ideg1][inod1][0], ordent[ideg2][inod2][0] + ordent[ideg1][inod1][0])
                   * binom(ordent[ideg1][inod1][1], ordent[ideg2][inod2][1] + ordent[ideg1][inod1][1])
                   * binom(ordent[ideg1][inod1][2], ordent[ideg2][inod2][2] + ordent[ideg1][inod1][2]));
        #endif
        double mulfac = binom2/(double)binom1;
      #endif

      coef_out(inods,0) += coef_inp1(inod1,0)*coef_inp2(inod2,0)*mulfac;


      #ifndef CT_LOOPS
      if(ndiff == 0) continue;
      #else
      if(ndiff == 0) CT_CONTINUE(inod2);
      #endif

      for(int ii = 0; ii < gdim; ii++)
        (*d1coef_out)(inods,ii) += (   coef_inp1(inod1,0 )*(*d1coef_inp2)(inod2,ii)
                                   + (*d1coef_inp1)(inod1,ii)*  coef_inp2(inod2,0 ) )
                                   *mulfac;

      #ifndef CT_LOOPS
      if(ndiff == 1) continue;
      #else
      if(ndiff == 1) CT_CONTINUE(inod2);
      #endif

      for(int ii = 0; ii < gdim; ii++)
        for(int jj = ii; jj < gdim; jj++)
          (*d2coef_out)(inods,sym2idx(ii,jj)) += ( (*d1coef_inp1)(inod1,ii)*(*d1coef_inp2)(inod2,jj)
                                                 + (*d1coef_inp1)(inod1,jj)*(*d1coef_inp2)(inod2,ii)
                                                 + (*d2coef_inp1)(inod1,sym2idx(ii,jj))
                                                    *coef_inp2(inod2,0)
                                                 + (*d2coef_inp2)(inod2,sym2idx(ii,jj))
                                                    *coef_inp1(inod1,0) )
                                                 *mulfac;

  #ifndef CT_LOOPS
    }
  }
  #else
    }CT_FOR1(inod2);
  }CT_FOR1(inod1);
  #endif

  return;
}


} // namespace
#endif