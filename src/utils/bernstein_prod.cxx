//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "bernstein_prod.hxx"

#include "../types.hxx"
#include "aux_misc.hxx"
#include "../fact_mulidx.hxx"
#include "../utils/binom.hxx"
#include "../linalg/symidx.hxx"
#include "../utils/CT_loop.hxx"



namespace Metris{

// Compute the square of a polynomial written in Bernstein basis 
// i.e. get its coeffs in the deg 2ideg basis
// Derivatives (with gdim entries) can be provided, up to second, returning as many. 
template<int gdim, int tdim, int ideg>
void square_bernstein(const dblAr2 &coef_inp, dblAr2 &coef_out, 
                      std::initializer_list<const dblAr2*> dcoef_inp, 
                      std::initializer_list<      dblAr2*> dcoef_out){
  METRIS_ASSERT(dcoef_inp.size() == dcoef_out.size());

  //static bool iinit = false;
  //if(!iinit){
  //  printf("instantiate sqr %d %d %d\n",gdim,tdim,ideg);
  //  iinit = true;
  //}

  int ndiff = dcoef_inp.size();

  constexpr int nhess = (gdim*(gdim+1))/2;
  constexpr int nnode = getnnode(tdim,  ideg);
  constexpr int nnod2 = getnnode(tdim,2*ideg);
  coef_out.set_n(nnod2);
  coef_out.fill(0);

  constexpr int fact_d = ifact<ideg>();
  constexpr int fact_2d = ifact<2*ideg>();

  constexpr auto ordent = ORDELT(tdim);

  //printf("# debug start tdim ideg %d %d \n",tdim,ideg);
  //printf(" inp = \n");
  //coef_inp.print();

  const dblAr2 *d1coef_inp,*d2coef_inp;
        dblAr2 *d1coef_out,*d2coef_out;
  if(ndiff > 0){
    auto iter_inp = dcoef_inp.begin();
    auto iter_out = dcoef_out.begin();
    d1coef_inp = *iter_inp;
    d1coef_out = *iter_out;
    d1coef_out->fill(0);

    if(ndiff > 1){
      iter_inp++;
      iter_out++;
      d2coef_inp = *iter_inp;
      d2coef_out = *iter_out;
      d2coef_out->fill(0);
    }
  }

  #ifdef CT_LOOPS
  CT_FOR0_EXC(0,nnode,inod1){
  #else
  for(int inod1 = 0; inod1 < nnode; inod1++){
  #endif
    #ifndef USE_BINOM
    #ifndef CT_LOOPS
    double fact_idx1 = (double) fact_mulidx::get(tdim,ideg,inod1);
    #else
    constexpr INT_FACT fact_idx1 = ifact<ordent[ideg][inod1][0]>()
                                 * ifact<ordent[ideg][inod1][1]>()
                                 * ifact<ordent[ideg][inod1][2]>();
    #endif
    #endif

    #ifdef CT_LOOPS
    CT_FOR0_EXC(inod1,nnode,inod2){
    #else
    for(int inod2 = inod1; inod2 < nnode; inod2++){
    #endif
    //for(int inod2 = inod1; inod2 < nnode; inod2++){
      #ifndef USE_BINOM
      #ifndef CT_LOOPS
      double fact_idx2 = (double) fact_mulidx::get(tdim,ideg,inod2);
      #else
      constexpr INT_FACT fact_idx2 = ifact<ordent[ideg][inod2][0]>()
                              * ifact<ordent[ideg][inod2][1]>()
                              * ifact<ordent[ideg][inod2][2]>();
      #endif
      #endif

      #ifdef CT_LOOPS
      constexpr
      #endif
      int inods = tdim == 1 ? mul2nod(ordent[ideg][inod1][0] + ordent[ideg][inod2][0],
                                      ordent[ideg][inod1][1] + ordent[ideg][inod2][1])
                : tdim == 2 ? mul2nod(ordent[ideg][inod1][0] + ordent[ideg][inod2][0],
                                      ordent[ideg][inod1][1] + ordent[ideg][inod2][1],
                                      ordent[ideg][inod1][2] + ordent[ideg][inod2][2])
                            : mul2nod(ordent[ideg][inod1][0] + ordent[ideg][inod2][0],
                                      ordent[ideg][inod1][1] + ordent[ideg][inod2][1],
                                      ordent[ideg][inod1][2] + ordent[ideg][inod2][2],
                                      ordent[ideg][inod1][3] + ordent[ideg][inod2][3]); 
      #ifndef USE_BINOM
        #ifndef CT_LOOPS
        double fact_idxs = (double) fact_mulidx::get(tdim,2*ideg,inods);
        #else
        constexpr INT_FACT fact_idxs = ifact<ordent[2*ideg][inods][0]>()
                                * ifact<ordent[2*ideg][inods][1]>()
                                * ifact<ordent[2*ideg][inods][2]>();
        #endif
        double mulfac = (1 + (inod2 != inod1))*(fact_d*fact_d*fact_idxs)
                                              /(double)(fact_idx1*fact_idx2*fact_2d);
      #else
        #ifndef CT_LOOPS
        int binom2 = binom(ordent[ideg][inod1][0], ordent[ideg][inod2][0] + ordent[ideg][inod1][0])
                   * binom(ordent[ideg][inod1][1], ordent[ideg][inod2][1] + ordent[ideg][inod1][1])
                   * binom(ordent[ideg][inod1][2], ordent[ideg][inod2][2] + ordent[ideg][inod1][2]);
        #else
        constexpr
        int binom2 = tdim == 1 ? ct_binom.get[ordent[ideg][inod1][0]+ordent[ideg][inod2][0]][ordent[ideg][inod1][0]]
                               * ct_binom.get[ordent[ideg][inod1][1]+ordent[ideg][inod2][1]][ordent[ideg][inod1][1]]
                   : tdim == 2 ? ct_binom.get[ordent[ideg][inod1][0]+ordent[ideg][inod2][0]][ordent[ideg][inod1][0]]
                               * ct_binom.get[ordent[ideg][inod1][1]+ordent[ideg][inod2][1]][ordent[ideg][inod1][1]]
                               * ct_binom.get[ordent[ideg][inod1][2]+ordent[ideg][inod2][2]][ordent[ideg][inod1][2]]
                   :             ct_binom.get[ordent[ideg][inod1][0]+ordent[ideg][inod2][0]][ordent[ideg][inod1][0]]
                               * ct_binom.get[ordent[ideg][inod1][1]+ordent[ideg][inod2][1]][ordent[ideg][inod1][1]]
                               * ct_binom.get[ordent[ideg][inod1][2]+ordent[ideg][inod2][2]][ordent[ideg][inod1][2]]
                               * ct_binom.get[ordent[ideg][inod1][3]+ordent[ideg][inod2][3]][ordent[ideg][inod1][3]];

        #endif
      double mulfac = (1 + (inod2 != inod1))*(fact_d*fact_d)*binom2/((double)fact_2d);
      #endif
      coef_out(inods,0) += coef_inp(inod1,0)*coef_inp(inod2,0)
                           *mulfac;

      #ifndef CT_LOOPS
      if(ndiff == 0) continue;
      #else
      if(ndiff == 0) CT_CONTINUE(inod2);
      #endif

      for(int ii = 0; ii < gdim; ii++)
        (*d1coef_out)(inods,ii) += (   coef_inp(inod1,0 )*(*d1coef_inp)(inod2,ii)
                                   + (*d1coef_inp)(inod1,ii)*  coef_inp(inod2,0 ) )
                                   *mulfac;

      #ifndef CT_LOOPS
      if(ndiff == 1) continue;
      #else
      if(ndiff == 1) CT_CONTINUE(inod2);
      #endif

      for(int ii = 0; ii < gdim; ii++)
        for(int jj = ii; jj < gdim; jj++)
          (*d2coef_out)(inods,sym2idx(ii,jj)) += ( (*d1coef_inp)(inod1,ii)*(*d1coef_inp)(inod2,jj)
                                               + (*d1coef_inp)(inod1,jj)*(*d1coef_inp)(inod2,ii)
                                               + (*d2coef_inp)(inod1,sym2idx(ii,jj))
                                                  *coef_inp(inod2,0)
                                               + (*d2coef_inp)(inod2,sym2idx(ii,jj))
                                                  *coef_inp(inod1,0) )
                                               *mulfac;
  #ifndef CT_LOOPS
    }
  }
  #else
    }CT_FOR1(inod2);
  }CT_FOR1(inod1);
  #endif

  //printf(" out = \n");
  //coef_out.print();
  return;
}


#define BOOST_PP_LOCAL_MACRO(n) \
template \
void square_bernstein<1,1,n>(const dblAr2 &coef_inp, dblAr2 &coef_out,\
                      std::initializer_list<const dblAr2*> dcoef_inp = {},\
                      std::initializer_list<      dblAr2*> dcoef_out = {});\
template \
void square_bernstein<2,2,n>(const dblAr2 &coef_inp, dblAr2 &coef_out,\
                      std::initializer_list<const dblAr2*> dcoef_inp = {},\
                      std::initializer_list<      dblAr2*> dcoef_out = {});\
template \
void square_bernstein<3,3,n>(const dblAr2 &coef_inp, dblAr2 &coef_out,\
                      std::initializer_list<const dblAr2*> dcoef_inp = {},\
                      std::initializer_list<      dblAr2*> dcoef_out = {});
#define BOOST_PP_LOCAL_LIMITS     (1, 3*METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

} // namespace
