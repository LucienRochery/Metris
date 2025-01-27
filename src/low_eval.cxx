//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_eval.hxx"
#include <cmath>

#include "ho_constants.hxx"
#include "metris_constants.hxx"
#include "codegen_lagrange.hxx"

#include <boost/preprocessor/iteration/local.hpp>


namespace Metris{


/*
Szfld was added to handle both metrics and points. 
In practice, it's not so much a variable as a paremeter either 2,3 or 6 (also reals).
Templated it in order to improve loop unrolling/vectorization
*/
template <int szfld, int ideg>
void eval1(const dblAr2 & __restrict__  rfld,  
           const int * __restrict__  lfld,
           FEBasis ibasis, DifVar idif1, DifVar idif2,
           const double * __restrict__  bary, 
           double * __restrict__  eval, 
           double * __restrict__  jmat, 
           double * __restrict__  hmat){

  METRIS_ASSERT(idif2 == DifVar::None || idif2 == DifVar::Bary);
  METRIS_ASSERT(idif1 == DifVar::None || idif1 == DifVar::Bary);
  METRIS_ASSERT(ibasis == FEBasis::Bezier || ibasis == FEBasis::Lagrange);

  if constexpr (ideg == 1){
    for(int i=0;i<szfld;i++){
      eval[i] = bary[0]*rfld[lfld[0]][i]
              + bary[1]*rfld[lfld[1]][i];
    }
    if(idif1 == DifVar::Bary){
      for(int i=0;i<szfld;i++){
        jmat[0*szfld+i] = rfld[lfld[1]][i] - rfld[lfld[0]][i];
      }
    }
    if(idif2 == DifVar::Bary){
      for(int i=0;i<szfld;i++){
        hmat[0*szfld+i] = 0.0;
      }
    }
    return;
  } 

  if(ibasis == FEBasis::Bezier){
    eval1_bezier<szfld,ideg>(rfld, lfld, idif1, idif2, bary, eval, jmat, hmat);
  }else if(ibasis == FEBasis::Lagrange){
    METRIS_ASSERT(idif2 == DifVar::None);
    eval1_lagrange<szfld,ideg>(rfld, lfld, idif1, bary, eval, jmat);
  }
  return;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void eval1<1,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);\
template void eval1<2,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);\
template void eval1<4,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);\
template void eval1<3,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);\
template void eval1<6,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG_JACOBIAN)
#include BOOST_PP_LOCAL_ITERATE()



template <int szfld, int ideg>
void eval2(const dblAr2 & __restrict__  rfld,  
          const int * __restrict__  lfld,
          FEBasis ibasis, DifVar idif1, DifVar idif2,
          const double * __restrict__  bary, 
          double * __restrict__  eval, 
          double * __restrict__  jmat, 
          double * __restrict__  hmat){

  METRIS_ASSERT(idif2 == DifVar::None || idif2 == DifVar::Bary);
  METRIS_ASSERT(idif1 == DifVar::None || idif1 == DifVar::Bary);
  METRIS_ASSERT(ibasis == FEBasis::Bezier || ibasis == FEBasis::Lagrange);

  if constexpr (ideg == 1){
    for(int i=0;i<szfld;i++){
      eval[i] = bary[0]*rfld[lfld[0]][i]
              + bary[1]*rfld[lfld[1]][i]
              + bary[2]*rfld[lfld[2]][i];
    }
    if(idif1 == DifVar::Bary){
      for(int i=0;i<szfld;i++){
        jmat[0*szfld+i] = rfld[lfld[1]][i] - rfld[lfld[0]][i];
        jmat[1*szfld+i] = rfld[lfld[2]][i] - rfld[lfld[0]][i];
      }
    }
    if(idif2 == DifVar::Bary){
      for(int i=0;i<szfld;i++){
        hmat[0*szfld+i] = 0.0;
        hmat[1*szfld+i] = 0.0;
        hmat[2*szfld+i] = 0.0;
      }
    }
    return;
  } 

  if(ibasis == FEBasis::Bezier){
    eval2_bezier<szfld,ideg>(rfld, lfld, idif1, idif2, bary, eval, jmat, hmat);
  }else if(ibasis == FEBasis::Lagrange){
    METRIS_ASSERT(idif2 == DifVar::None);

    #ifndef NDEBUG
    try{
    #endif
    eval2_lagrange<szfld,ideg>(rfld, lfld, idif1, bary, eval, jmat); 
    #ifndef NDEBUG
    }catch(const MetrisExcept& e){
      printf("eval2_lagrange exception szfld %d ideg %d \n",szfld,ideg);
      intAr1(facnpps[ideg],lfld).print();
      throw(e);
    }
    #endif
  }

  return;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template void eval2<1,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);\
template void eval2<2,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);\
template void eval2<4,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);\
template void eval2<3,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);\
template void eval2<6,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG_JACOBIAN)
#include BOOST_PP_LOCAL_ITERATE()


template <int szfld, int ideg>
void eval3(const dblAr2 & __restrict__  rfld,  
           const int * __restrict__  lfld,
           FEBasis ibasis, DifVar idif1, DifVar idif2,
           const double * __restrict__  bary, 
           double * __restrict__  eval, 
           double * __restrict__  jmat, 
           double * __restrict__  hmat){

  METRIS_ASSERT(!(idif1 == DifVar::None && idif2 == DifVar::Bary));
  METRIS_ASSERT(idif1 == DifVar::None || idif1 == DifVar::Bary);
  METRIS_ASSERT(idif2 == DifVar::None || idif2 == DifVar::Bary);
  METRIS_ASSERT(ibasis == FEBasis::Bezier || ibasis == FEBasis::Lagrange);
  
  if constexpr (ideg == 1){
    
    for(int i=0;i<szfld;i++){
      eval[i] = bary[0]*rfld[lfld[0]][i]
              + bary[1]*rfld[lfld[1]][i]
              + bary[2]*rfld[lfld[2]][i]
              + bary[3]*rfld[lfld[3]][i];
    }
    if(idif1 == DifVar::Bary){
      
      for(int i=0;i<szfld;i++){
        jmat[0*szfld+i] = rfld[lfld[1]][i] - rfld[lfld[0]][i];
        jmat[1*szfld+i] = rfld[lfld[2]][i] - rfld[lfld[0]][i];
        jmat[2*szfld+i] = rfld[lfld[3]][i] - rfld[lfld[0]][i];
      }
    }
    if(idif2 == DifVar::Bary){
      
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

  if(ibasis == FEBasis::Bezier){
    eval3_bezier<szfld,ideg>(rfld, lfld, idif1,idif2, bary, eval, jmat, hmat);
  }else if(ibasis == FEBasis::Lagrange){
    eval3_lagrange<szfld,ideg>(rfld, lfld, idif1,idif2, bary, eval, jmat,hmat); 
  }
  return;
}



#define BOOST_PP_LOCAL_MACRO(n)\
template void eval3<1,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);\
template void eval3<3,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);\
template void eval3<6,n>(const dblAr2 & __restrict__  rfld,\
                         const int * __restrict__  lfld,\
                         FEBasis ibasis, DifVar idif1, DifVar idif2,\
                         const double * __restrict__  bary,\
                         double * __restrict__  eval,\
                         double * __restrict__  jmat,\
                         double * __restrict__  hmat);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG_JACOBIAN)
#include BOOST_PP_LOCAL_ITERATE()


// De Casteljau evaluation is based on the following recursive relation:
// \sum_{a \in I^d } P_aB_a = \sum_{1 <= i <= n+1} \xi_i \sum_{b \in I^{d-1}} B_b P_{e^i + b}
// where e^i is the (n+1) tuple of general term (\delta_{ij})_j. 
// In other words, "a P^d element is a P^1 element of P^{d-1} elements". 
// This leads to overlap in access, e.g. P^3 triangle will be interested in P_{111} several times.
// De Casteljau would be at the level of eval3, not here. 

// It is worth trying to reduce random accesses by walking the control point array linearly and 
// spending some more on computations of the \xi_i^\alpha_i
// This low-level routine offers the kernel for that. 
template <int ideg, int idim>
double eval_bezierfunc(const int * __restrict__ idx, 
  const double  * __restrict__  bary, int ider, double * __restrict__  dbez){
    double eval;
    for(int ii = 0; ii < idim + 1; ii++){
      if(idx[ii] < 0) {
        eval = 0;
        if(ider > 0) 
          for(int jj = 0; jj < idim; jj++) dbez[jj] = 0;
        return eval;
      }
    }
    if constexpr (idim == 1){
      eval = cbzedg.s[ideg][mul2nod(idx[0],idx[1])];
      if(idx[0] < 0 || idx[1] < 0){
        eval = 0;
        METRIS_ASSERT(ider == 0); // Not implemented yet, just fill with 0
      }
    }else if(idim == 2){
      eval = cbzfac.s[ideg][mul2nod(idx[0],idx[1],idx[2])];
      if(idx[0] < 0 || idx[1] < 0 || idx[2] < 0){
        eval = 0;
        METRIS_ASSERT(ider == 0); // Not implemented yet, just fill with 0
      }
    }else if(idim == 3){
      eval = cbztet.s[ideg][mul2nod(idx[0],idx[1],idx[2],idx[3])];
      if(idx[0] < 0 || idx[1] < 0 || idx[2] < 0 || idx[3] < 0){
        eval = 0;
        METRIS_ASSERT(ider == 0); // Not implemented yet, just fill with 0
      }
    }else{
      printf("## DIMENSION > 3 NOT SUPPORTED \n");
    }

    //double coef = eval; // Used for derivative

    for(int i=0;i<idim+1;i++){
      eval *= pow(bary[i],idx[i]);
    }

    if(ider <= 0) return eval;

    if constexpr(ideg == 0){
      for(int jj = 0; jj < idim; jj++) dbez[jj] = 0;
    } else{

      int idx2[idim+1];
      for(int ii = 0; ii < idim + 1; ii++) idx2[ii] = idx[ii];
      double dd[idim+1];
      for(int ii = 0; ii < idim + 1; ii++){
        idx2[ii] -= 1;
        dd[ii] = idim*eval_bezierfunc<ideg-1,idim>(idx,bary,-1,NULL);
        idx2[ii] += 1;
      }
      for(int ii = 0; ii < idim; ii++){
        dbez[ii] = dd[ii+1] - dd[0];
      }

    }

    ////  Multiplying by idx[i] / bary[i] is not stable.
    //double d1 = idx[0] > 0 ? coef*pow(bary[0],idx[0]-1)*idx[0] : 0; // pow(0.0,-1) to be avoided..
    //for(int i=1;i< idim+1;i++){
    //  d1 *= pow(bary[i],idx[i]);
    //}

    //for(int i =1; i < idim +1 ;i++){
    //  dbez[i-1] = idx[i] > 0 ? coef*pow(bary[i],idx[i]-1)*idx[i] : 0;
    //  for(int j=0;j<idim+1;j++){
    //    if(i==j) continue;
    //    dbez[i-1] *= pow(bary[j],idx[j]);
    //  }
    //  dbez[i-1] -= d1;
    //}
 
    return eval;
}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double eval_bezierfunc< n ,1>(const int * __restrict__ idx, \
   const double  * __restrict__  bary, int ider, double * __restrict__  dbez);\
template double eval_bezierfunc< n ,2>(const int * __restrict__ idx, \
   const double  * __restrict__  bary, int ider, double * __restrict__  dbez);\
template double eval_bezierfunc< n ,3>(const int * __restrict__ idx, \
   const double  * __restrict__  bary, int ider, double * __restrict__  dbez);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG_JACOBIAN)
#include BOOST_PP_LOCAL_ITERATE()


//     -- Representation using roots. I thought this might be fast but it's the opposite. 
//     The main idea is to exploit the fact any Lagrange (with uniform nodes) writes as a product of d roots. 
//     This should have been faster than doing Bez2Lag and computing Bez on the side. 
//     It probably is, too. 
//     Right now, this is between 10 and 15 times slower than writing explicitly, even with -Ofast or -O3. 
//     At least for P2 which is all indirections and no computations. 
//     --- Initial comments:
//     Notice that the Lagrange of multi-index \alpha writes 
//       C \prod_{i=1,n}\prod_{0<=j<\alpha_i} (\xi_i - j/d)
//     Proof: Cardinal of product set is indeed d as there are alpha_i 
//     terms strictly under any alpha_i and at or above 0, and that  
//     |\alpha|_1 = d. The polynomial is thus of the right degree. 
//     For any \beta \ne \alpha there is i such that \beta_i < \alpha_i,
//     thus L(\beta/d) = 0. QED
//     This costs d products to evaluate versus ~d^n if we use a matrix 
//     to a reference basis (Bernstein or otherwise)

//     -- ider = get derivative, jmat stores it
//     Panzoult derivatives: d_i (physical) = d_{i+1} - d_1 (bary)
template <int ideg, int idim>
double eval_lagrangefunc(const int * __restrict__ idx, 
  const double * __restrict__ bary,int ider,double * __restrict__ dlag){
    double fac,eval,dd,d1,eval_lagrange;
    double ev[1+idim];
    int i,j,k;
    if(idim != 2 && idim != 3){
        printf("Only dimensions 2 and 3 supported, %d supplied\n",idim);
        exit(1);
    }

    if(ider <= 0){
        fac = 1.0;
        eval = 1.0;
        dd = 1.0/ideg;
        for(i = 0;i<idim+1;i++){
            for(j=0;j<idx[i];j++){
                eval *= (bary[i] - j*dd);
                fac  *= (idx[i]-j)*dd;
            }
        }
        return eval / fac;
    }

//     -- Derivatives in the sense of Panzoult: d_i (phy) = d_(i+1) - d_1 (bar)
//     This convention is so unnatural it cannot hurt to recall it every time. 

//     -- Compute eval, but separate contributions from each barycentric crd
//     This will minimize recomputations later on (when diff)
    fac = 1.0;
    dd  = 1.0/ideg;
    eval = 1.0;
    for (i = 0; i<(idim+1); i++){
        ev[i] = 1.0 ; //! separate the contributions to minimize recomputations
        for( j = 0; j < idx[i]; j++){
            ev[i] *=(bary[i] - j*dd);
            fac   *=(idx[i] - j)*dd;
        }
        eval *= ev[i];
    }
    eval_lagrange = eval / fac;

//     -- Now compute diff w.r.t 1 (bary) for Panzoult's d_(i+1) - d_1.
//      - Compute d_1 of factor depending on \xi_1
//      The factor writes (\xi_1 - i_1/d)(\xi_1 - i_2/d)... 
//      its derivative is sum_i \prod_{j\ne i} (\xi_1 - j/d)
    d1 = 0.0;
    for(i = 0;i<idx[0]; i++){ //! the factor we exclude
        eval = 1.0;
        for(j = 0; j < idx[0] ; j++){
            if(j == i) continue ;  // the factors remaining
            eval *=(bary[0]-j*dd);  // their contribution
        }
        d1 += eval ; //! the i-th term of the sum
    }
//      - Finalize by multiplying by the constant (w.r.t. \xi_1) factor
//      We can shave off one product by not *fac'ing this yet. 
    if(idim == 3){
      d1 *= ev[1]*ev[2]*ev[3];
    } else{
      d1 *= ev[1]*ev[2];
    }
    
//   -- Do the same thing with the others. 
    for(i = 1; i < idim+1; i++){   // this is truly idim, not idim + 1 
        dlag[i-1] = 0.0;
        for (j = 0; j < idx[i]; j++){
            eval = 1.0;
            for (k = 0; k < idx[i]; k++){
                if(k == j) continue;
                eval *= (bary[i]-k*dd);
            }
            dlag[i-1] += eval;
        }
//       -- Multiply by constant factors
        for(j = 0; j < idim+1; j++){
            if(j == i) continue;
            dlag[i-1] *= ev[j];
        }
//       -- Final factor and -d1, recall not multiplied by fac
        dlag[i-1] = (dlag[i-1] - d1) / fac;
    }

    return  eval_lagrange;
}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double eval_lagrangefunc< n ,2>(const int * __restrict__ idx, \
    const double * __restrict__ bary,int ider,double * __restrict__ dlag);\
template double eval_lagrangefunc< n ,3>(const int * __restrict__ idx, \
    const double * __restrict__ bary,int ider,double * __restrict__ dlag);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG_JACOBIAN)
#include BOOST_PP_LOCAL_ITERATE()


} // End namespace
