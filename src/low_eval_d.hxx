//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __LOW_EVAL_D__
#define __LOW_EVAL_D__

#include "types.hxx"
#include "ho_constants.hxx"
#include "utils/aux_misc.hxx"
#include "low_eval_d_SurrealS.hxx"
#include "low_eval_d_bezier.hxx"
#include "low_eval_d_helper.hxx"

#include "../SANS/Surreal/SurrealS.h"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
#include <boost/hana.hpp> 


namespace Metris{


/*
~~ low_eval_d.hxx: polynomial evaluation routines diff'd wrt poly coeffs ~~ 

~~ Main functions:
 - eval3_d: idem as eval3 with added ivar template param. Diff eval3 wrt poly coeff of 
 rank ivar. 
 - eval3_d_SurrealS0: Not to be called directly. Expects a tuple containing double* and Surreal<double,...>* 
 types determining the variable. eval3_d wraps around this for ease of use. 

~~ Inputs:
 - inherited from eval3: rfld,lfld,bary
 - ivar (template param): rank of the DoF in the element, comprised in [[0,getnnod3(ideg)-1]].
 - dfld (optional argument): derivatives of rfld[ivar][:] (default Id mat i.e. rfld[ivar] itself)
    dfld(i,j) = d_i rfld(idof,j) 

~~ Outputs
Derivatives w.r.t. coeff are D_i, w.r.t. bary are d_i
 - deval: szfld*szfld array. deval(i,j) = D_i eval[j]
 - djmat: szfld*szfld*tdim3 array. djmat(i,j)[k] = D_i d_j eval[k]
 - dhmat: szfld*[(szfld+1)*(szfld+2)/2]*tdim3 array. dhmat(i,j)[k] = D_i s_j eval[k]
 where s_j is a D^2 (bary) op ordered 1,1 1,2 etc. see low_eval.hxx. 

See low_eval.hxx header for other parameters and outputs. 
*/


/*
eval3_d_SurrealS0 is not designed to be called directly: call instead eval3_d. 

T designed to be a heterogeneous Boost::hana container such as a hana::tuple
Each type in the tuple can be a szfld allocated array. 
This tuple may hold SANS::SurrealS<szfld>* and no other n. 
Ideally, it should also hold exactly 1 SANS::SurrealS<szfld>*, doing otherwise
would be wasteful. 

This container can only be accessed by compile-time constants, which is why
every loop is replaced by a hana::while_ loop. 

Global indexing via lfld is no longer supported. This is because rfld is accessed
(in eval3) at lfld[i], but this cannot be made compile-time known (it depends on the mesh...)
As such, rfld is assumed ordered as given by ordedg, ordfac, ordtet. 
*/



/*
This is the wrapper that should be called. 
It sets up the tuple, gathers the derivatives.
- deval: derivatives of eval w.r.t. DoFs
  deval[N*i + j] with 0 <= i < szfld, N = szfld*szfld, 0 <= j < N
                 is the derivative w.r.t. to the j-th global (sized szfld) DoF
                 of the i-th evaluated coordinat
- djmat: similarly, deval[N*szfld*i + N*j + k] = d_k jmat(i,j)
- dhmat: idem, hmat (4D)

The plan is to constitute a pack of types using a hana::tuple
For instance hana::tuple<double, double, SANS::SurrealS, double>
Now, we'd want to populate it with values from rfld, or from our allocated SANS::SurrealS types. 
But we can't, because a hana::tuple is a purely compile-time construct. 
So we'll have to unpack its type template parameter into an std::tuple which can be 
modified at runtime! 
This is what we'll pass to eval3_d_SurrealS0. 

- hana offers hana::replicate to build a "constexpr dynamically sized" tuple
- we can "modify" types of this hana::tuple with the replace_at_c helper 
-> all to (double*) then ivar to (SANS::SurrealS<N,double>*)
- Problem: hana::tuple is fully constexpr so the double* etc values cannot be set
-> convert to std::tuple and pass that to the low-level routine
- Problem: std::tuple does not offer [] so we'll have to use std::get... which is not 
 overloadable to a basic dblAr2... 
-> wrap std::tuple in basic struct overloading []
- semi-Problem: [] argument must be compile-time known type; we can't deduce
 non-type template parameters... so we'll pass a hana::integral_constant
 and parameter deduction will catch its "value" (non-type parameter)
-> overload [] in MeshArrays to handle integral constants (simple cast to int)
*/

template <int szfld, int ideg, int ivar, int nvar = szfld>
void eval3_d(const dblAr2 & __restrict__ rfld,  
             const int * __restrict__ lfld,
             FEBasis ibasis, DifVar idif1, DifVar idif2, 
             const double * __restrict__ bary,
             double * __restrict__  eval,
             double * __restrict__  jmat,
             double * __restrict__  hmat,
             double * __restrict__  deval,
             double * __restrict__  djmat,
             double * __restrict__  dhmat,  
             const double * __restrict__ dfld = NULL);

template <int szfld, int ideg, int ivar, int nvar = szfld>
void eval2_d(const dblAr2 & __restrict__ rfld,  
             const int * __restrict__ lfld,
             FEBasis ibasis, DifVar idif1, DifVar idif2, 
             const double * __restrict__ bary,
             double * __restrict__  eval,
             double * __restrict__  jmat,
             double * __restrict__  hmat,
             double * __restrict__  deval,
             double * __restrict__  djmat,
             double * __restrict__  dhmat,  
             const double * __restrict__ dfld = NULL);






/* Secondary routines */
template <int szfld, int tdim, int ideg, int ivar, int nvar = szfld>
void eval_d_direct(const dblAr2 & __restrict__ rfld,  
                   const int * __restrict__ lfld,
                   FEBasis ibasis, DifVar idif1, DifVar idif2, 
                   const double * __restrict__ bary,  
                   double * __restrict__  eval,
                   double * __restrict__  jmat,
                   double * __restrict__  hmat,
                   double * __restrict__  deval,
                   double * __restrict__  djmat,
                   double * __restrict__  dhmat,
                   const double * __restrict__  dfld = NULL);
//// Attention: pointers to arrays; this is to handle NULLs (idif1/2 == 0).
//template <typename T, int szfld, int tdim, int ideg, int nvar = szfld>
//void eval_d_SurrealS0(const       T& __restrict__  rfld,  
//                      FEBasis ibasis, DifVar idif1, DifVar idif2, 
//                      const double * __restrict__  bary, 
//                      SANS::DLA::VectorS<                  szfld,SANS::SurrealS<nvar,double>> *eval, 
//                      SANS::DLA::MatrixS< tdim,            szfld,SANS::SurrealS<nvar,double>> *jmat, 
//                      SANS::DLA::MatrixS<(tdim*(tdim+1))/2,szfld,SANS::SurrealS<nvar,double>> *hmat);









/*
Main eval3_d routine. Calls either eval3_d_SurrealS or eval3_d_direct 
depending on arguments. Direct version can handle both Bézier and Lagrange. 
SANS::SurrealS can only handle Bézier (TODO: implement Lagrange evals with SANS::SurrealS)
*/
/* --- dfld
dfld of size (szfld,szfld) specifies the derivatives of the input field. 
This holds the derivatives of the ivar-th coeff in rfld w.r.t. considered
variable. 

Example, when evaluating metric at xi using metrics at control points, 
metrics at control points depend on the variable in some manner 
(determined by back mesh). 

We restrict ourselves to the case where only the ivar-th rfld entry (e.g. metric)
is non-constant, though we could imagine cases where all polynomial coeffs
depend on a common variable. 
*/
template <int szfld, int ideg, int ivar, int nvar>
void eval3_d(const dblAr2 & __restrict__ rfld,  
             const int * __restrict__ lfld,
             FEBasis ibasis, DifVar idif1, DifVar idif2, 
             const double * __restrict__ bary,
             double * __restrict__  eval,
             double * __restrict__  jmat,
             double * __restrict__  hmat,
             double * __restrict__  deval,
             double * __restrict__  djmat,
             double * __restrict__  dhmat,  
             const double * __restrict__ dfld){

  METRIS_ASSERT_MSG(!(dfld == NULL && nvar != szfld),
                     "If nvar != szfld, Jacobian matrix must be specified");

  // Only case where SANS::SurrealS benchmarked faster than direct method
  // But also the only one with a Hessian implementation...
  if(idif2 == DifVar::Bary || (ideg == 2 && idif1 == DifVar::Bary && ibasis == FEBasis::Bezier)){
    METRIS_ASSERT_MSG(ibasis == FEBasis::Bezier,
      "## EITHER IMPLEMENT HESSIAN IN DIRECT OR IMPLEMENT LAGRANGE IN SANS::SurrealS\n");
    eval_d_SurrealS<szfld, 3, ideg,  ivar, nvar>
                 (rfld,lfld,ibasis, idif1, idif2,bary,eval,jmat,hmat,deval,djmat,dhmat, dfld);
  }else{
    eval_d_direct<szfld, 3, ideg,  ivar, nvar>
                 (rfld,lfld,ibasis, idif1, idif2,bary,eval,jmat,hmat,deval,djmat,dhmat, dfld);
  }
}



template <int szfld, int ideg, int ivar, int nvar>
void eval2_d(const dblAr2 & __restrict__ rfld,  
             const int * __restrict__ lfld,
             FEBasis ibasis, DifVar idif1, DifVar idif2, 
             const double * __restrict__ bary,
             double * __restrict__  eval,
             double * __restrict__  jmat,
             double * __restrict__  hmat,
             double * __restrict__  deval,
             double * __restrict__  djmat,
             double * __restrict__  dhmat,  
             const double * __restrict__ dfld){

  METRIS_ASSERT_MSG(!(dfld == NULL && nvar != szfld),
                     "If nvar != szfld, Jacobian matrix must be specified");

  // Only case where SANS::SurrealS benchmarked faster than direct method
  // But also the only one with a Hessian implementation...
  if(idif2 == DifVar::Bary || (ideg == 2 && idif1 == DifVar::Bary && ibasis == FEBasis::Bezier)){
    METRIS_ASSERT_MSG(ibasis == FEBasis::Bezier,
      "## EITHER IMPLEMENT HESSIAN IN DIRECT OR IMPLEMENT LAGRANGE IN SANS::SurrealS\n");
    eval_d_SurrealS<szfld,2,  ideg,  ivar, nvar>
                 (rfld,lfld,ibasis, idif1, idif2,bary,eval,jmat,hmat,deval,djmat,dhmat, dfld);
  }else{
    eval_d_direct<szfld, 2, ideg,  ivar, nvar>
                 (rfld,lfld,ibasis, idif1, idif2,bary,eval,jmat,hmat,deval,djmat,dhmat, dfld);
  }
}


/* --------------------------------------------------------------------
  Secondary functions
   -------------------------------------------------------------------- */


/* Exploit F_K = sum f_\alpha P_\alpha with f either Bernstein or Lagrange
   and fill in derivatives directly. If we denote J_idof the Jacobian of 
   rfld[idof] w.r.t. variables, then the final Jacobian is simply 
   fun * J_idof where fun is either the Bernstein or Lagrange ivar-th 
   basis function. 
   Similar reasoning for jmat and hmat (not implemented). */
template <int szfld, int tdim, int ideg,  int ivar, int nvar>
void eval_d_direct(const dblAr2 & __restrict__ rfld,  
             const int * __restrict__ lfld,
             FEBasis ibasis, DifVar idif1, DifVar idif2, 
             const double * __restrict__ bary,  
             double * __restrict__   eval,
             double * __restrict__   jmat,
             double * __restrict__   hmat,
             double * __restrict__  deval,
             double * __restrict__  djmat,
             [[maybe_unused]] double * __restrict__  dhmat,
             const double * __restrict__  dfld){

  if constexpr(tdim == 3){
    eval3<szfld,ideg>(rfld,lfld,ibasis,idif1,idif2,bary,eval,jmat,hmat);
  }else{
    eval2<szfld,ideg>(rfld,lfld,ibasis,idif1,idif2,bary,eval,jmat,hmat);
  }

  auto ordent = ORDELT(tdim);

  double fun, dfun[tdim];
  int ideriv = idif1 == DifVar::None ? 0 : 1;
  if(ibasis == FEBasis::Bezier){
    fun = eval_bezierfunc<ideg,tdim>(ordent[ideg][ivar],bary,ideriv,dfun);
  }else{
    fun = eval_lagrangefunc<ideg,tdim>(ordent[ideg][ivar],bary,ideriv,dfun);
  } 

  //-- Evaluation and derivatives
  if(dfld == NULL){
    /* In this case, the Jacobian matrix (of ivar(:) wrt variables)
       is the identity. The gradient w.r.t. idof is simply 
       the basis function times (1, ..., 1). */
    for(int idof = 0; idof < nvar; idof++){
      for(int icmp = 0; icmp < szfld; icmp++){
        deval[szfld*idof + icmp] = 0;
      }
      deval[szfld*idof + idof] = fun;
    }
  }else{
    /* In this case, we may have nvar != szfld 
       The Jacobian of F_K w.r.t. idof is held in dfld:
       dfld(i,j) = d_i rfld(idof,j)  */
    for(int idof = 0; idof < nvar; idof++){
      for(int icmp = 0; icmp < szfld; icmp++){
        deval[szfld*idof + icmp] = fun*dfld[idof*szfld + icmp];
      }
    }
  }
  //-


  //-- Jacobian matrix and derivatives 
  /* Denote D_\alpha^i the i-th component (a tensor) of the derivative w.r.t. 
     P_\alpha. Then:
     D_\alpha^i J_{jk} = D_\alpha^i (d_k J_K^j) 
                       = (i==j) d_k B_\alpha   */
  if(idif1 == DifVar::Bary){
    if(dfld == NULL){
      for(int ii = 0; ii < nvar*tdim*szfld; ii++){
        djmat[ii] = 0;
      }
      for(int idof = 0; idof < nvar; idof++){
        for(int itdim = 0; itdim < tdim; itdim++){
          djmat[tdim*szfld*idof + itdim*szfld + idof] = dfun[itdim];
        }
      }
    }else{
      for(int idof = 0; idof < nvar; idof++){
        for(int itdim = 0; itdim < tdim; itdim++){
          for(int icmp = 0; icmp < szfld; icmp++){
            djmat[idof*tdim*szfld + itdim*szfld + icmp] = dfun[itdim]*dfld[idof*szfld + icmp];
          }
        }
      }
    }
  }
  //-

  //-- Hessian (bary) and derivatives
  // not implemented
  if(idif2 == DifVar::Bary){
    METRIS_THROW_MSG(TODOExcept(),"Unsupported diff2 in eval_d_direct");
    //for(int i=0; i<szfld; i++){
    //  for(int j = 0; j < 3; j++){
    //    dhmat[0*szfld*3 + i*3 + j] = shmat[0*szfld + i].deriv(j);
    //    dhmat[1*szfld*3 + i*3 + j] = shmat[1*szfld + i].deriv(j);
    //    dhmat[2*szfld*3 + i*3 + j] = shmat[2*szfld + i].deriv(j);
    //    dhmat[3*szfld*3 + i*3 + j] = shmat[3*szfld + i].deriv(j);
    //    dhmat[4*szfld*3 + i*3 + j] = shmat[4*szfld + i].deriv(j);
    //    dhmat[5*szfld*3 + i*3 + j] = shmat[5*szfld + i].deriv(j);
    //  }
    //}
    //}
  }
  //-
}












} // End namespace




#endif