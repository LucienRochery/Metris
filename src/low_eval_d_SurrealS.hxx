//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __LOW_EVAL_D_SURREALS__
#define __LOW_EVAL_D_SURREALS__

//#include "../SANS/Surreal/SurrealS.h"
//#include "ho_constants.hxx"
#include "low_eval_d.hxx"
//#include "types.hxx"



namespace Metris{




/*
--------------------------
--------------------------
--------------------------
--------------------------
  This is stupid: second derivatives are 0...
--------------------------
--------------------------
--------------------------
--------------------------

  Second derivatives wrt ivar1, ivar2. Potentially ivar1 = ivar2. 
  nvar effective (nvareff) = nvar if ivar1 == ivar2, 
                            2nvar otherwise
  - deval size nvareff * szfld
  - djmat size nvareff * szfld*tdim 
  - dhmat size nvareff * szfld*(tdim*(tdim+1))/2

  - d2eval size neffsym * etc. 
  neffsym is the size of a symmetric matrix with n = nvareff ; ordered as in sym2idx 

  Input (can be NULL):
  - dfld size nvareff * szfld. 

  Note the variables are ordered as 1 ...  nvar -> ivar1
                             nvar + 1 ... 2nvar -> ivar2
template <int szfld, int tdim, int ideg,  int ivar1, int ivar2, int nvar>
void eval_d2_SurrealS(const dblAr2 & __restrict__ rfld,  
                      const int * __restrict__ lfld,
                      FEBasis ibasis, DifVar idif1, DifVar idif2, 
                      const double * __restrict__ bary,
                      double * __restrict__  eval,
                      double * __restrict__  jmat,
                      double * __restrict__  hmat,
                      double * __restrict__  deval,
                      double * __restrict__  djmat,
                      double * __restrict__  dhmat, 
                      double * __restrict__  d2eval,
                      double * __restrict__  d2jmat,
                      double * __restrict__  d2hmat,  
                      const double * __restrict__ dfld){

  constexpr auto entnpps = ENTNPPS(tdim);

  // -- For clarity. 
  // Number of (r|l)fld entries
  constexpr int nrfld = entnpps[ideg];
  static_assert(nrfld < 31);
  // -
  static_assert(ivar1 >= 0 && ivar1 < entnpps[ideg]);
  static_assert(ivar2 >= 0 && ivar2 < entnpps[ideg]);

  auto rfld_0 = hana::replicate<hana::tuple_tag>((const double *) 1, hana::size_c<nrfld>);

  constexpr int nmem = 1 + (ivar1 != ivar2);
  constexpr int nvareff = nmem * nvar;

  // If ivar1 == ivar2, there's only nvar variables 
  // Else, there's 2nvar AND each dof is still a SurrealS<SurrealS>>. 
  typedef SANS::SurrealS<nvareff,SANS::SurrealS<nvareff,double>> d2real;


  //typedef typename typename std::conditional<ivar1 == ivar2, 
  //                                  d2real,
  //                                  d1real>::type direal;


  auto rfld_1 = to_std_tuple(replace_at_c<ivar1>(rfld_0,(d2real*) 0));
  auto rfld_2 = ivar1 != ivar2 ? 
                to_std_tuple(replace_at_c<ivar2>(rfld_1,(d2real*) 0)) 
              : to_std_tuple(replace_at_c<ivar1>(rfld_0,(d2real*) 0)) ;

  // Store the dof as either d1real smem[2][szfld]
  // or d2real smem[1][szfld]
  d2real smem[nmem][szfld];

  constexpr std::array<int,2> lvar = {ivar1,ivar2};

  if(dfld == NULL){
    assert(nvar == szfld);

    CT_FOR0_INC(0,1,imem){
      constexpr int ivar = lvar[imem];
      for(int icmp = 0; icmp < szfld; icmp++){
        smem[imem][icmp].value().value() = rfld[lfld[ivar]][icmp];
        for(int j = 0; j < nvar; j++){
          //smem[imem][icmp].deriv(j).value() = 0;
          smem[imem][icmp].value().deriv(j) = 0;
          for(int i = 0; i < nvar; i++){
            smem[imem][icmp].deriv(i).deriv(j) = 0;
            smem[imem][icmp].deriv(j).deriv(i) = 0;
          }
        }
        // See SurrealS(1,...,4)_btest.cpp in SANS ; only need to set value().deriv()
        //smem[imem][icmp].deriv(icmp).value() = 1;
        smem[imem][icmp].value().deriv(icmp) = 1;
      }
    }CT_FOR1(imem);
    //if constexpr(ivar1 != ivar2){
    //  for(int icmp = 0; icmp < szfld; icmp++){
    //    smem[0][icmp].value() = rfld[lfld[ivar1]][icmp];
    //    smem[1][icmp].value() = rfld[lfld[ivar2]][icmp];
    //    for(int j = 0; j < nvar; j++){
    //      smem[0][icmp].deriv(j) = 0;
    //      smem[1][icmp].deriv(j) = 0;
    //    }
    //    smem[0][icmp].deriv(icmp) = 1;
    //    smem[1][icmp].deriv(icmp) = 1;
    //  }
    //}else{
    //  for(int icmp = 0; icmp < szfld; icmp++){
    //    smem[icmp].value().value() = rfld[lfld[ivar1]][icmp];
    //    for(int j = 0; j < nvar; j++){
    //      smem[icmp].deriv(j).value() = 0;
    //      smem[icmp].value().deriv(j) = 0;
    //      for(int i = 0; i < nvar; i++){
    //        smem[icmp].deriv(i).deriv(j) = 0;
    //        smem[icmp].deriv(j).deriv(i) = 0;
    //      }
    //    }
    //    smem[icmp].value().deriv(icmp) = 1;
    //    smem[icmp].deriv(icmp).value() = 1;
    //  }
    //}
  }else{
    METRIS_THROW_MSG(TODOExcept(), "Implement d2_SurrealS with dfld");
    //for(int icmp = 0; icmp < szfld; icmp++){
    //  smem[icmp].value() = rfld[lfld[ivar]][icmp];
    //  for(int idof = 0; idof < nvar; idof++){
    //    smem[icmp].deriv(idof) = dfld[idof*szfld + icmp];
    //  }
    //}
  }

  tuple_wrapper w_op(rfld_1);

  // Populate w_op tuple with appropriate types. 
  // The idea is to get something like a tuple<double*, SANS::SurrealS*, double*, double*...>
  // Storing, for the double*'s, the non-dof rfld entries and, for the one SANS::SurrealS*, 
  // the DoF SANS::SurrealS array. 
  hana::while_(hana::less.than(hana::int_c<nrfld>), 0_c, [&](auto i_c){
    constexpr int i = i_c;
    if constexpr(i == ivar1){
      w_op[i_c] = smem[0];
    }else if(i == ivar2){
      w_op[i_c] = smem[1];
    }else{
      w_op[i_c] = rfld[lfld[i]];
    }
  return i_c+1_c;});  
//  SANS::SurrealS<nvar,double> seval[szfld];
  constexpr int nnsym = (tdim*(tdim+1))/2;

  SANS::DLA::VectorS<      szfld,d2real> seval;
  SANS::DLA::MatrixS<tdim ,szfld,d2real> sjmat;
  SANS::DLA::MatrixS<nnsym,szfld,d2real> shmat;


  //hana::while_(hana::less.than(hana::int_c<nrfld>), 0_c, [&](auto i_c){
  //  constexpr int i = i_c;
  //  printf("Debug rfld[%d] = \n",i);
  //  std::cout<<"0:"<<w_op[i_c][0]<<std::endl;
  //  std::cout<<"1:"<<w_op[i_c][1]<<std::endl;
  //  std::cout<<"2:"<<w_op[i_c][2]<<std::endl;
  //return i_c+1_c;});  


  eval_d_SurrealS0<decltype(w_op),szfld,tdim,ideg,nvareff>
                     (w_op,ibasis,idif1,idif2,bary,&seval,&sjmat,&shmat);

  for(int icmp=0; icmp<szfld; icmp++){
    eval[icmp] = seval[icmp].value().value();

    for(int idof = 0; idof < nvareff; idof++){
      //deval[szfld*idof + icmp] = seval[icmp].deriv(idof).value();
      deval[szfld*idof + icmp] = seval[icmp].value().deriv(idof);
    }

    if(idif1 == DifVar::Bary){
      for(int itdim = 0; itdim < tdim; itdim++){
        jmat[itdim*szfld + icmp] = sjmat(itdim,icmp).value();
        for(int idof = 0; idof < nvareff; idof++){
          djmat[idof*tdim*szfld + itdim*szfld + icmp] = sjmat(itdim,icmp).deriv(idof);
        } 
      }
    }

    if(idif2 == DifVar::Bary){
      for(int ii = 0; ii < nnsym; ii++) hmat[ii*szfld + icmp] = shmat(ii,icmp).value();
      //hmat[0*szfld + icmp] = shmat(0,icmp).value();
      //hmat[1*szfld + icmp] = shmat(1,icmp).value();
      //hmat[2*szfld + icmp] = shmat(2,icmp).value();
      //hmat[3*szfld + icmp] = shmat(3,icmp).value();
      //hmat[4*szfld + icmp] = shmat(4,icmp).value();
      //hmat[5*szfld + icmp] = shmat(5,icmp).value();
      for(int idof = 0; idof < tdim; idof++){
        for(int ii = 0; ii < nnsym; ii++) dhmat[idof*nnsym*szfld + 0*szfld + icmp] = shmat(ii,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 0*szfld + icmp] = shmat(0,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 1*szfld + icmp] = shmat(1,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 2*szfld + icmp] = shmat(2,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 3*szfld + icmp] = shmat(3,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 4*szfld + icmp] = shmat(4,icmp).deriv(idof);
        //dhmat[idof*sdim3*szfld + 5*szfld + icmp] = shmat(5,icmp).deriv(idof);
      }
    }

  }
}
*/












} // End namespace




#endif