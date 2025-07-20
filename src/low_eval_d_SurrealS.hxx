//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_EVAL_D_SURREALS__
#define __LOW_EVAL_D_SURREALS__


#include "../SANS/Surreal/SurrealS.h"
#include "ho_constants.hxx"
#include "types.hxx"
#include "low_eval_d_bezier.hxx"
#include "low_eval_d_helper.hxx"
#include "utils/CT_loop.hxx"
#include <boost/hana.hpp> 


namespace Metris{

template <typename T, int szfld, int tdim, int ideg,  int nvar>
void eval_d_SurrealS0(const T& __restrict__  rfld,   
                      FEBasis ibasis, DifVar idif1, DifVar idif2, 
                      const double * __restrict__  bary, 
                      SANS::DLA::VectorS<     szfld,SANS::SurrealS<nvar,double>> *eval, 
                      SANS::DLA::MatrixS<tdim,szfld,SANS::SurrealS<nvar,double>> *jmat, 
                      SANS::DLA::MatrixS<(tdim*(tdim+1))/2,szfld,SANS::SurrealS<nvar,double>> *hmat){
  if constexpr (ideg == 1){
    
    for(int icmp = 0; icmp < szfld; icmp++){
      (*eval)[icmp] = bary[0]*rfld[boost::hana::int_c<0>][icmp];
      // Note bounds included
      CT_FOR0_INC(1,tdim,ii){  
        // Macro creates c_ii
        (*eval)[icmp] += bary[ii]*rfld[c_ii][icmp];
      }CT_FOR1(ii);
    }

    if (idif1 == DifVar::Bary){
      
      for(int i = 0; i < szfld; i++){
        // Note upper bound excluded
        CT_FOR0_INC(1,tdim,jj){
          (*jmat)(jj-1,i) = rfld[c_jj][i] - rfld[boost::hana::int_c<0>][i];

        }CT_FOR1(jj);
      }
    }
    if (idif2 == DifVar::Bary){
      constexpr int nnsym = tdim*(tdim+1)/2;
      
      for(int i = 0; i < szfld; i++){
        for(int jj = 0; jj < nnsym; jj++) (*hmat)(jj,i) = 0;
      }
    }
    return;
  } 

  if(ibasis == FEBasis::Bezier){
    if constexpr(tdim == 2){
      eval2_d_bezier<T,szfld,ideg,nvar>(rfld, idif1, idif2, bary, eval, jmat, hmat);
    }else{
      eval3_d_bezier<T,szfld,ideg,nvar>(rfld, idif1, idif2, bary, eval, jmat, hmat);
    }
  }else{
    METRIS_THROW_MSG(TODOExcept(),"Implement eval3_d_LAGRANGE");
  } 
  return;
}


template <int szfld, int tdim, int ideg,  int nvar>
void eval_d_SurrealS0_simple(MeshArray2D<SANS::SurrealS<nvar,double>> &rfld,   
                             FEBasis ibasis, DifVar idif1, DifVar idif2, 
                             const double * __restrict__  bary, 
                             SANS::DLA::VectorS<     szfld,SANS::SurrealS<nvar,double>> *eval, 
                             SANS::DLA::MatrixS<tdim,szfld,SANS::SurrealS<nvar,double>> *jmat, 
                             SANS::DLA::MatrixS<(tdim*(tdim+1))/2,szfld,SANS::SurrealS<nvar,double>> *hmat){
  //printf("Debug surreal0 bary = %f %f %f %f field:\n",bary[0]
  //  ,bary[1],bary[2],bary[3]);
  //std::cout<<"0: "<<rfld[boost::hana::int_c<0>][0]<<" "<<rfld[boost::hana::int_c<0>][1]<<" "
  //                <<rfld[boost::hana::int_c<0>][2]<<std::endl;
  //std::cout<<"1: "<<rfld[boost::hana::int_c<1>][0]<<" "<<rfld[boost::hana::int_c<1>][1]<<" "
  //                <<rfld[boost::hana::int_c<1>][2]<<std::endl;
  //std::cout<<"2: "<<rfld[2_c][0]<<" "<<rfld[2_c][1]<<" "
  //                <<rfld[2_c][2]<<std::endl;
  //std::cout<<"3: "<<rfld[3_c][0]<<" "<<rfld[3_c][1]<<" "
  //                <<rfld[3_c][2]<<std::endl;
  if constexpr (ideg == 1){
    //hana::while_(hana::less.than(szfld_c), boost::hana::int_c<0>, [&](auto i_c){
    //  constexpr int i = i_c;
    
    for(int icmp = 0; icmp < szfld; icmp++){
      (*eval)[icmp] = bary[0]*rfld[boost::hana::int_c<0>][icmp];
      // Note bounds included
      for(int ii = 1; ii <= tdim; ii++){
        (*eval)[icmp] += bary[ii]*rfld[ii][icmp];
      }
      //CT_FOR0_INC(1,tdim,ii){  
      //  // Macro creates c_ii
      //  (*eval)[icmp] += bary[ii]*rfld[c_ii][icmp];
      //}CT_FOR1(ii);
      //(*eval)[icmp] = bary[0]*rfld[boost::hana::int_c<0>][icmp]
      //              + bary[1]*rfld[boost::hana::int_c<1>][icmp]
      //              + bary[2]*rfld[2_c][icmp]
      //              + bary[3]*rfld[3_c][icmp];
    }
    //return i_c+boost::hana::int_c<1>;});

    if (idif1 == DifVar::Bary){
      
      for(int i = 0; i < szfld; i++){
        // Note upper bound excluded
        for(int jj = 1;jj <= tdim; jj++){
          (*jmat)(jj-1,i) = rfld[jj][i] - rfld[0][i];
        }
        //CT_FOR0_INC(1,tdim,jj){
        //  (*jmat)(jj-1,i) = rfld[c_jj][i] - rfld[boost::hana::int_c<0>][i];
        //}CT_FOR1(jj);
        //(*jmat)(0,i) = rfld[boost::hana::int_c<1>][i] - rfld[boost::hana::int_c<0>][i];
        //(*jmat)(1,i) = rfld[2_c][i] - rfld[boost::hana::int_c<0>][i];
        //(*jmat)(2,i) = rfld[3_c][i] - rfld[boost::hana::int_c<0>][i];
      }
    }
    if (idif2 == DifVar::Bary){
      constexpr int nnsym = tdim*(tdim+1)/2;
      
      for(int i = 0; i < szfld; i++){
        for(int jj = 0; jj < nnsym; jj++) (*hmat)(jj,i) = 0;
        //(*hmat)(0,i) = 0;
        //(*hmat)(1,i) = 0;
        //(*hmat)(2,i) = 0;
        //(*hmat)(3,i) = 0;
        //(*hmat)(4,i) = 0;
        //(*hmat)(5,i) = 0;
      }
    }
    return;
  } 

  //std::cout<<"Debug simple:\n";
  //constexpr int npps = tdim == 2 ? getnnod2(ideg) : getnnod3(ideg);
  //CT_FOR0_INC(0,npps-1,ii){
  //  //std::cout<<"rfld["<<ii<<"] = "<<rfld[c_ii]<<"\n";
  //  std::cout<<"rfld["<<ii<<"] = "<<(*rfld[c_ii])<<"\n";
  //}CT_FOR1(ii);

  if(ibasis == FEBasis::Bezier){
    if constexpr(tdim == 2){
      eval2_d_bezier<MeshArray2D<SANS::SurrealS<nvar,double>>,szfld,ideg,nvar>(rfld, idif1, idif2, bary, eval, jmat, hmat);
    }else{
      eval3_d_bezier<MeshArray2D<SANS::SurrealS<nvar,double>>,szfld,ideg,nvar>(rfld, idif1, idif2, bary, eval, jmat, hmat);
    }
  }else{
    METRIS_THROW_MSG(TODOExcept(),"Implement eval3_d_LAGRANGE");
  } 
  return;
}







// SANS::SurrealS-based version can only handle Bézier for now. 
// Slower in all cases except P2 + Jacobian. 
template <int szfld, int tdim, int ideg, int ivar, int nvar = szfld>
void eval_d_SurrealS(const dblAr2 & __restrict__ rfld,  
                     const int * __restrict__ lfld,
                     FEBasis ibasis, DifVar idif1, DifVar idif2, 
                     const double * __restrict__ bary,
                     double * __restrict__  eval,
                     double * __restrict__  jmat,
                     double * __restrict__  hmat,
                     double * __restrict__  deval,
                     double * __restrict__  djmat,
                     double * __restrict__  dhmat,  
                     const double * __restrict__ dfld = NULL){


  // -- For clarity. 
  // Number of (r|l)fld entries
  constexpr int nrfld = getnnode(tdim,ideg);
  static_assert(nrfld < 31);
  // -

  static_assert(ivar >= 0 && ivar < nrfld);

  auto rfld_0 = hana::replicate<hana::tuple_tag>((const double *) 1, hana::size_c<nrfld>);

  auto rfld_1 = to_std_tuple(replace_at_c<ivar>(rfld_0,(SANS::SurrealS<nvar,double>*) 0));
  //// Store non-dof rfld entries
  //double rmem[szfld*(nrfld-1)];
  // Store the dof as SANS::SurrealS
  SANS::SurrealS<nvar,double> smem[szfld];

  if(dfld == NULL){
    METRIS_ASSERT(nvar == szfld);
    for(int icmp = 0; icmp < szfld; icmp++){
      smem[icmp].value() = rfld[lfld[ivar]][icmp];
      for(int j = 0; j < nvar; j++){
        smem[icmp].deriv(j) = 0;
      }
      smem[icmp].deriv(icmp) = 1;
    }
  }else{
    for(int icmp = 0; icmp < szfld; icmp++){
      smem[icmp].value() = rfld[lfld[ivar]][icmp];
      for(int idof = 0; idof < nvar; idof++){
        smem[icmp].deriv(idof) = dfld[idof*szfld + icmp];
      }
    }
  }

  //int imem = 0;
  //for(int ipoly = 0; ipoly < entnpps[ideg]; ipoly++){
  //  if(ipoly == ivar) continue;
  //  for(int icmp = 0; icmp < szfld; icmp++){
  //    rmem[imem*szfld + icmp] = rfld[lfld[ipoly]][icmp];
  //  }
  //  imem++;
  //}

  tuple_wrapper w_op(rfld_1);

  // Populate w_op tuple with appropriate types. 
  // The idea is to get something like a tuple<double*, SANS::SurrealS*, double*, double*...>
  // Storing, for the double*'s, the non-dof rfld entries and, for the one SANS::SurrealS*, 
  // the DoF SANS::SurrealS array. 
  //imem = 0;
  hana::while_(hana::less.than(hana::int_c<nrfld>), boost::hana::int_c<0>, [&](auto i_c){
    constexpr int i = i_c;
    if constexpr(i == ivar){
      w_op[i_c] = smem;
    }else{
//      w_op[i_c] = &rmem[szfld*imem];
      w_op[i_c] = rfld[lfld[i]];
//      imem++;
    }
  return i_c+boost::hana::int_c<1>;});  
//  SANS::SurrealS<nvar,double> seval[szfld];
  constexpr int nnsym = (tdim*(tdim+1))/2;

  SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> seval;
  SANS::DLA::MatrixS<tdim ,szfld,SANS::SurrealS<nvar,double>> sjmat;
//  SANS::SurrealS<nvar,double> sjmat[tdim3*szfld];
  SANS::DLA::MatrixS<nnsym,szfld,SANS::SurrealS<nvar,double>> shmat;
//  SANS::SurrealS<nvar,double> shmat[sdim3*szfld];


  //hana::while_(hana::less.than(hana::int_c<nrfld>), boost::hana::int_c<0>, [&](auto i_c){
  //  constexpr int i = i_c;
  //  printf("Debug rfld[%d] = \n",i);
  //  std::cout<<"0:"<<w_op[i_c][0]<<std::endl;
  //  std::cout<<"1:"<<w_op[i_c][1]<<std::endl;
  //  std::cout<<"2:"<<w_op[i_c][2]<<std::endl;
  //return i_c+boost::hana::int_c<1>;});  


  eval_d_SurrealS0<decltype(w_op),szfld,tdim,ideg>
                     (w_op,ibasis,idif1,idif2,bary,&seval,&sjmat,&shmat);

  for(int icmp=0; icmp<szfld; icmp++){
    eval[icmp] = seval[icmp].value();

    for(int idof = 0; idof < nvar; idof++){
      deval[szfld*idof + icmp] = seval[icmp].deriv(idof);
    }

    if(idif1 == DifVar::Bary){
      for(int itdim = 0; itdim < tdim; itdim++){
        jmat[itdim*szfld + icmp] = sjmat(itdim,icmp).value();
        for(int idof = 0; idof < nvar; idof++){
          djmat[idof*tdim*szfld + itdim*szfld + icmp] = sjmat(itdim,icmp).deriv(idof);
        } 
      }
    }

    if(idif2 == DifVar::Bary){
      for(int ii = 0; ii < nnsym; ii++){
        hmat[ii*szfld + icmp] = shmat(ii,icmp).value();
      }
      for(int idof = 0; idof < tdim; idof++){
        for(int ii = 0; ii < nnsym; ii++){
          dhmat[idof*nnsym*szfld + 0*szfld + icmp] = shmat(ii,icmp).deriv(idof);
        }
      }
    }// if idif2

  }
}



// SANS::SurrealS-based version can only handle Bézier for now. 
// Slower in all cases except P2 + Jacobian. 
template <int szfld, int tdim, int ideg, int nvar>
void eval_d_SurrealS(const dblAr2 & __restrict__ rfld,  
                     const int * __restrict__ lfld,
                     FEBasis ibasis, DifVar idif1, DifVar idif2, 
                     const double * __restrict__ bary,
                     int ivar_, 
                     SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>>* eval,
                     SANS::DLA::MatrixS<tdim ,szfld,SANS::SurrealS<nvar,double>>* jmat,
                     SANS::DLA::MatrixS<(tdim*(tdim+1))/2,szfld,SANS::SurrealS<nvar,double>>* hmat,
                     const SANS::DLA::MatrixS<szfld,nvar,double>* dfld = NULL){


  // -- For clarity. 
  // Number of (r|l)fld entries
  constexpr int nrfld = getnnode(tdim,ideg);
  static_assert(nrfld < 31);
  METRIS_ASSERT(ivar_ >= 0 && ivar_ < nrfld);


  auto rfld_0 = hana::replicate<hana::tuple_tag>((const double *) 1, hana::size_c<nrfld>);

  CT_FOR0_EXC(0, nrfld, ivar){if(ivar == ivar_){

    auto rfld_1 = to_std_tuple(replace_at_c<ivar>(rfld_0,(SANS::SurrealS<nvar,double>*) 0));
    // Store the dof as SANS::SurrealS
    SANS::SurrealS<nvar,double> smem[szfld];

    if(dfld == NULL){
      METRIS_ASSERT(nvar == szfld);
      // This is required to silence compiler warnings (iteration 2 invokes undefined behavior)
      if constexpr (nvar == szfld){
        for(int icmp = 0; icmp < szfld; icmp++){
          smem[icmp].value() = rfld[lfld[ivar]][icmp];
          for(int j = 0; j < nvar; j++){
            smem[icmp].deriv(j) = 0;
          }
          smem[icmp].deriv(icmp) = 1;
        }
      }
    }else{
      for(int icmp = 0; icmp < szfld; icmp++){
        smem[icmp].value() = rfld[lfld[ivar]][icmp];
        for(int idof = 0; idof < nvar; idof++){
          smem[icmp].deriv(idof) = (*dfld)(idof,icmp);
        }
      }
    }

    tuple_wrapper w_op(rfld_1);

    // Populate w_op tuple with appropriate types. 
    // The idea is to get something like a tuple<double*, SANS::SurrealS*, double*, double*...>
    // Storing, for the double*'s, the non-dof rfld entries and, for the one SANS::SurrealS*, 
    // the DoF SANS::SurrealS array. 
    // Replace c_ii by ii in this loop to bathe in compiler errors
    CT_FOR0_EXC(0, nrfld, ii){
      if constexpr(ii == ivar){
        w_op[c_ii] = smem;
      }else{
        w_op[c_ii] = rfld[lfld[ii]];
      }
    }CT_FOR1(ii); 

    eval_d_SurrealS0<decltype(w_op),szfld,tdim,ideg>
                       (w_op,ibasis,idif1,idif2,bary,eval,jmat,hmat);

  }}CT_FOR1(ivar);
}



// Unit tested: does not work currently (May 2025)
// Same as eval_d_SurrealS 
// but instead of computing szfld, compute size 1, then notice whole vector is repeated
template <int szfld, int tdim, int ideg,  int ivar, int nvar>
void eval_d_SurrealS_bcast(const dblAr2 & __restrict__ rfld,  
                     const int * __restrict__ lfld,
                     FEBasis ibasis, DifVar idif1, DifVar idif2, 
                     const double * __restrict__ bary,
                     double * __restrict__  eval,
                     double * __restrict__  jmat,
                     double * __restrict__  hmat,
                     double * __restrict__  deval,
                     double * __restrict__  djmat,
                     double * __restrict__  dhmat,  
                     const double * __restrict__ dfld = NULL){


  // -- For clarity. 
  // Number of (r|l)fld entries
  constexpr int nrfld = getnnode(tdim,ideg);
  static_assert(nrfld < 31);
  // -

  static_assert(ivar >= 0 && ivar < nrfld);

  auto rfld_0 = hana::replicate<hana::tuple_tag>((const double *) 1, hana::size_c<nrfld>);

  auto rfld_1 = to_std_tuple(replace_at_c<ivar>(rfld_0,(SANS::SurrealS<nvar,double>*) 0));
  // Store non-dof rfld entries
//  double rmem[szfld*(nrfld-1)];
  // Store the dof as SANS::SurrealS
  SANS::SurrealS<nvar,double> smem[szfld];

  if(dfld == NULL){
    METRIS_ASSERT(nvar == 1);
    for(int icmp = 0; icmp < szfld; icmp++){
      smem[icmp].value() = rfld[lfld[ivar]][icmp];
      smem[icmp].deriv(0) = 0;
    }
    smem[0].deriv(0) = 1; // They will all be the same
  }else{
    for(int icmp = 0; icmp < szfld; icmp++){
      smem[icmp].value() = rfld[lfld[ivar]][icmp];
      for(int idof = 0; idof < nvar; idof++){
        smem[icmp].deriv(idof) = dfld[idof*szfld + icmp];
      }
    }
  }

  //int imem = 0;
  //for(int ipoly = 0; ipoly < entnpps[ideg]; ipoly++){
  //  if(ipoly == ivar) continue;
  //  for(int icmp = 0; icmp < szfld; icmp++){
  //    rmem[imem*szfld + icmp] = rfld[lfld[ipoly]][icmp];
  //  }
  //  imem++;
  //}

  tuple_wrapper w_op(rfld_1);

  // Populate w_op tuple with appropriate types. 
  // The idea is to get something like a tuple<double*, SANS::SurrealS*, double*, double*...>
  // Storing, for the double*'s, the non-dof rfld entries and, for the one SANS::SurrealS*, 
  // the DoF SANS::SurrealS array. 
  //imem = 0;
  hana::while_(hana::less.than(hana::int_c<nrfld>), boost::hana::int_c<0>, [&](auto i_c){
    constexpr int i = i_c;
    if constexpr(i == ivar){
      w_op[i_c] = smem;
    }else{
//      w_op[i_c] = &rmem[szfld*imem];
      w_op[i_c] = rfld[lfld[i]];
//      imem++;
    }
  return i_c+boost::hana::int_c<1>;});  
//  SANS::SurrealS<nvar,double> seval[szfld];
  constexpr int nnsym = (tdim*(tdim+1))/2;

  SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> seval;
  SANS::DLA::MatrixS<tdim ,szfld,SANS::SurrealS<nvar,double>> sjmat;
//  SANS::SurrealS<nvar,double> sjmat[tdim3*szfld];
  SANS::DLA::MatrixS<nnsym,szfld,SANS::SurrealS<nvar,double>> shmat;
//  SANS::SurrealS<nvar,double> shmat[sdim3*szfld];


  //hana::while_(hana::less.than(hana::int_c<nrfld>), boost::hana::int_c<0>, [&](auto i_c){
  //  constexpr int i = i_c;
  //  printf("Debug rfld[%d] = \n",i);
  //  std::cout<<"0:"<<w_op[i_c][0]<<std::endl;
  //  std::cout<<"1:"<<w_op[i_c][1]<<std::endl;
  //  std::cout<<"2:"<<w_op[i_c][2]<<std::endl;
  //return i_c+boost::hana::int_c<1>;});  


  eval_d_SurrealS0<decltype(w_op),szfld,tdim,ideg>
                     (w_op,ibasis,idif1,idif2,bary,&seval,&sjmat,&shmat);

  //std::cout<<"Debug seval = \n"<<seval<<"\n";

  //std::cout<<"Debug sjmat = \n"<<sjmat<<"\n";

  //std::cout<<"debug ideg = "<<ideg<< " bary = "; 
  //dblAr1(tdim+1,bary).print();



  int nvareff = dfld == NULL ? szfld : nvar;

  for(int icmp=0; icmp<szfld; icmp++){
    eval[icmp] = seval[icmp].value();

    if(dfld != NULL){
      for(int idof = 0; idof < nvareff; idof++){
        deval[szfld*idof + icmp] = seval[icmp].deriv(idof);
      }
      if(idif1 == DifVar::Bary){
        for(int itdim = 0; itdim < tdim; itdim++){
          jmat[itdim*szfld + icmp] = sjmat(itdim,icmp).value();
          for(int idof = 0; idof < nvareff; idof++){
            djmat[idof*tdim*szfld + itdim*szfld + icmp] = sjmat(itdim,icmp).deriv(idof);
          } 
        }
      }
    }else{
      double de = seval[0].deriv(0);
      for(int idof = 0; idof < nvareff; idof++){
        deval[szfld*idof + icmp] = de * (idof == icmp);
      }
      if(idif1 == DifVar::Bary){
        double dj = sjmat(0,0).deriv(0);
        for(int itdim = 0; itdim < tdim; itdim++){
          jmat[itdim*szfld + icmp] = sjmat(itdim,icmp).value();
          for(int idof = 0; idof < nvareff; idof++){
            djmat[idof*tdim*szfld + itdim*szfld + icmp] = dj * (idof == icmp);
          } 
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
        for(int ii = 0; ii < nnsym; ii++) dhmat[idof*nnsym*szfld + 0*szfld + icmp] = shmat(ii,icmp).deriv(0);
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







// This version just makes one large SurrealS buffer instead of trying to limit it to 1 variable
template <int szfld, int tdim, int ideg, int nvar>
void eval_d_SurrealS_simple(const dblAr2 & __restrict__ rfld,  
                            const int * __restrict__ lfld,
                            FEBasis ibasis, DifVar idif1, DifVar idif2, 
                            const double * __restrict__ bary,
                            int ivar, 
                            double * __restrict__  eval,
                            double * __restrict__  jmat,
                            double * __restrict__  hmat,
                            double * __restrict__  deval,
                            double * __restrict__  djmat,
                            double * __restrict__  dhmat,  
                            const double * __restrict__ dfld = NULL){


  // -- For clarity. 
  // Number of (r|l)fld entries
  constexpr int nrfld = getnnode(tdim,ideg);
  // -


  SANS::SurrealS<nvar,double> buffer[nrfld][szfld];
  MeshArray2D<SANS::SurrealS<nvar,double>> rfllS(nrfld,szfld,buffer[0]);

  if(dfld == NULL){
    METRIS_ASSERT(nvar == szfld);
  }
  for(int ii = 0; ii < nrfld; ii++){
    for(int jj = 0; jj < szfld; jj++){
      rfllS[ii][jj].value() = rfld[lfld[ii]][jj];
      for(int kk = 0; kk < nvar; kk++){
        rfllS[ii][jj].deriv(kk) = 0;
      }
      if(ii == ivar) rfllS[ii][jj].deriv(jj) = 1;
    }
  }
  if(dfld != NULL){
    for(int icmp = 0; icmp < szfld; icmp++){
      rfllS[ivar][icmp].value() = rfld[lfld[ivar]][icmp];
      for(int idof = 0; idof < nvar; idof++){
        rfllS[ivar][icmp].deriv(idof) = dfld[idof*szfld + icmp];
      }
    }
  }




  constexpr int nnsym = (tdim*(tdim+1))/2;

  SANS::DLA::VectorS<      szfld,SANS::SurrealS<nvar,double>> seval;
  SANS::DLA::MatrixS<tdim ,szfld,SANS::SurrealS<nvar,double>> sjmat;
  SANS::DLA::MatrixS<nnsym,szfld,SANS::SurrealS<nvar,double>> shmat;

  eval_d_SurrealS0_simple<szfld,tdim,ideg>
                           (rfllS,ibasis,idif1,idif2,bary,&seval,&sjmat,&shmat);

  for(int icmp=0; icmp<szfld; icmp++){
    eval[icmp] = seval[icmp].value();

    for(int idof = 0; idof < nvar; idof++){
      deval[szfld*idof + icmp] = seval[icmp].deriv(idof);
    }

    if(idif1 == DifVar::Bary){
      for(int itdim = 0; itdim < tdim; itdim++){
        jmat[itdim*szfld + icmp] = sjmat(itdim,icmp).value();
        for(int idof = 0; idof < nvar; idof++){
          djmat[idof*tdim*szfld + itdim*szfld + icmp] = sjmat(itdim,icmp).deriv(idof);
        } 
      }
    }

    if(idif2 == DifVar::Bary){
      for(int ii = 0; ii < nnsym; ii++) hmat[ii*szfld + icmp] = shmat(ii,icmp).value();
      for(int idof = 0; idof < tdim; idof++){
        for(int ii = 0; ii < nnsym; ii++) dhmat[idof*nnsym*szfld + 0*szfld + icmp] = shmat(ii,icmp).deriv(idof);
      }
    }

  }
}




} // End namespace
#endif
