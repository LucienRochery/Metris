//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __CODEGEN_LAGRANGE__
#define __CODEGEN_LAGRANGE__
#include "types.hxx"
namespace Metris{


//------------------- Edges --------------------

template <int szfld, int ideg> 
struct eval1_lagrange{
  eval1_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1,
                 const double * __restrict__  bary,
                 double * __restrict__  eval, double * __restrict__  jmat);
};

template <int szfld>
struct eval1_lagrange<szfld,1>{
  eval1_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1,
                 const double * __restrict__  bary,
                 double * __restrict__  eval, double * __restrict__  jmat){
  for(int j=0; j < szfld; j++){
    eval[j] = rfld[lfld[0]][j]*(bary[0])
            + rfld[lfld[1]][j]*(bary[1]);
  }
  if (idif1 != DifVar::Bary) return;
  for(int j=0; j < szfld; j++){
    jmat[0*szfld+j] =  rfld[lfld[0]][j]*(-1.0)
                    +  rfld[lfld[1]][j]*(1.0);
  }
}
};
template <int szfld>
struct eval1_lagrange<szfld,2>{
  eval1_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1,
                 const double * __restrict__  bary,
                 double * __restrict__  eval, double * __restrict__  jmat){
  for(int j=0; j < szfld; j++){
    eval[j] = rfld[lfld[0]][j]*(bary[0]*( 2.0*bary[0]-1.0))
            + rfld[lfld[1]][j]*(bary[1]*( 2.0*bary[1]-1.0))
            + rfld[lfld[2]][j]*(4.0*bary[1]*bary[0]);
  }
  if (idif1 != DifVar::Bary) return;
  for(int j=0; j < szfld; j++){
    jmat[0*szfld+j] =  rfld[lfld[0]][j]*( -4.0*bary[0]+1.0)
                    +  rfld[lfld[1]][j]*( 4.0*bary[1]-1.0)
                    +  rfld[lfld[2]][j]*( -4.0*bary[1]+4.0*bary[0]);
  }
}
};
template <int szfld>
struct eval1_lagrange<szfld,3>{
  eval1_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1,
                 const double * __restrict__  bary,
                 double * __restrict__  eval, double * __restrict__  jmat){
  for(int j=0; j < szfld; j++){
    eval[j] = rfld[lfld[0]][j]*((1.0/2.0)*( 3.0*bary[0]-2.0)*bary[0]*( 3.0*bary[0]-1.0))
            + rfld[lfld[1]][j]*((1.0/2.0)*( 3.0*bary[1]-2.0)*( 3.0*bary[1]-1.0)*bary[1])
            + rfld[lfld[2]][j]*((9.0/2.0)*( 3.0*bary[0]-1.0)*bary[0]*bary[1])
            + rfld[lfld[3]][j]*((9.0/2.0)*( 3.0*bary[1]-1.0)*bary[0]*bary[1]);
  }
  if (idif1 != DifVar::Bary) return;
  for(int j=0; j < szfld; j++){
    jmat[0*szfld+j] =  rfld[lfld[0]][j]*( 9.0*bary[0]+-(27.0/2.0)*(bary[0]*bary[0])-1.0)
                    +  rfld[lfld[1]][j]*( -9.0*bary[1]+(27.0/2.0)*(bary[1]*bary[1])+1.0)
                    +  rfld[lfld[2]][j]*( (27.0/2.0)*(bary[0]*bary[0])+-(9.0/2.0)*( 6.0*bary[0]-1.0)*bary[1]+-(9.0/2.0)*bary[0])
                    +  rfld[lfld[3]][j]*( -(27.0/2.0)*(bary[1]*bary[1])+-(9.0/2.0)*bary[0]+(9.0/2.0)*( 6.0*bary[0]+1.0)*bary[1]);
  }
}
};

//------------------- Triangles --------------------

template <int szfld, int ideg> 
struct eval2_lagrange{
  eval2_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1,
                 const double * __restrict__  bary,
                 double * __restrict__  eval, double * __restrict__  jmat);
};

template <int szfld>
struct eval2_lagrange<szfld,1>{
  eval2_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1,
                 const double * __restrict__  bary,
                 double * __restrict__  eval, double * __restrict__  jmat){
  for(int j=0; j < szfld; j++){
    eval[j] = rfld[lfld[0]][j]*(bary[0])
            + rfld[lfld[1]][j]*(bary[1])
            + rfld[lfld[2]][j]*(bary[2]);
  }
  if (idif1 != DifVar::Bary) return;
  for(int j=0; j < szfld; j++){
    jmat[0*szfld+j] =  rfld[lfld[0]][j]*(-1.0)
                    +  rfld[lfld[1]][j]*(1.0);
  }
  for(int j=0; j < szfld; j++){
    jmat[1*szfld+j] =  rfld[lfld[0]][j]*(-1.0)
                    +  rfld[lfld[2]][j]*(1.0);
  }
}
};
template <int szfld>
struct eval2_lagrange<szfld,2>{
  eval2_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1,
                 const double * __restrict__  bary,
                 double * __restrict__  eval, double * __restrict__  jmat){
  for(int j=0; j < szfld; j++){
    eval[j] = rfld[lfld[0]][j]*(bary[0]*( 2.0*bary[0]-1.0))
            + rfld[lfld[1]][j]*(( 2.0*bary[1]-1.0)*bary[1])
            + rfld[lfld[2]][j]*(bary[2]*( 2.0*bary[2]-1.0))
            + rfld[lfld[3]][j]*(4.0*bary[1]*bary[2])
            + rfld[lfld[4]][j]*(4.0*bary[2]*bary[0])
            + rfld[lfld[5]][j]*(4.0*bary[1]*bary[0]);
  }
  if (idif1 != DifVar::Bary) return;
  for(int j=0; j < szfld; j++){
    jmat[0*szfld+j] =  rfld[lfld[0]][j]*( -4.0*bary[0]+1.0)
                    +  rfld[lfld[1]][j]*( 4.0*bary[1]-1.0)
                    +  rfld[lfld[3]][j]*(4.0*bary[2])
                    +  rfld[lfld[4]][j]*(-4.0*bary[2])
                    +  rfld[lfld[5]][j]*( -4.0*bary[1]+4.0*bary[0]);
  }
  for(int j=0; j < szfld; j++){
    jmat[1*szfld+j] =  rfld[lfld[0]][j]*( -4.0*bary[0]+1.0)
                    +  rfld[lfld[2]][j]*( 4.0*bary[2]-1.0)
                    +  rfld[lfld[3]][j]*(4.0*bary[1])
                    +  rfld[lfld[4]][j]*( -4.0*bary[2]+4.0*bary[0])
                    +  rfld[lfld[5]][j]*(-4.0*bary[1]);
  }
}
};
template <int szfld>
struct eval2_lagrange<szfld,3>{
  eval2_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1,
                 const double * __restrict__  bary,
                 double * __restrict__  eval, double * __restrict__  jmat){
  for(int j=0; j < szfld; j++){
    eval[j] = rfld[lfld[0]][j]*((1.0/2.0)*( 3.0*bary[0]-1.0)*( 3.0*bary[0]-2.0)*bary[0])
            + rfld[lfld[1]][j]*((1.0/2.0)*( 3.0*bary[1]-2.0)*bary[1]*( 3.0*bary[1]-1.0))
            + rfld[lfld[2]][j]*((1.0/2.0)*bary[2]*( 3.0*bary[2]-1.0)*( 3.0*bary[2]-2.0))
            + rfld[lfld[3]][j]*((9.0/2.0)*( 3.0*bary[1]-1.0)*bary[2]*bary[1])
            + rfld[lfld[4]][j]*((9.0/2.0)*( 3.0*bary[2]-1.0)*bary[1]*bary[2])
            + rfld[lfld[5]][j]*((9.0/2.0)*( 3.0*bary[2]-1.0)*bary[0]*bary[2])
            + rfld[lfld[6]][j]*((9.0/2.0)*bary[2]*bary[0]*( 3.0*bary[0]-1.0))
            + rfld[lfld[7]][j]*((9.0/2.0)*bary[0]*bary[1]*( 3.0*bary[0]-1.0))
            + rfld[lfld[8]][j]*((9.0/2.0)*bary[0]*( 3.0*bary[1]-1.0)*bary[1])
            + rfld[lfld[9]][j]*(27.0*bary[2]*bary[0]*bary[1]);
  }
  if (idif1 != DifVar::Bary) return;
  for(int j=0; j < szfld; j++){
    jmat[0*szfld+j] =  rfld[lfld[0]][j]*( -(27.0/2.0)*(bary[0]*bary[0])+9.0*bary[0]-1.0)
                    +  rfld[lfld[1]][j]*( -9.0*bary[1]+(27.0/2.0)*(bary[1]*bary[1])+1.0)
                    +  rfld[lfld[3]][j]*((9.0/2.0)*bary[2]*( 6.0*bary[1]-1.0))
                    +  rfld[lfld[4]][j]*((9.0/2.0)*( 3.0*bary[2]-1.0)*bary[2])
                    +  rfld[lfld[5]][j]*(-(9.0/2.0)*( 3.0*bary[2]-1.0)*bary[2])
                    +  rfld[lfld[6]][j]*(-(9.0/2.0)*( 6.0*bary[0]-1.0)*bary[2])
                    +  rfld[lfld[7]][j]*( -(9.0/2.0)*bary[0]+(27.0/2.0)*(bary[0]*bary[0])+-(9.0/2.0)*( 6.0*bary[0]-1.0)*bary[1])
                    +  rfld[lfld[8]][j]*( -(27.0/2.0)*(bary[1]*bary[1])+-(9.0/2.0)*bary[0]+(9.0/2.0)*( 6.0*bary[0]+1.0)*bary[1])
                    +  rfld[lfld[9]][j]*(27.0*bary[2]*( bary[0]-bary[1]));
  }
  for(int j=0; j < szfld; j++){
    jmat[1*szfld+j] =  rfld[lfld[0]][j]*( -(27.0/2.0)*(bary[0]*bary[0])+9.0*bary[0]-1.0)
                    +  rfld[lfld[2]][j]*( (27.0/2.0)*(bary[2]*bary[2])+-9.0*bary[2]+1.0)
                    +  rfld[lfld[3]][j]*((9.0/2.0)*( 3.0*bary[1]-1.0)*bary[1])
                    +  rfld[lfld[4]][j]*((9.0/2.0)*( 6.0*bary[2]-1.0)*bary[1])
                    +  rfld[lfld[5]][j]*( -(27.0/2.0)*(bary[2]*bary[2])+-(9.0/2.0)*bary[0]+(9.0/2.0)*( 6.0*bary[0]+1.0)*bary[2])
                    +  rfld[lfld[6]][j]*( -(9.0/2.0)*bary[0]*( 6.0*bary[2]+1.0)+(9.0/2.0)*bary[2]+(27.0/2.0)*(bary[0]*bary[0]))
                    +  rfld[lfld[7]][j]*(-(9.0/2.0)*( 6.0*bary[0]-1.0)*bary[1])
                    +  rfld[lfld[8]][j]*(-(9.0/2.0)*( 3.0*bary[1]-1.0)*bary[1])
                    +  rfld[lfld[9]][j]*(-27.0*( bary[2]-bary[0])*bary[1]);
  }
}
};

//------------------- Tetrahedra --------------------

template <int szfld, int ideg> 
struct eval3_lagrange{
  eval3_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1, DifVar idif2,
                 const double * __restrict__  bary,
                 double * __restrict__  eval,
                 double * __restrict__  jmat,
                 double * __restrict__  hmat);
};

template <int szfld>
struct eval3_lagrange<szfld,1>{
  eval3_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1, DifVar idif2,
                 const double * __restrict__  bary,
                 double * __restrict__  eval,
                 double * __restrict__  jmat,
                 double * __restrict__  hmat){
  for(int j=0; j < szfld; j++){
    eval[j] = rfld[lfld[0]][j]*(bary[0])
            + rfld[lfld[1]][j]*(bary[1])
            + rfld[lfld[2]][j]*(bary[2])
            + rfld[lfld[3]][j]*(bary[3]);
  }
  if (idif1 == DifVar::Bary){
    for(int j=0; j < szfld; j++){
      jmat[0*szfld+j] =  rfld[lfld[0]][j]*(-1.0)
                      +  rfld[lfld[1]][j]*(1.0);
    }
    for(int j=0; j < szfld; j++){
      jmat[1*szfld+j] =  rfld[lfld[0]][j]*(-1.0)
                      +  rfld[lfld[2]][j]*(1.0);
    }
    for(int j=0; j < szfld; j++){
      jmat[2*szfld+j] =  rfld[lfld[0]][j]*(-1.0)
                      +  rfld[lfld[3]][j]*(1.0);
    }
  }
  if (idif2 == DifVar::Bary){
    for(int j=0; j < szfld; j++){
      hmat[0*szfld+j] = 0;
      hmat[1*szfld+j] = 0;
      hmat[2*szfld+j] = 0;
      hmat[3*szfld+j] = 0;
      hmat[4*szfld+j] = 0;
      hmat[5*szfld+j] = 0;
    }
  }
}
};
template <int szfld>
struct eval3_lagrange<szfld,2>{
  eval3_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1, DifVar idif2,
                 const double * __restrict__  bary,
                 double * __restrict__  eval,
                 double * __restrict__  jmat,
                 double * __restrict__  hmat){
  for(int j=0; j < szfld; j++){
    eval[j] = rfld[lfld[0]][j]*(( 2.0*bary[0]-1.0)*bary[0])
            + rfld[lfld[1]][j]*(( 2.0*bary[1]-1.0)*bary[1])
            + rfld[lfld[2]][j]*(bary[2]*( 2.0*bary[2]-1.0))
            + rfld[lfld[3]][j]*(bary[3]*( 2.0*bary[3]-1.0))
            + rfld[lfld[4]][j]*(4.0*bary[0]*bary[1])
            + rfld[lfld[5]][j]*(4.0*bary[2]*bary[1])
            + rfld[lfld[6]][j]*(4.0*bary[2]*bary[0])
            + rfld[lfld[7]][j]*(4.0*bary[0]*bary[3])
            + rfld[lfld[8]][j]*(4.0*bary[3]*bary[1])
            + rfld[lfld[9]][j]*(4.0*bary[2]*bary[3]);
  }
  if (idif1 == DifVar::Bary){
    for(int j=0; j < szfld; j++){
      jmat[0*szfld+j] =  rfld[lfld[0]][j]*( -4.0*bary[0]+1.0)
                      +  rfld[lfld[1]][j]*( 4.0*bary[1]-1.0)
                      +  rfld[lfld[4]][j]*( 4.0*bary[0]+-4.0*bary[1])
                      +  rfld[lfld[5]][j]*(4.0*bary[2])
                      +  rfld[lfld[6]][j]*(-4.0*bary[2])
                      +  rfld[lfld[7]][j]*(-4.0*bary[3])
                      +  rfld[lfld[8]][j]*(4.0*bary[3]);
    }
    for(int j=0; j < szfld; j++){
      jmat[1*szfld+j] =  rfld[lfld[0]][j]*( -4.0*bary[0]+1.0)
                      +  rfld[lfld[2]][j]*( 4.0*bary[2]-1.0)
                      +  rfld[lfld[4]][j]*(-4.0*bary[1])
                      +  rfld[lfld[5]][j]*(4.0*bary[1])
                      +  rfld[lfld[6]][j]*( -4.0*bary[2]+4.0*bary[0])
                      +  rfld[lfld[7]][j]*(-4.0*bary[3])
                      +  rfld[lfld[9]][j]*(4.0*bary[3]);
    }
    for(int j=0; j < szfld; j++){
      jmat[2*szfld+j] =  rfld[lfld[0]][j]*( -4.0*bary[0]+1.0)
                      +  rfld[lfld[3]][j]*( 4.0*bary[3]-1.0)
                      +  rfld[lfld[4]][j]*(-4.0*bary[1])
                      +  rfld[lfld[6]][j]*(-4.0*bary[2])
                      +  rfld[lfld[7]][j]*( 4.0*bary[0]+-4.0*bary[3])
                      +  rfld[lfld[8]][j]*(4.0*bary[1])
                      +  rfld[lfld[9]][j]*(4.0*bary[2]);
    }
  }
  if (idif2 == DifVar::Bary){
    for(int j=0; j < szfld; j++){
      hmat[0*szfld+j] =  rfld[lfld[0]][j]*(4.0)
                      +  rfld[lfld[1]][j]*(4.0)
                      +  rfld[lfld[4]][j]*(-8.0);
    }
    for(int j=0; j < szfld; j++){
      hmat[1*szfld+j] =  rfld[lfld[0]][j]*(4.0)
                      +  rfld[lfld[4]][j]*(-4.0)
                      +  rfld[lfld[5]][j]*(4.0)
                      +  rfld[lfld[6]][j]*(-4.0);
    }
    for(int j=0; j < szfld; j++){
      hmat[2*szfld+j] =  rfld[lfld[0]][j]*(4.0)
                      +  rfld[lfld[2]][j]*(4.0)
                      +  rfld[lfld[6]][j]*(-8.0);
    }
    for(int j=0; j < szfld; j++){
      hmat[3*szfld+j] =  rfld[lfld[0]][j]*(4.0)
                      +  rfld[lfld[4]][j]*(-4.0)
                      +  rfld[lfld[7]][j]*(-4.0)
                      +  rfld[lfld[8]][j]*(4.0);
    }
    for(int j=0; j < szfld; j++){
      hmat[4*szfld+j] =  rfld[lfld[0]][j]*(4.0)
                      +  rfld[lfld[6]][j]*(-4.0)
                      +  rfld[lfld[7]][j]*(-4.0)
                      +  rfld[lfld[9]][j]*(4.0);
    }
    for(int j=0; j < szfld; j++){
      hmat[5*szfld+j] =  rfld[lfld[0]][j]*(4.0)
                      +  rfld[lfld[3]][j]*(4.0)
                      +  rfld[lfld[7]][j]*(-8.0);
    }
  }
}
};
template <int szfld>
struct eval3_lagrange<szfld,3>{
  eval3_lagrange(const dblAr2 & __restrict__  rfld,
                 const int * __restrict__  lfld,
                 DifVar idif1, DifVar idif2,
                 const double * __restrict__  bary,
                 double * __restrict__  eval,
                 double * __restrict__  jmat,
                 double * __restrict__  hmat){
  for(int j=0; j < szfld; j++){
    eval[j] = rfld[lfld[0]][j]*((1.0/2.0)*( 3.0*bary[0]-2.0)*bary[0]*( 3.0*bary[0]-1.0))
            + rfld[lfld[1]][j]*((1.0/2.0)*( 3.0*bary[1]-1.0)*bary[1]*( 3.0*bary[1]-2.0))
            + rfld[lfld[2]][j]*((1.0/2.0)*( 3.0*bary[2]-1.0)*bary[2]*( 3.0*bary[2]-2.0))
            + rfld[lfld[3]][j]*((1.0/2.0)*bary[3]*( 3.0*bary[3]-2.0)*( 3.0*bary[3]-1.0))
            + rfld[lfld[4]][j]*((9.0/2.0)*bary[0]*bary[1]*( 3.0*bary[0]-1.0))
            + rfld[lfld[5]][j]*((9.0/2.0)*bary[0]*bary[1]*( 3.0*bary[1]-1.0))
            + rfld[lfld[6]][j]*((9.0/2.0)*( 3.0*bary[1]-1.0)*bary[2]*bary[1])
            + rfld[lfld[7]][j]*((9.0/2.0)*bary[1]*bary[2]*( 3.0*bary[2]-1.0))
            + rfld[lfld[8]][j]*((9.0/2.0)*bary[0]*( 3.0*bary[2]-1.0)*bary[2])
            + rfld[lfld[9]][j]*((9.0/2.0)*bary[2]*( 3.0*bary[0]-1.0)*bary[0])
            + rfld[lfld[10]][j]*((9.0/2.0)*( 3.0*bary[0]-1.0)*bary[0]*bary[3])
            + rfld[lfld[11]][j]*((9.0/2.0)*( 3.0*bary[3]-1.0)*bary[0]*bary[3])
            + rfld[lfld[12]][j]*((9.0/2.0)*( 3.0*bary[1]-1.0)*bary[3]*bary[1])
            + rfld[lfld[13]][j]*((9.0/2.0)*( 3.0*bary[3]-1.0)*bary[3]*bary[1])
            + rfld[lfld[14]][j]*((9.0/2.0)*( 3.0*bary[2]-1.0)*bary[2]*bary[3])
            + rfld[lfld[15]][j]*((9.0/2.0)*( 3.0*bary[3]-1.0)*bary[2]*bary[3])
            + rfld[lfld[16]][j]*(27.0*bary[3]*bary[1]*bary[2])
            + rfld[lfld[17]][j]*(27.0*bary[0]*bary[2]*bary[3])
            + rfld[lfld[18]][j]*(27.0*bary[1]*bary[0]*bary[3])
            + rfld[lfld[19]][j]*(27.0*bary[0]*bary[1]*bary[2]);
  }
  if (idif1 == DifVar::Bary){
    for(int j=0; j < szfld; j++){
      jmat[0*szfld+j] =  rfld[lfld[0]][j]*( 9.0*bary[0]+-(27.0/2.0)*(bary[0]*bary[0])-1.0)
                      +  rfld[lfld[1]][j]*( (27.0/2.0)*(bary[1]*bary[1])+-9.0*bary[1]+1.0)
                      +  rfld[lfld[4]][j]*( -(9.0/2.0)*bary[0]+-(9.0/2.0)*( 6.0*bary[0]-1.0)*bary[1]+(27.0/2.0)*(bary[0]*bary[0]))
                      +  rfld[lfld[5]][j]*( -(9.0/2.0)*bary[0]+-(27.0/2.0)*(bary[1]*bary[1])+(9.0/2.0)*( 6.0*bary[0]+1.0)*bary[1])
                      +  rfld[lfld[6]][j]*((9.0/2.0)*bary[2]*( 6.0*bary[1]-1.0))
                      +  rfld[lfld[7]][j]*((9.0/2.0)*bary[2]*( 3.0*bary[2]-1.0))
                      +  rfld[lfld[8]][j]*(-(9.0/2.0)*( 3.0*bary[2]-1.0)*bary[2])
                      +  rfld[lfld[9]][j]*(-(9.0/2.0)*bary[2]*( 6.0*bary[0]-1.0))
                      +  rfld[lfld[10]][j]*(-(9.0/2.0)*( 6.0*bary[0]-1.0)*bary[3])
                      +  rfld[lfld[11]][j]*(-(9.0/2.0)*( 3.0*bary[3]-1.0)*bary[3])
                      +  rfld[lfld[12]][j]*((9.0/2.0)*bary[3]*( 6.0*bary[1]-1.0))
                      +  rfld[lfld[13]][j]*((9.0/2.0)*( 3.0*bary[3]-1.0)*bary[3])
                      +  rfld[lfld[16]][j]*(27.0*bary[3]*bary[2])
                      +  rfld[lfld[17]][j]*(-27.0*bary[2]*bary[3])
                      +  rfld[lfld[18]][j]*(-27.0*( bary[1]-bary[0])*bary[3])
                      +  rfld[lfld[19]][j]*(27.0*( bary[0]-bary[1])*bary[2]);
    }
    for(int j=0; j < szfld; j++){
      jmat[1*szfld+j] =  rfld[lfld[0]][j]*( 9.0*bary[0]+-(27.0/2.0)*(bary[0]*bary[0])-1.0)
                      +  rfld[lfld[2]][j]*( (27.0/2.0)*(bary[2]*bary[2])+-9.0*bary[2]+1.0)
                      +  rfld[lfld[4]][j]*(-(9.0/2.0)*( 6.0*bary[0]-1.0)*bary[1])
                      +  rfld[lfld[5]][j]*(-(9.0/2.0)*bary[1]*( 3.0*bary[1]-1.0))
                      +  rfld[lfld[6]][j]*((9.0/2.0)*( 3.0*bary[1]-1.0)*bary[1])
                      +  rfld[lfld[7]][j]*((9.0/2.0)*( 6.0*bary[2]-1.0)*bary[1])
                      +  rfld[lfld[8]][j]*( -(9.0/2.0)*bary[0]+(9.0/2.0)*( 6.0*bary[0]+1.0)*bary[2]+-(27.0/2.0)*(bary[2]*bary[2]))
                      +  rfld[lfld[9]][j]*( (27.0/2.0)*(bary[0]*bary[0])+(9.0/2.0)*bary[2]+-(9.0/2.0)*bary[0]*( 6.0*bary[2]+1.0))
                      +  rfld[lfld[10]][j]*(-(9.0/2.0)*( 6.0*bary[0]-1.0)*bary[3])
                      +  rfld[lfld[11]][j]*(-(9.0/2.0)*( 3.0*bary[3]-1.0)*bary[3])
                      +  rfld[lfld[14]][j]*((9.0/2.0)*bary[3]*( 6.0*bary[2]-1.0))
                      +  rfld[lfld[15]][j]*((9.0/2.0)*( 3.0*bary[3]-1.0)*bary[3])
                      +  rfld[lfld[16]][j]*(27.0*bary[3]*bary[1])
                      +  rfld[lfld[17]][j]*(27.0*( bary[0]-bary[2])*bary[3])
                      +  rfld[lfld[18]][j]*(-27.0*bary[1]*bary[3])
                      +  rfld[lfld[19]][j]*(27.0*bary[1]*( bary[0]-bary[2]));
    }
    for(int j=0; j < szfld; j++){
      jmat[2*szfld+j] =  rfld[lfld[0]][j]*( 9.0*bary[0]+-(27.0/2.0)*(bary[0]*bary[0])-1.0)
                      +  rfld[lfld[3]][j]*( -9.0*bary[3]+(27.0/2.0)*(bary[3]*bary[3])+1.0)
                      +  rfld[lfld[4]][j]*(-(9.0/2.0)*( 6.0*bary[0]-1.0)*bary[1])
                      +  rfld[lfld[5]][j]*(-(9.0/2.0)*bary[1]*( 3.0*bary[1]-1.0))
                      +  rfld[lfld[8]][j]*(-(9.0/2.0)*( 3.0*bary[2]-1.0)*bary[2])
                      +  rfld[lfld[9]][j]*(-(9.0/2.0)*bary[2]*( 6.0*bary[0]-1.0))
                      +  rfld[lfld[10]][j]*( (27.0/2.0)*(bary[0]*bary[0])+-(9.0/2.0)*bary[0]+-(9.0/2.0)*( 6.0*bary[0]-1.0)*bary[3])
                      +  rfld[lfld[11]][j]*( (9.0/2.0)*( 6.0*bary[0]+1.0)*bary[3]+-(27.0/2.0)*(bary[3]*bary[3])+-(9.0/2.0)*bary[0])
                      +  rfld[lfld[12]][j]*((9.0/2.0)*( 3.0*bary[1]-1.0)*bary[1])
                      +  rfld[lfld[13]][j]*((9.0/2.0)*bary[1]*( 6.0*bary[3]-1.0))
                      +  rfld[lfld[14]][j]*((9.0/2.0)*( 3.0*bary[2]-1.0)*bary[2])
                      +  rfld[lfld[15]][j]*((9.0/2.0)*bary[2]*( 6.0*bary[3]-1.0))
                      +  rfld[lfld[16]][j]*(27.0*bary[1]*bary[2])
                      +  rfld[lfld[17]][j]*(27.0*bary[2]*( bary[0]-bary[3]))
                      +  rfld[lfld[18]][j]*(27.0*( bary[0]-bary[3])*bary[1])
                      +  rfld[lfld[19]][j]*(-27.0*bary[1]*bary[2]);
    }
  }
  if (idif2 == DifVar::Bary){
    for(int j=0; j < szfld; j++){
      hmat[0*szfld+j] =  rfld[lfld[0]][j]*( 27.0*bary[0]-9.0)
                      +  rfld[lfld[1]][j]*( 27.0*bary[1]-9.0)
                      +  rfld[lfld[4]][j]*( -54.0*bary[0]+27.0*bary[1]+9.0)
                      +  rfld[lfld[5]][j]*( 27.0*bary[0]+-54.0*bary[1]+9.0)
                      +  rfld[lfld[6]][j]*(27.0*bary[2])
                      +  rfld[lfld[9]][j]*(27.0*bary[2])
                      +  rfld[lfld[10]][j]*(27.0*bary[3])
                      +  rfld[lfld[12]][j]*(27.0*bary[3])
                      +  rfld[lfld[18]][j]*(-54.0*bary[3])
                      +  rfld[lfld[19]][j]*(-54.0*bary[2]);
    }
    for(int j=0; j < szfld; j++){
      hmat[1*szfld+j] =  rfld[lfld[0]][j]*( 27.0*bary[0]-9.0)
                      +  rfld[lfld[4]][j]*( -27.0*bary[0]+27.0*bary[1]+(9.0/2.0))
                      +  rfld[lfld[5]][j]*( -27.0*bary[1]+(9.0/2.0))
                      +  rfld[lfld[6]][j]*( 27.0*bary[1]-(9.0/2.0))
                      +  rfld[lfld[7]][j]*( 27.0*bary[2]-(9.0/2.0))
                      +  rfld[lfld[8]][j]*( -27.0*bary[2]+(9.0/2.0))
                      +  rfld[lfld[9]][j]*( 27.0*bary[2]+-27.0*bary[0]+(9.0/2.0))
                      +  rfld[lfld[10]][j]*(27.0*bary[3])
                      +  rfld[lfld[16]][j]*(27.0*bary[3])
                      +  rfld[lfld[17]][j]*(-27.0*bary[3])
                      +  rfld[lfld[18]][j]*(-27.0*bary[3])
                      +  rfld[lfld[19]][j]*( 27.0*bary[0]+-27.0*bary[1]+-27.0*bary[2]);
    }
    for(int j=0; j < szfld; j++){
      hmat[2*szfld+j] =  rfld[lfld[0]][j]*( 27.0*bary[0]-9.0)
                      +  rfld[lfld[2]][j]*( 27.0*bary[2]-9.0)
                      +  rfld[lfld[4]][j]*(27.0*bary[1])
                      +  rfld[lfld[7]][j]*(27.0*bary[1])
                      +  rfld[lfld[8]][j]*( 27.0*bary[0]+-54.0*bary[2]+9.0)
                      +  rfld[lfld[9]][j]*( 27.0*bary[2]+-54.0*bary[0]+9.0)
                      +  rfld[lfld[10]][j]*(27.0*bary[3])
                      +  rfld[lfld[14]][j]*(27.0*bary[3])
                      +  rfld[lfld[17]][j]*(-54.0*bary[3])
                      +  rfld[lfld[19]][j]*(-54.0*bary[1]);
    }
    for(int j=0; j < szfld; j++){
      hmat[3*szfld+j] =  rfld[lfld[0]][j]*( 27.0*bary[0]-9.0)
                      +  rfld[lfld[4]][j]*( -27.0*bary[0]+27.0*bary[1]+(9.0/2.0))
                      +  rfld[lfld[5]][j]*( -27.0*bary[1]+(9.0/2.0))
                      +  rfld[lfld[9]][j]*(27.0*bary[2])
                      +  rfld[lfld[10]][j]*( -27.0*bary[0]+27.0*bary[3]+(9.0/2.0))
                      +  rfld[lfld[11]][j]*( -27.0*bary[3]+(9.0/2.0))
                      +  rfld[lfld[12]][j]*( 27.0*bary[1]-(9.0/2.0))
                      +  rfld[lfld[13]][j]*( 27.0*bary[3]-(9.0/2.0))
                      +  rfld[lfld[16]][j]*(27.0*bary[2])
                      +  rfld[lfld[17]][j]*(-27.0*bary[2])
                      +  rfld[lfld[18]][j]*( -27.0*bary[1]+27.0*bary[0]+-27.0*bary[3])
                      +  rfld[lfld[19]][j]*(-27.0*bary[2]);
    }
    for(int j=0; j < szfld; j++){
      hmat[4*szfld+j] =  rfld[lfld[0]][j]*( 27.0*bary[0]-9.0)
                      +  rfld[lfld[4]][j]*(27.0*bary[1])
                      +  rfld[lfld[8]][j]*( -27.0*bary[2]+(9.0/2.0))
                      +  rfld[lfld[9]][j]*( 27.0*bary[2]+-27.0*bary[0]+(9.0/2.0))
                      +  rfld[lfld[10]][j]*( -27.0*bary[0]+27.0*bary[3]+(9.0/2.0))
                      +  rfld[lfld[11]][j]*( -27.0*bary[3]+(9.0/2.0))
                      +  rfld[lfld[14]][j]*( 27.0*bary[2]-(9.0/2.0))
                      +  rfld[lfld[15]][j]*( 27.0*bary[3]-(9.0/2.0))
                      +  rfld[lfld[16]][j]*(27.0*bary[1])
                      +  rfld[lfld[17]][j]*( 27.0*bary[0]+-27.0*bary[2]+-27.0*bary[3])
                      +  rfld[lfld[18]][j]*(-27.0*bary[1])
                      +  rfld[lfld[19]][j]*(-27.0*bary[1]);
    }
    for(int j=0; j < szfld; j++){
      hmat[5*szfld+j] =  rfld[lfld[0]][j]*( 27.0*bary[0]-9.0)
                      +  rfld[lfld[3]][j]*( 27.0*bary[3]-9.0)
                      +  rfld[lfld[4]][j]*(27.0*bary[1])
                      +  rfld[lfld[9]][j]*(27.0*bary[2])
                      +  rfld[lfld[10]][j]*( -54.0*bary[0]+27.0*bary[3]+9.0)
                      +  rfld[lfld[11]][j]*( 27.0*bary[0]+-54.0*bary[3]+9.0)
                      +  rfld[lfld[13]][j]*(27.0*bary[1])
                      +  rfld[lfld[15]][j]*(27.0*bary[2])
                      +  rfld[lfld[17]][j]*(-54.0*bary[2])
                      +  rfld[lfld[18]][j]*(-54.0*bary[1]);
    }
  }
}
};

}//End namespace


#endif