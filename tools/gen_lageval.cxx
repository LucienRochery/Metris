//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../src/ho_constants.hxx"
#include "../src/metris_constants.hxx"

#include <string>
#include <cmath>
#include <sstream>
#include <fstream>
#include <ginac/ginac.h>

void gen_lageval_alldim();
void eval_lagrange3_symb(int i0,int j0, int k0, int l0,
  GiNaC::ex &eval, GiNaC::ex dlag[]);
void eval_lagrange2_symb(int i0,int j0, int k0, 
  GiNaC::ex &eval, GiNaC::ex dlag[]);
void eval_lagrange1_symb(int i0,int j0, 
  GiNaC::ex &eval, GiNaC::ex dlag[]);


using namespace Metris;



int main(int argc, char **argv){

  printf("Maximum degree = %d \n",METRIS_MAX_DEG);

  printf("------- Generate codegen_lagrange.hxx\n");
  gen_lageval_alldim();

  return 0;
}




void eval_lagrange3_symb(int i0,int j0, int k0, int l0,
  GiNaC::ex &eval, GiNaC::ex dlag[], GiNaC::ex d2lag[]){
  int fac;

  GiNaC::ex sub;
  GiNaC::ex ev[4];
  GiNaC::symbol bary[4] = {GiNaC::symbol("bary[0]"),
  GiNaC::symbol("bary[1]"),
  GiNaC::symbol("bary[2]"),
  GiNaC::symbol("bary[3]")};

  int ideg = i0 + j0 + k0 + l0;

  fac  = 1;
  eval = 1;
  for(int j=0;j<i0;j++){
    eval *= (bary[0]*ideg - j);
    fac  *= (i0 - j);
  }
  for(int j=0;j<j0;j++){
    eval *= (bary[1]*ideg - j);
    fac  *= (j0-j);
  }
  for(int j=0;j<k0;j++){
    eval *= (bary[2]*ideg - j);
    fac  *= (k0-j);
  }
  for(int j=0;j<l0;j++){
    eval *= (bary[3]*ideg - j);
    fac  *= (l0-j);
  }
  eval = eval / fac;


/*
  //   -- Derivatives in the sense of Panzoult: d_i (phy) = d_(i+1) - d_1 (bar)
  //   This convention is so unnatural it cannot hurt to recall it every time. 

  //   -- Compute eval, but separate contributions from each barycentric crd
  //   This will minimize recomputations later on (when diff)
  //  fac = 1.0;
  ev[0] = 1 ; //! separate the contributions to minimize recomputations
  for(int j = 0; j < i0; j++){
    ev[0] *=(bary[0]*ideg - j);
  }

  ev[1] = 1 ; //! separate the contributions to minimize recomputations
  for(int j = 0; j < j0; j++){
    ev[1] *=(bary[1]*ideg - j);
  }

  ev[2] = 1 ; //! separate the contributions to minimize recomputations
  for(int j = 0; j < k0; j++){
    ev[2] *=(bary[2]*ideg - j);
  }

  ev[3] = 1 ; //! separate the contributions to minimize recomputations
  for(int j = 0; j < l0; j++){
    ev[3] *=(bary[3]*ideg - j);
  }



  //  -- Compute d1, all others will be d_i - d1
  GiNaC::ex d1 = 0;
  for (int j = 0; j < i0; j++){
    sub = 1;
    for (int k = 0; k < i0; k++){
      if(k == j) continue;
      sub *= (bary[0]*ideg-k);
    }
    d1 += sub;
  }
  //   -- Multiply by constant factors
  d1 *= ev[1]*ev[2]*ev[3];
  //   -- Final factor and -d1, recall not multiplied by fac
  d1 = d1 * ideg;// / fac;


  //   -- Do the same thing with the others. 
  dlag[0] = 0;
  for (int j = 0; j < j0; j++){
    sub = 1;
    for (int k = 0; k < j0; k++){
      if(k == j) continue;
      sub *= (bary[1]*ideg-k);
    }
    dlag[0] += sub;
  }
  //   -- Multiply by constant factors
  dlag[0] *= ideg*ev[0]*ev[2]*ev[3];
  //   -- Final factor and -d1, recall not multiplied by fac
  dlag[0] = (dlag[0] - d1) / fac;


  dlag[1] = 0;
  for (int j = 0; j < k0; j++){
    sub = 1;
    for (int k = 0; k < k0; k++){
      if(k == j) continue;
      sub *= (bary[2]*ideg-k);
    }
    dlag[1] += sub;
  }
  //   -- Multiply by constant factors
  dlag[1] *= ideg*ev[0]*ev[1]*ev[3];
  //   -- Final factor and -d1, recall not multiplied by fac
  dlag[1] = (dlag[1] - d1) / fac;



  dlag[2] = 0;
  for (int j = 0; j < l0; j++){
    sub = 1;
    for (int k = 0; k < l0; k++){
      if(k == j) continue;
      sub *= (bary[3]*ideg-k);
    }
    dlag[2] += sub;
  }
  //   -- Multiply by constant factors
  dlag[2] *= ideg*ev[0]*ev[1]*ev[2];
  //   -- Final factor and -d1, recall not multiplied by fac
  dlag[2] = (dlag[2] - d1) / fac;

//    std::cout<<"\n\n";
//  a  std::cout<<i0<<j0<<k0<<l0<<" eval = "<<eval<<std::endl;
    //std::cout<<"end dlag1="<<dlag[0]<<std::endl;
    //std::cout<<"end dlag2="<<dlag[1]<<std::endl;
    //std::cout<<"end dlag3="<<dlag[2]<<std::endl;
//
    //std::cout<<"\n";
    //std::cout<<"collected"<<GiNaC::collect_common_factors(eval)<<std::endl;
    //std::cout<<"collected"<<GiNaC::collect_common_factors(dlag[0])<<std::endl;
    //std::cout<<"collected"<<GiNaC::collect_common_factors(dlag[1])<<std::endl;
    //std::cout<<"collected"<<GiNaC::collect_common_factors(dlag[2])<<std::endl;
//
    //std::cout<<"\n";
    //std::cout<<"factor "<<GiNaC::factor(eval)<<std::endl;
    //std::cout<<"factor "<<GiNaC::factor(dlag[0])<<std::endl;
    //std::cout<<"factor "<<GiNaC::factor(dlag[1])<<std::endl;
    //std::cout<<"factor "<<GiNaC::factor(dlag[2])<<std::endl;
*/
  auto d1 = eval.diff(bary[0]);
  dlag[0] = GiNaC::factor(eval.diff(bary[1]) - d1);
  dlag[1] = GiNaC::factor(eval.diff(bary[2]) - d1);
  dlag[2] = GiNaC::factor(eval.diff(bary[3]) - d1);

  auto tmp = GiNaC::factor(dlag[0].diff(bary[0]));
  // d_11
  d2lag[0] = dlag[0].diff(bary[1]) - tmp;
  // d_12
  d2lag[1] = dlag[0].diff(bary[2]) - tmp;
  // d_13
  d2lag[3] = dlag[0].diff(bary[3]) - tmp;

  tmp = GiNaC::factor(dlag[1].diff(bary[0]));
  // d_22
  d2lag[2] = dlag[1].diff(bary[2]) - tmp;
  // d_23
  d2lag[4] = dlag[1].diff(bary[3]) - tmp;

  //d_33
  d2lag[5] = dlag[2].diff(bary[3]) - dlag[2].diff(bary[0]);


  eval    = GiNaC::factor(eval) ;
  for(int ii = 0; ii < 6; ii++) d2lag[ii] = GiNaC::factor(d2lag[ii]);

  //std::cout<<"Debug d1lag[0] = "<<dlag[0]<<"\n";
  //std::cout<<"Debug d2lag[0] = "<<d2lag[0]<<"\n";

}


void eval_lagrange2_symb(int i0,int j0, int k0,
 GiNaC::ex &eval, GiNaC::ex dlag[])
{
  int fac;

  GiNaC::ex sub;
  GiNaC::ex ev[3];
  GiNaC::symbol bary[3] = {GiNaC::symbol("bary[0]"),
  GiNaC::symbol("bary[1]"),
  GiNaC::symbol("bary[2]")};

  int ideg = i0 + j0 + k0;

  fac  = 1;
  eval = 1;
  for(int j=0;j<i0;j++){
    eval *= (bary[0]*ideg - j);
    fac  *= (i0 - j);
  }
  for(int j=0;j<j0;j++){
    eval *= (bary[1]*ideg - j);
    fac  *= (j0-j);
  }
  for(int j=0;j<k0;j++){
    eval *= (bary[2]*ideg - j);
    fac  *= (k0-j);
  }
  eval = eval / fac;


//     -- Derivatives in the sense of Panzoult: d_i (phy) = d_(i+1) - d_1 (bar)
//     This convention is so unnatural it cannot hurt to recall it every time. 

//     -- Compute eval, but separate contributions from each barycentric crd
//     This will minimize recomputations later on (when diff)
//    fac = 1.0;
  ev[0] = 1 ; //! separate the contributions to minimize recomputations
  for(int j = 0; j < i0; j++){
    ev[0] *= (bary[0]*ideg - j);
  }

  ev[1] = 1 ; //! separate the contributions to minimize recomputations
  for(int j = 0; j < j0; j++){
    ev[1] *= (bary[1]*ideg - j);
  }

  ev[2] = 1 ; //! separate the contributions to minimize recomputations
  for(int j = 0; j < k0; j++){
    ev[2] *= (bary[2]*ideg - j);
  }


//   -- Compute d1, all others will be d_i - d1
  GiNaC::ex d1 = 0;
  for (int j = 0; j < i0; j++){
    sub = 1;
    for (int k = 0; k < i0; k++){
      if(k == j) continue;
      sub *= (bary[0]*ideg-k);
    }
    d1 += sub;
  }
//   -- Multiply by constant factors
  d1 *= ev[1]*ev[2];
//   -- Final factor and -d1, recall not multiplied by fac
  d1 = d1 * ideg;// / fac;


//   -- Do the same thing with the others. 
  dlag[0] = 0;
  for (int j = 0; j < j0; j++){
    sub = 1;
    for (int k = 0; k < j0; k++){
      if(k == j) continue;
      sub *= (bary[1]*ideg-k);
    }
    dlag[0] += sub;
  }
//   -- Multiply by constant factors
  dlag[0] *= ideg*ev[0]*ev[2];
//   -- Final factor and -d1, recall not multiplied by fac
  dlag[0] = (dlag[0] - d1) / fac;



  dlag[1] = 0;
  for (int j = 0; j < k0; j++){
    sub = 1;
    for (int k = 0; k < k0; k++){
      if(k == j) continue;
      sub *= (bary[2]*ideg-k);
    }
    dlag[1] += sub;
  }
//   -- Multiply by constant factors
  dlag[1] *= ideg*ev[0]*ev[1];
//   -- Final factor and -d1, recall not multiplied by fac
  dlag[1] = (dlag[1] - d1) / fac;



  eval    = GiNaC::factor(eval) ;
  dlag[0] = GiNaC::factor(dlag[0]);
  dlag[1] = GiNaC::factor(dlag[1]);

}


void eval_lagrange1_symb(int i0,int j0, 
 GiNaC::ex &eval, GiNaC::ex dlag[])
{
  int fac;

  GiNaC::ex sub;
  GiNaC::ex ev[2];
  GiNaC::symbol bary[2] = {GiNaC::symbol("bary[0]"),
  GiNaC::symbol("bary[1]")};

  int ideg = i0 + j0;

  fac  = 1;
  eval = 1;
  for(int j=0;j<i0;j++){
    eval *= (bary[0]*ideg - j);
    fac  *= (i0 - j);
  }
  for(int j=0;j<j0;j++){
    eval *= (bary[1]*ideg - j);
    fac  *= (j0-j);
  }
  eval = eval / fac;


//     -- Derivatives in the sense of Panzoult: d_i (phy) = d_(i+1) - d_1 (bar)
//     This convention is so unnatural it cannot hurt to recall it every time. 

//     -- Compute eval, but separate contributions from each barycentric crd
//     This will minimize recomputations later on (when diff)
//    fac = 1.0;
  ev[0] = 1 ; //! separate the contributions to minimize recomputations
  for(int j = 0; j < i0; j++){
    ev[0] *=(bary[0]*ideg - j);
  }

  ev[1] = 1 ; //! separate the contributions to minimize recomputations
  for(int j = 0; j < j0; j++){
    ev[1] *=(bary[1]*ideg - j);
  }

//   -- Compute d1, all others will be d_i - d1
  GiNaC::ex d1 = 0;
  for (int j = 0; j < i0; j++){
    sub = 1;
    for (int k = 0; k < i0; k++){
      if(k == j) continue;
      sub *= (bary[0]*ideg-k);
    }
    d1 += sub;
  }
//   -- Multiply by constant factors
  d1 *= ev[1];
//   -- Final factor and -d1, recall not multiplied by fac
  d1 = d1 * ideg;// / fac;


//   -- Do the same thing with the others. 
  dlag[0] = 0;
  for (int j = 0; j < j0; j++){
    sub = 1;
    for (int k = 0; k < j0; k++){
      if(k == j) continue;
      sub *= (bary[1]*ideg-k);
    }
    dlag[0] += sub;
  }
//   -- Multiply by constant factors
  dlag[0] *= ideg*ev[0];
//   -- Final factor and -d1, recall not multiplied by fac
  dlag[0] = (dlag[0] - d1) / fac;



  eval    = GiNaC::factor(eval) ;
  dlag[0] = GiNaC::factor(dlag[0]);

}

void gen_lageval_alldim(){
  constexpr int max_eval_deg = METRIS_MAX_DEG_JACOBIAN;
  std::ostringstream str3,str2,str1;
  str3 << GiNaC::csrc;
  str2 << GiNaC::csrc;
  str1 << GiNaC::csrc;

  GiNaC::ex eval[tetnpps[max_eval_deg]];
  GiNaC::ex dlag[tetnpps[max_eval_deg]][3], d2lag[tetnpps[max_eval_deg]][6];
  
 	// ---------------------- Tetras
  for(int ideg = 1; ideg <= max_eval_deg; ideg++){
    printf("Tetrahedra degree %d\n",ideg);

    for(int ifunc = 0; ifunc < tetnpps[ideg]; ifunc++){
      eval_lagrange3_symb(ordtet.s[ideg][ifunc][0],ordtet.s[ideg][ifunc][1],
                          ordtet.s[ideg][ifunc][2],ordtet.s[ideg][ifunc][3],
                          eval[ifunc], dlag[ifunc], d2lag[ifunc]);	
    }

    str3 << "\ntemplate <int szfld>\n";
    str3 << "struct eval3_lagrange<szfld,"+std::to_string(ideg)+">{\n";
    str3 << "  eval3_lagrange(const dblAr2 & __restrict__  rfld,\n"
    << "                 const int * __restrict__  lfld,\n"
    << "                 DifVar idif1, DifVar idif2,\n"
    << "                 const double * __restrict__  bary,\n"
    << "                 double * __restrict__  eval,\n"
    << "                 double * __restrict__  jmat,\n"
    << "                 double * __restrict__  hmat){\n";

    str3 << "  for(int j=0; j < szfld; j++){\n";
    str3 << "    eval[j] = rfld[lfld[0]][j]*("<<eval[0]<<")";
    for(int i = 1; i < tetnpps[ideg]; i++){
      str3 << "\n            + rfld[lfld["+std::to_string(i)+"]][j]*("<<
      eval[i]<<")";
    }
    str3 << ";\n  }\n";

    str3<<"  if (idif1 == DifVar::Bary){\n";

    //  --------------- Jacobian
    int ifirst;
    for(int ider = 0; ider < 3; ider ++){
      ifirst = 1;
      str3 << "    for(int j=0; j < szfld; j++){\n";
      for(int i = 0; i < tetnpps[ideg]; i++){
        if(dlag[i][ider] != 0){
          if(ifirst == 1){
            ifirst = 0;
            str3 << "      jmat["+std::to_string(ider)+"*szfld+j] =  rfld[lfld["+std::to_string(i)+"]][j]*("<<dlag[i][ider]<<")";
          }else{
            str3 << "\n                      +  rfld[lfld["+std::to_string(i)+"]][j]*("<<dlag[i][ider]<<")";
          }
        }
      }
      str3 << ";\n    }\n";
    }
    str3<<"  }\n";

    str3<<"  if (idif2 == DifVar::Bary){\n";

    // ---------------- Second derivatives
    if(ideg > 1){
      for(int ider = 0; ider < 6; ider ++){
        ifirst = 1;
        str3 << "    for(int j=0; j < szfld; j++){\n";
        for(int i = 0; i < tetnpps[ideg]; i++){
          if(d2lag[i][ider] != 0){
            if(ifirst == 1){
              ifirst = 0;
              str3 << "      hmat["+std::to_string(ider)+"*szfld+j] =  rfld[lfld["+std::to_string(i)+"]][j]*("<<d2lag[i][ider]<<")";
            }else{
              str3 << "\n                      +  rfld[lfld["+std::to_string(i)+"]][j]*("<<d2lag[i][ider]<<")";
            }
          }
        }
        str3 << ";\n    }\n";
      }
    }else{
      str3 << "    for(int j=0; j < szfld; j++){\n";
      for(int ider = 0; ider < 6; ider++){
        str3<<"      hmat["+std::to_string(ider)+"*szfld+j] = 0;\n";
      }
      str3 << "    }\n";
    }
    str3<<"  }\n";


    str3 << "}\n"; // Closing the constructor
 		str3<<"};"; //closing the struct.

 	}


 	// ---------------------- Triangles
 	for(int ideg = 1; ideg <= max_eval_deg; ideg++){
    printf("Triangles degree %d\n",ideg);


    for(int ifunc = 0; ifunc < facnpps[ideg]; ifunc++){
      eval_lagrange2_symb(ordfac.s[ideg][ifunc][0],ordfac.s[ideg][ifunc][1],
                          ordfac.s[ideg][ifunc][2],
                          eval[ifunc], dlag[ifunc]);	
    }

    str2 << "\ntemplate <int szfld>\n";
    str2 << "struct eval2_lagrange<szfld,"+std::to_string(ideg)+">{\n";
    str2 << "  eval2_lagrange(const dblAr2 & __restrict__  rfld,\n"
    << "                 const int * __restrict__  lfld,\n"
    << "                 DifVar idif1,\n"
    << "                 const double * __restrict__  bary,\n"
    << "                 double * __restrict__  eval, double * __restrict__  jmat){\n";

    str2 << "  for(int j=0; j < szfld; j++){\n";
    str2 << "    eval[j] = rfld[lfld[0]][j]*("<<eval[0]<<")";
    for(int i = 1; i < facnpps[ideg]; i++){
     str2 << "\n            + rfld[lfld["+std::to_string(i)+"]][j]*("<<
           eval[i]<<")";
   }
   str2 << ";\n  }\n";

   str2<<"  if (idif1 != DifVar::Bary) return;\n";

//  --------------- Jacobian
   int ifirst;

   for(int ider = 0; ider < 2; ider ++){

    ifirst = 1;
    str2 << "  for(int j=0; j < szfld; j++){\n";
    for(int i = 0; i < facnpps[ideg]; i++){
      if(dlag[i][ider] != 0){
        if(ifirst == 1){
          ifirst = 0;
          str2 << "    jmat["+std::to_string(ider)+"*szfld+j] =  rfld[lfld["+std::to_string(i)+"]][j]*("<<dlag[i][ider]<<")";
        }else{
          str2 << "\n                    +  rfld[lfld["+std::to_string(i)+"]][j]*("<<dlag[i][ider]<<")";
        }
      }
    }
    str2 << ";\n  }\n";
  }
  str2 << "}\n";
 		str2<<"};"; //closing the struct.
 	}

 	// ---------------------- Edges
 	for(int ideg = 1; ideg <= max_eval_deg; ideg++){
    printf("Edges degree %d\n",ideg);


    for(int ifunc = 0; ifunc < edgnpps[ideg]; ifunc++){
      eval_lagrange1_symb(ordedg.s[ideg][ifunc][0],ordedg.s[ideg][ifunc][1],
                          eval[ifunc], dlag[ifunc]);	
    }

    str1 << "\ntemplate <int szfld>\n";
    str1 << "struct eval1_lagrange<szfld,"+std::to_string(ideg)+">{\n";
    str1 << "  eval1_lagrange(const dblAr2 & __restrict__  rfld,\n"
    << "                 const int * __restrict__  lfld,\n"
    << "                 DifVar idif1,\n"
    << "                 const double * __restrict__  bary,\n"
    << "                 double * __restrict__  eval, double * __restrict__  jmat){\n";

    str1 << "  for(int j=0; j < szfld; j++){\n";
    str1 << "    eval[j] = rfld[lfld[0]][j]*("<<eval[0]<<")";
    for(int i = 1; i < edgnpps[ideg]; i++){
      str1 << "\n            + rfld[lfld["+std::to_string(i)+"]][j]*("<<
        eval[i]<<")";
    }
    str1 << ";\n  }\n";

    str1<<"  if (idif1 != DifVar::Bary) return;\n";

    //  --------------- Jacobian
    int ifirst;

    for(int ider = 0; ider < 1; ider ++){

      ifirst = 1;
      str1 << "  for(int j=0; j < szfld; j++){\n";
      for(int i = 0; i < edgnpps[ideg]; i++){
        if(dlag[i][ider] != 0){
          if(ifirst == 1){
            ifirst = 0;
            str1 << "    jmat["+std::to_string(ider)+"*szfld+j] =  rfld[lfld["+std::to_string(i)+"]][j]*("<<dlag[i][ider]<<")";
          }else{
            str1 << "\n                    +  rfld[lfld["+std::to_string(i)+"]][j]*("<<dlag[i][ider]<<")";
          }
        }
      }
      str1 << ";\n  }\n";
    }
    str1 << "}\n";
 		str1<<"};"; //closing the struct.
 	}

 	std::ofstream fc;
 	fc.open("codegen_lagrange.hxx", std::ios::out);

  fc << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  fc << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  fc << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  fc << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  fc << "#ifndef __CODEGEN_LAGRANGE__\n#define __CODEGEN_LAGRANGE__\n";
  fc << "#include \"types.hxx\"\n";
  fc << "#include \"metris_constants.hxx\"\n";
  fc << "namespace Metris{\n";
  fc << "\n\n//------------------- Edges --------------------\n\n";
  fc << "template <int szfld, int ideg> \n";
  fc << "struct eval1_lagrange{\n"
  << "  eval1_lagrange(const dblAr2 & __restrict__  rfld,\n"
  << "                 const int * __restrict__  lfld,\n"
  << "                 DifVar idif1,\n"
  << "                 const double * __restrict__  bary,\n"
  << "                 double * __restrict__  eval, double * __restrict__  jmat);\n"
  << "};\n";
  fc << str1.str();
  fc << "\n\n//------------------- Triangles --------------------\n\n";
  fc << "template <int szfld, int ideg> \n";
  fc << "struct eval2_lagrange{\n"
  << "  eval2_lagrange(const dblAr2 & __restrict__  rfld,\n"
  << "                 const int * __restrict__  lfld,\n"
  << "                 DifVar idif1,\n"
  << "                 const double * __restrict__  bary,\n"
  << "                 double * __restrict__  eval, double * __restrict__  jmat);\n"
  << "};\n";
  fc << str2.str();
  fc << "\n\n//------------------- Tetrahedra --------------------\n\n";
  fc << "template <int szfld, int ideg> \n";
  fc << "struct eval3_lagrange{\n"
  << "  eval3_lagrange(const dblAr2 & __restrict__  rfld,\n"
  << "                 const int * __restrict__  lfld,\n"
  << "                 DifVar idif1, DifVar idif2,\n"
  << "                 const double * __restrict__  bary,\n"
  << "                 double * __restrict__  eval,\n"
  << "                 double * __restrict__  jmat,\n"
  << "                 double * __restrict__  hmat);\n"
  << "};\n";
  fc << str3.str();
  fc << "\n\n}//End namespace\n\n";
  fc << "\n#endif";
  fc.close();
//  fclose(fh);
}

