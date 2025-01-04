//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


// Tool to generate source files to compute control coefficients.
// "Hand-written" routines turn out to be about 3x as quick, at least in P2, as 
// even compile-time routines, where the only runtime call is the determinant computation.
// This assumes an input Bézier mesh. 
#include "gen_ccoef.hxx"

#include "../src/ho_constants.hxx"
#include <cmath>

#include <string>
#include <sstream>
#include <fstream>


using namespace Metris;

void gen_ccoeff3_2();



static int ifactorial_func(int i){
	if(i == 0) return 1;
	return ifactorial_func(i-1)*i;
}

//template <int i1, int i2, int i3, int i4>
//void gen_idx_ccoef(char* str);
void simpfrac(int x, int y, int *xs, int *ys){
	int z = sqrt( x < y ? x : y );
	int i = z;
	*xs = x;
	*ys = y;
	while(i > 1){
		if(*xs % i == 0 && *ys % i == 0){
			*xs /= i;
			*ys /= i;
		}
		i--;
	}
}


int main(int argc, char** argv){

  printf("Maximum degree = %d \n",METRIS_MAX_DEG);


  printf("------- Generate codegen_ccoeff.hxx\n");
  std::ostringstream str;
  str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  str << "#ifndef __CODEGEN_CCOEF__\n";
  str << "#define __CODEGEN_CCOEF__\n";
  str << "\n";
  str << "#include \"types.hxx\"\n";
  str << "namespace Metris{\n\n";

  str << "template<int ideg>\n";
  str << "void ccoef_genbez3(const intAr2 & __restrict__ tet2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef);\n";
  str << "template<int ideg>\n";
  str << "void ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef);\n";
  str << "}//End namespace\n";
  str << "#endif\n";
  std::ofstream f;
  char fname[64]; snprintf(fname,64,"codegen_ccoef.hxx");
  f.open(fname, std::ios::out);
  f << str.str();
  f.close();

  str.str("");
  str.clear();
  printf("------- Generate codegen_ccoeff_d.hxx\n");
  str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  str << "#ifndef __CODEGEN_CCOEF_D__\n";
  str << "#define __CODEGEN_CCOEF_D__\n";
  str << "\n";
  str << "#include \"types.hxx\"\n";
  str << "namespace Metris{\n\n";


  str << "template<int ideg> void d_ccoef_genbez2"
      << "(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,\n" 
      << "                                        int ielem, double*__restrict__ ccoef, \n"
      << "                                        int icoor, dblAr2& __restrict__ d_ccoef);\n";

  str << "template<int ideg> void d_pt_ccoef_genbez2"
  << "(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,\n" 
  << "                                        int ielem,\n"
  << "                                        int icoor,\n"
  << "                                        int inode,\n"
  << "                                        double*__restrict__ ccoef,\n"
  << "                                        double*__restrict__ d_ccoef);\n";

  str << "template<int ideg, typename T>\n";
  str << "void d_ccoef_genbez3(const intAr2 & __restrict__ tet2poi,\n";
  str << "                     const dblAr2& __restrict__ coord,\n";
  str << "                     int ielem,\n";
  str << "                     int icoor,\n";
  str << "                     T* __restrict__ ccoef,\n";
  str << "                     dblAr2& __restrict__ d_ccoef);\n\n";


  str << "}//End namespace\n";
  str << "#endif\n";
  snprintf(fname,64,"codegen_ccoef_d.hxx");
  f.open(fname, std::ios::out);
  f << str.str();
  f.close();


  printf("------- Generate codegen_ccoeff3<d>.cxx\n");
  gen_ccoeff3();

  printf("------- Generate codegen_ccoeff2<d>.cxx\n");
  gen_ccoeff2();


  printf("------- Generate codegen_ccoeff2<d>_d.cxx\n");
  gen_ccoeff2_d();

  printf("------- Generate codegen_ccoeff2<d>_d_pt.cxx\n");
  gen_ccoef2_d_pt();


  printf("------- Generate codegen_ccoeff3<d>_d_pt.cxx\n");
  gen_ccoef3_d();

  return 0;
}



//void gen_ccoeff3_2(){
//
//  if(METRIS_MAX_DEG > 9){
//    printf("Welcome to 2050, we hope the flying cars are comfortable and powered by nuclear fusion.\n Change this routine, e.g. adding underscores between indices, as they are no longer single-digit.\n");
//    exit(1);
//  }
//
//  int max_jaco_deg = METRIS_MAX_DEG;
//
//  printf("-- Tetrahedra up to degree %d \n",METRIS_MAX_DEG);
//  for(int ideg = 2; ideg <= max_jaco_deg; ideg++){
//    std::ostringstream funargs;
//    std::ostringstream str;
//
//
//    for(int irnk1 = 0; irnk1 < tetnpps[ideg]; irnk1++){
//      int i = ordtet.s[ideg][irnk1][0];
//      int j = ordtet.s[ideg][irnk1][1];
//      int k = ordtet.s[ideg][irnk1][2];
//      int l = ordtet.s[ideg][irnk1][3];
//      char i_s[16]; snprintf(i_s,2,"%1d",i);
//      char j_s[16]; snprintf(j_s,2,"%1d",j);
//      char k_s[16]; snprintf(k_s,2,"%1d",k);
//      char l_s[16]; snprintf(l_s,2,"%1d",l);
//      funargs<<"double* __restrict__ cpoi"<<i_s<<j_s<<k_s<<l_s;
//      if(irnk1 < tetnpps[ideg]-1) funargs<<",";
//      if(irnk1 == 2 || irnk1 >  2 && (irnk1-2)%4 == 0) funargs<<"\n";
//    }
//
//
//    str << "#include \"codegen_ccoef.hxx\"\n\n";
//
//    str << "double det3_vdif(const double* x1,const double* x2
//                  ,const double* y1,const double* y2
//                  ,const double* z1,const double* z2);\n\n";
//
//
//    int npp_c = tetnpps[3*(ideg-1)];
//
//    int uuu1 = pow(ifactorial_func(ideg),3);
//    int lll1 = ifactorial_func(3*(ideg-1));
//
//    str << "void ccoef_genbez3_"<<ideg<<"("<<funargs.str()<< "){\n";
//    for(int irnkcc = 0; irnkcc < npp_c; irnkcc++){
//      int ifirst = 1;
//      int i1 = ordtet.s[3*(ideg-1)][irnkcc][0]; // constexpr
//      int i2 = ordtet.s[3*(ideg-1)][irnkcc][1]; // constexpr
//      int i3 = ordtet.s[3*(ideg-1)][irnkcc][2]; // constexpr
//      int i4 = ordtet.s[3*(ideg-1)][irnkcc][3]; // constexpr
//      int uuu2 = ifactorial_func(i1)*ifactorial_func(i2)*ifactorial_func(i3)*ifactorial_func(i4);
//
//      for(int j1 = 0; j1 <= ideg - 1 && j1 <= i1; j1 ++){
//        for(int j2 = 0; j2 <= ideg - 1 - j1 && j2 <= i2; j2 ++){
//          for(int j3 = 0; j3 <= ideg - 1 - j1 - j2 && j3 <= i3; j3 ++){
//            if(ideg-1-j1-j2-j3 > i4) continue;
//            int j4 = ideg-1-j1-j2-j3;
//
//            /*-----------------------------------K LOOP-------------------------------------------*/
//            for(int k1 = 0; k1 <= ideg - 1 && k1 <= i1 - j1 ; k1++){
//              for(int k2 = 0; k2 <= ideg - 1 - k1 && k2 <= i2 - j2 ; k2++){
//                for(int k3 = 0; k3 <= ideg - 1 - k1 - k2 && k3 <= i3 - j3 ; k3++){
//                  if(ideg-1-k1-k2-k3 > i4-j4) continue;
//                  int k4 = ideg - 1 - k1 - k2 - k3;
//                  int l1 = i1 - j1 - k1;
//                  int l2 = i2 - j2 - k2;
//                  int l3 = i3 - j3 - k3;
//                  int l4 = i4 - j4 - k4;
//
//                  int cc2 = ifactorial_func(j1)*ifactorial_func(j2)*ifactorial_func(j3)*ifactorial_func(j4);
//                  int cc3 = ifactorial_func(k1)*ifactorial_func(k2)*ifactorial_func(k3)*ifactorial_func(k4);
//                  int cc4 = ifactorial_func(l1)*ifactorial_func(l2)*ifactorial_func(l3)*ifactorial_func(l4);
//                  int lll2 = cc2*cc3*cc4;
//
//
//
//                  //int irnk1_1 = mul2nod(j1,j2+1,j3,j4);
//                  //int irnk1_2 = mul2nod(j1+1,j2,j3,j4);
//                  char j1_1_s[16]; snprintf(j1_1_s,2,"%1d",j1  );
//                  char j2_1_s[16]; snprintf(j2_1_s,2,"%1d",j2+1);
//                  char j3_1_s[16]; snprintf(j3_1_s,2,"%1d",j3  );
//                  char j4_1_s[16]; snprintf(j4_1_s,2,"%1d",j4  );
//
//                  char j1_2_s[16]; snprintf(j1_2_s,2,"%1d",j1+1);
//                  char j2_2_s[16]; snprintf(j2_2_s,2,"%1d",j2  );
//                  char j3_2_s[16]; snprintf(j3_2_s,2,"%1d",j3  );
//                  char j4_2_s[16]; snprintf(j4_2_s,2,"%1d",j4  );
//
//
//                  //int irnk2_1 = mul2nod(k1,k2,k3+1,k4);
//                  //int irnk2_2 = mul2nod(k1+1,k2,k3,k4);
//                  char k1_1_s[16]; snprintf(k1_1_s,2,"%1d",k1  );
//                  char k2_1_s[16]; snprintf(k2_1_s,2,"%1d",k2  );
//                  char k3_1_s[16]; snprintf(k3_1_s,2,"%1d",k3+1);
//                  char k4_1_s[16]; snprintf(k4_1_s,2,"%1d",k4  );
//
//                  char k1_2_s[16]; snprintf(k1_2_s,2,"%1d",k1+1);
//                  char k2_2_s[16]; snprintf(k2_2_s,2,"%1d",k2  );
//                  char k3_2_s[16]; snprintf(k3_2_s,2,"%1d",k3  );
//                  char k4_2_s[16]; snprintf(k4_2_s,2,"%1d",k4  );
//
//
//                  //int irnk3_1 = mul2nod(l1,l2,l3,l4+1);
//                  //int irnk3_2 = mul2nod(l1+1,l2,l3,l4);
//                  char l1_1_s[16]; snprintf(l1_1_s,2,"%1d",l1  );
//                  char l2_1_s[16]; snprintf(l2_1_s,2,"%1d",l2  );
//                  char l3_1_s[16]; snprintf(l3_1_s,2,"%1d",l3  );
//                  char l4_1_s[16]; snprintf(l4_1_s,2,"%1d",l4+1);
//
//                  char l1_2_s[16]; snprintf(l1_2_s,2,"%1d",l1+1);
//                  char l2_2_s[16]; snprintf(l2_2_s,2,"%1d",l2  );
//                  char l3_2_s[16]; snprintf(l3_2_s,2,"%1d",l3  );
//                  char l4_2_s[16]; snprintf(l4_2_s,2,"%1d",l4  );
//
//
//                  int up,lo;
//                  simpfrac(uuu1*uuu2,lll1*lll2,&up,&lo);
//
//                  char irnkcc_s[16]; snprintf(irnkcc_s,4,"%3d",irnkcc);
//                  char up_s[16]; snprintf(up_s,4,"%3d",up);
//                  char lo_s[16]; snprintf(lo_s,4,"%3d",lo);
//                  //char irnk1_1_s[16]; snprintf(irnk1_1_s,5,"%4d",irnk1_1);
//                  //char irnk1_2_s[16]; snprintf(irnk1_2_s,5,"%4d",irnk1_2);
//                  //char irnk2_1_s[16]; snprintf(irnk2_1_s,5,"%4d",irnk2_1);
//                  //char irnk2_2_s[16]; snprintf(irnk2_2_s,5,"%4d",irnk2_2);
//                  //char irnk3_1_s[16]; snprintf(irnk3_1_s,5,"%4d",irnk3_1);
//                  //char irnk3_2_s[16]; snprintf(irnk3_2_s,5,"%4d",irnk3_2);
//                  if(ifirst == 1){
//                    str << "\n  ccoef["<<irnkcc_s<<"] = "<<up_s<<"*det3_vdif(cpoi"<<j1_1_s<<j2_1_s<<j3_1_s<<j4_1_s<<
//                                                                           ",cpoi"<<j1_2_s<<j2_2_s<<j3_2_s<<j4_2_s<<"\n";
//                    ifirst = 0;
//                  }else{
//                    str << "\n             + "<<up_s<<"*det3_vdif(cpoi"<<j1_1_s<<j2_1_s<<j3_1_s<<j4_1_s<<
//                                                                ",cpoi"<<j1_2_s<<j2_2_s<<j3_2_s<<j4_2_s<<"\n";
//                  }
//                  str << "                            ,cpoi"<<k1_1_s<<k2_1_s<<k3_1_s<<k4_1_s<<
//                                                     ",cpoi"<<k1_2_s<<k2_2_s<<k3_2_s<<k4_2_s<<"\n";
//                  str << "                            ,cpoi"<<l1_1_s<<l2_1_s<<l3_1_s<<l4_1_s<<
//                                                     ",cpoi"<<l1_2_s<<l2_2_s<<l3_2_s<<l4_2_s<<")/"<<lo_s;
//                }
//              }
//            }
//            /*-----------------------------------K LOOP-------------------------------------------*/
//
//          }
//        }
//      }
//      str << ";\n";
//    }
//
//    str << "}\n";
//    str << "\n\n";
//    std::ofstream f;
//    char fname[64]; 
//    snprintf(fname,64,"codegen_ccoef3.%02d.cxx",ideg);
//    f.open(fname, std::ios::out);
//    f << str.str();
//    f.close();
//  }
//
//}



void gen_ccoeff3(){

  int max_jaco_deg = METRIS_MAX_DEG;

  printf("-- Tetrahedra up to degree %d \n",METRIS_MAX_DEG);
  for(int ideg = 2; ideg <= max_jaco_deg; ideg++){
    int npp_c = tetnpps[3*(ideg-1)];

    std::ostringstream str;

    str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
    str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
    str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
    str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

    str << "#include \"codegen_ccoef.hxx\"\n";
    str << "#include \"types.hxx\"\n\n";
    str << "namespace Metris{\n\n";

    str << "template<int ideg>\nvoid ccoef_genbez3(const intAr2 & __restrict__ tet2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){}\n\n";
    str << "double det3_vdif(const double* x1,const double* x2\n";
    str << "                ,const double* y1,const double* y2\n";
    str << "                ,const double* z1,const double* z2);\n\n";




    int uuu1 = pow(ifactorial_func(ideg),3);
    int lll1 = ifactorial_func(3*(ideg-1));

    str << "template<> void ccoef_genbez3<"<<ideg<<">(const intAr2 & __restrict__ tet2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){\n";
    for(int irnkcc = 0; irnkcc < npp_c; irnkcc++){
      int ifirst = 1;
      int i1 = ordtet.s[3*(ideg-1)][irnkcc][0]; // constexpr
      int i2 = ordtet.s[3*(ideg-1)][irnkcc][1]; // constexpr
      int i3 = ordtet.s[3*(ideg-1)][irnkcc][2]; // constexpr
      int i4 = ordtet.s[3*(ideg-1)][irnkcc][3]; // constexpr
      int uuu2 = ifactorial_func(i1)*ifactorial_func(i2)*ifactorial_func(i3)*ifactorial_func(i4);

      for(int j1 = 0; j1 <= ideg - 1 && j1 <= i1; j1 ++){
        for(int j2 = 0; j2 <= ideg - 1 - j1 && j2 <= i2; j2 ++){
          for(int j3 = 0; j3 <= ideg - 1 - j1 - j2 && j3 <= i3; j3 ++){
            if(ideg-1-j1-j2-j3 > i4) continue;
            int j4 = ideg-1-j1-j2-j3;

            /*-----------------------------------K LOOP-------------------------------------------*/
            for(int k1 = 0; k1 <= ideg - 1 && k1 <= i1 - j1 ; k1++){
              for(int k2 = 0; k2 <= ideg - 1 - k1 && k2 <= i2 - j2 ; k2++){
                for(int k3 = 0; k3 <= ideg - 1 - k1 - k2 && k3 <= i3 - j3 ; k3++){
                  if(ideg-1-k1-k2-k3 > i4-j4) continue;
                  int k4 = ideg - 1 - k1 - k2 - k3;
                  int l1 = i1 - j1 - k1;
                  int l2 = i2 - j2 - k2;
                  int l3 = i3 - j3 - k3;
                  int l4 = i4 - j4 - k4;

                  int cc2 = ifactorial_func(j1)*ifactorial_func(j2)*ifactorial_func(j3)*ifactorial_func(j4);
                  int cc3 = ifactorial_func(k1)*ifactorial_func(k2)*ifactorial_func(k3)*ifactorial_func(k4);
                  int cc4 = ifactorial_func(l1)*ifactorial_func(l2)*ifactorial_func(l3)*ifactorial_func(l4);
                  int lll2 = cc2*cc3*cc4;

                  int irnk1_1 = mul2nod(j1,j2+1,j3,j4);
                  int irnk1_2 = mul2nod(j1+1,j2,j3,j4);

                  int irnk2_1 = mul2nod(k1,k2,k3+1,k4);
                  int irnk2_2 = mul2nod(k1+1,k2,k3,k4);

                  int irnk3_1 = mul2nod(l1,l2,l3,l4+1);
                  int irnk3_2 = mul2nod(l1+1,l2,l3,l4);
                  int up,lo;
                  simpfrac(uuu1*uuu2,lll1*lll2,&up,&lo);

                  char irnkcc_s[16]; snprintf(irnkcc_s,4,"%3d",irnkcc);
                  char up_s[16]; snprintf(up_s,4,"%3d",up);
                  char lo_s[16]; snprintf(lo_s,4,"%3d",lo);
                  char irnk1_1_s[16]; snprintf(irnk1_1_s,5,"%4d",irnk1_1);
                  char irnk1_2_s[16]; snprintf(irnk1_2_s,5,"%4d",irnk1_2);
                  char irnk2_1_s[16]; snprintf(irnk2_1_s,5,"%4d",irnk2_1);
                  char irnk2_2_s[16]; snprintf(irnk2_2_s,5,"%4d",irnk2_2);
                  char irnk3_1_s[16]; snprintf(irnk3_1_s,5,"%4d",irnk3_1);
                  char irnk3_2_s[16]; snprintf(irnk3_2_s,5,"%4d",irnk3_2);
                  if(ifirst == 1){
                    str << "\n  ccoef["<<irnkcc_s<<"] = "<<up_s<<"*det3_vdif(coord[tet2poi[ielem]["<<irnk1_1_s<<"]],coord[tet2poi[ielem]["<<irnk1_2_s<<"]]\n";
                    ifirst = 0;
                  }else{
                    str << "\n             + "<<up_s<<"*det3_vdif(coord[tet2poi[ielem]["<<irnk1_1_s<<"]],coord[tet2poi[ielem]["<<irnk1_2_s<<"]]\n";
                  }
                  str << "                            ,coord[tet2poi[ielem]["<<irnk2_1_s<<"]],coord[tet2poi[ielem]["<<irnk2_2_s<<"]]\n";
                  str << "                            ,coord[tet2poi[ielem]["<<irnk3_1_s<<"]],coord[tet2poi[ielem]["<<irnk3_2_s<<"]])/"<<lo_s;
                }
              }
            }
            /*-----------------------------------K LOOP-------------------------------------------*/

          }
        }
      }
      str << ";\n";
    }

    str << "}\n";
    str << "\n\n";
    str << "} // End namespace\n";
    std::ofstream f;
    char fname[64]; 
    snprintf(fname,64,"codegen_ccoef3.%02d.cxx",ideg);
    f.open(fname, std::ios::out);
    f << str.str();
    f.close();
  }

}




void gen_ccoeff2(){

  int max_jaco_deg = METRIS_MAX_DEG;

  printf("-- Triangles up to degree %d \n",METRIS_MAX_DEG);
  for(int ideg = 2; ideg <= max_jaco_deg; ideg++){
    std::ostringstream str;

    str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
    str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
    str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
    str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";
    
    str << "#include \"codegen_ccoef.hxx\"\n\n";
    str << "#include \"types.hxx\"\n\n";
    str << "namespace Metris{\n\n";

    str << "template<int ideg>\nvoid ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){}\n\n";
    str << "double det2_vdif(const double* x1,const double* x2\n";
    str << "                ,const double* y1,const double* y2);\n\n";

                  
    int npp_c = facnpps[2*(ideg-1)];

    int uuu1 = pow(ifactorial_func(ideg),2);
    int lll1 = ifactorial_func(2*(ideg-1));

    str << "template<> void ccoef_genbez2<"<<ideg<<">(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){\n";
    for(int irnkcc = 0; irnkcc < npp_c; irnkcc++){
      int ifirst = 1;
      int i1 = ordfac.s[2*(ideg-1)][irnkcc][0]; // constexpr
      int i2 = ordfac.s[2*(ideg-1)][irnkcc][1]; // constexpr
      int i3 = ordfac.s[2*(ideg-1)][irnkcc][2]; // constexpr
      int uuu2 = ifactorial_func(i1)*ifactorial_func(i2)*ifactorial_func(i3);

      for(int j1 = 0; j1 <= ideg - 1 && j1 <= i1; j1 ++){
        for(int j2 = 0; j2 <= ideg - 1 - j1 && j2 <= i2; j2 ++){
          if(ideg-1-j1-j2 > i3) continue;
          int j3 = ideg-1-j1-j2;

          int k1 = i1 - j1; 
          int k2 = i2 - j2;
          int k3 = i3 - j3;

          /*-----------------------------------K LOOP-------------------------------------------*/

          int cc2 = ifactorial_func(j1)*ifactorial_func(j2)*ifactorial_func(j3);
          int cc3 = ifactorial_func(k1)*ifactorial_func(k2)*ifactorial_func(k3);
          int lll2 = cc2*cc3;

          int irnk1_1 = mul2nod(j1,j2+1,j3);
          int irnk1_2 = mul2nod(j1+1,j2,j3);

          int irnk2_1 = mul2nod(k1,k2,k3+1);
          int irnk2_2 = mul2nod(k1+1,k2,k3);

          int up,lo;
          simpfrac(uuu1*uuu2,lll1*lll2,&up,&lo);

          char irnkcc_s[16]; snprintf(irnkcc_s,4,"%3d",irnkcc);
          char up_s[16]; snprintf(up_s,4,"%3d",up);
          char lo_s[16]; snprintf(lo_s,4,"%3d",lo);
          char irnk1_1_s[16]; snprintf(irnk1_1_s,5,"%4d",irnk1_1);
          char irnk1_2_s[16]; snprintf(irnk1_2_s,5,"%4d",irnk1_2);
          char irnk2_1_s[16]; snprintf(irnk2_1_s,5,"%4d",irnk2_1);
          char irnk2_2_s[16]; snprintf(irnk2_2_s,5,"%4d",irnk2_2);
          if(ifirst == 1){
            str << "\n  ccoef["<<irnkcc_s<<"] = "<<up_s<<"*det2_vdif(coord[fac2poi[ielem]["<<irnk1_1_s<<"]],coord[fac2poi[ielem]["<<irnk1_2_s<<"]]\n";
            ifirst = 0;
          }else{
            str << "\n             + "<<up_s<<"*det2_vdif(coord[fac2poi[ielem]["<<irnk1_1_s<<"]],coord[fac2poi[ielem]["<<irnk1_2_s<<"]]\n";
          }
          str << "                            ,coord[fac2poi[ielem]["<<irnk2_1_s<<"]],coord[fac2poi[ielem]["<<irnk2_2_s<<"]])/"<<lo_s;
//          str << "                            ,coord[fac2poi[ielem]["<<irnk3_1_s<<"]],coord[fac2poi[ielem]["<<irnk3_2_s<<"]])/"<<lo_s;
          /*-----------------------------------K LOOP-------------------------------------------*/

        }
      }
      str << ";\n";
    }

    str << "}\n";
    str << "\n\n";
    str << "} // End namespace\n";
    std::ofstream f;
    char fname[64]; 
    snprintf(fname,64,"codegen_ccoef2.%02d.cxx",ideg);
    f.open(fname, std::ios::out);
    f << str.str();
    f.close();
  }

}





