//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "Mesh/MeshBase.hxx"

#include "ho_constants.hxx"
#include "types.hxx"
#include <array>

#include <cmath>

#include <string>
#include <sstream>
#include <fstream>

namespace Metris{

/* 
Bézier coefficients of the Jacobian determinant
*/

static int ifactorial_func(int i){
	if(i == 0) return 1;
	return ifactorial_func(i-1)*i;
}


//template <int i1, int i2, int i3, int i4>
//void gen_idx_ccoef(char* str);
static void simpfrac(int x, int y, int *xs, int *ys){
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

void getccoef_map(const MeshBase &msh, int ientt, double *ccoef, dblAr2 &d_ccoef) {
   int ideg = 2;

  std::ostringstream str;
  str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  str << "#include \"codegen_ccoef.hxx\"\n\n";
  str << "#include \"types.hxx\"\n\n";
  str << "#include \"low_geo.hxx\"\n\n";
  str << "namespace Metris{\n\n";

  str << "template<int ideg>\nvoid d_ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef, dblAr2& __restrict__ d_ccoef){}\n\n";
  str << "double det2_vdif(const double* x1,const double* x2\n";
  str << "                ,const double* y1,const double* y2);\n\n";
  str << "double vdiff_perp(const double* a,const double* b);\n\n";

  std::map<std::pair<int, int>, std::string> d_ccoef_map; // Declare d_ccoef_map

  int npp_c = facnpps[2*(ideg-1)];

  int uuu1 = pow(ifactorial_func(ideg),2);
  int lll1 = ifactorial_func(2*(ideg-1));

  str << "template<> void d_ccoef_genbez2<"<<ideg<<">(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef, dblAr2& __restrict__ d_ccoef){\n";
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

        int cc2 = ifactorial_func(j1)*ifactorial_func(j2)*ifactorial_func(j3);
        int cc3 = ifactorial_func(k1)*ifactorial_func(k2)*ifactorial_func(k3);
        int lll2 = cc2*cc3;

        int irnk1_1 = mul2nod(j1+1,j2,j3); // gam1  =  P_{\alpha  + e^1} in other notation
        int irnk1_2 = mul2nod(j1,j2,j3+1); // gam2 ''  P_{\alpha  + e^3}      ''

        int irnk2_1 = mul2nod(k1,k2+1,k3); // gam2 ''  P_{\alpha' + e^2}      ''
        int irnk2_2 = mul2nod(k1,k2,k3+1); // gam4 ''  P_{\alpha' + e^3}      ''
      

        int up,lo;
        simpfrac(uuu1*uuu2,lll1*lll2,&up,&lo);

        char irnkcc_s[16];  snprintf(irnkcc_s,4,"%3d",irnkcc);
        char up_s[16];      snprintf(up_s,4,"%3d",up);
        char lo_s[16];      snprintf(lo_s,4,"%3d",lo);
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
        str << "                            ,coord[fac2poi[ielem]["<<irnk2_1_s<<"]],coord[fac2poi[ielem]["<<irnk2_2_s<<"]])/"<<lo_s<<";\n\n";


        std::pair<int, int> key1(irnkcc, irnk1_1);
        std::pair<int, int> key2(irnkcc, irnk1_2);
        std::pair<int, int> key3(irnkcc, irnk2_1);
        std::pair<int, int> key4(irnkcc, irnk2_2);

        d_ccoef_map[key1] += " + vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]])";
        d_ccoef_map[key2] += " + vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]])";
        d_ccoef_map[key3] += " + vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]])";
        d_ccoef_map[key4] += " + vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]])";
    }
}
}
// NEED TO BE CAREFUL OF THE CONSTANT FACTORS IN THE END
    for (const auto& entry : d_ccoef_map) {
        const auto& key = entry.first;
        const std::string& value = entry.second.substr(3); // removes front " + "
        str << "\n d_ccoef[" << key.first << "][" << key.second << "] = " << value << ";\n";
    }

    str << "}\n";
    str << "\n\n";
    str << "} // End namespace\n";

    std::ofstream f;
    char fname[64];
    snprintf(fname, 64, "codegen_linear_comps.%02d.cxx", ideg);
    f.open(fname, std::ios::out);
    f << str.str();
    f.close();
}

void rev_getccoef_map(const MeshBase &msh, int ientt, double *ccoef, dblAr2 &d_ccoef){
  int ideg = 2;

  std::ostringstream str;
  str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  str << "#include \"codegen_ccoef.hxx\"\n\n";
  str << "#include \"types.hxx\"\n\n";
  str << "#include \"low_geo.hxx\"\n\n";
  str << "namespace Metris{\n\n";

  str << "template<int ideg>\nvoid rev_d_ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef, dblAr2& __restrict__ d_ccoef){}\n\n";
  str << "double det2_vdif(const double* x1,const double* x2\n";
  str << "                ,const double* y1,const double* y2);\n\n";
  str << "double vdiff_perp(const double* a,const double* b);\n\n";

  std::map<std::pair<int, int>, std::string> d_ccoef_map; // Declare d_ccoef_map

  int npp_c = facnpps[2*(ideg-1)];

  int uuu1 = pow(ifactorial_func(ideg),2);
  int lll1 = ifactorial_func(2*(ideg-1));

  str << "template<> void rev_d_ccoef_genbez2<"<<ideg<<">(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef, dblAr2& __restrict__ d_ccoef){\n";
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

        int cc2 = ifactorial_func(j1)*ifactorial_func(j2)*ifactorial_func(j3);
        int cc3 = ifactorial_func(k1)*ifactorial_func(k2)*ifactorial_func(k3);
        int lll2 = cc2*cc3;

        int irnk1_1 = mul2nod(j1,j2+1,j3);
        int irnk1_2 = mul2nod(j1+1,j2,j3);

        int irnk2_1 = mul2nod(k1,k2,k3+1);
        int irnk2_2 = mul2nod(k1+1,k2,k3);

        int up,lo;
        simpfrac(uuu1*uuu2,lll1*lll2,&up,&lo);

        char irnkcc_s[16];  snprintf(irnkcc_s,4,"%3d",irnkcc);
        char up_s[16];      snprintf(up_s,4,"%3d",up);
        char lo_s[16];      snprintf(lo_s,4,"%3d",lo);
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
        str << "                            ,coord[fac2poi[ielem]["<<irnk2_1_s<<"]],coord[fac2poi[ielem]["<<irnk2_2_s<<"]])/"<<lo_s<<";\n\n";


        std::pair<int, int> key1(irnkcc, irnk1_1);
        std::pair<int, int> key2(irnkcc, irnk1_2);
        std::pair<int, int> key3(irnkcc, irnk2_1);
        std::pair<int, int> key4(irnkcc, irnk2_2);

        d_ccoef_map[key1] += " + vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]])";
        d_ccoef_map[key2] += " + vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]])";
        d_ccoef_map[key3] += " + vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]])";
        d_ccoef_map[key4] += " + vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]])";
      }
    }
  }
  // NEED TO BE CAREFUL OF THE CONSTANT FACTORS IN THE END
  for (const auto& entry : d_ccoef_map) {
    const auto& key = entry.first;
    const std::string& value = entry.second.substr(3); // takes away leading " + "
    str << "\n  d_ccoef[" << key.first << "][" << key.second << "] = " << value << ";\n";
  }

  str << "}\n";
  str << "\n\n";
  str << "} // End namespace\n";

  std::ofstream f;
  char fname[64];
  snprintf(fname, 64, "rev_codegen_linear_comps.%02d.cxx", ideg);
  f.open(fname, std::ios::out);
  f << str.str();
  f.close();
}

// template<int gdim, int ideg>
// void getccoef_P2(const MeshBase &msh, int ientt, double *ccoef, dblAr2 &d_ccoef, int idx) {
//    int ideg = 2;

//   // printf("-- Triangles up to degree %d \n",METRIS_MAX_DEG);
//   // for(int ideg = 2; ideg <= max_jaco_deg; ideg++){
//   std::ostringstream str;

//   str << "#include "codegen_ccoef.hxx"\n\n";
//   str << "#include "types.hxx"\n\n";
//   str << "namespace Metris{\n\n";

//   str << "template<int ideg>\nvoid ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){}\n\n";
//   str << "double det2_vdif(const double* x1,const double* x2\n";
//   str << "                ,const double* y1,const double* y2);\n\n";

                
//   int npp_c = facnpps[2*(ideg-1)];

//   int uuu1 = pow(ifactorial_func(ideg),2);
//   int lll1 = ifactorial_func(2*(ideg-1));

//   str << "template<> void ccoef_genbez2<"<<ideg<<">(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){\n";
//   for(int irnkcc = 0; irnkcc < npp_c; irnkcc++){
//     int ifirst = 1;
//     int i1 = ordfac.s[2*(ideg-1)][irnkcc][0]; // constexpr
//     int i2 = ordfac.s[2*(ideg-1)][irnkcc][1]; // constexpr
//     int i3 = ordfac.s[2*(ideg-1)][irnkcc][2]; // constexpr
//     int uuu2 = ifactorial_func(i1)*ifactorial_func(i2)*ifactorial_func(i3);
//     printf("--------------------------------------------------------------\n");
//     printf("BETA = (%d%d%d)\n", i1, i2, i3);
// //                         beta        = (i1, i2, i3), 
// //                         alpha       = (j1, j2, j3),
// //                         alpha_prime = (k1, k2, k3) 
//     for(int j1 = 0; j1 <= ideg - 1 && j1 <= i1; j1 ++){
//       for(int j2 = 0; j2 <= ideg - 1 - j1 && j2 <= i2; j2 ++){
//         if(ideg-1-j1-j2 > i3) continue;
//         int j3 = ideg-1-j1-j2;
      
//         int k1 = i1 - j1; 
//         int k2 = i2 - j2;
//         int k3 = i3 - j3;

//         printf("a = (%d%d%d) \n", j1, j2, j3);
//         printf("a' = (%d%d%d)\n", k1, k2, k3);

//         int cc2 = ifactorial_func(j1)*ifactorial_func(j2)*ifactorial_func(j3);
//         int cc3 = ifactorial_func(k1)*ifactorial_func(k2)*ifactorial_func(k3);
//         int lll2 = cc2*cc3;

        
// //         terms in the matrix: w/ index notation (1-u-v,u,v) instead of thesis: (u,v,1-u-v) 
// //         need to be careful of this notation in the end

//         // int dalpha_e1      = vdiff(gam3,gam4);
//         // int dalpha_e3      = -dalpha_e1; // or ? -dalpha_e1
//         // int dalphaprime_e2 = vdiff(gam2,gam1);
//         // int dalphaprime_e3 = -dalphaprime_e2;

//         // int irnk1_1 = mul2nod(j1,j2+1,j3); // gam1  =  P_{\alpha  + e^1} in other notation
//         // int irnk1_2 = mul2nod(j1+1,j2,j3); // gam2 ''  P_{\alpha  + e^3}      ''

//         // int irnk2_1 = mul2nod(k1,k2,k3+1); // gam2 ''  P_{\alpha' + e^2}      ''
//         // int irnk2_2 = mul2nod(k1+1,k2,k3); // gam4 ''  P_{\alpha' + e^3}      ''

//         int irnk1_1 = mul2nod(j1+1,j2,j3); // gam1  =  P_{\alpha  + e^1} in other notation
//         int irnk1_2 = mul2nod(j1,j2,j3+1); // gam2 ''  P_{\alpha  + e^3}      ''

//         int irnk2_1 = mul2nod(k1,k2+1,k3); // gam2 ''  P_{\alpha' + e^2}      ''
//         int irnk2_2 = mul2nod(k1,k2,k3+1); // gam4 ''  P_{\alpha' + e^3}      ''

//         printf("irnkX_X = (%d%d%d%d)\n", irnk1_1, irnk1_2, irnk2_1, irnk2_2);

//         // double dalpha_e1      = vdiff_perp(mul2nod(k1, k2+1, k3),mul2nod(k1, k2, k3+1));  // (P_{\alpha' + e^2} - P_{\alpha' + e^3})^{\perp}
//         // double dalpha_e3      = vdiff_perp(mul2nod(k1, k2, k3+1),mul2nod(k1, k2+1, k3));  // (P_{\alpha' + e^3} - P_{\alpha' + e^2})^{\perp}
//         // double dalphaprime_e2 = vdiff_perp(mul2nod(j1, j2, j3+1),mul2nod(j1+1, j2, j3));  // (P_{\alpha  + e^3} - P_{\alpha  + e^1})^{\perp}
//         // double dalphaprime_e3 = vdiff_perp(mul2nod(j1+1, j2, j3),mul2nod(j1, j2, j3+1));  // (P_{\alpha  + e^1} - P_{\alpha  + e^3})^{\perp}

//         // double d_gam1 = (irnk1_1 == irnk2_1) ? (dalpha_e1 + dalphaprime_e2) : 
//         //                 (irnk1_1 == irnk2_2) ? (dalpha_e1 + dalphaprime_e3) : 
//         //                 dalpha_e1;

//         // double d_gam2 = (irnk1_2 == irnk2_1) ? (dalpha_e3 + dalphaprime_e2) : 
//         //                 (irnk1_2 == irnk2_2) ? (dalpha_e3 + dalphaprime_e3) :
//         //                 dalpha_e3;

//         // REMINDER: include proper terms from up & lo in derivative expressions

//         int up,lo;
//         simpfrac(uuu1*uuu2,lll1*lll2,&up,&lo);

//         char irnkcc_s[16];  snprintf(irnkcc_s,4,"%3d",irnkcc);
//         char up_s[16];      snprintf(up_s,4,"%3d",up);
//         char lo_s[16];      snprintf(lo_s,4,"%3d",lo);
//         char irnk1_1_s[16]; snprintf(irnk1_1_s,5,"%4d",irnk1_1);
//         char irnk1_2_s[16]; snprintf(irnk1_2_s,5,"%4d",irnk1_2);
//         char irnk2_1_s[16]; snprintf(irnk2_1_s,5,"%4d",irnk2_1);
//         char irnk2_2_s[16]; snprintf(irnk2_2_s,5,"%4d",irnk2_2);

//       (irnk1_1 == irnk2_1) ? str << "\n d_ccoef[["<<irnkcc_s<<"]["<<irnk1_1_s<<"]] = vdiff_perp(coord[fac2poi[ielem]]["<<irnk2_1_s<<"], coord[fac2poi[ielem]]["<<irnk2_2_s<<"])+vdiff_perp(coord[fac2poi[ielem]["<<irnk1_2_s<<"], coord[fac2poi[ielem]]["<<irnk1_1_s<<"]);" :  
//       (irnk1_1 == irnk2_2) ? str << "\n d_ccoef[[" << irnkcc_s << "][" << irnk1_2_s << "]] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_1_s << "], coord[fac2poi[ielem]][" << irnk2_2_s << "]) + vdiff_perp(coord[fac2poi[ielem]][" << irnk1_1_s << "], coord[fac2poi[ielem]][" << irnk1_2_s << "]);" : 
//       str << "\n d_ccoef[[" << irnkcc_s << "][" << irnk1_1_s << "]] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_1_s << "], coord[fac2poi[ielem]][" << irnk2_2_s << "]);" ;
      
    
//       (irnk1_2 == irnk2_1) ? str << "\n d_ccoef[[" << irnkcc_s << "][" << irnk2_1_s << "]] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_2_s << "], coord[fac2poi[ielem]][" << irnk2_1_s << "]) + vdiff_perp(coord[fac2poi[ielem]][" << irnk1_1_s << "], coord[fac2poi[ielem]][" << irnk1_2_s << "]);" :
//       (irnk1_2 == irnk2_2) ? str << "\n d_ccoef[[" << irnkcc_s << "][" << irnk2_2_s << "]] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_2_s << "], coord[fac2poi[ielem]][" << irnk2_1_s << "]) + vdiff_perp(coord[fac2poi[ielem]][" << irnk1_2_s << "], coord[fac2poi[ielem]][" << irnk1_1_s << "]);" :
//       str << "\n d_ccoef[[" << irnkcc_s << "][" << irnk1_2_s << "]] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_2_s << "], coord[fac2poi[ielem]][" << irnk2_1_s << "]);";
          
                
//       }
      
//     }
//     str << ";\n";
//   }
//   str << "}\n";
//   str << "\n\n";
//   str << "} // End namespace\n";
//   std::ofstream f;
//   char fname[64]; 
//   snprintf(fname,64,"codegen_linear_comps.%02d.cxx",ideg);
//   f.open(fname, std::ios::out);
//   f << str.str();
//   f.close(); 
// }


// void getccoef_naive(const MeshBase &msh, int ientt, double *ccoef, dblAr2 &d_ccoef, int idx) {
//    int ideg = 2;

//   std::ostringstream str;

//   str << "#include "codegen_ccoef.hxx"\n\n";
//   str << "#include "types.hxx"\n\n";
//   str << "namespace Metris{\n\n";

//   str << "template<int ideg>\nvoid ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){}\n\n";
//   str << "double det2_vdif(const double* x1,const double* x2\n";
//   str << "                ,const double* y1,const double* y2);\n\n";

//   std::map<std::pair<int, int>, std::string> d_ccoef_map; // Declare d_ccoef_map

//   int npp_c = facnpps[2*(ideg-1)];

//   int uuu1 = pow(ifactorial_func(ideg),2);
//   int lll1 = ifactorial_func(2*(ideg-1));

//   str << "template<> void d_ccoef_genbez2<"<<ideg<<">(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){\n";
//   for(int irnkcc = 0; irnkcc < npp_c; irnkcc++){
//     int ifirst = 1;
//     int i1 = ordfac.s[2*(ideg-1)][irnkcc][0]; // constexpr
//     int i2 = ordfac.s[2*(ideg-1)][irnkcc][1]; // constexpr
//     int i3 = ordfac.s[2*(ideg-1)][irnkcc][2]; // constexpr
//     int uuu2 = ifactorial_func(i1)*ifactorial_func(i2)*ifactorial_func(i3);
//     printf("--------------------------------------------------------------\n");
//     printf("BETA = (%d%d%d)\n", i1, i2, i3);
// //                         beta        = (i1, i2, i3), 
// //                         alpha       = (j1, j2, j3),
// //                         alpha_prime = (k1, k2, k3) 
//     for(int j1 = 0; j1 <= ideg - 1 && j1 <= i1; j1 ++){
//       for(int j2 = 0; j2 <= ideg - 1 - j1 && j2 <= i2; j2 ++){
//         if(ideg-1-j1-j2 > i3) continue;
//         int j3 = ideg-1-j1-j2;
      
//         int k1 = i1 - j1; 
//         int k2 = i2 - j2;
//         int k3 = i3 - j3;

//         printf("a = (%d%d%d) \n", j1, j2, j3);
//         printf("a' = (%d%d%d)\n", k1, k2, k3);

//         int cc2 = ifactorial_func(j1)*ifactorial_func(j2)*ifactorial_func(j3);
//         int cc3 = ifactorial_func(k1)*ifactorial_func(k2)*ifactorial_func(k3);
//         int lll2 = cc2*cc3;

//         int irnk1_1 = mul2nod(j1+1,j2,j3); // gam1  =  P_{\alpha  + e^1} in other notation
//         int irnk1_2 = mul2nod(j1,j2,j3+1); // gam2 ''  P_{\alpha  + e^3}      ''

//         int irnk2_1 = mul2nod(k1,k2+1,k3); // gam2 ''  P_{\alpha' + e^2}      ''
//         int irnk2_2 = mul2nod(k1,k2,k3+1); // gam4 ''  P_{\alpha' + e^3}      ''

//         printf("irnkX_X = (%d%d%d%d)\n", irnk1_1, irnk1_2, irnk2_1, irnk2_2);

//         int up,lo;
//         simpfrac(uuu1*uuu2,lll1*lll2,&up,&lo);

//         char irnkcc_s[16];  snprintf(irnkcc_s,4,"%3d",irnkcc);
//         char up_s[16];      snprintf(up_s,4,"%3d",up);
//         char lo_s[16];      snprintf(lo_s,4,"%3d",lo);
//         char irnk1_1_s[16]; snprintf(irnk1_1_s,5,"%4d",irnk1_1);
//         char irnk1_2_s[16]; snprintf(irnk1_2_s,5,"%4d",irnk1_2);
//         char irnk2_1_s[16]; snprintf(irnk2_1_s,5,"%4d",irnk2_1);
//         char irnk2_2_s[16]; snprintf(irnk2_2_s,5,"%4d",irnk2_2);

//       // these will compute overlaps, and the commented lines fix things 
//       str << "\n d_ccoef[" << irnkcc_s << "][" << irnk1_1_s << "] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_1_s << "], coord[fac2poi[ielem]][" << irnk2_2_s << "]);" ;
//       str << "\n d_ccoef[" << irnkcc_s << "][" << irnk1_2_s << "] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_2_s << "], coord[fac2poi[ielem]][" << irnk2_1_s << "]);" ;
//       str << "\n d_ccoef[" << irnkcc_s << "][" << irnk2_1_s << "] = vdiff_perp(coord[fac2poi[ielem]][" << irnk1_2_s << "], coord[fac2poi[ielem]][" << irnk1_1_s << "]);" ;
//       str << "\n d_ccoef[" << irnkcc_s << "][" << irnk2_2_s << "] = vdiff_perp(coord[fac2poi[ielem]][" << irnk1_1_s << "], coord[fac2poi[ielem]][" << irnk1_2_s << "]);" ;
    
//       // (irnk1_1 == irnk2_1) ? str << "\n d_ccoef[["<<irnkcc_s<<"]["<<irnk1_1_s<<"]] = vdiff_perp(coord[fac2poi[ielem]]["<<irnk2_1_s<<"], coord[fac2poi[ielem]]["<<irnk2_2_s<<"])+vdiff_perp(coord[fac2poi[ielem]["<<irnk1_2_s<<"], coord[fac2poi[ielem]]["<<irnk1_1_s<<"]);" :  
//       // (irnk1_1 == irnk2_2) ? str << "\n d_ccoef[[" << irnkcc_s << "][" << irnk1_2_s << "]] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_1_s << "], coord[fac2poi[ielem]][" << irnk2_2_s << "]) + vdiff_perp(coord[fac2poi[ielem]][" << irnk1_1_s << "], coord[fac2poi[ielem]][" << irnk1_2_s << "]);" : 
//       // str << "\n d_ccoef[[" << irnkcc_s << "][" << irnk1_1_s << "]] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_1_s << "], coord[fac2poi[ielem]][" << irnk2_2_s << "]);" ;
    
//       // (irnk1_2 == irnk2_1) ? str << "\n d_ccoef[[" << irnkcc_s << "][" << irnk2_1_s << "]] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_2_s << "], coord[fac2poi[ielem]][" << irnk2_1_s << "]) + vdiff_perp(coord[fac2poi[ielem]][" << irnk1_1_s << "], coord[fac2poi[ielem]][" << irnk1_2_s << "]);" :
//       // (irnk1_2 == irnk2_2) ? str << "\n d_ccoef[[" << irnkcc_s << "][" << irnk2_2_s << "]] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_2_s << "], coord[fac2poi[ielem]][" << irnk2_1_s << "]) + vdiff_perp(coord[fac2poi[ielem]][" << irnk1_2_s << "], coord[fac2poi[ielem]][" << irnk1_1_s << "]);" :
//       // str << "\n d_ccoef[[" << irnkcc_s << "][" << irnk1_2_s << "]] = vdiff_perp(coord[fac2poi[ielem]][" << irnk2_2_s << "], coord[fac2poi[ielem]][" << irnk2_1_s << "]);";
          
                
//       }
      
//     }
//     str << ";\n";
//   }
//   str << "}\n";
//   str << "\n\n";
//   str << "} // End namespace\n";
//   std::ofstream f;
//   char fname[64]; 
//   snprintf(fname,64,"codegen_linear_comps.%02d.cxx",ideg);
//   f.open(fname, std::ios::out);
//   f << str.str();
//   f.close(); 
// }            
}