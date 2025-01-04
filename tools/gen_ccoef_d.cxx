//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "gen_ccoef.hxx"
#include "../src/ho_constants.hxx"
#include "../src/aux_exceptions.hxx"

#include <sstream>
#include <fstream>

using namespace Metris;

static int ifactorial_func(int i){
  if(i == 0) return 1;
  return ifactorial_func(i-1)*i;
}

void gen_ccoeff2_d(){
  int max_jaco_deg = METRIS_MAX_DEG;
  printf("-- Triangles ccoef_d up to degree %d \n",METRIS_MAX_DEG);
  for(int ideg = 2; ideg <= max_jaco_deg; ideg++){
    getccoef2_map_coord(ideg);
    if(ideg > 2){
      printf("Degrees > 2 not tested yet \n");
      break;
    }
  }
}

void getccoef2_map(int ideg){

  // Other degrees not tested yet
  METRIS_ENFORCE(ideg == 2);

  std::ostringstream str;
  str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  str << "#include \"codegen_ccoef_d.hxx\"\n\n";
  str << "#include \"types.hxx\"\n\n";
  //str << "#include \"low_geo.hxx\"\n\n";
  str << "namespace Metris{\n\n";

  str << "template<int ideg>\nvoid d_ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double*__restrict__ ccoef, dblAr2& __restrict__ d_ccoef){}\n\n";
  str << "double det2_vdif(const double* x1,const double* x2\n";
  str << "                ,const double* y1,const double* y2);\n\n";
  str << "double vdiff_perp(const double* a,const double* b);\n\n";

  std::map<std::pair<int, int>, std::string> d_ccoef_map; // Declare d_ccoef_map

  int npp_c = facnpps[2*(ideg-1)];

  int uuu1 = pow(ifactorial_func(ideg),2);
  int lll1 = ifactorial_func(2*(ideg-1));

  str << "template<> void d_ccoef_genbez2<"<<ideg<<">(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double*__restrict__ ccoef, dblAr2& __restrict__ d_ccoef){\n";
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

        d_ccoef_map[key1] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]])/"+std::string(lo_s)+"";
        d_ccoef_map[key2] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]])/"+std::string(lo_s)+"";
        d_ccoef_map[key3] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]])/"+std::string(lo_s)+"";
        d_ccoef_map[key4] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]])/"+std::string(lo_s)+"";
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
  snprintf(fname, 64, "codegen_ccoef2.%02d_d.cxx", ideg);
  f.open(fname, std::ios::out);
  f << str.str();
  f.close();
}

void rev_getccoef2_map(int ideg){
  // Other degrees not tested yet
  METRIS_ENFORCE(ideg == 2);

  std::ostringstream str;
  str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  str << "#include \"codegen_ccoef_d.hxx\"\n\n";
  str << "#include \"types.hxx\"\n\n";
  //str << "#include \"low_geo.hxx\"\n\n";
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

        d_ccoef_map[key1] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]])";
        d_ccoef_map[key2] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]])";
        d_ccoef_map[key3] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]])";
        d_ccoef_map[key4] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]])";
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
  //snprintf(fname, 64, "rev_codegen_linear_comps.%02d.cxx", ideg);
  snprintf(fname, 64, "codegen_ccoef2.%02d_d.cxx", ideg);
  f.open(fname, std::ios::out);
  f << str.str();
  f.close();
}

void getccoef2_map_coord(int ideg){
  // Other degrees not tested yet
  METRIS_ENFORCE(ideg == 2);

  std::ostringstream str;
  str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  str << "#include \"codegen_ccoef_d.hxx\"\n\n";
  str << "#include \"types.hxx\"\n\n";
  //str << "#include \"low_geo.hxx\"\n\n";
  str << "namespace Metris{\n\n";

  str << "template<int ideg> void d_ccoef_genbez2"
      << "(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,\n" 
      << "                                        int ielem, double*__restrict__ ccoef, \n"
      << "                                        int icoor, dblAr2& __restrict__ d_ccoef){}\n";

  str << "double det2_vdif(const double* x1,const double* x2\n";
  str << "                ,const double* y1,const double* y2);\n\n";
  str << "double vdiff_perp_x(const double* a,const double* b);\n\n";
  str << "double vdiff_perp_y(const double* a,const double* b);\n\n";

  std::map<std::pair<int, int>, std::string> d_ccoef_map_x; // Declare d_ccoef_map_x
  std::map<std::pair<int, int>, std::string> d_ccoef_map_y;

  int npp_c = facnpps[2*(ideg-1)];

  int uuu1 = pow(ifactorial_func(ideg),2);
  int lll1 = ifactorial_func(2*(ideg-1));

  str << "template<> void d_ccoef_genbez2<"<<ideg<<">"
      << "(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,\n" 
      << "                                   int ielem, double*__restrict__ ccoef, \n"
      << "                                   int icoor, dblAr2& __restrict__ d_ccoef){\n";
  str << "  METRIS_ASSERT(icoor == 0 || icoor == 1);\n";
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
          str << "                            ,coord[fac2poi[ielem]["<<irnk2_1_s<<"]],coord[fac2poi[ielem]["<<irnk2_2_s<<"]])/"<<lo_s;


        std::pair<int, int> key1(irnkcc, irnk1_1);
        std::pair<int, int> key2(irnkcc, irnk1_2);
        std::pair<int, int> key3(irnkcc, irnk2_1);
        std::pair<int, int> key4(irnkcc, irnk2_2);
        
        d_ccoef_map_x[key1] += " + vdiff_perp_x(coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]])";
        d_ccoef_map_x[key2] += " + vdiff_perp_x(coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]])";
        d_ccoef_map_x[key3] += " + vdiff_perp_x(coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]])";
        d_ccoef_map_x[key4] += " + vdiff_perp_x(coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]])";

        d_ccoef_map_y[key1] += " + vdiff_perp_y(coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]])";
        d_ccoef_map_y[key2] += " + vdiff_perp_y(coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]])";
        d_ccoef_map_y[key3] += " + vdiff_perp_y(coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]])";
        d_ccoef_map_y[key4] += " + vdiff_perp_y(coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]])";
      }
    }
    str << ";\n";
  }
  str << "\n";
  // NEED TO BE CAREFUL OF THE CONSTANT FACTORS IN THE END
  str<<"  if(icoor == 0){\n";
  for (const auto& entry : d_ccoef_map_x) {
    const auto& key = entry.first;
    int factor;
    factor = (key.first == 0 || key.first == 1 || key.first == 2) ? 4 : 2;
    char factor_s[16];  snprintf(factor_s,4,"%3d",factor);
    const std::string& value = entry.second.substr(3); // removes front " + "
    str << "\n    d_ccoef[" << key.first << "][" << key.second << "] = " << factor_s << "*(" << value << ");\n";
  }
  str<<"\n  }else{\n";
  for (const auto& entry : d_ccoef_map_y) {
    const auto& key = entry.first;
    int factor;
    factor = (key.first == 0 || key.first == 1 || key.first == 2) ? 4 : 2;
    char factor_s[16];  snprintf(factor_s,4,"%3d",factor);
    const std::string& value = entry.second.substr(3); // removes front " + "
    str << "\n    d_ccoef[" << key.first << "][" << key.second << "] = " << factor_s << "*(" << value << ");\n";
  }
  str<<"  }\n";

  str << "}\n";
  str << "\n\n";
  str << "} // End namespace\n";

  std::ofstream f;
  char fname[64];
  snprintf(fname, 64, "codegen_ccoef2.%02d_d.cxx", ideg);
  f.open(fname, std::ios::out);
  f << str.str();
  f.close();
}

void gen_ccoef2_d_pt(){
    int max_jaco_deg = METRIS_MAX_DEG;
    printf("-- Triangles ccoef_d up to degree %d \n",METRIS_MAX_DEG);
    for(int ideg = 2; ideg <= max_jaco_deg; ideg++){
        get_point_derivatives(ideg);
        if(ideg > 2){
            printf("Degrees > 2 not tested yet \n");
            break;
        }
    } 
}


void get_point_derivatives(int ideg){
  // Other degrees not tested yet
  METRIS_ENFORCE(ideg == 2);

  std::ostringstream str;
  str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  str << "#include \"codegen_ccoef_d.hxx\"\n\n";
  str << "#include \"types.hxx\"\n\n";
  //str << "#include \"low_geo.hxx\"\n\n";
  str << "namespace Metris{\n\n";

  str << "template<int ideg> void d_pt_ccoef_genbez2"
  << "(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,\n" 
  << "                                        int ielem,\n"
  << "                                        int icoor,\n"
  << "                                        int inode,\n"
  << "                                        double*__restrict__ ccoef,\n"
  << "                                        double*__restrict__ d_ccoef){}\n\n";

  str << "double det2_vdif(const double* x1,const double* x2\n";
  str << "                ,const double* y1,const double* y2);\n\n";

  str << "double* vdiff_perp(const double* a,const double* b);\n\n";
 

  std::map<std::pair<int, int>, std::string> d_ccoef_map; 

  int npp_c = facnpps[2*(ideg-1)];

  int uuu1 = pow(ifactorial_func(ideg),2);
  int lll1 = ifactorial_func(2*(ideg-1));

  str << "template<> void d_pt_ccoef_genbez2<"<<ideg<<">"
      << "(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,\n" 
      << "                                        int ielem,\n"
      << "                                        int icoor,\n"
      << "                                        int inode,\n"
      << "                                        double*__restrict__ ccoef,\n"
      << "                                        double*__restrict__ d_ccoef){\n";

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
          str << "                            ,coord[fac2poi[ielem]["<<irnk2_1_s<<"]],coord[fac2poi[ielem]["<<irnk2_2_s<<"]])/"<<lo_s;


        std::pair<int, int> key1(irnkcc, irnk1_1);
        std::pair<int, int> key2(irnkcc, irnk1_2);
        std::pair<int, int> key3(irnkcc, irnk2_1);
        std::pair<int, int> key4(irnkcc, irnk2_2);
        
        d_ccoef_map[key1] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]])[icoor] / " + std::string(lo_s);
        d_ccoef_map[key2] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]])[icoor] / " + std::string(lo_s);
        d_ccoef_map[key3] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]])[icoor] / " + std::string(lo_s);
        d_ccoef_map[key4] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]])[icoor] / " + std::string(lo_s);

      }
    }
    str << ";\n";
  }
  str << "\n";
  str << "    for(int i = 0; i < " << npp_c << "; i++) {\n";
  str << "      d_ccoef[i] = 0;\n";
  str << "    }\n";

  for (int inode = 0; inode < 6; inode++){ 
    char inode_s[16];  snprintf(inode_s,4,"%3d",inode);
    if (inode==0){
      str<<"    if(inode  == "<<inode_s<<"){\n";
    }else{
      str<<"    else if(inode  == "<<inode_s<<"){\n";
    }
    for (const auto& entry : d_ccoef_map) {
        const auto& key = entry.first;
        if (key.first==inode)
        {
            const std::string& value = entry.second.substr(3); // removes front " + "
            str << "\n      d_ccoef[" << key.second << "] = " << value << ";";
        }
    }
    str<<"\n    }\n";
  }
  
  str << "\n  } \n \n";
  str << "} // End namespace\n";

  std::ofstream f;
  char fname[64];
  snprintf(fname, 64, "codegen_ccoef2.%02d_d_pt.cxx", ideg);
  f.open(fname, std::ios::out);
  f << str.str();
  f.close();
}


void gen_ccoef3_d(){
    int max_jaco_deg = METRIS_MAX_DEG;
    printf("-- Triangles ccoef_d up to degree %d \n",METRIS_MAX_DEG);
    for(int ideg = 2; ideg <= max_jaco_deg; ideg++){
        get_ccoeff3d(ideg);
        if(ideg > 2){
            printf("Degrees > 2 not tested yet \n");
            break;
        }
    } 
}

void get_ccoeff3d(int ideg){

  int max_jaco_deg = METRIS_MAX_DEG;
  // Other degrees not tested yet
  METRIS_ENFORCE(ideg == 2);

  printf("-- Tetrahedra up to degree %d \n",METRIS_MAX_DEG);
  for(int ideg = 2; ideg <= max_jaco_deg; ideg++){
    int npp_c = tetnpps[3*(ideg-1)];

    std::ostringstream str;

    str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
    str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
    str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
    str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

    str << "#include <src/codegen_ccoef_d.hxx>\n";
    str << "#include <src/types.hxx>\n\n";
    //str << "#include <src/low_geo.hxx>\n\n";

    str << "namespace Metris{\n\n";

    str << "double det3_vdif(const double* x1,const double* x2\n";
    str << "                ,const double* y1,const double* y2\n";
    str << "                ,const double* z1,const double* z2);\n\n";

    str << "double* vdiff(const double* a,const double* b);\n\n";
    str << "double* vproduct(const double* a, const double* b);\n\n";


    std::map<std::pair<int, int>, std::string> d_ccoef_map; // Declare d_ccoef_map_x

    int uuu1 = pow(ifactorial_func(ideg),3);
    int lll1 = ifactorial_func(3*(ideg-1));

    str << "template<int ideg, typename T>\n";
    str << "void d_ccoef_genbez3(const intAr2 & __restrict__ tet2poi,\n";
    str << "                     const dblAr2& __restrict__ coord,\n";
    str << "                     int ielem,\n";
    str << "                     int icoor,\n";
    str << "                     T* __restrict__ ccoef,\n";
    str << "                     dblAr2& __restrict__ d_ccoef){\n\n";

    str << "  METRIS_ASSERT(ideg==2);\n";

    
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

                  std::pair<int, int> key1(irnkcc, irnk1_1);
                  std::pair<int, int> key2(irnkcc, irnk1_2);

                  std::pair<int, int> key3(irnkcc, irnk2_1);
                  std::pair<int, int> key4(irnkcc, irnk2_2);

                  std::pair<int, int> key5(irnkcc, irnk3_1);
                  std::pair<int, int> key6(irnkcc, irnk3_2);

    
                  d_ccoef_map[key1] += " + " + std::string(up_s) + "*vproduct(vdiff(coord[tet2poi[ielem][" +  std::string(irnk2_1_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk2_2_s) + "]]), vdiff(coord[tet2poi[ielem][" +  std::string(irnk3_1_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk3_2_s) + "]]))[icoor] / "+ std::string(lo_s);
                  d_ccoef_map[key2] += " - " + std::string(up_s) + "*vproduct(vdiff(coord[tet2poi[ielem][" +  std::string(irnk2_2_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk2_1_s) + "]]), vdiff(coord[tet2poi[ielem][" +  std::string(irnk3_2_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk3_1_s) + "]]))[icoor] / "+ std::string(lo_s);

                  d_ccoef_map[key3] += " - " + std::string(up_s) + "*vproduct(vdiff(coord[tet2poi[ielem][" +  std::string(irnk1_2_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk1_1_s) + "]]), vdiff(coord[tet2poi[ielem][" +  std::string(irnk3_2_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk3_1_s) + "]]))[icoor] / "+ std::string(lo_s);
                  d_ccoef_map[key4] += " + " + std::string(up_s) + "*vproduct(vdiff(coord[tet2poi[ielem][" +  std::string(irnk1_1_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk1_2_s) + "]]), vdiff(coord[tet2poi[ielem][" +  std::string(irnk3_1_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk3_2_s) + "]]))[icoor] / "+ std::string(lo_s);

                  d_ccoef_map[key5] += " + " + std::string(up_s) + "*vproduct(vdiff(coord[tet2poi[ielem][" +  std::string(irnk1_1_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk1_2_s) + "]]), vdiff(coord[tet2poi[ielem][" +  std::string(irnk2_1_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk2_2_s) + "]]))[icoor] / "+ std::string(lo_s);
                  d_ccoef_map[key6] += " - " + std::string(up_s) + "*vproduct(vdiff(coord[tet2poi[ielem][" +  std::string(irnk1_2_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk1_1_s) + "]]), vdiff(coord[tet2poi[ielem][" +  std::string(irnk2_2_s)+ "]], coord[tet2poi[ielem]["  + std::string(irnk2_1_s) + "]]))[icoor] / "+ std::string(lo_s);

                }
              }
            }
            /*-----------------------------------K LOOP-------------------------------------------*/
          }
        }
      }
      str << ";\n";
    }
    str << "\n";

    str << "  for(int i = 0; i < " << npp_c << "; i++) {\n";
    str << "    for(int j = 0; j < " << 10 << "; j++) {\n";
    str << "      d_ccoef[i][j] = 0;\n";
    str << "    }\n  }";

    for (const auto& entry : d_ccoef_map) {
      const auto& key = entry.first;
      str << "\n    d_ccoef[" << key.first << "][" << key.second << "] = " << entry.second << ";\n";
    }    

    str << "}\n";

    str << "template ";
    str << "void d_ccoef_genbez3<2,double>(const intAr2 & __restrict__ tet2poi,\n";
    str << "                     const dblAr2& __restrict__ coord,\n";
    str << "                     int ielem,\n";
    str << "                     int icoor,\n";
    str << "                     double* __restrict__ ccoef,\n";
    str << "                     dblAr2& __restrict__ d_ccoef);\n\n";

    //str << "template ";
    //str << "void d_ccoef_genbez3<2,float8>(const intAr2 & __restrict__ tet2poi,\n";
    //str << "                     const dblAr2& __restrict__ coord,\n";
    //str << "                     int ielem,\n";
    //str << "                     int icoor,\n";
    //str << "                     float8* __restrict__ ccoef,\n";
    //str << "                     dblAr2& __restrict__ d_ccoef);\n\n";


    str << "\n\n";
    str << "} // End namespace\n";
    std::ofstream f;
    char fname[64]; 
    snprintf(fname,64,"codegen_ccoef3.%02d_d.cxx",ideg);
    f.open(fname, std::ios::out);
    f << str.str();
    f.close();
  }

}

