//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "gen_ccoef.hxx"
#include "../src/ho_constants.hxx"
#include "../src/aux_exceptions.hxx"
#include <../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h>

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


void insert_map2(std::map<std::pair<int, int>, std::vector<SANS::DLA::VectorS<4, int>>> &d_ccoef_map,
            std::pair<int,int> &key, int up, int lo, int irnk1, int irnk2){

  // Might be necessary to call oneself if an entry is modified (could then cancel/combine with new entries)
  int new_entry[4];

  //printf(" - Insert into key %d %d : %d %d %d %d \n",key.first,key.second,up,lo,irnk1,irnk2);

  for(auto& entry : d_ccoef_map[key]){
    if(entry[2] < 0) continue;

    int up2 = entry[0];
    int lo2 = entry[1];
    if(up*lo2 == up2*lo && up != up2){
      printf("## ERROR SIMPLIFYING FRACTIONS\n");
      exit(1);
    }
    
    //if(key.first == 5 && key.second == 5) printf("Debug irnk1 irnk2 %d %d up %d lo %d entry[2] %d entry[3] %d \n",
    //  irnk1,irnk2,up,lo,entry[2],entry[3]);
    // Cases where there is simplification    
    // Recall the terms are irnk1 - irnk2, or entry[2] - entry[3]
    if(entry[3] == irnk1 && up*lo2 == up2*lo){
      entry[3] = irnk2;
      goto makenew;
    }
    if(entry[2] == irnk2 && up*lo2 == up2*lo){
      entry[2] = irnk1;
      goto makenew;
    }
    // Sum the terms
    if(entry[2] == irnk1 && entry[3] == irnk2){
      int upn = up*lo2 + up2*lo;
      int lon = lo2*lo;
      int upn2, lon2;
      simpfrac(upn,lon,&upn2,&lon2);
      entry[0] = upn;
      entry[1] = lon;
      goto makenew;
    }
    continue;

    makenew:
    new_entry[0] = entry[0];
    new_entry[1] = entry[1];
    new_entry[2] = entry[2];
    new_entry[3] = entry[3];
    // Kill this entry 
    entry[2] = -1;
    entry[3] = -1;
    //printf("Debug killed vector is length %zu \n",d_ccoef_map[key].size());
    //for(int ii = 0; ii < d_ccoef_map[key].size(); ii++){
    //  printf(" %d : %d %d %d %d \n",ii,d_ccoef_map[key][ii][0],d_ccoef_map[key][ii][1],d_ccoef_map[key][ii][2],d_ccoef_map[key][ii][3]);
    //}
    insert_map2(d_ccoef_map,key,new_entry[0],new_entry[1],new_entry[2],new_entry[3]);
    return;
  }


  d_ccoef_map[key].push_back(SANS::DLA::VectorS<4, int>{up,lo,irnk1,irnk2});
  return;

}


void getccoef2_map_coord(int ideg){
  // Other degrees not tested yet
  if(ideg != 2){
    printf("##Degrees > 2 not tested yet: write unit test for getccoef2_d\n");
  }

  std::ostringstream str;
  str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  str << "#include \"codegen_ccoef_d.hxx\"\n";
  str << "#include \"types.hxx\"\n\n\n";
  //str << "#include \"low_geo.hxx\"\n\n";
  str << "namespace Metris{\n\n\n";

  str << "template<int ideg> void d_ccoef_genbez2"
      << "(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,\n" 
      << "                                        int ielem, int icoor,\n"
      << "                                        dblAr2& __restrict__ d_ccoef){}\n";


  //std::map<std::pair<int, int>, std::string> d_ccoef_map_x; // Declare d_ccoef_map_x
  //std::map<std::pair<int, int>, std::string> d_ccoef_map_y;

  // Store up, lo, irnk1, irnk2 s.t. term is up*(coord(irnk1,1-icoor) - coord(irnk2,1-icoor))/lo
  // Store all entries to be summed for a given pair idof icoef
  std::map<std::pair<int, int>, std::vector<SANS::DLA::VectorS<4, int>>> d_ccoef_map;

  int npp_c = facnpps[2*(ideg-1)];

  int uuu1 = pow(ifactorial_func(ideg),2);
  int lll1 = ifactorial_func(2*(ideg-1));

  str << "template<> void d_ccoef_genbez2<"<<ideg<<">"
      << "(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,\n" 
      << "                                   int ielem, int icoor,\n"
      << "                                   dblAr2& __restrict__ d_ccoef){\n";
  str << "  METRIS_ASSERT(icoor == 0 || icoor == 1);\n";
  str << "  const int sg = 1 - 2*(icoor%2);\n";
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
             
        char up_s[16] = ""; if(up != 1) snprintf(up_s,5,"%3d*",up);
        char lo_s[16] = ""; if(lo != 1) snprintf(lo_s,5,"/%3d",lo);

        char irnk1_1_s[16]; snprintf(irnk1_1_s,5,"%4d",irnk1_1);
        char irnk1_2_s[16]; snprintf(irnk1_2_s,5,"%4d",irnk1_2);
        char irnk2_1_s[16]; snprintf(irnk2_1_s,5,"%4d",irnk2_1);
        char irnk2_2_s[16]; snprintf(irnk2_2_s,5,"%4d",irnk2_2);
        //if(ifirst == 1){
        //  str << "\n  ccoef["<<irnkcc_s<<"] = "<<up_s<<"det2_vdif(coord[fac2poi(ielem,"<<irnk1_1_s<<")],coord[fac2poi(ielem,"<<irnk1_2_s<<")]\n";
        //  ifirst = 0;
        //}else{
        //  str << "\n             + "<<up_s<<"det2_vdif(coord[fac2poi(ielem,"<<irnk1_1_s<<")],coord[fac2poi(ielem,"<<irnk1_2_s<<")]\n";
        //}
        //str << "                            ,coord[fac2poi(ielem,"<<irnk2_1_s<<")],coord[fac2poi(ielem,"<<irnk2_2_s<<")])"<<lo_s;


        std::pair<int, int> key1(irnkcc, irnk1_1);
        std::pair<int, int> key2(irnkcc, irnk1_2);
        std::pair<int, int> key3(irnkcc, irnk2_1);
        std::pair<int, int> key4(irnkcc, irnk2_2);
        
        std::string up_ss = std::string(up_s);
        std::string lo_ss = std::string(lo_s);


        //d_ccoef_map_x[key1] += " + " + up_ss + "(coord(fac2poi(ielem," + std::string(irnk2_1_s) + "),1-icoor) - coord(fac2poi(ielem," + std::string(irnk2_2_s) + "),1-icoor))" + lo_ss;
        //d_ccoef_map_x[key2] += " + " + up_ss + "(coord(fac2poi(ielem," + std::string(irnk2_2_s) + "),1-icoor) - coord(fac2poi(ielem," + std::string(irnk2_1_s) + "),1-icoor))" + lo_ss;
        //d_ccoef_map_x[key3] += " + " + up_ss + "(coord(fac2poi(ielem," + std::string(irnk1_2_s) + "),1-icoor) - coord(fac2poi(ielem," + std::string(irnk1_1_s) + "),1-icoor))" + lo_ss;
        //d_ccoef_map_x[key4] += " + " + up_ss + "(coord(fac2poi(ielem," + std::string(irnk1_1_s) + "),1-icoor) - coord(fac2poi(ielem," + std::string(irnk1_2_s) + "),1-icoor))" + lo_ss;

        insert_map2(d_ccoef_map,key1, up, lo, irnk2_1, irnk2_2);
        insert_map2(d_ccoef_map,key2, up, lo, irnk2_2, irnk2_1);
        insert_map2(d_ccoef_map,key3, up, lo, irnk1_2, irnk1_1);
        insert_map2(d_ccoef_map,key4, up, lo, irnk1_1, irnk1_2);

        //d_ccoef_map_x[key1] += " + vdiff_perp_x(coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]])";
        //d_ccoef_map_x[key2] += " + vdiff_perp_x(coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]])";
        //d_ccoef_map_x[key3] += " + vdiff_perp_x(coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]])";
        //d_ccoef_map_x[key4] += " + vdiff_perp_x(coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]])";

        //d_ccoef_map_y[key1] += " + vdiff_perp_y(coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]])";
        //d_ccoef_map_y[key2] += " + vdiff_perp_y(coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]])";
        //d_ccoef_map_y[key3] += " + vdiff_perp_y(coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]])";
        //d_ccoef_map_y[key4] += " + vdiff_perp_y(coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]])";
      }
    }
    //str << ";\n";
  }
  //str << "\n";

  for (const auto& entry : d_ccoef_map) {
    const auto& key = entry.first;
    std::string value = "";
    for(const auto &term : d_ccoef_map[key]){
      //printf("debug key %d %d term %d %d %d %d\n",key.first,key.second,term[0],term[1],term[2],term[3]);
      if(term[2] < 0) continue;
      int up = term[0];
      int lo = term[1];
      std::string up_s = "";
      std::string lo_s = "";
      if(up != 1) up_s = std::to_string(up)+"*";
      if(lo != 1) lo_s = "/"+std::to_string(lo);
      int irnk1 = term[2];
      int irnk2 = term[3];
      value += " + sg*" + up_s + "(coord(fac2poi(ielem," + std::to_string(irnk1) + "),1-icoor) - coord(fac2poi(ielem," + std::to_string(irnk2) + "),1-icoor))" + lo_s;
    }
    //std::cout<<"value = "<<value<<"\n";
    //fflush(stdout);
    value = value.substr(3); // removes front " + "
    //str << "\n  d_ccoef(" << key.first << "," << key.second << ") = " << factor_s << "*(" << value << ");\n";
    str << "\n  d_ccoef(" << key.first << "," << key.second << ") = " << value << ";\n";
  }

  //// NEED TO BE CAREFUL OF THE CONSTANT FACTORS IN THE END
  ////str<<"  if(icoor == 0){\n";
  //for (const auto& entry : d_ccoef_map_x) {
  //  const auto& key = entry.first;
  //  //int factor;
  //  //factor = (key.first == 0 || key.first == 1 || key.first == 2) ? 4 : 2;
  //  //char factor_s[16];  snprintf(factor_s,4,"%3d",factor);
  //  const std::string& value = entry.second.substr(3); // removes front " + "
  //  //str << "\n  d_ccoef(" << key.first << "," << key.second << ") = " << factor_s << "*(" << value << ");\n";
  //  str << "\n  d_ccoef(" << key.first << "," << key.second << ") = " << value << ";\n";
  //}
  //str<<"\n  }else{\n";
  //for (const auto& entry : d_ccoef_map_y) {
  //  const auto& key = entry.first;
  //  int factor;
  //  factor = (key.first == 0 || key.first == 1 || key.first == 2) ? 4 : 2;
  //  char factor_s[16];  snprintf(factor_s,4,"%3d",factor);
  //  const std::string& value = entry.second.substr(3); // removes front " + "
  //  str << "\n    d_ccoef[" << key.first << "][" << key.second << "] = " << factor_s << "*(" << value << ");\n";
  //}
  //str<<"  }\n";

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
  if(ideg != 2){
    printf("##Degrees > 2 not tested yet: write unit test for getccoef2_d\n");
  }

  std::ostringstream str;
  str << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  str << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  str << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  str << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";

  str << "#include \"codegen_ccoef_d.hxx\"\n\n";
  str << "#include \"types.hxx\"\n\n";
  str << "namespace Metris{\n\n";

  str << "template<int ideg>\n"
      << "void d_pt_ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,\n" 
      << "                        int ielem, int inode,\n"
      << "                        dblAr2&__restrict__ d_ccoef){}\n\n";

  str << "void vdiff_perp(const double* a, const double* b, int up, int lo, double *res);\n\n";
  str << "void vdiff_perp_sum(const double* a, const double* b, int up, int lo, double *res);\n\n";
 

  std::map<std::pair<int, int>, std::string> d_ccoef_map_old; 
  std::map<std::pair<int, int>, std::vector<SANS::DLA::VectorS<4, int>>> d_ccoef_map;

  int nnodj = facnpps[2*(ideg-1)];
  int nnode = facnpps[ideg];

  int uuu1 = pow(ifactorial_func(ideg),2);
  int lll1 = ifactorial_func(2*(ideg-1));

  str << "template<> void d_pt_ccoef_genbez2<"<<ideg<<">"
      << "(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,\n" 
      << "                                        int ielem, int inode,\n"
      << "                                        dblAr2&__restrict__ d_ccoef){\n";

  for(int irnkcc = 0; irnkcc < nnodj; irnkcc++){
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
        //if(ifirst == 1){
        //  str << "\n  ccoef["<<irnkcc_s<<"] = "<<up_s<<"*det2_vdif(coord[fac2poi[ielem]["<<irnk1_1_s<<"]],coord[fac2poi[ielem]["<<irnk1_2_s<<"]]\n";
        //  ifirst = 0;
        //}else{
        //  str << "\n             + "<<up_s<<"*det2_vdif(coord[fac2poi[ielem]["<<irnk1_1_s<<"]],coord[fac2poi[ielem]["<<irnk1_2_s<<"]]\n";
        //}
        //str << "                            ,coord[fac2poi[ielem]["<<irnk2_1_s<<"]],coord[fac2poi[ielem]["<<irnk2_2_s<<"]])/"<<lo_s;


        std::pair<int, int> key1(irnkcc, irnk1_1);
        std::pair<int, int> key2(irnkcc, irnk1_2);
        std::pair<int, int> key3(irnkcc, irnk2_1);
        std::pair<int, int> key4(irnkcc, irnk2_2);

        insert_map2(d_ccoef_map,key1,up,lo, irnk2_1, irnk2_2);
        insert_map2(d_ccoef_map,key2,up,lo, irnk2_2, irnk2_1);
        insert_map2(d_ccoef_map,key3,up,lo, irnk1_2, irnk1_1);
        insert_map2(d_ccoef_map,key4,up,lo, irnk1_1, irnk1_2);
        
        //d_ccoef_map_old[key1] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]])[icoor] / " + std::string(lo_s);
        //d_ccoef_map_old[key2] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk2_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk2_1_s) + "]])[icoor] / " + std::string(lo_s);
        //d_ccoef_map_old[key3] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]])[icoor] / " + std::string(lo_s);
        //d_ccoef_map_old[key4] += " + " + std::string(up_s) + "*vdiff_perp(coord[fac2poi[ielem][" + std::string(irnk1_1_s) + "]], coord[fac2poi[ielem][" + std::string(irnk1_2_s) + "]])[icoor] / " + std::string(lo_s);

      }
    }
    //str << ";\n";
  }
  str << "\n";
  str << "  d_ccoef.fill("<<std::to_string(nnodj)<<",2,0);\n\n";
  //str << "    for(int ii = 0; ii < " << npp_c << "; ii++) {\n";
  //str << "      d_ccoef[ii] = 0;\n";
  //str << "    }\n";

  //for (int inode = 0; inode < 6; inode++){ 
  //  char inode_s[16];  snprintf(inode_s,4,"%3d",inode);
  //  if (inode==0){
  //    str<<"    if(inode  == "<<inode_s<<"){\n";
  //  }else{
  //    str<<"    else if(inode  == "<<inode_s<<"){\n";
  //  }
  //  for (const auto& entry : d_ccoef_map_old) {
  //      const auto& key = entry.first;
  //      if (key.first==inode)
  //      {
  //          const std::string& value = entry.second.substr(3); // removes front " + "
  //          str << "\n      d_ccoef[" << key.second << "] = " << value << ";";
  //      }
  //  }
  //  str<<"\n    }\n";
  //}

  //str<<"\n\n\n";

  std::ostringstream str_nodes[nnode];
  std::string funnames[2] = {"    vdiff_perp", "vdiff_perp_sum"};

  for (const auto& entry : d_ccoef_map) {
    const auto& key = entry.first;
    int icoef = key.first;
    int inode = key.second; 
    int ifirst = 0;
    std::string termstr = "";
    for(const auto &term : d_ccoef_map[key]){
      if(term[2] < 0) continue;

      std::string funname = funnames[ifirst];
      ifirst = 1;
      //printf("debug key %d %d term %d %d %d %d\n",key.first,key.second,term[0],term[1],term[2],term[3]);

      int up = term[0];
      int lo = term[1];
      int irnk1 = term[2];
      int irnk2 = term[3];
      termstr += funname + "(coord[fac2poi(ielem," + std::to_string(irnk1) + ")], coord[fac2poi(ielem," + std::to_string(irnk2) + ")],"
                + std::to_string(up) + "," + std::to_string(lo) + ",d_ccoef[" + std::to_string(icoef) + "]);\n";
    }
    str_nodes[inode] << termstr;
  }

  for(int inode = 0; inode < nnode; inode++){
    char inode_s[16];  snprintf(inode_s,4,"%3d",inode);
    //std::cout<<"inode "<<inode<<" str = "<<str_nodes[inode].str()<<"\n";
    if (inode==0){
      str<<"  if(inode == "<<inode_s<<"){\n";
    }else{
      str<<"  else if(inode == "<<inode_s<<"){\n";
    }
    str<<str_nodes[inode].str();
    str<<"  }\n";
  }

  
  str << "\n} \n \n";
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
  if(ideg != 2){
    printf("## Degrees > 2 not tested yet: write unit test for getccoef2_d\n");
  }

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

    //str << "double* vdiff(const double* a,const double* b);\n\n";
    //str << "double* vproduct(const double* a, const double* b);\n\n";

    str << "static double* vdiff(const double* a, const double* b){";
    str << "  METRIS_THROW_MSG(TODOExcept(),\"Reimplement ccoef3_d\");";
    str << "}";
    str << "static double* vproduct(const double* a, const double* b){";
    str << "  METRIS_THROW_MSG(TODOExcept(),\"Reimplement ccoef3_d\");";
    str << "}";

    str << "template<int ideg>\n";
    str << "void d_ccoef_genbez3(const intAr2 & __restrict__ tet2poi,\n";
    str << "                     const dblAr2& __restrict__ coord,\n";
    str << "                     int ielem,\n";
    str << "                     int icoor,\n";
    str << "                     dblAr2& __restrict__ d_ccoef){}\n\n";


    std::map<std::pair<int, int>, std::string> d_ccoef_map; // Declare d_ccoef_map_x

    int uuu1 = pow(ifactorial_func(ideg),3);
    int lll1 = ifactorial_func(3*(ideg-1));

    str << "template<>\n";
    str << "void d_ccoef_genbez3<"<<ideg<<">(const intAr2 & __restrict__ tet2poi,\n";
    str << "                     const dblAr2& __restrict__ coord,\n";
    str << "                     int ielem,\n";
    str << "                     int icoor,\n";
    str << "                     dblAr2& __restrict__ d_ccoef){\n\n";

    
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
                  //if(ifirst == 1){
                  //  str << "\n  ccoef["<<irnkcc_s<<"] = "<<up_s<<"*det3_vdif(coord[tet2poi[ielem]["<<irnk1_1_s<<"]],coord[tet2poi[ielem]["<<irnk1_2_s<<"]]\n";
                  //  ifirst = 0;
                  //}else{
                  //  str << "\n             + "<<up_s<<"*det3_vdif(coord[tet2poi[ielem]["<<irnk1_1_s<<"]],coord[tet2poi[ielem]["<<irnk1_2_s<<"]]\n";
                  //}
                  //str << "                            ,coord[tet2poi[ielem]["<<irnk2_1_s<<"]],coord[tet2poi[ielem]["<<irnk2_2_s<<"]]\n";
                  //str << "                            ,coord[tet2poi[ielem]["<<irnk3_1_s<<"]],coord[tet2poi[ielem]["<<irnk3_2_s<<"]])/"<<lo_s;

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
      //str << ";\n";
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

    //str << "template ";
    //str << "void d_ccoef_genbez3<2,double>(const intAr2 & __restrict__ tet2poi,\n";
    //str << "                     const dblAr2& __restrict__ coord,\n";
    //str << "                     int ielem,\n";
    //str << "                     int icoor,\n";
    //str << "                     double* __restrict__ ccoef,\n";
    //str << "                     dblAr2& __restrict__ d_ccoef);\n\n";

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

