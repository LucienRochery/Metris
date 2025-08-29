//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include "common_setup.hxx"

#include "../src/utils/tuple_hashtable.hxx"
#include "../src/utils/aux_misc.hxx"


#include <random>
#include <cmath>

namespace utf = boost::unit_test;

using namespace Metris;
typedef MetricFieldAnalytical MFT;

BOOST_AUTO_TEST_CASE(interperr) 
{


  std::vector<std::string> meshes = {"../cases/1200_p1.meshb"};

  for(auto mesh : meshes){

    cargHandler arg("-in " + mesh + " -anamet 1");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);

    printf("-- Start test mesh = %s\n",mesh.c_str());
    fflush(stdout);

    std::unordered_map<std::tuple<int, int>, int, tup2_hash::hash> edgHshTab_std;
    TupleHashTable<2,int> edgHshTab_mtr;

    const int tdim = msh.get_tdim();
    const int nentt = msh.nentt(tdim);
    const int nedgl = (tdim*(tdim+1))/2;
    const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                               tdim == 2 ? lnoed2[0] : lnoed3[0]);
    const intAr2 &ent2poi = msh.ent2poi(tdim);

    edgHshTab_std.reserve(2*nentt);
    edgHshTab_mtr.reserve(2*nentt);

    int nhash = edgHshTab_mtr.get_nhash();
    printf("nhash = %d \n",nhash);
    fflush(stdout);

    for(int ientt = 0; ientt < nentt; ientt++){
      for(int ied = 0; ied < nedgl; ied++){
        int i1 = ent2poi(ientt, lnoed(ied,0));
        int i2 = ent2poi(ientt, lnoed(ied,1));
        auto key = stup2(i1, i2);
        auto tt = edgHshTab_std.find(key);
        if(tt != edgHshTab_std.end()) continue;
        edgHshTab_std.insert({key, ientt});
      }
    }
    double t0_std = get_cpu_time();
    int ninser_std = 0;
    for(int ientt = 0; ientt < nentt; ientt++){
      for(int ied = 0; ied < nedgl; ied++){
        int i1 = ent2poi(ientt, lnoed(ied,0));
        int i2 = ent2poi(ientt, lnoed(ied,1));
        auto key = stup2(i1, i2);
        auto tt = edgHshTab_std.find(key);
        if(tt != edgHshTab_std.end()) continue;
        edgHshTab_std.insert({key, ientt});
        ninser_std++;
      }
    }
    double t1_std = get_cpu_time();


    for(int ientt = 0; ientt < nentt; ientt++){
      for(int ied = 0; ied < nedgl; ied++){
        int i1 = ent2poi(ientt, lnoed(ied,0));
        int i2 = ent2poi(ientt, lnoed(ied,1));
        uint32_t key[2] = {(uint32_t) i1, (uint32_t) i2};
        sortupto8_dec(key,2);
        int ifnd = edgHshTab_mtr.find(key);
        if(ifnd >= 0) continue;
        edgHshTab_mtr.insert(key, ientt);
      }
    }
    double t0_mtr = get_cpu_time();
    int ninser_mtr = 0;
    for(int ientt = 0; ientt < nentt; ientt++){
      for(int ied = 0; ied < nedgl; ied++){
        int i1 = ent2poi(ientt, lnoed(ied,0));
        int i2 = ent2poi(ientt, lnoed(ied,1));
        uint32_t key[2] = {(uint32_t) i1, (uint32_t) i2};
        sortupto8_dec(key,2);
        int ifnd = edgHshTab_mtr.find(key);
        if(ifnd >= 0) continue;
        edgHshTab_mtr.insert(key, ientt);
        ninser_mtr++;
      }
    }
    double t1_mtr = get_cpu_time();

    printf("Dummy %d %d \n",ninser_std, ninser_mtr);

    printf("Time std = %e Metris = %e factor %ex\n", t1_std - t0_std, t1_mtr - t0_mtr,
      (t1_mtr - t0_mtr) / (t1_std - t0_std));

    int ncol_min, ncol_max, nempty;
    double ncol_avg;
    edgHshTab_mtr.stat(&ncol_min, &ncol_max, &ncol_avg, &nempty);
    nhash = edgHshTab_mtr.get_nhash();
    printf("Metris hash table ncol min %d max %d avg %e nempty %d/%d = %e%%",
          ncol_min, ncol_max, ncol_avg, nempty,nhash,
          nempty/(double)nhash *100);

  }


}