//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include <bunit/common_setup.hxx>

#include <boost/timer/progress_display.hpp>

#include "../src/ho_constants.hxx"
#include "../src/io_libmeshb.hxx"
#include "../src/utils/aux_misc.hxx"
#include "../src/quality/low_metqua.hxx"
#include "../SANS/Surreal/SurrealS.h"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
#include "../src/low_geo/misc.hxx"
#include "../src/Adaptation/low_cavqual.cxx"
#include <boost/hana.hpp> 

using namespace Metris;

typedef MetricFieldAnalytical MFT;


BOOST_AUTO_TEST_CASE(test_eval3) 
{

  // bool is whether straight
  std::vector<std::pair<std::string,double>> meshes = {
     {METRIS_CASES_DIR "/unit/3D/cube/iso.p1.16k", 1.0}
    ,{METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k.meshb",0.21}
    }; 


  CavOprOpt opts;
  CavOprInfo info;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = false;
  opts.dryrun = false;
  //opts.allow_remove_points = icollapse; 
  opts.allow_remove_points = true; 
  opts.allow_remove_points_superdim = true; // For boundary
  opts.qmax_nec = -1;
  opts.qmax_suf = -1;
  opts.qmax_iff = -1;

  const int ithrd1 = 0;
  const int ithrd2 = 1;

  for(auto cases : meshes)
  { 
    std::string mesh = cases.first;
    double sclmet = cases.second;
    cargHandler arg("-in " + mesh + "  -anamet 1 -sclmet " + std::to_string(sclmet) + " -verb 2 -vdepth 1 -prefix tmp/");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g); 
    run.degElevate();


    std::cout<<"\n\n------------------------------------------------\n";
    std::cout<<"Mesh "<<mesh<<"\n";
    const int tdim = msh.get_tdim();

    MeshStat stat_ini;
    run.statMesh(tdim,&stat_ini);

    dblAr2 coor0(msh.npoin,msh.idim);
    msh.coord.copyTo(coor0);

    // Apply regular smoothing, then get stats
    msh.param->opt_niter = 20;
    msh.param->iverb = 0;

    double t0_1 = get_wall_time();
    run.optimMesh();
    double t1_1 = get_wall_time();

    std::cout<<"\n\n------------------------------------------------\n";
    printf("---- Post smoothing stat\n");
    MeshStat stat_qua;
    run.statMesh(tdim,&stat_qua);

    // Reset, then apply length based smoothing
    coor0.copyTo(msh.coord);
    msh.param->iverb = 0;
    MshCavity cav(100,100,1);
    int iopen;
    double t0_2 = get_wall_time();
    for(int niter = 0; niter < 10; niter++){
      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        if(msh.poi2bpo[ipoin] >= 0) continue;
        if(msh.poi2ent(ipoin,0) < 0) continue;

        cav.reset();
        int ierro = ball(msh,ipoin,cav.lcedg,cav.lcfac,cav.lctet,&iopen, false, ithrd1);
        METRIS_ENFORCE(ierro == 0);

        cav.ipins = ipoin;
        int iseed = msh.poi2ent(ipoin,0);
        ierro = movePointCavLen(msh,cav,tdim,iseed,3,ithrd1);
        METRIS_ENFORCE(ierro == 0);
      }
    }
    double t1_2 = get_wall_time();

    std::cout<<"\n\n------------------------------------------------\n";
    printf("---- Post len based stat\n");
    MeshStat stat_len;
    run.statMesh(tdim,&stat_len);


    printf("\n\n Final time smoothing %f len based %f\n",
           t1_1 - t0_1, t1_2 - t0_2);

    stat_ini.print("ini");
    stat_qua.print("qua");
    stat_len.print("len");

  }// fortest cases

}// end boost test case











