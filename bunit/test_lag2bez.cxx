//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include "common_setup.hxx"

#include "../src/msh_lag2bez.hxx"
#include "../src/mprintf.hxx"
#include "../src/utils/CT_loop.hxx"
#include "../src/aux_timer.hxx"

#include <cmath>

namespace utf = boost::unit_test;

using namespace Metris;
typedef MetricFieldAnalytical MFT;

// 1.0e-8 relative error -> 1.0e-6% utf::tolerance()
BOOST_AUTO_TEST_CASE(lag2bez)  // , * utf::tolerance(double(1.0e-6)) 
{//METRIS_MAX_DEG


  //#if METRIS_MAX_DEG >= 4
  //  std::vector<std::string> meshes = {"../cases/1200_p2.meshb","../cases/1200_p3.meshb","../cases/1200_p4.meshb",
  //                                     "../cases/curved_p2.meshb","../cases/curved_p3.meshb","../cases/curved_p4.meshb"};
  //#elif METRIS_MAX_DEG >= 3
  //  std::vector<std::string> meshes = {"../cases/1200_p2.meshb","../cases/1200_p3.meshb",
  //                                     "../cases/curved_p2.meshb","../cases/curved_p3.meshb"};
  //#else
  //  std::vector<std::string> meshes = {"../cases/1200_p2.meshb",
  //                                     "../cases/curved_p2.meshb"};
  //#endif
  std::vector<std::string> meshes = {"../cases/1200_p2.meshb",
                                     "../cases/curved_p2.meshb"};
  //std::vector<std::string> meshes = {"cases/1200_p2.meshb","cases/1200_p3.meshb"};

  double tol = 1.0e-14;
  int nconv = 200;
  double dummy;
  for(auto s : meshes)
  {


    for(int ideg = 2; ideg <= METRIS_MAX_DEG;  ideg++){
      printf("-- Reading mesh %s and elev to degree %d \n",s.c_str(),ideg);

      cargHandler arg("-in " + s + " -anamet 1 -vdepth 0 -verb 0 -t" + std::to_string(ideg));
      MetrisRunner run(arg.c, arg.v);
      Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);


      dblAr2 coor0(msh.npoin, msh.idim);

      msh.coord.copyTo(coor0);

      double t0 = get_wall_time();
      for(int ii = 0; ii < nconv; ii++){
        if(msh.getBasis() == FEBasis::Lagrange){
          msh.setBasis(FEBasis::Bezier);
          dummy += msh.coord(0,0);
          msh.setBasis(FEBasis::Lagrange);
          dummy += msh.coord(0,0);
        }else{
          msh.setBasis(FEBasis::Lagrange);
          dummy += msh.coord(0,0);
          msh.setBasis(FEBasis::Bezier);
          dummy += msh.coord(0,0);
        }
      }
      double t1 = get_wall_time();

      int npcomp = 0;
      double erri = -1.0e30;
      double err2 = 0.0;
      for(int ipoin=0; ipoin<msh.npoin;ipoin++){
        if(msh.poi2ent(ipoin,0) < 0) continue;
        npcomp++;
        err2 += geterrl2<3>(coor0[ipoin],msh.coord[ipoin]);
        for(int i=0; i < 3;i++){
          double err0 = abs(msh.coord(ipoin,i) - coor0(ipoin,i));
          erri  = erri > err0 ? erri : err0;
          #ifndef NDEBUG
          if(erri >= tol){
            printf("## LARGE ERROR FOR POINT %d \n",ipoin);
          }
          #endif
        }
      }
      err2 /= npcomp*3;
      err2 = sqrt(err2);
      #ifndef NDEBUG
        BOOST_TEST(erri < tol);
        BOOST_TEST(err2 < tol);
      #endif
      printf(" %d double conversion abs coord error: inf = %23.16e, l2 = %23.16e\n",
              nconv,erri,err2);
      printf(" Time = %f s = %f / conv = %e /conv.elt \n", t1-t0, (t1-t0)/nconv, (t1-t0)/((double)nconv * msh.nentt(msh.get_tdim())) );
    }
  }// for s meshes

  printf("Dummy %f \n",dummy);

}