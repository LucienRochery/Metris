//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include "common_setup.hxx"

#include <random>
#include "../src/adapt/low_delaunay.hxx"
#include "../src/utils/mprintf.hxx"
#include "../src/adapt/msh_insert.hxx"

#include <boost/hana.hpp> 
namespace hana = boost::hana;
using namespace hana::literals;

namespace utf = boost::unit_test;

using namespace Metris;

typedef MetricFieldAnalytical MFT;

// Test indelsphere works in surface case 
BOOST_AUTO_TEST_CASE(test_indelsphere) 
{

  // bool is whether straight
  std::vector<std::string> meshes = {
    METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
  }; 


  const int nsamp = 100;
  dblAr2 bary(nsamp,3);
  dblAr2 bary_out(nsamp,3);
  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(0);
  for(int ii = 0; ii < nsamp; ii++){
    double sum = 0;
    do{
      for(int jj = 0; jj < 3; jj++){
        bary(ii,jj) = unif(rng);
        sum += bary(ii,jj);
      }
    }while(abs(sum) < 1.0e-16);

    for(int jj = 0; jj < 3; jj++){
      bary(ii,jj) /= sum;
    }
  }
  std::uniform_real_distribution<double> unif_out(-100.0,100.0);
  std::default_random_engine rng_out(0);
  for(int ii = 0; ii < nsamp; ii++){
    double sum = 0;
    do{
      for(int jj = 0; jj < 3; jj++){
       bary_out(ii,jj) = unif_out(rng_out);
       sum += bary_out(ii,jj);
      }
    }while(abs(sum) < 1.0e-16);

    for(int jj = 0; jj < 3; jj++){
      bary_out(ii,jj) /= sum;
    }
  }



  for(std::string s : meshes)
  { 
    cargHandler arg("-in " + s + " -anamet 1 -sclmet 0.1 -vdepth 0 -verb 0 -opt-niter 0 -adp-opt-niter 0 -adapt 1");
    MetrisRunner run2D(arg.c, arg.v);
    Mesh<MFT> &msh2D = *((Mesh<MFT>*) run2D.msh_g);


    std::cout<<"\n\n------------------------------------------------\n";
    std::cout<<"Mesh "<<s<<"\n";
    std::cout<<"------------------------------------------------\n";


    // We can add a dimension to coord, it shouldn't change anything as 
    // long as all we don't flatten (IO only)
    // indelsphere is also not looking at msh.idim. We can fool it easily
    msh2D.coord.set_stride(3);
    for(int ipoin = 0; ipoin < msh2D.npoin; ipoin++){
      msh2D.coord(ipoin,2) = 0;
    }

    msh2D.idim = 3;
    msh2D.set_nelem(0);
    writeMesh(s+".3D",msh2D);


    cargHandler arg3D("-in " + s + ".3D -anamet 1 -sclmet 0.1 -vdepth 0 -verb 0 -opt-niter 0 -adp-opt-niter 0 -adapt 1");
    MetrisRunner run3D(arg3D.c, arg3D.v);
    Mesh<MFT> &msh3D = *((Mesh<MFT>*) run3D.msh_g);


    std::cout<<"\n\n------------------------------------------------\n";
    std::cout<<"Setup done\n";
    std::cout<<"------------------------------------------------\n";


    INCVDEPTH(msh2D.param);

    double hmet = 1;
    for(int iscal = 0; iscal < 10; iscal++){
      hmet /= 2;
      double metl[6] = {hmet, 0, hmet, 0, 0, hmet};

      printf("-- Met scaling %f \n",hmet);

      for(int itest = 0; itest <= 1; itest++){
        dblAr2& bary1 = itest == 0 ? bary : bary_out;
        printf(" - TEST %d: sample in and locate ",itest);
        if(itest == 0) printf("in      element\n");
        else           printf("outside element\n");
        for(int iface = 0; iface < msh2D.nface; iface++){
          if(isdeadent(iface,msh2D.fac2poi)) continue;

          for(int isamp = 0; isamp < nsamp; isamp++){
            double coop[3];
            eval2<2,1>(msh2D.coord,msh2D.fac2poi[iface],FEBasis::Lagrange,
                       DifVar::None, DifVar::None, bary1[isamp], coop, NULL, NULL);
            coop[2] = 0;

            bool indel2 = indelsphere<2,2>(msh2D,coop,metl,msh2D.fac2poi[iface]);
            bool indel3 = indelsphere<3,2>(msh2D,coop,metl,msh2D.fac2poi[iface]);
            bool indel4 = indelsphere<3,2>(msh3D,coop,metl,msh3D.fac2poi[iface]);

            BOOST_TEST(indel2 == indel3);
            BOOST_TEST(indel2 == indel4);
          }// for isamp
        }//for iface
      }// for itest

    }// for iscal


    #if 0
    std::cout<<"\n\n------------------------------------------------\n";
    std::cout<<"Test II: adapt in 2D and in 3D + z, check same results.\n";
    std::cout<<"------------------------------------------------\n";

    msh2D.coord.set_stride(2);
    msh2D.idim = 2;

    // Next adapt both 2D and 3D mesh and test 
    MeshStat stat2D_ini, stat2D, stat3D;
    int ninser;

    run2D.statMesh(&stat2D_ini);

    run2D.msh_g->param->iverb = 2;
    run2D.msh_g->param->ivdepth = 2;
    run3D.msh_g->param->iverb = 2;
    run3D.msh_g->param->ivdepth = 2;

    // Don't use adaptMesh because there are other dimension dependent things.
    // Insertion is also the only one that use Delaunay.
    //printf("-- Adapt 2D\n");
    insertLongEdges<MFT,2,1>(msh2D, &ninser, 0, 1, 2);
    //run2D.adaptMesh();

    //printf("-- Adapt 3D\n");
    insertLongEdges<MFT,3,1>(msh3D, &ninser, 0, 1, 2);
    //run3D.adaptMesh();


    run2D.statMesh(&stat2D);
    run3D.statMesh(&stat3D);

    stat2D_ini.print("2D initial");
    stat2D.print("2D");
    stat3D.print("3D");

    BOOST_TEST(abs(stat2D.pctunit - stat3D.pctunit) < 1.0e-2);
    #endif


  }// for s

}// end boost test case











