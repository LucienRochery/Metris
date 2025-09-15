//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_libkdtree 

#include "common_setup.hxx"
#include "../libs/libkdtree/kdtree++/kdtree.hpp"

#include <random>
#include <vector>

using namespace Metris;
typedef MetricFieldAnalytical MFT;

template<int gdim>
struct kdtreeNode_iso
{
  typedef double value_type;

  double operator[](size_t n) const{
    METRIS_ASSERT(n < gdim);
    return xyz[n];
  }

  double distance(const kdtreeNode_iso &node) const {
    return geterrl2<gdim>(xyz,node.xyz);
    // this is not correct   return sqrt( x*x+y*y+z*z);
    
    // this is what kdtree checks with find_within_range()
    // the "manhattan distance" from the search point.
    // effectively, distance is the maximum distance in any one dimension.
    //   return max(fabs(x),max(fabs(y),fabs(z)));

  }

  double xyz[gdim];
  int ipoin;
  size_t index;
};


template<int gdim>
struct kdtreeNode_aniso
{
  typedef double value_type;
  static constexpr int nnmet = (gdim*(gdim+1))/2;

  double operator[](size_t n) const{
    METRIS_ASSERT(n < gdim);
    return xyz[n];
  }

  double distance(const kdtreeNode_aniso &node) const {
    double metl[nnmet];
    for(int ii = 0; ii < nnmet; ii++) metl[ii] = 0.5*(lmet[ii] + node.lmet[ii]);
    getexpmet_inp<gdim,double>(metl);

    return getlenedgsq<gdim>(xyz,node.xyz,metl);
  }

  double xyz[gdim];
  double lmet[nnmet];
  int ipoin;
  size_t index;
};

BOOST_AUTO_TEST_CASE(test_libkdtree) 
{

  const int nsamp = 1e3;
  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(0);

  // anamet 2 is non-uniform with some curvature in both dims (though not similar)
  std::vector<std::string> meshes = {
     METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k -anamet 2"
    ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k -anamet 2"
  };
  double dum = 0;

  for(std::string mesh : meshes){

    cargHandler arg("-in " + mesh + " -verb 0");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
    msh.cleanup();

    CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
      constexpr int nnmet = (gdim*(gdim+1))/2;
      typedef KDTree::KDTree<gdim,kdtreeNode_iso<gdim>> t_KDTree_iso;
      typedef KDTree::KDTree<gdim,kdtreeNode_aniso<gdim>> t_KDTree_aniso;

      // Initialize kdtree with all points in the mesh
      double t0i_iso = get_cpu_time();
      t_KDTree_iso kdtree_iso;
      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        kdtreeNode_iso<gdim> node;
        node.ipoin = ipoin;
        for(int ii = 0; ii < gdim; ii++) node.xyz[ii] = msh.coord(ipoin,ii);
        kdtree_iso.insert(node);
      }
      double t1i_iso = get_cpu_time();
      kdtree_iso.optimise();
      double t2i_iso = get_cpu_time();
      int nkpsi1_iso = (int)(msh.npoin / (t1i_iso - t0i_iso) / 1000);
      int nkpsi2_iso = (int)(msh.npoin / (t2i_iso - t1i_iso) / 1000);
      printf("-- %dD   iso init times fill = %f = %dk/s optim = %f = %dk/s\n",
             gdim,t1i_iso - t0i_iso,nkpsi1_iso,t2i_iso - t1i_iso,nkpsi2_iso);

      // Generate samples once and for all.
      dblAr2 coord_samp(nsamp, gdim);
      dblAr2   met_samp(nsamp, nnmet);
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < gdim; ii++){
          double tt = unif(rng);
          double crd = msh.bb[ii][0]*tt + msh.bb[ii][1]*(1 - tt);
          coord_samp(isamp,ii) = crd;
        }
        msh.met.getMetFullinfo(AsDeg::P1,DifVar::None,MetSpace::Log,NULL,0,NULL,
                               coord_samp[isamp],met_samp[isamp],NULL);
      }

      // Now loop over points, do quadratic search for closest as reference
      for(int isamp = 0; isamp < nsamp; isamp++){
        double dminq = 1.0e30;
        int iminq = -1;
        for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
          double dist = geterrl2<gdim>(coord_samp[isamp],msh.coord[ipoin]);
          if(dist > dminq) continue;
          dminq = dist;
          iminq = ipoin;
        }

        kdtreeNode_iso<gdim> node;
        for(int ii = 0; ii < gdim; ii++) node.xyz[ii] = coord_samp(isamp,ii);

        std::pair<typename t_KDTree_iso::const_iterator,double> found 
          = kdtree_iso.find_nearest(node);
        BOOST_REQUIRE(found.first != kdtree_iso.end());
        BOOST_TEST(found.first->ipoin == iminq);

        if(found.first->ipoin != iminq){
          printf("kdtree found %d we found %d\n",found.first->ipoin,iminq);
          double dist1 = geterrl2<gdim>(coord_samp[isamp],coord_samp[iminq]);
          double dist2 = geterrl2<gdim>(coord_samp[isamp],found.first->xyz);
          printf("dist to our node %e to kdtree %e\n",dist1,dist2);
        }

      }



      printf("Debug first test \n");




      // Now make it aniso
      msh.met.setSpace(MetSpace::Log);
      double t0i_aniso = get_cpu_time();
      t_KDTree_aniso kdtree_aniso;
      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        kdtreeNode_aniso<gdim> node;
        node.ipoin = ipoin;
        for(int ii = 0; ii < gdim; ii++)  node.xyz[ii] = msh.coord(ipoin,ii);
        for(int ii = 0; ii < nnmet; ii++) node.lmet[ii] = msh.met(ipoin,ii);
        kdtree_aniso.insert(node);
      }
      double t1i_aniso = get_cpu_time();
      kdtree_aniso.optimise();
      double t2i_aniso = get_cpu_time();
      int nkpsi1_aniso = (int)(msh.npoin / (t1i_aniso - t0i_aniso) / 1000);
      int nkpsi2_aniso = (int)(msh.npoin / (t2i_aniso - t1i_aniso) / 1000);
      printf("-- %dD aniso init times fill = %f = %dk/s optim = %f = %dk/s\n",
             gdim,t1i_aniso - t0i_aniso,nkpsi1_aniso,t2i_aniso - t1i_aniso,nkpsi2_aniso);

      // Test 1 is figure out if any close points. We don't need to hit the target
      // but we need to know there exists at least one.
      double lentar = 1.0/sqrt(2);
      int nsucc1_tot = 0, nsucc2_tot = 0, nerro1_tot = 0, nerro2_tot = 0;
      for(int niter = 0; niter < 10; niter++){

        for(double rslop = 1; rslop < 10; rslop *= 1.1){
          int nsucc1 = 0, nsucc2 = 0, nerro1 = 0, nerro2 = 0;

          for(int isamp = 0; isamp < nsamp; isamp++){
            int nclose = 0;
            for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
              double metl[nnmet];
              for(int ii = 0; ii < nnmet; ii++) metl[ii] = 0.5*( met_samp(isamp,ii) 
                                                               + msh.met(ipoin,ii) );
              getexpmet_inp<gdim,double>(metl);
              double dist = getlenedgsq<gdim>(msh.coord[ipoin],coord_samp[isamp],metl);
              if(dist > lentar) continue;
              nclose++;
            }

            kdtreeNode_aniso<gdim> node;
            for(int ii = 0; ii < gdim; ii++) node.xyz[ii] = coord_samp(isamp,ii);
            for(int ii = 0; ii < nnmet; ii++) node.lmet[ii] = met_samp(isamp,ii);

            // Add some slop
            std::vector<kdtreeNode_aniso<gdim>> howClose;
            kdtree_aniso.find_within_range(node,lentar*rslop,
                    std::back_insert_iterator<std::vector<kdtreeNode_aniso<gdim>> >(howClose));

            //BOOST_TEST((nclose == 0 && howClose.size() == 0 
            //                   || nclose > 0 && howClose.size() > 0));
            if(nclose == 0){
              if(howClose.size() == 0) nsucc1++;
              else nerro1++;
            }else{
              //if(niter == 4 || niter == 5) printf("debug lentar %e close %d vs %d\n",
              //  lentar,nclose,(int)howClose.size());
              if(howClose.size() > 0) nsucc2++;
              else nerro2++;
            }
          }
          printf("Test 1 aniso len %e slop %e nerro1 %d nsucc1 %d nerro2 %d nsucc2 %d\n",
                 lentar,rslop,nerro1,nsucc1,nerro2,nsucc2);
          if(nerro2 == 0) break;
        }
        lentar/=5;
      }

      wait();

      // Now loop over points, do quadratic search for closest as well as kdtree search
      for(int isamp = 0; isamp < nsamp; isamp++){
        double dminq = 1.0e30;
        int iminq = -1;
        for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
          double metl[nnmet];
          for(int ii = 0; ii < nnmet; ii++) metl[ii] = 0.5*( met_samp(isamp,ii) 
                                                           + msh.met(ipoin,ii) );
          getexpmet_inp<gdim,double>(metl);
          double dist = getlenedgsq<gdim>(msh.coord[ipoin],coord_samp[isamp],metl);
          if(dist > dminq) continue;
          dminq = dist;
          iminq = ipoin;
        }

        kdtreeNode_aniso<gdim> node;
        for(int ii = 0; ii < gdim; ii++) node.xyz[ii] = coord_samp(isamp,ii);
        for(int ii = 0; ii < nnmet; ii++) node.lmet[ii] = met_samp(isamp,ii);

        std::pair<typename t_KDTree_aniso::const_iterator,double> found 
          = kdtree_aniso.find_nearest(node);
        BOOST_REQUIRE(found.first != kdtree_aniso.end());
        BOOST_TEST(found.first->ipoin == iminq);

        if(found.first->ipoin != iminq){
          printf("aniso kdtree found %d we found %d\n",found.first->ipoin,iminq);
          // Distance between our node and the sought one
          double metl[nnmet];
          for(int ii = 0; ii < nnmet; ii++) metl[ii] = 0.5*( met_samp(isamp,ii) 
                                                           + msh.met(iminq,ii) );
          getexpmet_inp<gdim,double>(metl);
          double dist1 = getlenedgsq<gdim>(msh.coord[iminq],coord_samp[isamp],metl);

          // Distance to kdtree's using node info
          for(int ii = 0; ii < nnmet; ii++) metl[ii] = 0.5*( met_samp(isamp,ii) 
                                                           + found.first->lmet[ii] );
          getexpmet_inp<gdim,double>(metl);
          double dist2 = getlenedgsq<gdim>(coord_samp[isamp],found.first->xyz,metl);
          printf("dist to our node %e to kdtree %e dbg min dst = %e\n",dist1,dist2,
            dminq);

        }
      }

      printf("Debug aniso test \n");
      wait();


      // Lastly run a benchmark
      double t0_iso = get_cpu_time();
      double dum = 0;
      for(int ipoi1 = 0; ipoi1 < msh.npoin; ipoi1++){
        kdtreeNode_iso<gdim> node;
        for(int ii = 0; ii < gdim; ii++) node.xyz[ii] = msh.coord(ipoi1,ii);

        std::pair<typename t_KDTree_iso::const_iterator,double> found 
          = kdtree_iso.find_nearest(node);
        dum += found.first->ipoin;
      }
      double t1_iso = get_cpu_time();

      int nkps_iso = (int)(msh.npoin / (t1_iso - t0_iso) / 1000);
      printf("-- Benchmark %dD   iso time = %f : %dk pt/s\n",gdim,t1_iso-t0_iso,nkps_iso);


      double t0_aniso = get_cpu_time();
      for(int ipoi1 = 0; ipoi1 < msh.npoin; ipoi1++){
        kdtreeNode_aniso<gdim> node;
        for(int ii = 0; ii < gdim; ii++) node.xyz[ii] = msh.coord(ipoi1,ii);
        for(int ii = 0; ii < nnmet; ii++) node.lmet[ii] = msh.met(ipoi1,ii);

        std::pair<typename t_KDTree_aniso::const_iterator,double> found 
          = kdtree_aniso.find_nearest(node);
        dum -= found.first->ipoin;
      }
      double t1_aniso = get_cpu_time();

      int nkps_aniso = (int)(msh.npoin / (t1_aniso - t0_aniso) / 1000);
      printf("-- Benchmark %dD aniso time = %f : %dk pt/s\n",gdim,t1_aniso-t0_aniso,nkps_aniso);

    }}CT_FOR1(gdim);
  }





}// test_DIRECT