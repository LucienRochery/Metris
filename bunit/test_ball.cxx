//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_ball

#include "common_setup.hxx"

using namespace Metris;

typedef MetricFieldAnalytical MFT;

// 1.0e-8 relative error -> 1.0e-6% utf::tolerance()
BOOST_AUTO_TEST_CASE(test_ball) 
{ 

  // bool is whether straight
  std::vector<std::string> meshes = {
   METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
  ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
  };

  for(std::string mesh : meshes){

    cargHandler arg("-in " + mesh + "  -anamet 1 -verb 0");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
    msh.cleanup();

    std::cout<<"-- Mesh "<<mesh<<" dim "<<msh.idim<<"\n";

    intAr1 lbedg0(10), lbfac0(10), lbtet0(10);
    intAr1 lbedg1(10), lbfac1(10), lbtet1(10);
    // Not going to be reset by ball in 3D
    lbtet0.set_n(0);
    lbtet1.set_n(0);
    MeshArray1D<intAr1> poi2edg(msh.npoin);
    MeshArray1D<intAr1> poi2fac(msh.npoin);
    MeshArray1D<intAr1> poi2tet(msh.npoin);
    auto get_poi2ent = [&](int tdim) -> MeshArray1D<intAr1>& {
           if(tdim == 1) return poi2edg;
      else if(tdim == 2) return poi2fac;
      else               return poi2tet;
    };

    auto get_lbent0 = [&](int tdim) -> intAr1& {
           if(tdim == 1) return lbedg0;
      else if(tdim == 2) return lbfac0;
      else               return lbtet0;
    };
    auto get_lbent1 = [&](int tdim) -> intAr1& {
           if(tdim == 1) return lbedg1;
      else if(tdim == 2) return lbfac1;
      else               return lbtet1;
    };

    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      if(msh.poi2ent(ipoin,0) < 0) continue;
      poi2edg[ipoin].allocate(4);
      poi2fac[ipoin].allocate(10);
      if(msh.get_tdim() >= 3) poi2tet[ipoin].allocate(100);
    }
    for(int tdim = 1; tdim <= msh.get_tdim(); tdim++){
      const intAr2& ent2poi = msh.ent2poi(tdim);
      MeshArray1D<intAr1>& poi2ent = get_poi2ent(tdim);
      for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
        if(isdeadent(ientt,ent2poi)) continue;
        for(int iver = 0; iver < tdim + 1; iver++){
          int ipoin = ent2poi(ientt, iver);
          poi2ent[ipoin].stack(ientt);
        }
      }// for ientt
    }// for tdim

    printf("-- Test 1: empty ball\n");

    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      if(msh.poi2ent(ipoin,0) < 0) continue;
      int iopen;

      int ierro = ball(msh, ipoin, lbedg0, lbfac0, lbtet0, &iopen, false, 0);
      BOOST_REQUIRE(ierro == 0);

      msh.tag[0]++;
      for(int iedge : lbedg0){
        BOOST_TEST(msh.edg2tag(0, iedge) < msh.tag[0]); // check no duplicates
        msh.edg2tag(0, iedge) = msh.tag[0];
      }
      for(int iface : lbfac0){
        BOOST_TEST(msh.fac2tag(0, iface) < msh.tag[0]); // check no duplicates
        msh.fac2tag(0, iface) = msh.tag[0];
      }
      for(int itetr : lbtet0){
        BOOST_TEST(msh.tet2tag(0, itetr) < msh.tag[0]); // check no duplicates
        msh.tet2tag(0, itetr) = msh.tag[0];
      }

      for(int tdim = 1; tdim <= msh.get_tdim(); tdim++){
        const intAr2& ent2tag = msh.ent2tag(tdim);
        const intAr1& poi2ent = get_poi2ent(tdim)[ipoin];
        const intAr1& lbent0  = get_lbent0(tdim);
        BOOST_TEST(poi2ent.get_n() == lbent0.get_n());
        for(int ientt : poi2ent){
          BOOST_TEST(ent2tag(0,ientt) == msh.tag[0]);
        }
      }// for tdim

    }// for ipoin


    printf("-- Test 2: append\n");

    // Next test ball's append feature
    // We'll be using edges for this.
    const int tdim = msh.get_tdim();
    const intAr2& ent2poi = msh.ent2poi(tdim);
    const auto lnoed = tdim == 1 ? lnoed1 : 
                       tdim == 2 ? lnoed2 : lnoed3;
    const int nedgl = (tdim*(tdim+1))/2;
    for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
      if(isdeadent(ientt,ent2poi)) continue;
      for(int ied = 0; ied < nedgl; ied++){
        int ipoi1 = ent2poi(ientt, lnoed[ied][0]);
        int ipoi2 = ent2poi(ientt, lnoed[ied][1]);
        lbedg0.set_n(0);
        lbfac0.set_n(0);
        lbtet0.set_n(0);
        //msh.param->iverb = 5;
        //msh.param->ivdepth = 5;
        //printf(" \n\n-- debug");
        int iopen, ierro;
        ierro = ball(msh, ipoi1, lbedg0, lbfac0, lbtet0, &iopen, true, 0);
        BOOST_REQUIRE(ierro == 0);
        ierro = ball(msh, ipoi2, lbedg0, lbfac0, lbtet0, &iopen, true, 0);
        BOOST_REQUIRE(ierro == 0);

        bool isucc = true;
        // Gather the common ball
        msh.tag[1]++;
        for(int tdim1 = 1; tdim1 <= msh.get_tdim(); tdim1++){
          intAr1& lbent1  = get_lbent1(tdim1);
          lbent1.set_n(0);
          for(int ipoin : {ipoi1, ipoi2}){
            const intAr1& poi2ent = get_poi2ent(tdim1)[ipoin];
            for(int ientt : poi2ent){
              if(msh.ent2tag(tdim1)(1, ientt) == msh.tag[1]) continue;
              msh.ent2tag(tdim1)(1, ientt) = msh.tag[1];
              lbent1.stack(ientt);
            }// for ientt
          }// for ipoin
        }// for tdim1

        // Check they coincide
        msh.tag[0]++;
        for(int tdim1 = 1; tdim1 <= msh.get_tdim(); tdim1++){
          const intAr1& lbent0 = get_lbent0(tdim1);
          const intAr1& lbent1 = get_lbent1(tdim1);
          BOOST_TEST(lbent1.get_n() == lbent0.get_n());
          if(!(lbent1.get_n() == lbent0.get_n())) isucc = false;
          for(int ientt : lbent0){
            BOOST_TEST(msh.ent2tag(tdim1)(0,ientt) < msh.tag[0]);
            msh.ent2tag(tdim1)(0,ientt) = msh.tag[0];
          }
          for(int ientt : lbent1){
            // Check was seen in lbent0.
            BOOST_TEST(msh.ent2tag(tdim1)(0,ientt) == msh.tag[0]);
            //printf("Debug tdim1 %d ientt %d tag %d other %d\n",tdim1,ientt,msh.ent2tag(tdim1)(0,ientt),msh.tag[0]);
            if(!(msh.ent2tag(tdim1)(0,ientt) == msh.tag[0])) isucc = false;
          }
        }// for tdim1


        if(!isucc){
          printf("ball found: \n");
          printf("lbedg0 = "); lbedg0.print();
          printf("lbfac0 = "); lbfac0.print();
          printf("lbtet0 = "); lbtet0.print();
          printf("quadratic:\n");
          printf("lbedg1 = "); lbedg1.print();
          printf("lbfac1 = "); lbfac1.print();
          printf("lbtet1 = "); lbtet1.print();
          printf("ipoi1 = %d ipoi2 = %d\n",ipoi1,ipoi2);

          if(msh.poi2bpo[ipoi1] >= 0){
            printf("ipoi1 bpo:\n");
            for(int ibpoi = msh.poi2bpo[ipoi1]; ibpoi >= 0;ibpoi = msh.bpo2ibi(ibpoi,3)){
              printf("%d : ",ibpoi);
              intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
            }
          }
          if(msh.poi2bpo[ipoi2] >= 0){
            printf("ipoi2 bpo:\n");
            for(int ibpoi = msh.poi2bpo[ipoi2]; ibpoi >= 0;ibpoi = msh.bpo2ibi(ibpoi,3)){
              printf("%d : ",ibpoi);
              intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
            }
          }
          writeMesh("test_ball",msh);
          wait();
        }
      }

    }// for ientt

    check_topo(msh,2);

    printf("-- Test 3: shell then 2 ball append\n");
    for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
      if(isdeadent(ientt,ent2poi)) continue;
      for(int ied = 0; ied < nedgl; ied++){
        int ipoi1 = ent2poi(ientt, lnoed[ied][0]);
        int ipoi2 = ent2poi(ientt, lnoed[ied][1]);
        lbedg0.set_n(0);
        lbfac0.set_n(0);
        lbtet0.set_n(0);

        int iopen, ierro;
        shell(msh, ipoi1, ipoi2, tdim, ientt, lbedg0, lbfac0, lbtet0, &iopen);

        //msh.param->iverb = 5;
        //msh.param->ivdepth = 5;
        //printf(" \n\n-- debug");
        ierro = ball(msh, ipoi1, lbedg0, lbfac0, lbtet0, &iopen, true, 0);
        BOOST_REQUIRE(ierro == 0);
        ierro = ball(msh, ipoi2, lbedg0, lbfac0, lbtet0, &iopen, true, 0);
        BOOST_REQUIRE(ierro == 0);

        bool isucc = true;
        // Gather the common ball
        msh.tag[1]++;
        for(int tdim1 = 1; tdim1 <= msh.get_tdim(); tdim1++){
          intAr1& lbent1  = get_lbent1(tdim1);
          lbent1.set_n(0);
          for(int ipoin : {ipoi1, ipoi2}){
            const intAr1& poi2ent = get_poi2ent(tdim1)[ipoin];
            for(int ientt : poi2ent){
              if(msh.ent2tag(tdim1)(1, ientt) == msh.tag[1]) continue;
              msh.ent2tag(tdim1)(1, ientt) = msh.tag[1];
              lbent1.stack(ientt);
            }// for ientt
          }// for ipoin
        }// for tdim1

        // Check they coincide
        msh.tag[0]++;
        for(int tdim1 = 1; tdim1 <= msh.get_tdim(); tdim1++){
          const intAr1& lbent0 = get_lbent0(tdim1);
          const intAr1& lbent1 = get_lbent1(tdim1);
          BOOST_TEST(lbent1.get_n() == lbent0.get_n());
          if(!(lbent1.get_n() == lbent0.get_n())) isucc = false;
          for(int ientt : lbent0){
            BOOST_TEST(msh.ent2tag(tdim1)(0,ientt) < msh.tag[0]);
            msh.ent2tag(tdim1)(0,ientt) = msh.tag[0];
          }
          for(int ientt : lbent1){
            // Check was seen in lbent0.
            BOOST_TEST(msh.ent2tag(tdim1)(0,ientt) == msh.tag[0]);
            //printf("Debug tdim1 %d ientt %d tag %d other %d\n",tdim1,ientt,msh.ent2tag(tdim1)(0,ientt),msh.tag[0]);
            if(!(msh.ent2tag(tdim1)(0,ientt) == msh.tag[0])) isucc = false;
          }
        }// for tdim1


        if(!isucc){
          printf("ball found: \n");
          printf("lbedg0 = "); lbedg0.print();
          printf("lbfac0 = "); lbfac0.print();
          printf("lbtet0 = "); lbtet0.print();
          printf("quadratic:\n");
          printf("lbedg1 = "); lbedg1.print();
          printf("lbfac1 = "); lbfac1.print();
          printf("lbtet1 = "); lbtet1.print();
          printf("ipoi1 = %d ipoi2 = %d\n",ipoi1,ipoi2);

          if(msh.poi2bpo[ipoi1] >= 0){
            printf("ipoi1 bpo:\n");
            for(int ibpoi = msh.poi2bpo[ipoi1]; ibpoi >= 0;ibpoi = msh.bpo2ibi(ibpoi,3)){
              printf("%d : ",ibpoi);
              intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
            }
          }
          if(msh.poi2bpo[ipoi2] >= 0){
            printf("ipoi2 bpo:\n");
            for(int ibpoi = msh.poi2bpo[ipoi2]; ibpoi >= 0;ibpoi = msh.bpo2ibi(ibpoi,3)){
              printf("%d : ",ibpoi);
              intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
            }
          }
          writeMesh("test_ball",msh);
          wait();
        }
      }

    }// for ientt
    check_topo(msh,2);


  }// for testcase

}