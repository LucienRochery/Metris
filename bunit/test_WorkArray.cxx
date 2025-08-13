//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_WorkArray

#include "common_setup.hxx"


using namespace Metris;

typedef MetricFieldAnalytical MFT;

BOOST_AUTO_TEST_CASE(test_WorkArray) 
{

  std::string mesh = METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k";


  cargHandler arg("-in " + mesh + "  -anamet 1 -verb 0");
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
  
  // Begin by checking all work arrays are unlocked if not 
  // within a routine that uses any. 
  const intAr1& msh_iwork_lock = msh.debug_get_iwork_lock();
  const intAr1& msh_rwork_lock = msh.debug_get_rwork_lock();
  for(int ilock : msh_iwork_lock) BOOST_REQUIRE(ilock == 0);
  for(int ilock : msh_rwork_lock) BOOST_REQUIRE(ilock == 0);

  
  #define TEST_WORK_TYPE1(T)\
  printf("-- Testing work array type %s\n", #T);\
  auto& msh_##T##work = msh.debug_get_##T##work();\
  int n##T##work0 = msh_##T##work.get_n();\
\
  int maxsize_##T##work = -1;\
  for(auto& T##work : msh_##T##work){\
    maxsize_##T##work = MAX(maxsize_##T##work, T##work.size());\
  }\
\
  BOOST_REQUIRE(n##T##work0 == 0 || maxsize_##T##work > 0);\
\
  if(maxsize_##T##work < 0){\
    printf("-- No " #T "work allocated by mesh yet, create dummy one\n");\
    auto T##work1 = msh.get_##T##work(1000);\
    maxsize_##T##work = 1000;\
    n##T##work0 = 1;\
  }/* exiting this scope, the 1000 slot array should be freed for use */\
\
\
  for(int ialloc = 0; ialloc < 100; ialloc++){\
    auto T##work1 = msh.get_##T##work(maxsize_##T##work - 1);\
    BOOST_REQUIRE(T##work1.get_n() == maxsize_##T##work - 1);\
    BOOST_REQUIRE(T##work1.get_array().get_n() == maxsize_##T##work - 1);\
  }\
  /* Since we stayed under the size of existing free arrays,*/\
  /* we expect this size to be unchanged.*/ \
  BOOST_REQUIRE(n##T##work0 == msh_##T##work.get_n());\
  printf("-- SUCCESS: 100 work arrays created and destroyed,  no memory allocated\n");\
  for(int ialloc = 0; ialloc < 10; ialloc++){\
    auto T##work1 = msh.get_##T##work(maxsize_##T##work + ialloc + 1);\
    BOOST_REQUIRE(T##work1.get_n() == maxsize_##T##work + ialloc + 1);\
    BOOST_REQUIRE(T##work1.get_array().get_n() == maxsize_##T##work + ialloc + 1);\
  }\
  BOOST_REQUIRE(msh_##T##work.get_n() == n##T##work0 + 10);\
  printf("-- SUCCESS: 10 work arrays created and destroyed of increasing sizes, 10 memory regions allocated\n");


  TEST_WORK_TYPE1(i);
  TEST_WORK_TYPE1(r);
  #undef TEST_WORK_TYPE1

  // Now we're going to tie work arrays to mesh entities.
  // We will then grow the number of entities to force a reallocation,
  // and make sure we find the work array reallocated.
  std::vector<std::tuple<MeshSize, const int&, const int&, std::function<void(int)>>> work_types = {
     {MeshSize::Point, msh.npoin, msh.mpoin, [&msh](int n){ msh.set_npoin(n);}}
    ,{MeshSize::Edge, msh.nedge, msh.medge, [&msh](int n){ msh.set_nedge(n);}}
    ,{MeshSize::Face, msh.nface, msh.mface, [&msh](int n){ msh.set_nface(n);}}
    ,{MeshSize::Tetra, msh.nelem, msh.melem, [&msh](int n){ msh.set_nelem(n);}}
    ,{MeshSize::BPoint, msh.nbpoi, msh.mbpoi, [&msh](int n){ msh.set_nbpoi(n);}}
  };

  std::map<MeshSize, std::string> work_type_names = {
     {MeshSize::Point, "Point"}
    ,{MeshSize::Edge, "Edge"}
    ,{MeshSize::Face, "Face"}
    ,{MeshSize::Tetra, "Tetra"}
    ,{MeshSize::BPoint, "BPoint"}
  };

  #define TEST_WORK_TYPE2(DATATYPE)\
  {\
    for(auto work_type : work_types){\
      const MeshSize itype = std::get<0>(work_type);\
      const int& nentt = std::get<1>(work_type);\
      const int& mentt = std::get<2>(work_type);\
      auto set_nentt = std::get<3>(work_type);\
\
      printf("-- Start test data " #DATATYPE " entity type %s.\n", work_type_names[itype].c_str());\
\
      auto itype_work = msh.get_work<DATATYPE>(itype);\
      BOOST_REQUIRE(itype_work.get_n() == nentt);\
      BOOST_REQUIRE(itype_work.size() >= mentt);\
\
      MeshArray1D<DATATYPE> save_itype_work(nentt);\
      int nent0 = nentt;\
      for(int ii = 0; ii < nentt; ii++){\
        itype_work[ii] = (DATATYPE) ii;\
        save_itype_work[ii] = itype_work[ii];\
      }\
\
      int ment0 = mentt;\
      set_nentt(mentt + 1);\
      BOOST_REQUIRE(mentt > ment0);\
      BOOST_REQUIRE(itype_work.size() >= mentt);\
      BOOST_REQUIRE(itype_work.get_n() == nentt);\
      for(int ii = 0; ii < nent0; ii++){\
        BOOST_REQUIRE(itype_work[ii] == save_itype_work[ii]);\
      }\
    }\
  }

  TEST_WORK_TYPE2(int);
  TEST_WORK_TYPE2(double);
  #undef TEST_WORK_TYPE2


}// end boost test case
