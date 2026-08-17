//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_WorkArray

#include "common_setup.hxx"

#include <functional>


using namespace Metris;

typedef MetricFieldAnalytical MFT;

template<typename T>
void test_workArray1(MeshBase& msh){

  printf("-- Testing work array type %s\n", std::is_same_v<T,int> ? "int" : "double");
  auto& msh_work = msh.debug_get_work<T>();
  auto& msh_work_lock = msh.debug_get_work_lock<T>();
  int nwork0 = msh_work.get_n();

  BOOST_REQUIRE(msh_work.get_n() == msh_work_lock.get_n());


  // With no work arrays currently in scope, ensure that nothing is locked.
  int nlocked_beg = 0;
  for(int ilock = 0; ilock < msh_work_lock.get_n(); ilock++){
    if(msh_work_lock[ilock] == 0) continue;
    nlocked_beg++;
    printf("## At start, work array %d/%d is locked type %d\n",ilock,msh_work_lock.get_n(),msh_work_lock[ilock]);
  }
  BOOST_REQUIRE(nlocked_beg == 0);


  int maxsize_work = -1;
  for(auto& work : msh_work){
    maxsize_work = MAX(maxsize_work, work.size());
  }

  BOOST_REQUIRE(nwork0 == 0 || maxsize_work > 0);

  if(maxsize_work < 0){
    printf("-- No work allocated by mesh yet, create dummy one\n");
    auto work1 = msh.get_work<T>(1000);
    maxsize_work = 1000;
    nwork0 = 1;
  }/* exiting this scope, the 1000 slot array should be freed for use */


  for(int ialloc = 0; ialloc < 100; ialloc++){
    auto work1 = msh.get_work<T>(maxsize_work - 1);
    BOOST_REQUIRE(work1.get_n() == maxsize_work - 1);
    BOOST_REQUIRE(work1.get_array().get_n() == maxsize_work - 1);
  }
  /* Since we stayed under the size of existing free arrays,*/
  /* we expect this size to be unchanged. */ 
  BOOST_REQUIRE(nwork0 == msh_work.get_n());
  printf("-- SUCCESS: 100 work arrays created and destroyed,  no memory allocated\n");
  for(int ialloc = 0; ialloc < 10; ialloc++){
    auto work1 = msh.get_work<T>(maxsize_work + ialloc + 1); 
    BOOST_REQUIRE(work1.get_n() == maxsize_work + ialloc + 1);
    BOOST_REQUIRE(work1.get_array().get_n() == maxsize_work + ialloc + 1);
  }
  BOOST_REQUIRE(msh_work.get_n() == nwork0 + 10);
  printf("-- SUCCESS: 10 work arrays created and destroyed of increasing sizes, 10 memory regions allocated\n");

  // With no work arrays currently in scope, ensure that nothing is locked.
  int nlocked_end = 0;
  for(int ilock = 0; ilock < msh_work_lock.get_n(); ilock++){
    if(msh_work_lock[ilock] == 0) continue;
    nlocked_end++;
    printf("## Work array %d is locked type %d\n",ilock,msh_work_lock[ilock]);
  }
  BOOST_REQUIRE(nlocked_end == 0);
  printf("-- SUCCESS: mesh has currently 0 locked arrays\n");
}

template void test_workArray1<int>(MeshBase& msh);
template void test_workArray1<double>(MeshBase& msh);



template<typename T>
void test_workArray2(MeshBase& msh){

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

  for(auto work_type : work_types){
    const MeshSize itype = std::get<0>(work_type);
    const int& nentt = std::get<1>(work_type);
    const int& mentt = std::get<2>(work_type);
    auto set_nentt = std::get<3>(work_type);

    printf("-- Start test data %s entity type %s.\n", std::is_same_v<T,int> ? "int" : "double", work_type_names[itype].c_str());

    auto itype_work = msh.get_work<T>(itype);
    BOOST_REQUIRE(itype_work.get_n() == nentt);
    BOOST_REQUIRE(itype_work.size() + 1 >= mentt + 1);

    MeshArray1D<T> save_itype_work(nentt);
    int nent0 = nentt;
    for(int ii = 0; ii < nentt; ii++){
      itype_work[ii] = (T) ii;
      save_itype_work[ii] = itype_work[ii];
    }

    int ment0 = mentt;
    set_nentt(mentt + 1);
    BOOST_REQUIRE(mentt > ment0);
    if(itype_work.size() + 2 < mentt + 2){
      printf("work type %s ment0 = %d mentt = %d, itype_work.size() = %d\n",work_type_names[itype].c_str(),ment0,mentt,itype_work.size());
    }
    BOOST_REQUIRE(itype_work.size() + 2 >= mentt + 2);
    BOOST_REQUIRE(itype_work.get_n() == nentt);
    for(int ii = 0; ii < nent0; ii++){
      //if(itype_work[ii] != save_itype_work[ii]){
      //  fmt::print("## ii = %d itype_work[ii] = {}, save_itype_work[ii] = {}\n", ii, itype_work[ii], save_itype_work[ii]);
      //}
      //fmt::print("## ii = %d itype_work[ii] = {}, save_itype_work[ii] = {}\n", ii, itype_work[ii], save_itype_work[ii]);
      METRIS_ASSERT_MSG(itype_work[ii] == save_itype_work[ii], "itype_work[ii] = {}, save_itype_work[ii] = {} ii = {} nent0 = {}", itype_work[ii], save_itype_work[ii], ii, nent0);
      BOOST_REQUIRE(itype_work[ii] == save_itype_work[ii]);
    }
  }
}
template void test_workArray2<int>(MeshBase& msh);
template void test_workArray2<double>(MeshBase& msh);


BOOST_AUTO_TEST_CASE(test_WorkArray) 
{

  std::string mesh = METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k";


  cargHandler arg("-in " + mesh + "  -anamet 1 -verb 0");
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
  
  // Begin by checking all work arrays are unlocked if not 
  // within a routine that uses any. 
  intAr1& msh_iwork_lock = msh.debug_get_iwork_lock();
  intAr1& msh_rwork_lock = msh.debug_get_rwork_lock();
  for(int ilock : msh_iwork_lock) BOOST_REQUIRE(ilock == 0);
  for(int ilock : msh_rwork_lock) BOOST_REQUIRE(ilock == 0);

  // Create a work array, then force a reallocation. Ensure the initial 
  // work array remains valid.

  #if 0
  intWrkAr1 idbg = msh.get_iwork(MeshSize::Point);
  printf("Debug in main routine address of array %p\n", &(idbg.get_array()));

  {
    printf("Restricted scope\n");
    intWrkAr1 idbg1 = msh.get_iwork(MeshSize::Point);
    printf("Debug in main routine address of array %p\n", &(idbg1.get_array()));
  }

  printf("After scope\n");
  intWrkAr1 idbg2 = msh.get_iwork(MeshSize::Point);
  printf("Debug in main routine address of array %p\n", &(idbg2.get_array()));

  printf("Double tests\n");
  intWrkAr1 rdbg = msh.get_iwork(MeshSize::Point);
  printf("Debug in main routine address of array %p\n", &(rdbg.get_array()));

  {
    printf("Restricted scope\n");
    intWrkAr1 rdbg1 = msh.get_iwork(MeshSize::Point);
    printf("Debug in main routine address of array %p\n", &(rdbg1.get_array()));
  }

  printf("After scope\n");
  intWrkAr1 rdbg2 = msh.get_iwork(MeshSize::Point);
  printf("Debug in main routine address of array %p\n", &(rdbg2.get_array()));

  exit(0);
  #endif

  #if 0
  intWrkAr1 idbg = msh.get_iwork(MeshSize::Point);
  printf("Debug 1 in main routine address of array %p\n", &(idbg.get_array()));
  msh.set_npoin(msh.npoin+1);
  printf("Debug 2 in main routine address of array %p\n", &(idbg.get_array()));

  // Force reallocation of tracking array.
  auto work_tracked = msh.debug_get_work_tracked<int>(MeshSize::Point);
  work_tracked.allocate(work_tracked.size() + 1);
  printf("Debug 3 in main routine address of array %p\n", &(idbg.get_array()));
  msh.set_npoin(msh.npoin+1);
  printf("Debug 4 in main routine address of array %p\n", &(idbg.get_array()));

  exit(0);

  #endif

  #if 0
  {
    MeshArray1D<intAr1>& msh_iwork = msh.debug_get_iwork();
    intWrkAr1 wrkar0 = msh.get_iwork(100);
    msh_iwork_lock.inc_n();
    msh_iwork.inc_n();
    BOOST_REQUIRE(wrkar0.size() == 100);
    BOOST_REQUIRE(wrkar0.get_n() == 100);

    for(int ilock = 0; ilock < msh_iwork_lock.get_n(); ilock++){
      printf("## At start 1, work array %d/%d is locked type %d\n",ilock,msh_iwork_lock.get_n(),msh_iwork_lock[ilock]);
    }

    msh_iwork_lock.set_n(msh_iwork_lock.get_n() - 1);
    msh_iwork.set_n(msh_iwork.get_n() - 1);
  }

  for(int ilock = 0; ilock < msh_iwork_lock.get_n(); ilock++){
    printf("## At start 2, work array %d/%d is locked type %d\n",ilock,msh_iwork_lock.get_n(),msh_iwork_lock[ilock]);
  }

  #endif

  test_workArray1<int>(msh);
  test_workArray1<double>(msh);


  test_workArray2<int>(msh);
  test_workArray2<double>(msh);



  // Test currently n_iwork_tracked and n_rwork_tracked are all 
  // zeroes, since no tracked work arrays currently alive. 
  const int* n_iwork_tracked = msh.debug_n_work_tracked<int>();
  const int* n_rwork_tracked = msh.debug_n_work_tracked<double>();
  for(int ii = 0; ii < (int) MeshSize::nTrackedType; ii++){
    BOOST_REQUIRE(n_iwork_tracked[ii] == 0);
    BOOST_REQUIRE(n_rwork_tracked[ii] == 0);
  }

  // Create two tracked arrays and check the corresponding
  // entry is two.
  {
    auto iwork1 = msh.get_iwork(MeshSize::Point);
    auto iwork2 = msh.get_iwork(MeshSize::Point);
    BOOST_REQUIRE(iwork1.get_n() == msh.npoin);
    BOOST_REQUIRE(iwork2.get_n() == msh.npoin);

    // Check the tracked work arrays are now two.
    BOOST_REQUIRE(n_iwork_tracked[(int) MeshSize::Point] == 2);
  }
  BOOST_REQUIRE(n_iwork_tracked[(int) MeshSize::Point] == 0);

}// end boost test case
