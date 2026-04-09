//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Mesh/MeshBase.hxx"
#include "fmt/format.h"
#include "../utils/fmt_formatters.hxx"

namespace Metris{

MeshArray1D<intWrkAr1*>& MeshBase::get_iwork_tracked(MeshSize itype){
  switch(itype){
    case(MeshSize::Point):  return iwork_Point;
    case(MeshSize::Edge):   return iwork_Edge;
    case(MeshSize::Face):   return iwork_Face;
    case(MeshSize::Tetra):  return iwork_Tetra;
    case(MeshSize::BPoint): return iwork_BPoint;
    default: METRIS_THROW_MSG("TODO: implement work type {}",(int)itype)
  }
}
MeshArray1D<dblWrkAr1*>& MeshBase::get_rwork_tracked(MeshSize itype){
  switch(itype){
    case(MeshSize::Point):  return rwork_Point;
    case(MeshSize::Edge):   return rwork_Edge;
    case(MeshSize::Face):   return rwork_Face;
    case(MeshSize::Tetra):  return rwork_Tetra;
    case(MeshSize::BPoint): return rwork_BPoint;
    default: METRIS_THROW_MSG("TODO: implement work type {}",(int)itype)
  }
}


template<typename T>
void MeshBase::free_work(int ii){
  lwork_lock_map<T>()[ii] = 0;
}
template void MeshBase::free_work<int>(int);
template void MeshBase::free_work<double>(int);

template<typename T>
void MeshBase::free_work(int ii, WorkArray1D<T>* obj){
  free_work<T>(ii);

  // This ampersand is load bearing:
  // auto& is a dangerous return type that doesn't actually return
  // a reference if auto myret = function_that_returns_auto&().
  auto& work_tracked = get_work_tracked<T>(obj->itype);
  //if(obj->itype == MeshSize::Point && std::is_same_v<T,int>){
  //  fmt::print("Debug: free_work<int> Point addr {} iref_tracked {}\n",(void*) obj, obj->iref_tracked);
  //  fmt::print("Debug 0: work_tracked = {}\n",work_tracked);
  //}

  // Check we were storing the pointer at the right position
  METRIS_ASSERT(work_tracked[obj->iref_tracked] == obj);
  work_tracked[obj->iref_tracked] = NULL;

  //if(obj->itype == MeshSize::Point && std::is_same_v<T,int>){
  //  fmt::print("Debug 1: work_tracked = {}\n",work_tracked);
  //  fmt::print("         ntracked = {}\n",n_work_tracked<T>()[(int) obj->itype]);
  //}

  // See if we can reset the tracked work array.
  n_work_tracked<T>()[(int) obj->itype]--;
  METRIS_ASSERT(n_work_tracked<T>()[(int) obj->itype] >= 0);
  if(n_work_tracked<T>()[(int) obj->itype] == 0){
    // Check all currenty elements are empty:
    #ifndef NDEBUG
    int ii = 0;
    int ntracked = 0;
    for(auto itracked : work_tracked){
      ii++;
      if(itracked == NULL) continue;
      fmt::print("## ii = {}/{} itracked = {} type {} data {}\n", ii-1, work_tracked.get_n(), (void*)itracked, (int) obj->itype, std::is_same_v<T,int> ? "int" : "double");
      ntracked++;
    }
    METRIS_ASSERT(ntracked == 0);
    #endif
    work_tracked.set_n(0);
  }

}
template void MeshBase::free_work<int>(int, WorkArray1D<int>* obj);
template void MeshBase::free_work<double>(int, WorkArray1D<double>* obj);

// To be used both by get_iwork and get_tagarray
// Returns the index of the locked array
template<typename T>
inline int get_locked_array(int nn, MeshArray1D<T> &array_pool, intAr1 &array_locks){
  int npool = array_pool.get_n();
  int ibest = -1, nbest = -1;
  for(int iarray = 0; iarray < npool; iarray++){
    if(array_locks[iarray] > 0) continue;
    int array_size = array_pool[iarray].size();
    if(array_size < nn) continue;

    if(array_size < nbest || nbest < 0){
      ibest = iarray;
      nbest = array_size;
    }
  }

  return ibest;
}


template<typename T>
WorkArray1D<T> MeshBase::get_work(int nn){
  int iarray = get_locked_array<MeshArray1D<T>>(nn, lwork_map<T>(), lwork_lock_map<T>());
  if(iarray < 0){
    iarray = lwork_map<T>().get_n();
    lwork_map<T>().inc_n();
    lwork_lock_map<T>().inc_n();
    lwork_map<T>()[iarray].allocate(nn);
  }
  lwork_map<T>()[iarray].set_n(nn);
  lwork_lock_map<T>()[iarray] = (int) MeshSize::Untracked;
  return WorkArray1D<T>(*this, iarray, lwork_map<T>()[iarray]);
}

template intWrkAr1 MeshBase::get_work<int   >(int nn);
template dblWrkAr1 MeshBase::get_work<double>(int nn);


template<typename T>
WorkArray1D<T> MeshBase::get_work(MeshSize itype){
  int nentt = -1, mentt = -1;
  get_nMeshSize(itype, &nentt, &mentt);

  int iarray = get_locked_array<MeshArray1D<T>>(mentt, lwork_map<T>(), lwork_lock_map<T>());
  if(iarray < 0){
    iarray = lwork_map<T>().get_n();
    lwork_map<T>().inc_n();
    lwork_lock_map<T>().inc_n();
    lwork_map<T>()[iarray].allocate(mentt);
  }
  lwork_map<T>()[iarray].set_n(nentt);
  lwork_lock_map<T>()[iarray] = (int) itype;

  auto& work_tracked = get_work_tracked<T>(itype);
  int ntracked = work_tracked.get_n();

  // With mandatory RVO since C++17, this here ret should be
  // the same object we find in the caller routine... so its address
  // is safe to use.
  WorkArray1D<T> ret(*this, iarray, lwork_map<T>()[iarray], itype, ntracked);
  work_tracked.stack(&ret);

  //if(itype == MeshSize::Point && std::is_same_v<T,int>) fmt::print("Debug type int Point track {} new ntrack {}\n",(void*) &ret, n_work_tracked<T>()[(int) itype]+1);
  n_work_tracked<T>()[(int) itype]++;

  return ret;
}

template intWrkAr1 MeshBase::get_work<int   >(MeshSize itype);
template dblWrkAr1 MeshBase::get_work<double>(MeshSize itype);


intWrkAr1 MeshBase::get_iwork(int nn){
  return get_work<int>(nn);
}
dblWrkAr1 MeshBase::get_rwork(int nn){
  return get_work<double>(nn);
}
intWrkAr1 MeshBase::get_iwork(MeshSize size){
  return get_work<int>(size);
}
dblWrkAr1 MeshBase::get_rwork(MeshSize size){
  return get_work<double>(size);
}


TagArray MeshBase::get_tagarray(MeshSize itype){
  METRIS_ASSERT((int) itype > 0);

  int nentt = -1, mentt = -1;
  get_nMeshSize(itype, &nentt, &mentt);

  int iarray = get_locked_array<intAr1>(mentt, tagarrs, tagarr_locks);

  if(iarray < 0){
    iarray = tagarrs.get_n();

    tagarrs.inc_n();
    tagarr_locks.inc_n();
    itags.stack(0);

    tagarrs[iarray].allocate(mentt);
    tagarrs[iarray].set_n(nentt);
    tagarrs[iarray].fill(0);

  }

  itags[iarray]++;
  tagarr_locks[iarray] = (int) itype;
  tagarrs[iarray].set_n(nentt); // could be originally larger than needed

  return TagArray(*this, iarray, tagarrs[iarray], itags[iarray]);
}

void MeshBase::reset_tags(){

  for (int ii = 0; ii < METRIS_MAXTAGS; ii++) tag[ii] = 1;

  tet2tag.fill(0);
  fac2tag.fill(0);
  edg2tag.fill(0);
  poi2tag.fill(0);
  bpo2tag.fill(0);
  dom2tag.fill(0);
  cfa2tag.fill(0);
  ced2tag.fill(0);
  cno2tag.fill(0);

}

void MeshBase::update_tracked_work_arrays(MeshSize itype, int mentt, int nentt){
  MeshArray1D<intWrkAr1*>& iwork_tracked = get_work_tracked<int>(itype);
  for(int ii = 0; ii < iwork_tracked.get_n(); ii++){
    intWrkAr1* wrkArMem = iwork_tracked[ii];
    if(wrkArMem == NULL) continue;
    // Reallocate our pooled array, then replace the Work's array
    // memory with this one.
    int ilock = wrkArMem->ilock; // index in the work array pool iwork
    iwork[ilock].allocate(mentt);
    iwork[ilock].set_n(nentt);

    wrkArMem->array = iwork[ilock]; // we are friend of WorkArray1D
  }

  MeshArray1D<dblWrkAr1*>& rwork_tracked = get_work_tracked<double>(itype);
  for(int ii = 0; ii < rwork_tracked.get_n(); ii++){
    dblWrkAr1* wrkArMem = rwork_tracked[ii];
    if(wrkArMem == NULL) continue;
    // Reallocate our pooled array, then replace the Work's array
    // memory with this one.
    int ilock = wrkArMem->ilock; // index in the work array pool iwork
    rwork[ilock].allocate(mentt);
    rwork[ilock].set_n(nentt);

    wrkArMem->array = rwork[ilock]; // we are friend of WorkArray1D
  }
}

}// namespace Metris