//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Mesh/MeshBase.hxx"


namespace Metris{
           
MeshArray1D<intAr1*>& MeshBase::get_iwork_tracked(MeshSize itype){
  switch(itype){
    case(MeshSize::Point):  return iwork_Point;
    case(MeshSize::Edge):   return iwork_Edge;
    case(MeshSize::Face):   return iwork_Face;
    case(MeshSize::Tetra):  return iwork_Tetra;
    case(MeshSize::BPoint): return iwork_BPoint;
    default: METRIS_THROW(TODOExcept())
  }
}
MeshArray1D<dblAr1*>& MeshBase::get_rwork_tracked(MeshSize itype){
  switch(itype){
    case(MeshSize::Point):  return rwork_Point;
    case(MeshSize::Edge):   return rwork_Edge;
    case(MeshSize::Face):   return rwork_Face;
    case(MeshSize::Tetra):  return rwork_Tetra;
    case(MeshSize::BPoint): return rwork_BPoint;
    default: METRIS_THROW(TODOExcept())
  }
}


template<typename T>
void MeshBase::free_work(int ii){
  static_assert(std::is_same<T,double>::value || std::is_same<T,int>::value);
  lwork_lock_map<T>()[ii] = 0;
}
template void MeshBase::free_work<int>(int);
template void MeshBase::free_work<double>(int);

template<typename T>
void MeshBase::free_work(int ii, WorkArray1D<T>* obj){
  free_work<T>(ii);

  auto work_tracked = get_work_tracked<T>(obj->itype);
  // Check we were storing the pointer at the right position
  METRIS_ASSERT(work_tracked[obj->iref_tracked] == &(obj->get_array()));
  work_tracked[obj->iref_tracked] = NULL;

  // See if we can reset the tracked work array.
  n_work_tracked<T>()[(int) obj->itype]--;
  METRIS_ASSERT(n_work_tracked<T>()[(int) obj->itype] >= 0);
  if(n_work_tracked<T>()[(int) obj->itype] == 0){
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

  WorkArray1D<T> ret(*this, iarray, lwork_map<T>()[iarray], itype, ntracked);
  work_tracked.stack(&(ret.get_array()));

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




}// namespace Metris