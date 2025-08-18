//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "WorkArray.hxx"
#include "MeshBase.hxx"

#include "fmt/format.h"

namespace Metris{

template<typename T>
WorkArray1D<T>::WorkArray1D(MeshBase& msh_, int ilock_, MeshArray1D<T>& array_) 
    : ilock(ilock_), msh(msh_), itype(MeshSize::Untracked), iref_tracked(-1) {
  array = array_;
}


template<typename T>
WorkArray1D<T>::WorkArray1D(MeshBase& msh_, int ilock_, MeshArray1D<T>& array_, MeshSize itype_, int iref_tracked) 
    : ilock(ilock_), msh(msh_), itype(itype_), iref_tracked(iref_tracked) {
  array = array_;
}


// Copy and move constructors called when returning a temp
template<typename T>
WorkArray1D<T>::WorkArray1D(const WorkArray1D& other) 
    : ilock(other.ilock), msh(other.msh), itype(other.itype), iref_tracked(other.iref_tracked) {

  printf("Debug copy constructor called\n");
  array = other.array;  // Copy the MeshArray1D
  
  // Update tracking pointer to point to our copy
  if(itype != MeshSize::Untracked) {
    auto work_tracked = msh.get_work_tracked<T>(itype);
    work_tracked[iref_tracked] = this;
  }
}

// Internal use:
template<typename T>
WorkArray1D<T>::WorkArray1D(WorkArray1D&& other) 
    : ilock(other.ilock), msh(other.msh), itype(other.itype), iref_tracked(other.iref_tracked) {
  array = std::move(other.array);  // Move the MeshArray1D

  printf("Debug move constructor called\n");
  // Update tracking pointer to point to our copy
  if(itype != MeshSize::Untracked) {
    auto work_tracked = msh.get_work_tracked<T>(itype);
    work_tracked[iref_tracked] = this;
  }
  
  // Invalidate the moved-from object
  other.itype = MeshSize::Untracked;
}

template<typename T>
WorkArray1D<T>& WorkArray1D<T>::operator=(const WorkArray1D& other) {
    printf("Debug copy assignment called\n");
    if (this != &other) {
        // Don't reassign msh, ilock as they should be the same
        itype = other.itype;
        iref_tracked = other.iref_tracked;
        array = other.array;
        
        // Update tracking pointer
        if(itype != MeshSize::Untracked) {
            auto work_tracked = msh.get_work_tracked<T>(itype);
            work_tracked[iref_tracked] = this;
        }
    }
    return *this;
}

template<typename T>
WorkArray1D<T>& WorkArray1D<T>::operator=(WorkArray1D&& other) {
    printf("Debug move assignment called\n");
    if (this != &other) {
        itype = other.itype;
        iref_tracked = other.iref_tracked;
        array = std::move(other.array);
        
        // Update tracking pointer
        if(itype != MeshSize::Untracked) {
            auto work_tracked = msh.get_work_tracked<T>(itype);
            work_tracked[iref_tracked] = this;
        }
        
        // Invalidate moved-from object
        other.itype = MeshSize::Untracked;
    }
    return *this;
}

template<typename T>
WorkArray1D<T>::~WorkArray1D(){
  if(itype == MeshSize::Untracked){
    msh.free_work<T>(ilock);
  }else{
    //if(itype == MeshSize::Point && std::is_same_v<T,int>) 
    //  fmt::print("Debug: ~WorkArray1D<int> Point calling free with addr {}\n",(void*) this);
    msh.free_work<T>(ilock,this);
  }
}

template class WorkArray1D<double>;
template class WorkArray1D<int>;

} // namespace Metris