//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "WorkArray.hxx"
#include "MeshBase.hxx"

namespace Metris{

template<typename T>
WorkArray1D<T>::WorkArray1D(MeshBase& msh_, int ilock_, MeshArray1D<T>& array_) 
    : ilock(ilock_), array(array_), msh(msh_) {}

template<typename T>
WorkArray1D<T>::~WorkArray1D(){
  msh.free_work<T>(ilock);
}

template class WorkArray1D<double>;
template class WorkArray1D<int>;

} // namespace Metris