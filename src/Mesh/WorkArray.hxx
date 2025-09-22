//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MESH_WORK_ARRAY__
#define __METRIS_MESH_WORK_ARRAY__

#include "MeshFwd.hxx"
#include "../types_arrays.hxx"


namespace Metris{


// For tag arrays. Only positive values can be used, 0 is "not tagged".
enum class MeshSize;


template<typename T>
class WorkArray1D{
public:
friend class MeshBase;
  WorkArray1D(MeshBase& msh_, int ilock_, MeshArray1D<T>& array_);
  WorkArray1D(MeshBase& msh_, int ilock_, MeshArray1D<T>& array_, MeshSize itype_, int iref_tracked_);
  ~WorkArray1D();

  ALWAYS_INLINE T &operator[](const int &ii){
    return array[ii];
  }
  ALWAYS_INLINE const T &operator[](const int &ii) const {
    return array[ii];
  }

  int size() const {return array.size();} // allocated memory
  int get_n() const {return array.get_n();} // bounds
  void set_n(int nn){array.set_n(nn);}

  void stack(T val){array.stack(val);}
  T pop(){return array.pop();}
  bool allocate(int m){return array.allocate(m);}

  MeshArray1D<T>& get_array(){return array;}
  const MeshArray1D<T>& get_array() const {return array;}


  // Internal use:
  WorkArray1D(WorkArray1D&& other);


private:
  int ilock;
  MeshArray1D<T> array;
  MeshBase& msh;
  MeshSize itype;
  int iref_tracked;
  bool imoved = false;
};


}// namespace Metris
#endif