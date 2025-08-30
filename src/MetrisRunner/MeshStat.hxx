//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MESH_STAT__
#define __METRIS_MESH_STAT__


#include "../utils/aux_MinMaxAvg.hxx"

#include "fmt/format.h"
#include "nlohmann/json_fwd.hpp"
#include "../types_arrays.hxx"
#include "../metris_constants.hxx"


namespace Metris{


class MeshStat{
public:
  MeshStat(){ reset(); }
  // Having an operator= means the compiler generates a copy cstr
  // But this will be deprecated and generates a warning. 
  MeshStat(const MeshStat& rhs);

  void reset();

  void setLength(int tdim, const dblAr1& rlened);
  void setQuality(int tdim, const dblAr1& rqualel, AsDeg asdeg);

  // Fixed dimension:
  dblAr1 pctunit;
  MeshArray1D<MinMaxAvg,int> len, quaP1, quaPk;

  int tdim_max;

  MeshStat& operator=(const MeshStat& rhs);
  bool operator==(const MeshStat& rhs);

  // STREAM_T can be FILE* or std::ostream&
  template<typename STREAM_T = FILE*>
  void print(std::string name = "", STREAM_T logfile = stdout) const;

  friend std::ostream& operator<<(std::ostream& _os, const MeshStat& stat){
    stat.print<std::ostream&>("", _os);
    return _os;   
  }

};


void to_json(nlohmann::json& jj, const MeshStat& stat);
void from_json(const nlohmann::json& jj, MeshStat& stat);




} // End namespace



#endif