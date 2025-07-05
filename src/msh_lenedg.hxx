//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php
#ifndef __METRIS_MSH_LENEDG__
#define __METRIS_MSH_LENEDG__


#include "types.hxx"


namespace Metris{

template<class MFT> class MeshMetric;

enum class LenTyp {GeoSiz, Quad, LogIntrp, BdryCor, MetCrv};

struct lenStat{
  double prop_unit;// proportion unit (0 to 1)
  // Qualities are between 0 (exactly unit) and 1
  // Quality of 1/sqrt(2) = qual of sqrt(2). 
  double qua_short;// 1 - len
  double qua_long; // 1 - 1 / len
  double qua_glo;  // the maximum of the two, for convenience.
};

// ilned returns edge vertices
// rlned its length in the metric field 
template<class MFT>
void getLengthEdges(MeshMetric<MFT> &msh, int tdim, int iref, intAr2 &ilned, dblAr1 &rlned, 
                    lenStat& stat, LenTyp itype = LenTyp::GeoSiz);


}// end namespace
#endif