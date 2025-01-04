//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php
#ifndef __METRIS_MSH_LENEDG__
#define __METRIS_MSH_LENEDG__


#include "types.hxx"


namespace Metris{

template<class MFT> class MeshMetric;

enum class LenTyp {GeoSiz, Quad};

// ilned returns edge vertices
// rlned its length in the metric field 
template<class MFT>
double getLengthEdges(MeshMetric<MFT> &msh, intAr2 &ilned, dblAr1 &rlned, 
                      LenTyp itype = LenTyp::GeoSiz);

// internal 
template<class MFT, int gdim, int ideg>
double getLengthEdges0(MeshMetric<MFT> &msh, intAr2 &ilned, dblAr1 &rlned, 
                       LenTyp itype = LenTyp::GeoSiz);


// ilned returns edge vertices
// rlned its length in the metric field 
template<class MFT>
double getLengthEdges_Bdry(MeshMetric<MFT> &msh, intAr2 &ilned, dblAr1 &rlned, 
                           LenTyp itype = LenTyp::GeoSiz);

// internal 
template<class MFT, int gdim, int ideg>
double getLengthEdges_Bdry0(MeshMetric<MFT> &msh, intAr2 &ilned, dblAr1 &rlned, 
                            LenTyp itype = LenTyp::GeoSiz);

}// end namespace
#endif