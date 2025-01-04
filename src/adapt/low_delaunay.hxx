//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_DELAUNAY__
#define __METRIS_LOW_DELAUNAY__

#include "../types.hxx"

namespace Metris{

template <int gdim>
bool indelsphere(const double *coop, const double *metl, 
                 const dblAr2 &coord, const int *ent2pol);

}//end namespace

#endif