//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php
#ifndef __METRIS_AUX_HISTOGRAM__
#define __METRIS_AUX_HISTOGRAM__

#include "types.hxx"

namespace Metris{
class MeshBase;

enum class IntrpTyp{Linear, Geometric};
void print_histogram(const MeshBase &msh, dblAr1 &values, IntrpTyp iinter, dblAr1 &bounds,
                     std::string symb, std::string name);

}//end namespace
#endif