//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_AUX_EGADSPRINTERR__
#define __METRIS_AUX_EGADSPRINTERR__

#include <string>

namespace Metris{

  void print_EGADS_error(std::string fname, int ierro);

}// end namespace
#endif