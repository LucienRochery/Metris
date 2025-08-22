//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "msh_anamet.hxx"
#include "MetrisRunner/MetrisParameters.hxx"

namespace Metris{

AnaMetCtx::AnaMetCtx(const MetrisParameters& param){
  dx = param.anamet_dx;
  dy = param.anamet_dy;
  dz = param.anamet_dz;
}

}// namespace Metris
