//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_LINALG_UTILS__
#define __METRIS_LINALG_UTILS__

#include "../SANS/Surreal/SurrealS.h"


namespace Metris{

template<int gdim, int nvar>
void getmet_dbl2SurS(const double*__restrict__ met, const double*__restrict__ dmet, 
							       SANS::SurrealS<nvar,double>*__restrict__ metS);

template<int gdim, int nvar>
void getmet_SurS2dbl(const SANS::SurrealS<nvar,double>*__restrict__ metS,
										 double*__restrict__ met, double*__restrict__ dmet);

}// End namsepace



#endif