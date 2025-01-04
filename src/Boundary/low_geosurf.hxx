//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


// Deprecated 

#ifndef __LOW_GEO_SURF__
#define __LOW_GEO_SURF__

#include "../aux_exceptions.hxx"


namespace Metris{



void projptsurf(MeshBase &msh, int ibpoi, double *coop, double tol = 0.1);

// To be used for entities not handled by CAD e.g. interior planes. 
void projptsurf_disc(MeshBase &msh, int ibpoi, double *coop);
//{
//	METRIS_THROW(TODOExcept());
//}

// Get normal at ibpoi using face. 
void bpo2CADnormal(MeshBase &msh, int ibpoi, double *du, double *dv, double *nrmal);






} // End namespace

#endif