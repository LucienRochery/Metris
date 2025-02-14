//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_MSH_COLLAPSE__
#define __METRIS_MSH_COLLAPSE__

#include "../Mesh/MeshFwd.hxx"
#include "../types.hxx"


namespace Metris{

class MshCavity;
struct CavWrkArrs;

enum insedgesuf_Errors {INS2D_NOERR = 0, 
                        INS2D_ERR_INTERPMETBACK = 1,
                        INS2D_ERR_EGEVALUATE = 2,
                        INS2D_ERR_INCCAV2D = 3,
                        INS2D_ERR_INCCAV2D2 = 4,
                        INS2D_ERR_SHORTEDG = 5,
                        INS2D_ERR_BDRYNOCORR = 6,
                        INS2D_ERR_INCCAVDEL = 7,
                        INS2D_ERR_NERROR = 8
                        };


// Collapse edge iedl of triangle iface
// bar1 is t along the edge with 1 if lnoed[iedl][0]
template<class MetricFieldType>
int insedgesurf(Mesh<MetricFieldType>& msh, int iface, int iedl, 
               double* coop, double bar1, 
               MshCavity &cav, CavWrkArrs &work, 
               intAr1 &lerro, int ithrd1 = 0, int ithrd2 = 1);


// Collapse edge iedl of triangle iface
// bar1 is t along the edge with 1 if lnoed[iedl][0]
template<class MetricFieldType>
int insfacsurf(Mesh<MetricFieldType>& msh, int iface, double* coop, 
               MshCavity &cav, CavWrkArrs &work, 
               intAr1 &lerro, int ithrd1 = 0, int ithrd2 = 1);



} // end namespace

#endif