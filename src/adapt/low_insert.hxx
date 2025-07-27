//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_LOW_INSERT__
#define __METRIS_LOW_INSERT__

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
                        INS2D_ERR_CAVITYOPERATOR = 8,
                        INS2D_ERR_MOVEPT = 9,
                        INS2D_ERR_SHORTCSTR = 10,
                        INS2D_ERR_BISECTION = 11,
                        INS2D_ERR_LENQUA = 12,
                        INS2D_ERR_NORDEV = 13,
                        INS2D_ERR_NERROR = 14
                        };


// Collapse edge iedl of triangle iface
// bar1 is t along the edge with 1 if lnoed[iedl][0]
template<class MFT>
int insertEdge(Mesh<MFT>& msh, int tdim, int ientt, int iedl, 
               double lenqua_short_max,
               bool icollapse,
               MshCavity &cav, CavWrkArrs &work, 
               intAr1 &lerro, int ithrd1, int ithrd2);


} // end namespace

#endif