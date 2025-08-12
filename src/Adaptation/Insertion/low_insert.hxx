//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_LOW_INSERT__
#define __METRIS_LOW_INSERT__

#include "../../Mesh/MeshFwd.hxx"
#include "../../types.hxx"


namespace Metris{

class MshCavity;
struct CavWrkArrs;
struct EdgeSeed;

enum insedgesuf_Errors {INS2D_NOERR = 0, 
                        INS2D_ERR_INTERPMETBACK1 = 1,
                        INS2D_ERR_INTERPMETBACK2 = 2,
                        INS2D_ERR_EGEVALUATE = 3,
                        INS2D_ERR_INCCAVVAL1 = 4,
                        INS2D_ERR_INCCAVVAL2 = 5,
                        INS2D_ERR_INCCAVVAL3 = 6,
                        INS2D_ERR_SHORTEDG = 7,
                        INS2D_ERR_BDRYNOCORR = 8,
                        INS2D_ERR_INCCAVDEL = 9,
                        INS2D_ERR_CAVITYOPERATOR = 10,
                        INS2D_ERR_MOVEPT = 11,
                        INS2D_ERR_SHORTCSTR = 12,
                        INS2D_ERR_BISECTION = 13,
                        INS2D_ERR_LENQUA = 14,
                        INS2D_ERR_NORDEV = 15,
                        INS2D_ERR_BISECLEN = 16,
                        INS2D_ERR_NOOPERATION = 17,
                        INS2D_ERR_NOOPERATION2 = 18,
                        INS2D_ERR_COLPDIM = 19,
                        INS2D_ERR_COLCORNER = 20,
                        INS2D_ERR_NERROR = 21
                        };


// Collapse edge iedl of triangle iface
// bar1 is t along the edge with 1 if lnoed[iedl][0]
template<class MFT>
int insertEdge(Mesh<MFT>& msh, 
               const EdgeSeed &insertionSeed,
               double lenqua_short_max,
               bool icollapse,
               MshCavity &cav, CavWrkArrs &work, 
               intAr1 &lerro, int ithrd1, int ithrd2);


} // end namespace

#endif