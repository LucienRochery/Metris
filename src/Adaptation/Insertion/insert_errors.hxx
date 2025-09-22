//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_INSERT_ERRORS__
#define __METRIS_INSERT_ERRORS__

namespace Metris{
enum insedgesuf_Errors {INS2D_NOERR = 0, 
                        INS2D_ERR_INTERPMETBACK1 = 1,
                        INS2D_ERR_INTERPMETBACK2 = 2,
                        INS2D_ERR_EGEVALUATE = 3,
                        INS2D_ERR_INCCAVVAL1 = 4,
                        INS2D_ERR_INCCAVVAL2 = 5,
                        INS2D_ERR_INCCAVVAL3 = 6,
                        INS2D_ERR_SHORTEDG1 = 7,
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
                        INS2D_ERR_LOCALIZATION = 21,
                        INS2D_ERR_STEINER0NORMAL = 22,
                        INS2D_ERR_STEINERINCCAV = 23,
                        INS2D_ERR_STEINERCAVOPR = 24,
                        INS2D_ERR_INTERPMETBACK3 = 25,
                        INS2D_ERR_MOVPTCAVLEN = 26,
                        INS2D_ERR_SHORTCSTR2 = 27,
                        INS2D_ERR_SHORTEDG2 = 28,
                        INS2D_ERR_SHORTEDG3 = 29,
                        INS2D_ERR_SHORTEDG4 = 30,
                        INS2D_ERR_SHORTEDG5 = 31,
                        INS2D_ERR_SHORTEDG6 = 32,
                        INS2D_ERR_INCCAVLEN = 33,
                        INS2D_ERR_COLREF = 34,
                        INS2D_ERR_NERROR = 35
                        };

} // end namespace Metris

#endif