//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "EGADSprinterr.hxx"
#include <cstdio>

namespace Metris{
std::string EG_err2str(int ierro){
  if(ierro == -37)
    return "EGADS_EXTRAPOL";
  if(ierro == -36)
    return "EGADS_EFFCTOBJ";
  if(ierro == -35)
    return "EGADS_UVMAP";
  if(ierro == -34)
    return "EGADS_SEQUERR";
  if(ierro == -33)
    return "EGADS_CNTXTHRD";
  if(ierro == -32)
    return "EGADS_READERR";
  if(ierro == -31)
    return "EGADS_TESSTATE";
  if(ierro == -30)
    return "EGADS_EXISTS";
  if(ierro == -29)
    return "EGADS_ATTRERR";
  if(ierro == -28)
    return "EGADS_TOPOCNT";
  if(ierro == -27)
    return "EGADS_OCSEGFLT";
  if(ierro == -26)
    return "EGADS_BADSCALE";
  if(ierro == -25)
    return "EGADS_NOTORTHO";
  if(ierro == -24)
    return "EGADS_DEGEN";
  if(ierro == -23)
    return "EGADS_CONSTERR";
  if(ierro == -22)
    return "EGADS_TOPOERR";
  if(ierro == -21)
    return "EGADS_GEOMERR";
  if(ierro == -20)
    return "EGADS_NOTBODY";
  if(ierro == -19)
    return "EGADS_WRITERR";
  if(ierro == -18)
    return "EGADS_NOTMODEL";
  if(ierro == -17)
    return "EGADS_NOLOAD";
  if(ierro == -16)
    return "EGADS_RANGERR";
  if(ierro == -15)
    return "EGADS_NOTGEOM";
  if(ierro == -14)
    return "EGADS_NOTTESS";
  if(ierro == -13)
    return "EGADS_EMPTY";
  if(ierro == -12)
    return "EGADS_NOTTOPO";
  if(ierro == -11)
    return "EGADS_REFERCE";
  if(ierro == -10)
    return "EGADS_NOTXFORM";
  if(ierro ==  -9)
    return "EGADS_NOTCNTX";
  if(ierro ==  -8)
    return "EGADS_MIXCNTX";
  if(ierro ==  -7)
    return "EGADS_NODATA";
  if(ierro ==  -6)
    return "EGADS_NONAME";
  if(ierro ==  -5)
    return "EGADS_INDEXERR";
  if(ierro ==  -4)
    return "EGADS_MALLOC";
  if(ierro ==  -3)
    return "EGADS_NOTOBJ";
  if(ierro ==  -2)
    return "EGADS_NULLOBJ";
  if(ierro ==  -1)
    return "EGADS_NOTFOUND";
  if(ierro ==   0)
    return "EGADS_SUCCESS";
  if(ierro ==   1)
    return "EGADS_OUTSIDE";
  return "EGADS_UNKNOWN_ERROR";
}
}// end namespace