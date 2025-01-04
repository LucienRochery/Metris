//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "aux_EGADSprinterr.hxx"
#include <cstdio>

namespace Metris{
void print_EGADS_error(std::string fname, int ierro){
  printf("## %s error (%d): ",fname.c_str(),ierro);
  if(ierro == -37)
    printf("EGADS_EXTRAPOL");
  if(ierro == -36)
    printf("EGADS_EFFCTOBJ");
  if(ierro == -35)
    printf("EGADS_UVMAP");
  if(ierro == -34)
    printf("EGADS_SEQUERR");
  if(ierro == -33)
    printf("EGADS_CNTXTHRD");
  if(ierro == -32)
    printf("EGADS_READERR");
  if(ierro == -31)
    printf("EGADS_TESSTATE");
  if(ierro == -30)
    printf("EGADS_EXISTS");
  if(ierro == -29)
    printf("EGADS_ATTRERR");
  if(ierro == -28)
    printf("EGADS_TOPOCNT");
  if(ierro == -27)
    printf("EGADS_OCSEGFLT");
  if(ierro == -26)
    printf("EGADS_BADSCALE");
  if(ierro == -25)
    printf("EGADS_NOTORTHO");
  if(ierro == -24)
    printf("EGADS_DEGEN");
  if(ierro == -23)
    printf("EGADS_CONSTERR");
  if(ierro == -22)
    printf("EGADS_TOPOERR");
  if(ierro == -21)
    printf("EGADS_GEOMERR");
  if(ierro == -20)
    printf("EGADS_NOTBODY");
  if(ierro == -19)
    printf("EGADS_WRITERR");
  if(ierro == -18)
    printf("EGADS_NOTMODEL");
  if(ierro == -17)
    printf("EGADS_NOLOAD");
  if(ierro == -16)
    printf("EGADS_RANGERR");
  if(ierro == -15)
    printf("EGADS_NOTGEOM");
  if(ierro == -14)
    printf("EGADS_NOTTESS");
  if(ierro == -13)
    printf("EGADS_EMPTY");
  if(ierro == -12)
    printf("EGADS_NOTTOPO");
  if(ierro == -11)
    printf("EGADS_REFERCE");
  if(ierro == -10)
    printf("EGADS_NOTXFORM");
  if(ierro ==  -9)
    printf("EGADS_NOTCNTX");
  if(ierro ==  -8)
    printf("EGADS_MIXCNTX");
  if(ierro ==  -7)
    printf("EGADS_NODATA");
  if(ierro ==  -6)
    printf("EGADS_NONAME");
  if(ierro ==  -5)
    printf("EGADS_INDEXERR");
  if(ierro ==  -4)
    printf("EGADS_MALLOC");
  if(ierro ==  -3)
    printf("EGADS_NOTOBJ");
  if(ierro ==  -2)
    printf("EGADS_NULLOBJ");
  if(ierro ==  -1)
    printf("EGADS_NOTFOUND");
  if(ierro ==   0)
    printf("EGADS_SUCCESS");
  if(ierro ==   1)
    printf("EGADS_OUTSIDE");
  printf("\n");
}
}// end namespace