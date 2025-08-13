//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "EdgeSeed.hxx"

#include "../../Mesh/MeshBase.hxx"
#include "../../cavity/msh_cavity.hxx"
#include "../../metris_constants.hxx"

#include "../../low_topo.hxx"
#include "../../utils/mprintf.hxx"

namespace Metris{

class MshCavity;

EdgeSeed::EdgeSeed(MeshBase& msh, MshCavity& cav_, int tdim_adp_, int tdim_ent, int ientt, int iedl) : tdim_adp(tdim_adp_), cav(cav_){
  
  obj = NULL;
  tdimp = -1;
  iseed = -1;
  iref  = -1;
  
  cav.reset();
  cav.lcedg.allocate(1);
  cav.lcfac.allocate(10);
  if(msh.get_tdim() >= 3) cav.lctet.allocate(10);

  const auto lnoed = tdim_ent == 1 ? lnoed1 : 
                     tdim_ent == 2 ? lnoed2 : lnoed3;

  ipedg[0] = msh.ent2poi(tdim_ent)(ientt,lnoed[iedl][0]);
  ipedg[1] = msh.ent2poi(tdim_ent)(ientt,lnoed[iedl][1]);

  int iopen;
  shell(msh,ipedg[0],ipedg[1],tdim_ent,ientt,cav.lcedg,cav.lcfac,cav.lctet,&iopen);


  tdimp = -1;
       if(cav.lcedg.get_n() > 0) tdimp = 1;
  else if(cav.lcfac.get_n() > 0) tdimp = 2;
  else                           tdimp = 3;

  if(tdimp == 1){
    int iedge = cav.lcedg[0];
    METRIS_ASSERT(iedge >= 0);
    iseed = iedge;
    iref = msh.edg2ref[iedge];
    if(msh.isboundary_edges() && msh.CAD()) 
      obj = msh.CAD.cad2edg[iref];
  }else if(tdimp == 2){
    int iface = cav.lcfac[0];
    METRIS_ASSERT(iface >= 0);
    iseed = iface;
    iref  = msh.fac2ref[iface];
    if(msh.isboundary_faces() && msh.CAD())
      obj = msh.CAD.cad2fac[iref];
  }else{
    iseed = ientt;
    iref = msh.tet2ref[ientt];
  }

  #ifndef NDEBUG
  bool ithrow = !(iseed >= 0 && iseed < msh.nentt(tdimp))
             ||  isdeadent(iseed,msh.ent2poi(tdimp))
             ||  iref < 0;
  if(ithrow){
    GETVDEPTH(msh.param);
    PRINTF("## EdgeSeed failed: tdimp = {} iseed = {} iref = {}\n",
           tdimp,iseed,iref);
    PRINTF("Called with tdim = {} ientt = {} iedl = {}\n",
           tdim_ent,ientt,iedl);
    METRIS_THROW(TopoExcept());
  }
  #endif
}



}// namespace Metris
