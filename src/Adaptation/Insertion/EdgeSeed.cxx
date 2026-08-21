//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "EdgeSeed.hxx"
#include "insert_errors.hxx"

#include "../../Mesh/MeshBase.hxx"
#include "../../utils/aux_misc.hxx"
#include "../../cavity/msh_cavity.hxx"
#include "../../metris_constants.hxx"

#include "../../low_topo.hxx"
#include "../../utils/mprintf.hxx"

namespace Metris{

class MshCavity;

EdgeSeed::EdgeSeed(MeshBase& msh, MshCavity& cav_, int tdim_adp_, int tdim_ent, int ientt, int iedl, bool icollapse, int ithrd) : tdim_adp(tdim_adp_), cav(cav_) {
  
  ierro = 0;
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

  // Dimension of the edge proper, as gathered from the shell. Necessary to ensure collapse compatibility.
  int shell_dim = -1;
  if (cav.lcedg.get_n() > 0)
      shell_dim = 1;
  else if (cav.lcfac.get_n() > 0)
      shell_dim = 2;
  else
      shell_dim = 3;

  if (icollapse) {
      int pdim0 = msh.getpoitdim(ipedg[0]);
      int pdim1 = msh.getpoitdim(ipedg[1]);
      // Collapse compatibility. If any of the boundary points is lower dimension than the edge,
      // then the points are either incompatible OR (if only one is boundary), we need to do this 
      // as a VERTEX collapse of the higher dimensional point.
      if (shell_dim > METRIS_MIN(pdim0, pdim1)) {
          ierro = INS2D_ERR_COLEDGDIM;
          return;
      }
      // Seed collapse cavity.
      for (int ii = 0; ii < 2; ii++)
          ball(msh, ipedg[ii], cav.lcedg, cav.lcfac, cav.lctet, &iopen, true, ithrd);
  }

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
  if(ithrow) METRIS_THROW_MSG(
              "## EdgeSeed failed: tdimp = {} iseed = {} iref = {}\n"
              "Called with tdim = {} ientt = {} iedl = {}\n",
              tdimp,iseed,iref,tdim_ent,ientt,iedl);
#endif
}



}// namespace Metris
