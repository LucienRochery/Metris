//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_Steiner.hxx"
#include "aux_insert.hxx"
#include "insert_errors.hxx"
#include "EdgeSeed.hxx"

#include "../low_cavqual.hxx"
#include "../low_increasecav.hxx"

#include "../../Mesh/Mesh.hxx"
#include "../../MetrisRunner/MetrisParameters.hxx"

#include "../../utils/mprintf.hxx"
#include "../../utils/fmt_formatters.hxx"

#include "../../Localization/msh_localization.hxx"



#include "../../cavity/msh_cavity.hxx"
#include "../../aux_topo.hxx"
#include "../../msh_structs.hxx"
#include "../../low_topo.hxx"
#include "../../low_geo/normal.hxx"
#include "../../low_geo/measure.hxx"
#include "../../low_geo/lenedg.hxx"
#include "../../low_geo/misc.hxx"
#include "../../io_libmeshb.hxx"
#include "../../linalg/det.hxx"

#include "../../msh_checktopo.hxx"

namespace Metris{

// Return 0 if done nothing, 1 if error, -1 if inserted new point
template<class MFT>
int insertSteiner(Mesh<MFT>& msh, 
                  const EdgeSeed &insertionSeed,
                  MshCavity &cav, CavWrkArrs &work, 
                  intAr1 &lcaverr, int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);
  const int ntry = 5;
  
  const int iseed = insertionSeed.iseed;
  const int tdimp = insertionSeed.tdimp;
  const int iref  = insertionSeed.iref ;

  if(msh.get_tdim() == tdimp) return 0;

  //fmt::print("## DEBUG SET MAX PRINTS\n");
  //msh.param->iverb = 10;
  //msh.param->ivdepth = 15;
  //iverb__ = 10;
  //ivdepth__ = 15;


  CPRINTF1("-- START insertSteiner iseed {} tdimp {} ipedg {} {}\n",
           iseed, tdimp, insertionSeed.ipedg[0], insertionSeed.ipedg[1]);


  int ierro = 0;
  CavOprInfo info;
  CavOprOpt opts;

  static int nwarnprt = 5;
  if(tdimp != 2){
    if(nwarnprt --> 0) fmt::print("## insertSteiner only implemented for tdim = 2");
    return 0;
  }
  

  int ipoi1 = insertionSeed.ipedg[0];
  int ipoi2 = insertionSeed.ipedg[1];


  // Get the seed's normal, as well as its neighbour.
  // Note getnorfacP1 gives outgoing normals.
  double noredg[3] = {0,0,0};
  double norfac[3];

  getnorfacP1(msh.fac2poi[iseed], msh.coord, norfac);
  for(int ii = 0; ii < msh.idim; ii++) noredg[ii] += norfac[ii];

  const int iedl = getedgent(msh, tdimp, iseed, 
                             ipoi1, ipoi2);
  int isee2 = msh.fac2fac(iseed,iedl);
  METRIS_ASSERT(isee2 >= 0);
  getnorfacP1(msh.fac2poi[isee2], msh.coord, norfac);
  for(int ii = 0; ii < msh.idim; ii++) noredg[ii] += norfac[ii];

  if(normalize_vec<3>(noredg)) return INS2D_ERR_STEINER0NORMAL;

  CPRINTF1(" - noredg {} {} {}\n",noredg[0],noredg[1],noredg[2]);

  // Now we extrude along noredg. 
  // For this, let's use the edge end metrics to get an average size
  
  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
  double h1 = getlenedg<3>(noredg, msh.met[ipoi1]);
  METRIS_ASSERT(h1 > 0);
  double h2 = getlenedg<3>(noredg, msh.met[ipoi2]);
  METRIS_ASSERT(h2 > 0);
  double hh = sqrt(h1*h2);

  CPRINTF1(" - using h1 = {:.2e} h2 = {:.2e} hh = {:.2e}\n",h1,h2,hh);

  // Now we have a characteristic metric size for the direction
  // we'll be putting a point in.

  // Get rid of cavity entities of dim above tdimp
  cav.lctet.set_n(0);
  int ielem = msh.fac2tet(iseed, 0);
  METRIS_ASSERT(ielem >= 0);
  METRIS_ASSERT_MSG(msh.fac2tet(iseed,1) < 0, "tdimp = "<<tdimp<<" iseed = "<<iseed<<" fac2tet = "<<msh.fac2tet(iseed,0)<<" "<<msh.fac2tet(iseed,1)); 
  int iref_sup = msh.tet2ref[ielem];
  cav.lctet.stack(ielem);

  // Save off the other cavity entities. This Steiner insertion
  // will not affect them, and the overarching insertion does need
  // them again.
  // Theoretically, the data is still there in the array
  // and we don't even need to save it off. 
  intWrkAr1 lcedg = msh.get_iwork(cav.lcedg.get_n());
  intWrkAr1 lcfac = msh.get_iwork(cav.lcfac.get_n());
  cav.lcedg.copyTo(lcedg.get_array());
  cav.lcfac.copyTo(lcfac.get_array());



  cav.ipins = msh.newpoitopo(tdimp+1, -1);

  double step = 1;
  double *uvsrf = NULL; // will be useful when we implement tdimp = 1
  double *algnd = NULL; // idem
  for(int itry = 1; itry <= ntry; itry++, step /= 2){
    INCVDEPTH(msh.param);


    for(int ii = 0; ii < msh.idim; ii++)
      msh.coord(cav.ipins,ii) = (msh.coord(ipoi1,ii) + msh.coord(ipoi2,ii))/2
                              - step*noredg[ii]/hh;
    CPRINTF2(" - try {}/{} step {:.2e} ipins coord {} loc seed {}\n",
             itry,ntry,step,dblAr1(msh.idim,msh.coord[cav.ipins]),ielem);

    if(DOPRINTS3()){
      int ipdbg = msh.newpoitopo(0, -1);
      msh.newbpotopo(ipdbg, 0, ipdbg);
      for(int ii = 0; ii < msh.idim; ii++)
        msh.coord(ipdbg,ii) = msh.coord(cav.ipins,ii);
      writeMesh("Steiner_locMesh_" + std::to_string(itry),msh);
      msh.met.writeMetricFile("Steiner_locMesh_" + std::to_string(itry));
      msh.killpoint(ipdbg);
    }

    double coopr[3], bary[4];
    ielem = cav.lctet[0];
    MSH_DIM_DEG0(msh){
      if constexpr (gdim == 3)
        ierro = locMesh<gdim,3,ideg>(msh, &ielem, msh.coord[cav.ipins], 
                                     tdimp, uvsrf,
                                     iref_sup, algnd, 
                                     coopr, bary);
      else ierro = 1;
    }MSH_DIM_DEG1();

    //printf("## DEBUG REMOVE THIS\n");
    //check_topo(msh,msh.nbpoi,msh.npoin-1,msh.nedge,msh.nface,msh.nelem,ithrd1);

    
    if(ierro == 0){
      CPRINTF1(" - locMesh succeeded on try {}\n",itry);
      break;
    }

    CPRINTF2(" # locMesh error {} final element {}\n",ierro,ielem);

  }

  if(ierro != 0){
    CPRINTF1("## END insertSteiner: all {} localizations failed\n",ntry);
    ierro = INS2D_ERR_LOCALIZATION;
    goto cleanup;
  }

  if(ielem != cav.lctet[0]) cav.lctet.stack(ielem);
  msh.poi2ent(cav.ipins, 0) = ielem;
  msh.poi2ent(cav.ipins, 1) = 3;

  ierro = msh.interpMetBack(cav.ipins, tdimp+1, ielem, iref_sup, algnd);
  if(ierro != 0){
    CPRINTF1("## END insertSteiner: interpMetBack ierro {}\n", ierro);
    ierro = INS2D_ERR_INTERPMETBACK3;
    goto cleanup;
  }

  //printf("## DEBUG REMOVE THIS\n");
  //check_topo(msh,msh.nbpoi,msh.npoin-1,msh.nedge,msh.nface,msh.nelem,ithrd1);

  ierro = increase_cavity_validity(msh, cav, ithrd1);
  if(ierro != 0){
    CPRINTF1("## END insertSteiner: increase_cavity_validity ierro {}\n", ierro);
    ierro = INS2D_ERR_STEINERINCCAV;
    goto cleanup;
  }
  //printf("## DEBUG REMOVE THIS\n");
  //check_topo(msh,msh.nbpoi,msh.npoin-1,msh.nedge,msh.nface,msh.nelem,ithrd1);

  // Only after increase_cavity_validity do we remove the 
  // boundary cavity. This ensures all elements between
  // original boundary elements and the Steiner point
  // are destroyed. 
  cav.lcedg.set_n(0);
  cav.lcfac.set_n(0);

  
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = false;
  opts.dryrun = false;
  opts.allow_remove_points = true; 
  opts.allow_remove_points_superdim = true; // For boundary
  opts.qmax_nec = -1;
  opts.qmax_suf = -1;
  opts.qmax_iff = -1;
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd1);
  }}CT_FOR1(ideg);

  //printf("## DEBUG REMOVE THIS\n");
  //check_topo(msh,msh.nbpoi,msh.npoin - (ierro > 0),msh.nedge,msh.nface,msh.nelem,ithrd1);

  if(ierro != 0){
    CPRINTF1("## END insertSteiner: increase_cavity_validity ierro {}\n", ierro);
    ierro = INS2D_ERR_STEINERCAVOPR;
    lcaverr[ierro - 1]++;
    goto cleanup;
  }

  if(info.done){
    CPRINTF1(" - cavity operator did operation\n");
    ierro = -1;
  }



  cleanup:
  lcedg.get_array().copyTo(cav.lcedg);
  lcfac.get_array().copyTo(cav.lcfac);
  if(ierro > 0) msh.killpoint(cav.ipins);
  cav.lctet.set_n(0);


  return ierro;
}



template int insertSteiner<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                         const EdgeSeed &insertionSeed,
                         MshCavity &cav, CavWrkArrs &work, 
                         intAr1 &lerro, int ithrd1, int ithrd2);
template int insertSteiner<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                         const EdgeSeed &insertionSeed,
                         MshCavity &cav, CavWrkArrs &work, 
                         intAr1 &lerro, int ithrd1, int ithrd2);




} // end namespace
