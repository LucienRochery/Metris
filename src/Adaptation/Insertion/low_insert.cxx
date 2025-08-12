//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_insert.hxx"
#include "aux_insert.hxx"
#include "seed_edge.hxx"

#include "../low_cavqual.hxx"
#include "../low_increasecav.hxx"

#include "../../Mesh/Mesh.hxx"
#include "../../MetrisRunner/MetrisParameters.hxx"

#include "../../utils/mprintf.hxx"
#include "../../cavity/msh_cavity.hxx"
#include "../../aux_topo.hxx"
#include "../../msh_structs.hxx"
#include "../../low_topo.hxx"
#include "../../low_geo/normal.hxx"
#include "../../low_geo/measure.hxx"
#include "../../io_libmeshb.hxx"
#include "../../linalg/det.hxx"
#include "../../low_geo/lenedg.hxx"

#include "../../msh_checktopo.hxx"

namespace Metris{

// Return 0 if done nothing, 1 if error, -1 if done swap
template<class MFT>
int insertEdge(Mesh<MFT>& msh, 
               const EdgeSeed &insertionSeed,
               double lenqua_short_max, // maximum quality (error) a new short edge can have
               bool icollapse,
               MshCavity &cav, CavWrkArrs &work, 
               intAr1 &lerro, int ithrd1, int ithrd2){

  int iverb0   = msh.param->iverb;
  int ivdepth0 = msh.param->ivdepth;
  int ierro_cavity = 0;

  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1   != ithrd2);

  CavOprOpt opts;
  CavOprInfo info;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = false;
  opts.dryrun = false;
  //opts.allow_remove_points = icollapse; 
  opts.allow_remove_points = true; 
  opts.allow_remove_points_superdim = true; // For boundary
  opts.qmax_nec = -1;
  opts.qmax_suf = -1;
  opts.qmax_iff = -1;
  bool idbg = false;
  cav.inewp = 1;

  int mgrow = 100;

  CPRINTF1("-- START insertEdge tdimp = {} iseed = {}\n",
           insertionSeed.tdimp,insertionSeed.iseed);

  cav.print(msh);
  CPRINTF1(" - cavity seed nedge {} nface {} ntetr {}\n",
           cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());
  METRIS_ASSERT(cav.lcedg.get_n() > 0 || cav.lcfac.get_n() > 0 || cav.lctet.get_n() > 0);



  int ierro = 0;
  bool irestart_cav;


  int nced0 = cav.lcedg.get_n();
  int ncfa0 = cav.lcfac.get_n();
  int ncte0 = cav.lctet.get_n();

  double algnd[3];

  // Create the point, set info for localization 
  cav.ipins = msh.newpoitopo(insertionSeed.tdimp, insertionSeed.iseed);
  int ibins = -1;
  if(msh.isboundary_tdim(insertionSeed.tdimp)) 
    ibins = msh.newbpotopo(cav.ipins,insertionSeed.tdimp,insertionSeed.iseed);
  
  if(msh.CAD()) METRIS_ASSERT(insertionSeed.obj != NULL 
                    || insertionSeed.tdimp == 2 && !msh.isboundary_faces() || insertionSeed.tdimp == 3);

  // Append edge ends balls to the shell:
  if(icollapse){
    for(int ii = 0; ii < 2; ii++){
      int ipoin = insertionSeed.ipedg[ii];
      int iopen;
      ball(msh, ipoin, cav.lcedg, cav.lcfac, cav.lctet, &iopen, true, ithrd1);
    }
  }


  // work for collrejcav_lenqua
  #ifndef NDEBUG
  static int nwarnprt = 0;
  if(nwarnprt++ < 10) printf("## WARNING REMOVE STATI FROM NOCOMP\n");
  #endif
  static std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp;


  CPRINTF1(" - create ipins {} tdim = {} seed {} ref {} icollapse {}\n",cav.ipins,insertionSeed.tdimp,insertionSeed.iseed, insertionSeed.iref,icollapse);

  bool imoved_point = false;

  ierro = aux_bisecPointLen(msh, insertionSeed, ibins, icollapse, cav);
  if(ierro != 0){
    CPRINTF1(" # Failed aux_bisecPointLen ierro = {}\n",ierro);
    goto cleanup;
  }
  // Seed the cavity properly
  ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);

  

  if(icollapse){
    ierro = collrejcav_lenqua(msh, cav, false, false, false, -1, nocomp, ithrd2);
    if(ierro > 0){
      ierro = INS2D_ERR_SHORTEDG;
      CPRINTF1(" # collrejcav_lenqua rejects cavity, try fix\n");
      CPRINTF1(" # reject cavity\n");

      // Skeptical but keeping it for now to keep collapses unchanged
      ierro = aux_movePointCav(msh, cav, insertionSeed.tdimp, insertionSeed.iseed, insertionSeed.iref, algnd);
      imoved_point = true;

      if(ierro != 0){
        CPRINTF1(" - Failed to move point in insertEdge\n");
        PRINTF("## DEBUG WAIT HERE\n");
        wait();
        goto cleanup;
      }
      ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
      ierro = collrejcav_lenqua(msh, cav, false, false, false, -1, nocomp, ithrd2);
      if(ierro > 0){
        CPRINTF1(" # collrejcav_lenqua rejects cavity after fix\n");
        goto cleanup;
      }
    }
    goto call_cavity;
  }

  // -- This section only if !icollapse

  ierro = setCavityInsertion(msh,cav,opts,mgrow,lenqua_short_max,nocomp,ithrd1,ithrd2);
  if(ierro != 0) goto cleanup;


call_cavity:

  // nordev is now checked in the cavity
  #if 0
  // Effects both insertions and collapses
  if(tdimp == 2 && msh.idim == 3){
    int iverb0 = msh.param->iverb;
    int ivdepth0 = msh.param->ivdepth;
    //if(msh.fac2ref[cav.lcfac[0]] == 4){
    //  msh.param->iverb = 5;
    //  msh.param->ivdepth = 10;
    //  idbg = true;
    //  printf("## DEBUG SET MAX PRINTS\n");
    //  writeMesh("debug_insert.meshb",msh);
    //}
    bool irej = rejcavnordev(msh,cav,ibins,ithrd1);
    msh.param->iverb = iverb0;
    msh.param->ivdepth = ivdepth0;
    //if(idbg && !irej){
    //  printf("## DEBUG WAIT HERE irej = {}\n",irej);
    //  wait();
    //}
    if(irej){
      ierro = INS2D_ERR_NORDEV;
      goto cleanup;
    }
  }
  #endif

  irestart_cav = false;
restart_cavity:
  ierro = 0;
  if(!irestart_cav) irestart_cav = true;

  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd1);
  }}CT_FOR1(ideg);
  ierro_cavity = ierro;


  //if(DOPRINTS1()){
  //  msh.param->iverb = iverb0;
  //  msh.param->ivdepth = ivdepth0;
  //  printf("## END OF OPERATION after cavity_operator wait\n");
  //  printf("lerro:");
  //  lerro.print();
  //  wait();
  //  if(ierro > 0){
  //    printf("Error {} wait \n",ierro);
  //    wait();
  //  }
  //}

  if(ierro == CAV_ERR_REMPT && !irestart_cav && !imoved_point){
    cav.lcedg.set_n(nced0);
    cav.lcfac.set_n(ncfa0);
    cav.lctet.set_n(ncte0);
    int ierr2 = aux_movePointCav(msh, cav, insertionSeed.tdimp, insertionSeed.iseed, insertionSeed.iref, algnd);
    imoved_point = true;
    writeMeshCavity("insert_cavity_fail_move0.meshb", 
                                    msh,cav);
    if(ierr2 != 0){
      ierro = INS2D_ERR_MOVEPT;
      goto cleanup;
    }
    ierro = increase_cavity_Delaunay(msh, cav, -1, ithrd1);
    if(ierro != 0){
      CPRINTF1(" - +cav error {}\n",ierro);
      ierro = INS2D_ERR_INCCAVDEL;
      goto cleanup;
    }
    ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
    if(ierro != 0){
      CPRINTF1(" - +cav error {}\n",ierro);
      ierro = INS2D_ERR_INCCAVVAL3;
      goto cleanup;
    }

    goto restart_cavity;
  }

  if(ierro > 0) lerro[ierro-1]++;


  if(ierro != 0) ierro = INS2D_ERR_CAVITYOPERATOR;

  if(info.done){

    if(idbg){
      PRINTF("## CAVITY SUCCESSFUL inserted ipoin {} \n",cav.ipins);
      writeMeshCavity("insert_cavity_success.meshb", msh, cav);
      writeMesh("insert_mesh_success.meshb", msh);
      wait();
    }
    CPRINTF1("-- END insertEdge ipins = {}  \n",cav.ipins);
    #ifndef NDEBUG
      if(DOPRINTS2()) writeMesh("debug_insert1.meshb",msh);
    #endif
    msh.param->iverb = iverb0;
    msh.param->ivdepth = ivdepth0;
    return -1; // Return did op
  }

  cleanup:
  msh.killpoint(cav.ipins);
  msh.param->iverb = iverb0;
  msh.param->ivdepth = ivdepth0;
  //if(DOPRINTS1() && ierro_cavity > 0){
  //  printf("## DEBUG IERRO = {} \n",ierro_cavity);
  //  wait();
  //}
  return ierro;
}



template int insertEdge<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                         const EdgeSeed &insertionSeed,
                         double lenqua_short_max, bool icollapse,
                         MshCavity &cav, CavWrkArrs &work, 
                         intAr1 &lerro, int ithrd1, int ithrd2);
template int insertEdge<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                         const EdgeSeed &insertionSeed,
                         double lenqua_short_max, bool icollapse,
                         MshCavity &cav, CavWrkArrs &work, 
                         intAr1 &lerro, int ithrd1, int ithrd2);




} // end namespace
