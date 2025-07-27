//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_insert.hxx"
#include "aux_insert.hxx"
#include "low_cavqual.hxx"
#include "low_increasecav.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../utils/mprintf.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../msh_structs.hxx"
#include "../low_topo.hxx"
#include "../low_geo/normal.hxx"
#include "../low_geo/measure.hxx"
#include "../io_libmeshb.hxx"
#include "../linalg/det.hxx"
#include "../low_geo/lenedg.hxx"

#include "../msh_checktopo.hxx"

namespace Metris{

// Return 0 if done nothing, 1 if error, -1 if done swap
// bar1 is t along the edge with 1 if lnoed[iedl][0]
template<class MFT>
int insertEdge(Mesh<MFT>& msh, 
               int tdim, int ientt, int iedl, 
               double lenqua_short_max, // maximum quality (error) a new short edge can have
               bool icollapse,
               MshCavity &cav, CavWrkArrs &work, 
               intAr1 &lerro, int ithrd1, int ithrd2){

  int iverb0   = msh.param->iverb;
  int ivdepth0 = msh.param->ivdepth;
  //if(icollapse){
  //  printf("## DEBUG SET MAX PRINTS \n");
  //  writeMesh("debug",msh);
  //  msh.param->iverb  = 5;
  //  msh.param->ivdepth= 5;
  //  if(lerro[19] > 0){
  //    printf("## DEBUG START WITH LERRO[19] = %d \n",lerro[19]);
  //    wait();
  //  }
  //}

  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1   != ithrd2);
  METRIS_ASSERT(!isdeadent(ientt,msh.ent2poi(tdim)));

  const int nedgl = (tdim*(tdim+1))/2;
  METRIS_ASSERT(iedl >= 0 && iedl < nedgl);

  int nentt = msh.nentt(tdim);

  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const auto lnoed = tdim == 1 ? lnoed1 : 
                     tdim == 2 ? lnoed2 : lnoed3;
  METRIS_ASSERT(ientt >= 0 && ientt < nentt && !isdeadent(ientt, ent2poi));

  //if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(), "Implement + tet nelem = "<<msh.nelem)

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

  int mgrow = 100;

  cav.reset();
  cav.lcedg.allocate(10);
  cav.lcfac.allocate(10);
  cav.lctet.allocate(10);
  cav.inewp = 1;

  // work for collrejcav_lenqua
  #ifndef NDEBUG
  static int nwarnprt = 0;
  if(nwarnprt++ < 10) printf("## WARNING REMOVE STATI FROM NOCOMP\n");
  #endif
  static std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp;


  int ierro = 0;
  bool irestart_cav;
  int ip1 = ent2poi(ientt,lnoed[iedl][0]);
  int ip2 = ent2poi(ientt,lnoed[iedl][1]);

  CPRINTF1("-- START insertEdge tdim = %d ientt = %d ied %d = %d %d\n",tdim,ientt,iedl,ip1,ip2);
  // The shell does not need pdim to gather elements: always use
  int iopen;
  shell(msh,ip1,ip2,tdim,ientt,cav.lcedg,cav.lcfac,cav.lctet,&iopen);
  CPRINTF1(" - cavity seed nedge %d nface %d ntetr %d\n",
           cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());
  if(DOPRINTS2()){
    cav.print(msh);
  }

  METRIS_ASSERT(cav.lcedg.get_n() > 0 || cav.lcfac.get_n() > 0 || cav.lctet.get_n() > 0);
  #ifndef NDEBUG
  if(msh.get_tdim() == 2) METRIS_ASSERT(cav.lcfac.get_n() > 0);
  if(msh.get_tdim() == 3) METRIS_ASSERT(cav.lctet.get_n() > 0);
  if(msh.param->dbgfull){
    for(int ientt : cav.lcent(msh.get_tdim())){
      METRIS_ASSERT_MSG(!isdeadent(ientt,msh.ent2poi(msh.get_tdim())),
        "Shell element dead");
    }
  }
  #endif  


  int nced0 = cav.lcedg.get_n();
  int ncfa0 = cav.lcfac.get_n();
  int ncte0 = cav.lctet.get_n();

  int tdimp = -1;
       if(nced0 > 0) tdimp = 1;
  else if(ncfa0 > 0) tdimp = 2;
  else               tdimp = 3;

  int ibins = -1;

  // Create the point, set info for localization 
  int iseed, iref;
  ego obj = NULL;
  double algnd[3];
  if(tdimp == 1){
    int iedge = cav.lcedg[0];
    METRIS_ASSERT(iedge >= 0);
    cav.ipins = msh.newpoitopo(1,iedge);
    ibins = msh.newbpotopo(cav.ipins,1,iedge);
    iseed = iedge;
    iref = msh.edg2ref[iedge];
    if(msh.CAD()) obj  = msh.CAD.cad2edg[iref];
  }else if(tdimp == 2){
    int iface = cav.lcfac[0];
    METRIS_ASSERT(iface >= 0);
    cav.ipins = msh.newpoitopo(2,iface);
    iseed = iface;
    iref  = msh.fac2ref[iface];
    if(msh.isboundary_faces()){
      ibins = msh.newbpotopo(cav.ipins,2,iface);
      if(msh.CAD()) obj = msh.CAD.cad2fac[iref];
    }
  }else{
    cav.ipins = msh.newpoitopo(3,ientt);
    iseed = ientt;
    iref = msh.tet2ref[ientt];
  }
  if(msh.CAD()) METRIS_ASSERT(obj != NULL 
                    || tdimp == 2 && !msh.isboundary_faces() || tdimp == 3);

  // The point is well seeded for ball now
  if(icollapse){
    for(int ii = 0; ii < 2; ii++){
      int ipoin = ent2poi(ientt, lnoed[iedl][ii]);
      ball(msh, ipoin, cav.lcedg, cav.lcfac, cav.lctet, &iopen, true, ithrd1);
    }
  }


  CPRINTF1(" - create ipins %d tdim = %d seed %d ref %d\n",cav.ipins,tdimp,iseed,iref);

  bool imoved_point = false;

  ierro = aux_bisecPointLen(msh, tdim, ientt, iedl, ibins, tdimp, iseed, iref, icollapse, cav);
  if(ierro != 0){
    CPRINTF1(" # Failed aux_bisecPointLen ierro = %d\n",ierro);
    goto cleanup;
  }


  if(icollapse){
    ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
    ierro = collrejcav_lenqua(msh, cav, false, false, false, -1, nocomp, ithrd2);
    if(ierro > 0){
      ierro = INS2D_ERR_SHORTEDG;
      CPRINTF1(" # collrejcav_lenqua rejects cavity, try fix\n");
      CPRINTF1(" # reject cavity\n");

      // Skeptical but keeping it for now to keep collapses unchanged
      ierro = aux_movePointCav(msh, cav, tdimp, iseed, iref, algnd);
      imoved_point = true;

      if(ierro != 0){
        ierro = INS2D_ERR_INTERPMETBACK;
        CPRINTF1(" - Failed to move point in insertEdge\n");
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

  // Effects both insertions and collapses
  if(tdimp == 2 && msh.idim == 3){
    if(rejcavnordev(msh,cav,ibins,ithrd1)) return INS2D_ERR_NORDEV;
  }

  irestart_cav = false;
restart_cavity:
  ierro = 0;
  if(!irestart_cav) irestart_cav = true;

  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd1);
  }}CT_FOR1(ideg);

  //if(DOPRINTS1()){
  //  msh.param->iverb = iverb0;
  //  msh.param->ivdepth = ivdepth0;
  //  printf("## END OF OPERATION after cavity_operator wait\n");
  //  printf("lerro:");
  //  lerro.print();
  //  wait();
  //  if(ierro > 0){
  //    printf("Error %d wait \n",ierro);
  //    wait();
  //  }
  //}

  if(ierro == CAV_ERR_REMPT && !irestart_cav && !imoved_point){
    cav.lcedg.set_n(nced0);
    cav.lcfac.set_n(ncfa0);
    cav.lctet.set_n(ncte0);
    int ierr2 = aux_movePointCav(msh, cav, tdimp, iseed, iref, algnd);
    imoved_point = true;
    writeMeshCavity("insert_cavity_fail_move0.meshb", 
                                    msh,cav);
    if(ierr2 != 0){
      ierro = INS2D_ERR_MOVEPT;
      goto cleanup;
    }
    ierro = increase_cavity_Delaunay(msh, cav, -1, ithrd1);
    if(ierro != 0){
      CPRINTF1(" - +cav error %d\n",ierro);
      ierro = INS2D_ERR_INCCAV2D;
      goto cleanup;
    }
    ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
    if(ierro != 0){
      CPRINTF1(" - +cav error %d\n",ierro);
      ierro = INS2D_ERR_INCCAV2D;
      goto cleanup;
    }

    goto restart_cavity;
  }

  if(ierro > 0) lerro[ierro-1]++;

  if(ierro != 0) ierro = INS2D_ERR_CAVITYOPERATOR;

  if(info.done){

    CPRINTF1("-- END insertEdge ipins = %d  \n",cav.ipins);
    #ifndef NDEBUG
      if(DOPRINTS2()) writeMesh("debug_insert1.meshb",msh);
    #endif
    msh.param->iverb = iverb0;
    msh.param->ivdepth = ivdepth0;
    return -1; // Return did op
  }

  cleanup:
  msh.killpoint(cav.ipins);
  if(DOPRINTS1()){
    msh.param->iverb = iverb0;
    msh.param->ivdepth = ivdepth0;
    //printf("## END OF OPERATION WAIT ierro = %d\n",ierro);
    //printf("lerro:");
    //lerro.print();
    //wait();
  }
  return ierro;
}



template int insertEdge<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                         int tdim, int ientt, int iedl, 
                         double lenqua_short_max, bool icollapse,
                         MshCavity &cav, CavWrkArrs &work, 
                         intAr1 &lerro, int ithrd1, int ithrd2);
template int insertEdge<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                         int tdim, int ientt, int iedl, 
                         double lenqua_short_max, bool icollapse,
                         MshCavity &cav, CavWrkArrs &work, 
                         intAr1 &lerro, int ithrd1, int ithrd2);




} // end namespace
