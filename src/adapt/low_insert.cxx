//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_insert.hxx"
#include "low_cavqual.hxx"
#include "low_increasecav.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../utils/mprintf.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../msh_structs.hxx"
#include "../low_topo.hxx"
#include "../low_normal.hxx"
#include "../low_geo.hxx"
#include "../io_libmeshb.hxx"
#include "../linalg/det.hxx"
#include "../low_lenedg.hxx"

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
  //  //wait();
  //  msh.param->iverb  = 5;
  //  msh.param->ivdepth= 5;
  //}

  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
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
  opts.allow_remove_points = icollapse; 
  opts.allow_remove_points_superdim = true; // For boundary
  opts.qmax_nec = -1;
  opts.qmax_suf = -1;
  opts.qmax_iff = -1;

  int mgrow = 100;
  int ngrow = 0;
  intWrkAr1 lrempoi = msh.get_iwork(10);

  cav.reset();
  cav.lcedg.allocate(10);
  cav.lcfac.allocate(10);
  cav.lctet.allocate(10);


  // work for collrejcav_lenqua
  #ifndef NDEBUG
  static int nwarnprt = 0;
  if(nwarnprt++ < 10) printf("## WARNING REMOVE STATI FROM NOCOMP\n");
  #endif
  static std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp;


  int ierro = 0, nprem;
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

  double bar1_min = 1.0e-6, bar1_max = 1 - 1.0e-6;

  constexpr int mnode = getnnode(1,METRIS_MAX_DEG);
  const int nnode = getnnode(1,msh.curdeg);
  int edg2pol[mnode], edg2po2[2] = {cav.ipins, -1};
  int idx[4] = {0};
  int idx1[2];
  for(int inoed = 0; inoed <= msh.curdeg; inoed++){
    idx[lnoed[iedl][0]] = msh.curdeg - inoed;
    idx[lnoed[iedl][1]] = inoed;
    idx1[0] = msh.curdeg - inoed;
    idx1[1] = inoed;
    int inode_sup = mul2nod(tdim,idx);
    int inode_sub = mul2nod(1,idx1);
    edg2pol[inode_sub] = ent2poi(ientt,inode_sup);
  }
  CPRINTF2(" - edg2pol = ");
  if(DOPRINTS2()) intAr1(nnode,edg2pol).print();


  // Get bar1 s.t. new edges are not short. There can be other short edges, but 
  // not from splitting the parent edge. 
  bool fnd_len = false;
  double bar1_opt = -1, err_opt = 1.0e30, bar1;
  for(int ntry_len = 0; ntry_len < 10; ntry_len++){
    INCVDEPTH(msh.param);
    bar1 = (bar1_min + bar1_max) / 2;
    double bar2[2] = {bar1, 1 - bar1};

    // Evaluate ipins on CAD or element, also get algnd for interpMetBack 
    if(ibins >= 0 && msh.CAD()){
      int ib[2];
      // Correct ibs : attach to ref or edge/face as needed
      for(int ii = 0; ii < 2; ii++){
        ib[ii] = msh.poi2ebp(edg2pol[ii],tdimp,iseed,iref);
        METRIS_ASSERT(ib[ii] >= 0);
      }

      for(int ii = 0; ii < 2; ii++) msh.bpo2rbi(ibins,ii) = 
          bar1*msh.bpo2rbi[ib[0]][ii] + (1.0 - bar1)*msh.bpo2rbi[ib[1]][ii];

      CPRINTF1(" - boundary point new t/(u,v) = %f %f\n",
               msh.bpo2rbi(ibins,0),msh.bpo2rbi(ibins,1));

      double result[18];
      METRIS_ASSERT(obj != NULL);
      ierro = EG_evaluate(obj, msh.bpo2rbi[ibins], result);
      if(ierro != 0){
        ierro = INS2D_ERR_EGEVALUATE; 
        goto cleanup;
      }
      for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins,ii) = result[ii];

      if(tdimp == 1){
        for(int ii = 0; ii < msh.idim; ii++) algnd[ii] = result[3+ii];
      }else{
        vecprod(&result[3], &result[6], algnd);
      }
    }else if(ibins >= 0 && !msh.CAD()){ 
      METRIS_ASSERT(tdimp <= 2);
      // No reevaluation, but initialize algnd to edge tangent 
      CPRINTF1(" - discrete algnd initialization tdimp %d \n",tdimp);
      if(tdimp == 1){
        // To compute at higher degree, copy more vertices into ip
        MSH_DIM_DEG0(msh){
          eval1<gdim,ideg>(msh.coord, edg2pol,
                           msh.getBasis(), DifVar::Bary, DifVar::None,
                           bar2, msh.coord[cav.ipins], algnd, NULL);
        }MSH_DIM_DEG1();
      }else{
        MSH_DIM_DEG0(msh){
          eval1<gdim,ideg>(msh.coord, edg2pol,
                           msh.getBasis(), DifVar::None, DifVar::None,
                           bar2, msh.coord[cav.ipins], NULL, NULL);
        }MSH_DIM_DEG1();
        if(msh.idim == 3) getnorfacP1(ent2poi[ientt],msh.coord,algnd);
      }
    }else{
      MSH_DIM_DEG0(msh){
        eval1<gdim,ideg>(msh.coord, edg2pol,
                         msh.getBasis(), DifVar::None, DifVar::None,
                         bar2, msh.coord[cav.ipins], NULL, NULL);
      }MSH_DIM_DEG1();
    }

    ierro = msh.interpMetBack(cav.ipins, tdimp, iseed, iref, algnd);
    if(ierro != 0){
      ierro = INS2D_ERR_INTERPMETBACK;
      goto cleanup;
    }

    //if(DOPRINTS3()){
    //  int ipnew = msh.newpoitopo(0);
    //  msh.newbpotopo(ipnew, 0, ipnew);
    //  const int nnmet = (msh.idim*(msh.idim+1))/2;
    //  for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipnew, ii) = msh.coord(cav.ipins, ii);
    //  for(int ii = 0; ii < nnmet; ii++) msh.met(ipnew, ii) = msh.met(cav.ipins, ii);
    //}


    double sz[2];
    edg2po2[1] = ip1;
    double len1 = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,edg2po2,sz)
                                : getlenedg_geosz<MFT,3,1>(msh,edg2po2,sz);
    edg2po2[1] = ip2;
    double len2 = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,edg2po2,sz)
                                : getlenedg_geosz<MFT,3,1>(msh,edg2po2,sz);

    CPRINTF1(" - %d bar1 = %e lens = %e %e valid %d %d dist %e %e sumlen %e\n",
              ntry_len,bar1,len1,len2,
              len1 > 1/sqrt(2), len2 > 1/sqrt(2),
              abs(len1-1), abs(len2-1),
              len1+len2);

    if(len1 > len2){
      // make len1 shorter by increasing bar1 (pulling ipins towards ip1)
      bar1_min = bar1;
    }else{
      bar1_max = bar1;
    }

    bool ivalid = !icollapse ? len1 >= 1/sqrt(2) && len2 >= 1/sqrt(2) : true;

    if(ivalid){
      fnd_len = true;
      // Once a viable is found, it is possible length distance to 1 will 
      // make a couple of iterates not viable, so it's important to keep a viable
      // bar1. 
      double err = abs(abs(1-len1) - abs(1-len2));
      if(err < err_opt){
        err_opt = err;
        bar1_opt = bar1;
      }
      CPRINTF1(" - config error %e \n",err);
      if(err < 1.0e-2) break;
    }
  }// for ntry_len

  //if(DOPRINTS3()){
  //  writeMesh("debug_bisection",msh);
  //  msh.met.writeMetricFile("debug_bisection");
  //  for(int ii = npoi0; ii < msh.npoin; ii++){
  //    msh.killpoint(ii);
  //  }
  //  msh.set_npoin(npoi0);
  //  msh.set_nbpoi(nbpo0);
  //}

  if(!fnd_len){
    ierro = INS2D_ERR_SHORTEDG;
    goto cleanup;
  }
  bar1 = bar1_opt;
  CPRINTF1(" - end bisection using bar1 = %f\n",bar1);


  if(icollapse){
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
      ierro = collrejcav_lenqua(msh, cav, false, false, false, -1, nocomp, ithrd2);
      if(ierro > 0){
        CPRINTF1(" # collrejcav_lenqua rejects cavity after fix\n");
        goto cleanup;
      }
    }
    goto call_cavity;
  }

  // This section only if !icollapse
  do{

    // If need to revert elements
    int nced1 = cav.lcedg.get_n();
    int ncfa1 = cav.lcfac.get_n();
    int ncte1 = cav.lctet.get_n();
    ierro = 0;

    if(DOPRINTS2()){
      std::string fname = "insert_cavity0."+std::to_string(ngrow);
      writeMeshCavity(fname,msh,cav);
      msh.met.writeMetricFile(fname);
    }
    CPRINTF1(" - step %d cavity nedge %d nface %d nelem %d\n",ngrow,
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());

    //if(opts.allow_remove_points_superdim && tdimp < msh.get_tdim()){
    //  // nprem < 0 is a rejection, but we don't care anymore as we use 
    //  // collrejcav_lenqua.
    //  nprem = increase_cavity_lenedg(msh,cav,opts,cav.ipins,ithrd1,ithrd2);
    //  CPRINTF1(" - +len cavity size %d nprem = %d\n", cav.lcfac.get_n(),nprem); 
    //  if(DOPRINTS2()){
    //    writeMeshCavity("insert_cavity1."+std::to_string(ngrow), 
    //                                msh,cav);
    //  }
    //}

    // -- 1 step Delaunay increase
    ierro = increase_cavity_Delaunay(msh, cav, 1, ithrd1);
    if(ierro != 0){
      CPRINTF1(" # +del error %d\n",ierro);
      ierro = INS2D_ERR_INCCAV2D;
      goto finish_grow_step;
    }
    CPRINTF1(" - +del nedge %d nface %d nelem %d\n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());
    if(DOPRINTS2())
      writeMeshCavity("insert_cavity1."+std::to_string(ngrow), 
                                  msh,cav);


    // -- increase for validity
    ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
    if(ierro != 0){
      CPRINTF1(" # +cav error %d\n",ierro);
      ierro = INS2D_ERR_INCCAV2D;
      goto finish_grow_step;
    }
    CPRINTF1(" - +cav nedge %d nface %d nelem %d\n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());
    if(DOPRINTS2())
      writeMeshCavity("insert_cavity2."+std::to_string(ngrow), 
                                  msh,cav);
 
    //ierro = collrejcav_lenqua(msh, cav, true, false, true, lenqua_short_max, nocomp, ithrd2);
    //if(ierro > 0){
    //  ierro = INS2D_ERR_SHORTEDG;
    //  CPRINTF1(" # collrejcav_lenqua rejects cavity, try fix\n");
    //  CPRINTF1(" # reject cavity\n");
    //  goto finish_grow_step;
    //}

    // Check if the cavity needs fixing.
    // This is only if points are going to be removed, and they have length to 
    // ipins too short. 

    check_cavity_rempoint(msh, cav, opts, lrempoi.get_array(), true, ithrd1);
    if(lrempoi.get_n() > 0){
      ierro = INS2D_ERR_SHORTEDG;
      CPRINTF1(" # error nrem point = %d\n",lrempoi.get_n());
      goto finish_grow_step;
    }

    finish_grow_step:
    if(ierro > 0){
      if(lrempoi.get_n() == 0){
        CPRINTF1(" # Unfixable cavity\n");
        ierro = INS2D_ERR_MOVEPT;
        goto cleanup;
      }

      // Now we need to remove all the newly added elements that contain 
      // one of the lrempoi.
      msh.tag[ithrd1]++;
      for(int ii = 0; ii < lrempoi.get_n(); ii++){
        int ipoin = lrempoi[ii];
        msh.poi2tag(ithrd1, ipoin) = msh.tag[ithrd1];
      }
      for(int tdimc = 1; tdimc <= msh.get_tdim(); tdimc++){
        intAr1 &lcent = cav.lcent(tdimc);
        const int ncen0 = tdimc == 1 ? nced1 : 
                          tdimc == 2 ? ncfa1 : ncte1;
        const int ncent = lcent.get_n();
        const intAr2& ent2poc = msh.ent2poi(tdimc);
        int nrem = 0;
        for(int ii = ncen0; ii < lcent.get_n();){
          INCVDEPTH(msh.param);
          int icent = lcent[ii];
          bool remelt = false;
          for(int iver = 0; iver < tdimc + 1; iver++){
            int ipoin = ent2poc(icent,iver);
            if(msh.poi2tag(ithrd1, ipoin) < msh.tag[ithrd1]) continue;
            remelt = true;
            break;
          }// for iver
          if(!remelt){
            ii++;
            continue;
          }
          CPRINTF1(" - remove %d from cavity dim %d\n",icent,tdimc);
          int icend = lcent.pop();
          // This can only happen if we're the last element. In that case we 
          // shrank the array and can quit. 
          if(icend == icent) break;
          // otherwise place last here.
          icent = icend;
          nrem++;
        }// for icent
        CPRINTF1(" - removed %d dim %d cavity elements\n",nrem,tdimc);
      }// for tdimc
    }// if ierro > 0
    
    ierro = 0;

    // Make sure not shrinking (would be a bug)
    METRIS_ASSERT(cav.lcedg.get_n() >= nced1);
    METRIS_ASSERT(cav.lcfac.get_n() >= ncfa1);
    METRIS_ASSERT(cav.lctet.get_n() >= ncte1);

    // Check if the cavity has grown; break if not
    bool igrow =  cav.lcedg.get_n() > nced1 
               || cav.lcfac.get_n() > ncfa1 
               || cav.lctet.get_n() > ncte1;
    if(!igrow) break;

  }while(ierro > 0 && ngrow++ <= mgrow);

  if(ierro > 0) goto cleanup;


call_cavity:

  irestart_cav = false;
restart_cavity:
  ierro = 0;
  if(!irestart_cav) irestart_cav = true;

  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd1);
  }}CT_FOR1(ideg);

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

  #if 0
  if(ierro == CAV_ERR_REMPT && irestart_cav && ncfa0 == 0){
    printf("## DEBUG DESPITE RESTART CAVITY DOES NOT WORK\n");
    writeMeshCavity("insert_cavity_fail1.meshb", 
                                    msh,cav);
    cav.lcedg.set_n(nced0);
    cav.lcfac.set_n(ncfa0);
    cav.lctet.set_n(ncte0);
    writeMeshCavity("insert_cavity_fail0.meshb", 
                                    msh,cav);

    int ierr2 = aux_movePointCav(msh, cav, tdimp, iseed, iref, algnd);
    writeMeshCavity("insert_cavity_fail_move0.meshb", 
                                    msh,cav);

    wait();
  }
  #endif

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
  #if 0
  if(ierro != INS2D_ERR_CAVITYOPERATOR && DOPRINTS1()){
    printf("## DEBUG WAIT HERE ierro = %d\n",ierro);
    msh.param->iverb = 0;
    msh.param->ivdepth = 0;
    if(ierro == INS2D_ERR_SHORTEDG){
      printf("Debug try many barys and compute length:\n");
      int edg2pol[2] = {ip1, ip2};
      int edg2po2[2] = {cav.ipins, -1};
      for(int itry = 0; itry < 20; itry++){
        bar1 = (itry + 1.0) / 21;
        if(itry == 19) bar1 = 1 - 0.547745;
        double bar2[2] = {bar1, 1 - bar1};
        METRIS_ENFORCE(msh.idim == 3);
        METRIS_ENFORCE(msh.curdeg == 1);
        eval1<3,1>(msh.coord, edg2pol, msh.getBasis(), 
                         DifVar::None, DifVar::None, 
                         bar2, msh.coord[cav.ipins], NULL, NULL);
        ierro = msh.interpMetBack(cav.ipins,tdimp,iseed,iref,algnd);
        METRIS_ENFORCE(ierro == 0);

        edg2po2[1] = ip1;
        double sz[2];
        double len1 = getlenedg_geosz<MFT,3,1>(msh,edg2po2,sz);
        edg2po2[1] = ip2;
        double len2 = getlenedg_geosz<MFT,3,1>(msh,edg2po2,sz);

        printf("bar1 = %e lens = %e %e valid %d %d \n",bar1,len1,len2,
          len1 > 1/sqrt(2), len2 > 1/sqrt(2));
      }

      // Try a bisection search, figure out number of iterations necessary
      printf("-- bissection\n");
      double bar1_min = 1.0e-6, bar1_max = 1 - 1.0e-6;
      for(int itry = 0; itry < 20; itry++){
        bar1 = (bar1_min + bar1_max)/2;
        double bar2[2] = {bar1, 1 - bar1};
        METRIS_ENFORCE(msh.idim == 3);
        METRIS_ENFORCE(msh.curdeg == 1);
        eval1<3,1>(msh.coord, edg2pol, msh.getBasis(), 
                         DifVar::None, DifVar::None, 
                         bar2, msh.coord[cav.ipins], NULL, NULL);
        ierro = msh.interpMetBack(cav.ipins,tdimp,iseed,iref,algnd);
        METRIS_ENFORCE(ierro == 0);

        edg2po2[1] = ip1;
        double sz[2];
        double len1 = getlenedg_geosz<MFT,3,1>(msh,edg2po2,sz);
        edg2po2[1] = ip2;
        double len2 = getlenedg_geosz<MFT,3,1>(msh,edg2po2,sz);

        printf("bar1 = %e lens = %e %e valid %d %d \n",bar1,len1,len2,
          len1 > 1/sqrt(2), len2 > 1/sqrt(2));

        if(len1 <= 1/sqrt(2)){
          bar1_max = bar1;
        }else if(len2 <= 1/sqrt(2)){
          bar1_min = bar1;
        }

      }
    }
    wait();
  }
  #endif
  if(DOPRINTS1()){
    msh.param->iverb = iverb0;
    msh.param->ivdepth = ivdepth0;
    printf("## END OF OPERATION WAIT\n");
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






template<class MFT>
int aux_movePointCav(Mesh<MFT>& msh, MshCavity &cav, 
                     int tdimp, int iseed, int iref, double *algnd){
  GETVDEPTH(msh.param);
  int ierro = 0;

  // Interior and non-surface case, most straightforward
  if(tdimp == msh.idim){
    const intAr2 &ent2poi = msh.ent2poi(msh.get_tdim());
    int tdim = tdimp;

    double bary[4];
    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim + 1);

    double newp[3] = {0,0,0};
    double eval[3];
    //double metl[6];
    //const int nnmet = (msh.idim*(msh.idim+1))/2;
    double meast = 0;
    const intAr1& lcent = cav.lcent(tdim);
    for(int ientt : lcent){
      double wt;
      bool iflat;
      CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
        // gdim == tdim here
        if constexpr(gdim == 2){
          eval2<gdim,1>(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                        DifVar::None,bary,eval,NULL,NULL);
          wt = getmeasentP1<gdim,2>(msh, ent2poi[ientt], algnd, &iflat);
        }else{
          eval3<gdim,1>(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                        DifVar::None,bary,eval,NULL,NULL);
          wt = getmeasentP1<gdim,3>(msh, ent2poi[ientt], algnd, &iflat);
        }
      }}CT_FOR1(gdim);
      if(iflat) continue;

      //msh.met.getMetBary(AsDeg::P1, DifVar::None, msh.met.getSpace(), 
      //                   ent2poi[ientt], msh.get_tdim(), bary, 
      //                   metl, NULL);

      // For simply barycentre, use meas0
      // To skew towards the largest elements use meas0*meas0
      for(int ii = 0; ii < msh.idim;ii++) newp[ii] += wt*eval[ii];
      meast += wt;
    }
    METRIS_ASSERT(meast > 0);
    for(int ii = 0; ii < msh.idim;ii++){
      newp[ii] /= meast;
      msh.coord(cav.ipins,ii) = newp[ii];
    }

  }else{
    // Boundary case. 
    return 0;
  }// if tdimp == msh.get_tdim()

  // reinterp metric. This is always interior case, no need for ref of bdry dir
  ierro = msh.interpMetBack(cav.ipins,tdimp,iseed,iref,algnd);
  if(ierro != 0){
    CPRINTF1(" - interpMetBack failed ierro = %d \n",ierro);
    ierro = INS2D_ERR_INTERPMETBACK;
  }
  
  return ierro;
}

















} // end namespace
