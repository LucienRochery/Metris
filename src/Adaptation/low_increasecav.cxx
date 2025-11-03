//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_increasecav.hxx"
#include "low_delaunay.hxx"
#include "Insertion/low_insert.hxx" // for error codes
#include "Insertion/aux_insert.hxx"
#include "Insertion/EdgeSeed.hxx"
#include "low_cavqual.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../low_geo/normal.hxx"
#include "../low_geo/measure.hxx"
#include "../low_geo/lenedg.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../Mesh/Mesh.hxx"
#include "../io_libmeshb.hxx"
#include "../smoothing/low_smoolen.hxx"

#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"

//#define NODELSURF

namespace Metris{



template<class MFT>
int setCavityInsertion2(Mesh<MFT>& msh, MshCavity &cav,
                        const EdgeSeed &insertionSeed,
                        int miter, int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);
  //const double dens_ideal = msh.get_tdim() == 2 ? pi / 4 : 0.54;
  int ierro = 0;
  double dens0, dens1;

  int ncent1_prev = cav.lcedg.get_n(), ncent1;
  int ncent2_prev = cav.lcfac.get_n(), ncent2;
  int ncent3_prev = cav.lctet.get_n(), ncent3;
  for(int niter = 0; niter < miter; niter++){

    int ncedg0 = cav.lcedg.get_n();
    int ncfac0 = cav.lcfac.get_n();
    int nctet0 = cav.lctet.get_n();

    ierro = movePointCavLen<MFT>(msh, cav, 5, ithrd1);
    if(ierro != 0) goto cleanup_loop;

    ierro = increase_cavity_Delaunay(msh, cav, insertionSeed.tdimp, -1, ithrd1);
    if(ierro != 0) goto cleanup_loop;

    ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
    if(ierro != 0) goto cleanup_loop;

    collrejcav_dens(msh,cav,&dens0,&dens1,ithrd1,ithrd2);
    ncent1 = cav.lcedg.get_n();
    ncent2 = cav.lcfac.get_n();
    ncent3 = cav.lctet.get_n();
    CPRINTF1(" - iter {} + del, dens0 {} dens1 {} ncent {} {} {}\n",niter,dens0,dens1,ncent1,ncent2,ncent3);

    if(dens1 < dens0){
      CPRINTF1(" # insertion is leading to lower density: {} -> {}, reject\n",dens0,dens1);
      goto cleanup_loop;
    }

    if(ncent1 == ncent1_prev && ncent2 == ncent2_prev && ncent3 == ncent3_prev){
      CPRINTF1(" - no change in cavity, stopping\n");
      return 0;
    }

    ncent1_prev = ncent1;
    ncent2_prev = ncent2;
    ncent3_prev = ncent3;

    continue;
    cleanup_loop:
    cav.lcedg.set_n(ncedg0);
    cav.lcfac.set_n(ncfac0);
    cav.lctet.set_n(nctet0);
    return ierro;
  }

  return 0;
}

template
int setCavityInsertion2<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, MshCavity &cav,
                        const EdgeSeed &insertionSeed,
                        int miter, int ithrd1, int ithrd2);
template
int setCavityInsertion2<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, MshCavity &cav,
                        const EdgeSeed &insertionSeed,
                        int miter, int ithrd1, int ithrd2);



template<class MFT>
int setCavityInsertion(Mesh<MFT>& msh, MshCavity &cav, const CavOprOpt &opts,
                       const EdgeSeed &insertionSeed, int mgrow, double lenqua_short_max,
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);
  int ierro;

  const int tdim  = insertionSeed.tdim_adp;
  const int tdimp = insertionSeed.tdimp;
  const int iseed = insertionSeed.iseed;
  const int iref  = insertionSeed.iref;
  const int nnmet = (msh.idim*(msh.idim+1))/2;

  bool filter_long = true;

  static int nwarnprt = 0;
  if(nwarnprt++ < 10 && !filter_long) printf("## FILTER_LONG SET TO TRUE\n");

  CPRINTF1("-- START setCavityInsertion tdim = {} mgrow = {}\n",tdim,mgrow);
  intWrkAr1 lrempoi = msh.get_iwork(10);
  lrempoi.set_n(0);

  // Try commenting this out since we now move the point.
  //// Check any close constrained points
  //ierro = aux_findCloseConstrained(msh, cav, ithrd1, ithrd2);
  //if(ierro > 0) return INS2D_ERR_SHORTCSTR;

  const int ibins = msh.poi2ebp(cav.ipins,tdimp,iseed,iref);

  int nprem;
  double coor0[3], met0[6], uv0[2];
  for(int ngrow = 0; ngrow < mgrow; ngrow++){
    INCVDEPTH(msh.param);

    // If need to revert elements
    int nced1 = cav.lcedg.get_n();
    int ncfa1 = cav.lcfac.get_n();
    int ncte1 = cav.lctet.get_n();
    for(int ii = 0; ii < msh.idim; ii++) coor0[ii] = msh.coord(cav.ipins, ii);
    for(int ii = 0; ii < nnmet   ; ii++) met0[ii] = msh.met(cav.ipins, ii);
    if(ibins >= 0) for(int ii = 0; ii < 2 ; ii++) uv0[ii] = msh.bpo2rbi(ibins, ii);

    ierro = 0;
    if(DOPRINTS2()){
      std::string fname = "insert_cavity0."+std::to_string(ngrow);
      writeMeshCavity(fname,msh,cav);
      msh.met.writeMetricFile(fname);
    }
    CPRINTF1(" - step {} cavity nedge {} nface {} nelem {}\n",ngrow,
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());

    ierro = movePointCavLen<MFT>(msh, cav, 5, ithrd1);
    if(DOPRINTS2()){
      writeMeshCavity("insert_cavity1."+std::to_string(ngrow),msh,cav);
      msh.met.writeMetricFile("insert_cavity1."+std::to_string(ngrow));
    }
    if(ierro > 0){
      CPRINTF1(" # movePointCavLen error {}\n",ierro);
      ierro = INS2D_ERR_MOVPTCAVLEN;
      goto finish_grow_step;
    }

    nprem = increase_cavity_lenedg(msh,cav,opts,cav.ipins,ithrd1,ithrd2);
    CPRINTF1(" - +remp nedge {} nface {} nelem {} nprem = {}\n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n(),nprem);
    if(DOPRINTS2()){
      writeMeshCavity("insert_cavity2."+std::to_string(ngrow),msh,cav);
      msh.met.writeMetricFile("insert_cavity2."+std::to_string(ngrow));
    }
    if(nprem < 0){
      ierro = -nprem;
      CPRINTF1(" # increase_cavity_lenedg error {}\n",nprem);
      ierro = INS2D_ERR_INCCAVLEN;
      goto finish_grow_step;
    }

    // -- 1 step Delaunay increase
    ierro = increase_cavity_Delaunay(msh, cav, tdim, 1, ithrd1);
    if(ierro != 0){
      CPRINTF1(" # +del error {}\n",ierro);
      ierro = INS2D_ERR_INCCAVDEL;
      goto finish_grow_step;
    }
    CPRINTF1(" - +del nedge {} nface {} nelem {}\n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());
    if(DOPRINTS2()){
      writeMeshCavity("insert_cavity3."+std::to_string(ngrow),msh,cav);
      msh.met.writeMetricFile("insert_cavity3."+std::to_string(ngrow));
    }


    // -- increase for validity
    ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
    if(ierro != 0){
      CPRINTF1(" # +cav error {}\n",ierro);
      ierro = INS2D_ERR_INCCAVVAL1;
      goto finish_grow_step;
    }
    CPRINTF1(" - +cav nedge {} nface {} nelem {}\n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());
    if(DOPRINTS2()){
      writeMeshCavity("insert_cavity4."+std::to_string(ngrow),msh,cav);
      msh.met.writeMetricFile("insert_cavity4."+std::to_string(ngrow));
    }

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
      ierro = INS2D_ERR_SHORTEDG1;
      CPRINTF1(" # error nrem point = {}\n",lrempoi.get_n());
      goto finish_grow_step;
    }

    finish_grow_step:
    if(ierro > 0){
      ierro = 0;


      for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins, ii) = coor0[ii];
      for(int ii = 0; ii < nnmet   ; ii++) msh.met(cav.ipins, ii)   = met0[ii];
      if(ibins >= 0) for(int ii = 0; ii < 2       ; ii++) msh.bpo2rbi(ibins, ii)   = uv0[ii];

      bool unfixable = false;
      if(lrempoi.get_n() > 0){
        // Now we need to remove all the newly added elements that contain
        // one of the lrempoi.
        CPRINTF2(" # Fix cavity, lrempoi = {}\n", lrempoi.get_n());
        msh.tag[ithrd1]++;
        for(int ii = 0; ii < lrempoi.get_n(); ii++){
          int ipoin = lrempoi[ii];
          msh.poi2tag(ithrd1, ipoin) = msh.tag[ithrd1];
        }
        for(int tdimc = 1; tdimc <= msh.get_tdim(); tdimc++){
          intAr1 &lcent = cav.lcent(tdimc);
          const int ncen0 = tdimc == 1 ? nced1 :
                            tdimc == 2 ? ncfa1 : ncte1;
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
            CPRINTF1(" - remove {} from cavity dim {}\n",icent,tdimc);
            int icend = lcent.pop();
            // This can only happen if we're the last element. In that case we
            // shrank the array and can quit.
            if(icend == icent) break;
            // otherwise place last here.
            icent = icend;
            nrem++;
          }// for icent
          CPRINTF1(" - removed {} dim {} cavity elements\n",nrem,tdimc);
        }// for tdimc

        // Try correcting cavity for validity then rechecking
        ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
        if(ierro != 0){
          CPRINTF1(" # +cav error after fix {}\n",ierro);
          unfixable = true;
          goto finish_correction;
        }

        check_cavity_rempoint(msh, cav, opts, lrempoi.get_array(), true, ithrd1);
        if(lrempoi.get_n() > 0){
          ierro = INS2D_ERR_SHORTEDG2;
          CPRINTF1(" # error nrem point = {} after fix\n",lrempoi.get_n());
          unfixable = true;
          goto finish_correction;
        }

      }

      finish_correction:
      if(lrempoi.get_n() == 0 || unfixable){
        CPRINTF1(" # Unfixable cavity: reset to: {} edges, {} faces, {} tetra and test\n",
                 nced1, ncfa1, ncte1);
        // The cavity can't be fixed to continue iterating. Simply stop it now.
        cav.lcedg.set_n(nced1);
        cav.lcfac.set_n(ncfa1);
        cav.lctet.set_n(ncte1);

        ierro = collrejcav_lenqua(msh, cav, filter_long, false, true, lenqua_short_max, nocomp, ithrd2);
        if(ierro > 0){
          CPRINTF1(" # collrejcav_lenqua rejects cavity\n");
          return INS2D_ERR_SHORTEDG3;
        }

        ierro = 0;
        break;
      }
      if(DOPRINTS2()) writeMeshCavity("insert_cavity4."+std::to_string(ngrow),msh,cav);
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

  }// for ngrow

  if(ierro > 0) return ierro;

  ierro = collrejcav_lenqua(msh, cav, filter_long, false, true, lenqua_short_max, nocomp, ithrd2);
  if(ierro > 0) return INS2D_ERR_LENQUA;

  return 0;
}

template
int setCavityInsertion<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, MshCavity &cav, const CavOprOpt &opts,
                       const EdgeSeed &insertionSeed, int mgrow, double lenqua_short_max,
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2);
template
int setCavityInsertion<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, MshCavity &cav, const CavOprOpt &opts,
                       const EdgeSeed &insertionSeed, int mgrow, double lenqua_short_max,
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2);


// Same as setCavityInsertion but try to cut down on time
// by moving the point only once.
template<class MFT>
int setCavityInsertion3(Mesh<MFT>& msh, MshCavity &cav, const CavOprOpt &opts,
                       const EdgeSeed &insertionSeed, int mgrow, double lenqua_short_max,
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);

  int ierro;

  const int tdim  = insertionSeed.tdim_adp;
  //const int tdimp = insertionSeed.tdimp;
  //const int iseed = insertionSeed.iseed;
  //const int iref  = insertionSeed.iref;
  //const int nnmet = (msh.idim*(msh.idim+1))/2;

  const bool filter_long = true;

  static int nwarnprt = 0;
  if(nwarnprt++ < 10 && !filter_long) printf("## FILTER_LONG SET TO TRUE\n");


  int nced1 = cav.lcedg.get_n();
  int ncfa1 = cav.lcfac.get_n();
  int ncte1 = cav.lctet.get_n();

  CPRINTF1("-- START setCavityInsertion tdim = {} mgrow = {}\n",tdim,mgrow);
  intWrkAr1 lrempoi = msh.get_iwork(10);
  lrempoi.set_n(0);

  // Try commenting this out since we now move the point.
  //// Check any close constrained points
  //ierro = aux_findCloseConstrained(msh, cav, ithrd1, ithrd2);
  //if(ierro > 0) return INS2D_ERR_SHORTCSTR;

  //const int ibins = msh.poi2ebp(cav.ipins,tdimp,iseed,iref);

  ierro = movePointCavLen<MFT>(msh, cav, 5, ithrd1);
  if(DOPRINTS2()){
    writeMeshCavity("insert_cavity0",msh,cav);
    msh.met.writeMetricFile("insert_cavity0");
  }
  if(ierro > 0){
    CPRINTF1(" # movePointCavLen error {}\n",ierro);
    return INS2D_ERR_MOVPTCAVLEN;
  }

  ierro = increase_cavity_Delaunay(msh,cav,tdim,5,ithrd1);
  if(DOPRINTS2()){
    writeMeshCavity("insert_cavity1",msh,cav);
  }
  if(ierro > 0){
    CPRINTF1(" # increase_cavity_Delaunay error {}\n",ierro);
    return INS2D_ERR_INCCAVDEL;
  }

  ierro = increase_cavity_validity(msh,cav,ithrd1);
  if(DOPRINTS2()){
    writeMeshCavity("insert_cavity2",msh,cav);
  }
  if(ierro > 0){
    CPRINTF1(" # increase_cavity_validity error {}\n",ierro);
    return INS2D_ERR_INCCAVVAL1;
  }

  // Check if the cavity needs fixing.
  check_cavity_rempoint(msh, cav, opts, lrempoi.get_array(), true, ithrd1);
  if(lrempoi.get_n() == 0) goto finish_cavity;

  // Now we need to remove all the newly added elements that contain
  // one of the lrempoi.
  CPRINTF2(" # Fix cavity, lrempoi = {}\n", lrempoi.get_n());
  msh.tag[ithrd1]++;
  for(int ii = 0; ii < lrempoi.get_n(); ii++){
    int ipoin = lrempoi[ii];
    msh.poi2tag(ithrd1, ipoin) = msh.tag[ithrd1];
  }
  for(int tdimc = 1; tdimc <= msh.get_tdim(); tdimc++){
    intAr1 &lcent = cav.lcent(tdimc);
    const int ncen0 = tdimc == 1 ? nced1 :
                      tdimc == 2 ? ncfa1 : ncte1;
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
      CPRINTF1(" - remove {} from cavity dim {}\n",icent,tdimc);
      int icend = lcent.pop();
      // This can only happen if we're the last element. In that case we
      // shrank the array and can quit.
      if(icend == icent) break;
      // otherwise place last here.
      icent = icend;
      nrem++;
    }// for icent
    CPRINTF1(" - removed {} dim {} cavity elements\n",nrem,tdimc);
  }// for tdimc

  // Try correcting cavity for validity then rechecking
  ierro = increase_cavity_validity(msh,cav,ithrd1);
  if(ierro != 0){
    CPRINTF1(" # +cav error after fix {}\n",ierro);
    return INS2D_ERR_INCCAVVAL2;
  }

  check_cavity_rempoint(msh, cav, opts, lrempoi.get_array(), true, ithrd1);
  if(lrempoi.get_n() > 0){
    CPRINTF1(" # error nrem point = {} after fix\n",lrempoi.get_n());
    return INS2D_ERR_SHORTEDG2;
  }


finish_cavity:

  ierro = collrejcav_lenqua(msh, cav, filter_long, false, true, lenqua_short_max, nocomp, ithrd2);
  if(ierro > 0) return INS2D_ERR_LENQUA;

  return 0;
}

template
int setCavityInsertion3<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, MshCavity &cav, const CavOprOpt &opts,
                       const EdgeSeed &insertionSeed, int mgrow, double lenqua_short_max,
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2);
template
int setCavityInsertion3<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, MshCavity &cav, const CavOprOpt &opts,
                       const EdgeSeed &insertionSeed, int mgrow, double lenqua_short_max,
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2);


template<class MFT>
int setCavityInsertion2(Mesh<MFT>& msh, MshCavity &cav, const CavOprOpt &opts,
                       int mgrow, double lenqua_short_max,
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);
  int ierro;

  intWrkAr1 lrempoi = msh.get_iwork(10);

  const int tdim = msh.get_tdim();

  // Check any close constrained points
  ierro = aux_findCloseConstrained(msh, cav, ithrd1, ithrd2);
  if(ierro > 0) return INS2D_ERR_SHORTCSTR2;

  for(int ngrow = 0; ngrow < mgrow; ngrow++){

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
    CPRINTF1(" - step {} cavity nedge {} nface {} nelem {}\n",ngrow,
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());

    int nprem = increase_cavity_lenedg(msh,cav,opts,cav.ipins,ithrd1,ithrd2);
    CPRINTF1(" - +remp nedge {} nface {} nelem {} nprem = {}\n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n(),nprem);
    if(DOPRINTS2()) writeMeshCavity("insert_cavity1."+std::to_string(ngrow),msh,cav);

    // -- 1 step Delaunay increase
    ierro = increase_cavity_Delaunay(msh, cav, tdim, 1, ithrd1);
    if(ierro != 0){
      CPRINTF1(" # +del error {}\n",ierro);
      ierro = INS2D_ERR_INCCAVDEL;
      goto finish_grow_step;
    }
    CPRINTF1(" - +del nedge {} nface {} nelem {}\n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());
    if(DOPRINTS2()) writeMeshCavity("insert_cavity2."+std::to_string(ngrow),msh,cav);


    // -- increase for validity
    ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
    if(ierro != 0){
      CPRINTF1(" # +cav error {}\n",ierro);
      ierro = INS2D_ERR_INCCAVVAL2;
      goto finish_grow_step;
    }
    CPRINTF1(" - +cav nedge {} nface {} nelem {}\n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());
    if(DOPRINTS2()) writeMeshCavity("insert_cavity3."+std::to_string(ngrow),msh,cav);

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
      ierro = INS2D_ERR_SHORTEDG4;
      CPRINTF1(" # error nrem point = {}\n",lrempoi.get_n());
      goto finish_grow_step;
    }

    finish_grow_step:
    if(ierro > 0){
      ierro = 0;
      if(lrempoi.get_n() == 0){
        CPRINTF1(" # Unfixable cavity: reset to: {} edges, {} faces, {} tetra and test\n",
                 nced1, ncfa1, ncte1);
        // The cavity can't be fixed to continue iterating. Simply stop it now.
        cav.lcedg.set_n(nced1);
        cav.lcfac.set_n(ncfa1);
        cav.lctet.set_n(ncte1);

        ierro = collrejcav_lenqua(msh, cav, true, false, true, lenqua_short_max, nocomp, ithrd2);
        if(ierro > 0) return INS2D_ERR_SHORTEDG5;

        ierro = 0;
        break;
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
          CPRINTF1(" - remove {} from cavity dim {}\n",icent,tdimc);
          int icend = lcent.pop();
          // This can only happen if we're the last element. In that case we
          // shrank the array and can quit.
          if(icend == icent) break;
          // otherwise place last here.
          icent = icend;
          nrem++;
        }// for icent
        CPRINTF1(" - removed {} dim {} cavity elements\n",nrem,tdimc);
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
    if(ierro > 0) break;

  }// for ngrow

  if(ierro > 0) return ierro;

  ierro = collrejcav_lenqua(msh, cav, true, false, true, lenqua_short_max, nocomp, ithrd2);
  if(ierro > 0) return INS2D_ERR_LENQUA;

  return 0;
}

template
int setCavityInsertion2<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, MshCavity &cav, const CavOprOpt &opts,
                       int mgrow, double lenqua_short_max,
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2);
template
int setCavityInsertion2<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, MshCavity &cav, const CavOprOpt &opts,
                       int mgrow, double lenqua_short_max,
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2);

// Cavity growth based on quality
template<class MFT>
int setCavityInsertionQuality(Mesh<MFT>& msh, MshCavity &cav, const CavOprOpt &opts,
                       const EdgeSeed &insertionSeed, int mgrow, double lenqua_short_max,
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);

  int ierro;

  const int tdim  = insertionSeed.tdim_adp;

  const bool filter_long = true;

  static int nwarnprt = 0;
  if(nwarnprt++ < 10 && !filter_long) printf("## FILTER_LONG SET TO TRUE\n");


  int nced1 = cav.lcedg.get_n();
  int ncfa1 = cav.lcfac.get_n();
  int ncte1 = cav.lctet.get_n();

  CPRINTF1("-- START setCavityInsertionQuality tdim = {} mgrow = {}\n",tdim,mgrow);
  intWrkAr1 lrempoi = msh.get_iwork(10);
  lrempoi.set_n(0);

  ierro = increase_cavity_quality(msh,cav,tdim,5,ithrd1);
  if(DOPRINTS2()){
    writeMeshCavity("insert_cavity1",msh,cav);
  }
  if(ierro > 0){
    CPRINTF1(" # increase_cavity_quality error {}\n",ierro);
    return INS2D_ERR_INCCAVDEL;
  }

  return 0;
}

template
int setCavityInsertionQuality<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, MshCavity &cav, const CavOprOpt &opts,
                              const EdgeSeed &insertionSeed, int mgrow, double lenqua_short_max,
                              std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                              int ithrd1, int ithrd2);
template
int setCavityInsertionQuality<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, MshCavity &cav, const CavOprOpt &opts,
                              const EdgeSeed &insertionSeed, int mgrow, double lenqua_short_max,
                              std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                              int ithrd1, int ithrd2);



// Check if any removed points; only those > 1/sqrt(2) from ipins if chklen
// This can possibly be reworked to be faster, for now we check everything every
// time, even though this is called in iterative cavity building.
template<class MFT>
void check_cavity_rempoint(MeshMetric<MFT> &msh, MshCavity &cav, const CavOprOpt &opts,
                           intAr1 &lrempoi, bool chklen, int ithrd1){
  GETVDEPTH(msh.param);

  lrempoi.set_n(0);

  for(int ientt : cav.lcedg) msh.edg2tag(ithrd1,ientt) = msh.tag[ithrd1];
  for(int ientt : cav.lcfac) msh.fac2tag(ithrd1,ientt) = msh.tag[ithrd1];
  for(int ientt : cav.lctet) msh.tet2tag(ithrd1,ientt) = msh.tag[ithrd1];

  // Points to be removed are those that are surrounded by only cavity elements.
  // Hence, loop over cavity elements and tag any points that belong to a
  // non-cavity neighbour.
  // Lastly, count untagged vertices.

  // If it belongs to any lower dim elements, that should be in the cavity.
  // It suffice there is one, as if it doesnt belong to all, it would be tagged.

  int tdimn = cav.lctet.get_n() > 0 ? 3
            : cav.lcfac.get_n() > 0 ? 2
                                    : 1;
  const intAr1&  lcent = cav.lcent(tdimn);
  const intAr2&  ent2ent = msh.ent2ent(tdimn);
  const intAr2&  ent2poi = msh.ent2poi(tdimn);
  const intAr2r& ent2tag = msh.ent2tag(tdimn);

  // ipins should always be seeded with a newbpotopo if it is going to be bdry
  const int pdim_ipins = msh.getpoitdim(cav.ipins);
  METRIS_ASSERT_MSG(pdim_ipins >= 0 && pdim_ipins <= msh.get_tdim(),
                    "invalid pdim_ipins = {}", pdim_ipins);

  // Tag points that won't be deleted: there is at least one elt outside
  // the cavity that has the point.
  for(int ientt : lcent){
    for(int ii = 0; ii < tdimn + 1; ii++){
      int ipoin = ent2poi(ientt,ii);
      // Cycle neighbours that have ii (i.e. all but ii-th neighbour)
      for(int jj = 0; jj < tdimn + 1; jj++){
        if(jj == ii) continue;
        int ient2 = ent2ent(ientt,jj);
        if(ient2 < 0) continue;
        // Tag point if the adjacent element is not in the cavity
        // This point is not set to be deleted.
        if(ent2tag(ithrd1,ient2) < msh.tag[ithrd1]){
          msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];
          CPRINTF2("  - not rem point {} \n", ipoin);
        }
      }
    }
  }

  // Go over elements, counting vertices that have not been tagged.
  for(int ientt : lcent){
    for(int ii = 0; ii < tdimn + 1; ii++){
      int ipoin = ent2poi(ientt,ii);
      if(ipoin == cav.ipins) continue;
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
      CPRINTF2("  - rem pt ? {} \n", ipoin);

      // Check the point dimension wrt to option allow_remove_points_superdim
      int pdim = msh.getpoitdim(ipoin);
      if(pdim > pdim_ipins && opts.allow_remove_points_superdim){
        CPRINTF1(" - point dim {} > {} = dim(ipins) "
                 "with allow_remove_points_superdim, skip check\n",
                 pdim, pdim_ipins);
        continue;
      }

      // point going to be deleted, but only if any existing lower dim entities
      // are also in the cavity.
      if(tdimn == 3){
        // If there is a face attached, check it is in the cavity.
        int iface = getpoifac(msh, ipoin);
        // If not, this point won't be removed. Continue.
        if(iface >= 0 && msh.fac2tag(ithrd1,iface) < msh.tag[ithrd1]) continue;
      }

      if(tdimn >= 2){
        // If there is an edge attached, check it is in the cavity.
        int iedge = getpoiedg(msh,ipoin);
        // If not, this point won't be removed. Continue.
        if(iedge >= 0 && msh.edg2tag(ithrd1,iedge) < msh.tag[ithrd1]) continue;
      }

      // If we're here, that means that there are either no attached lower dim
      // or there are and they are all in the cavity; indeed, assume there exist
      // at least one, and at least one not in the cav. Then the point is not
      // tagged. Then we wouldn't be here.

      CPRINTF1(" ## point {} will be removed \n",ipoin);
      // tag point so we don't check for it again
      msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];
      if(chklen){
        int edg2pol[2] = {cav.ipins, ipoin};
        double sz[2];
        double len = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,edg2pol,sz)
                                   : getlenedg_geosz<MFT,3,1>(msh,edg2pol,sz);
        CPRINTF1(" -> found len = {} >? 1/sqrt(2): {}\n",len,len*sqrt(2) > 1);
        if(len > 1.0/sqrt(2)) lrempoi.stack(ipoin);
      }else{
        lrempoi.stack(ipoin);
      }
    }
  }

  return;
}

template void check_cavity_rempoint<MetricFieldAnalytical>
  (MeshMetric<MetricFieldAnalytical> &msh, MshCavity &cav, const CavOprOpt &opts,
   intAr1 &lrempoi, bool chklen, int ithrd1);
template void check_cavity_rempoint<MetricFieldFE        >
  (MeshMetric<MetricFieldFE        > &msh, MshCavity &cav, const CavOprOpt &opts,
   intAr1 &lrempoi, bool chklen, int ithrd1);


// Increase for validity and Delaunay (if idelaunay == true) both.
// Argument ref2nordev is optional unless surface is involved. It need not be filled prior.
template<class MFT>
int increase_cavity(MeshMetric<MFT>& msh, MshCavity& cav,
                    bool idelaunay, int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(cav.ipins >= 0 && cav.ipins < msh.npoin);

  static int nwarnprt = 0;
  if(nwarnprt++ < 4){
    if(msh.param->iverb > 0) PRINTF("## Could move nordev checks to increase_cavity. For this, we need to precompute the ccos as in reconnect_faccav.\n");
  }

  #ifndef NDEBUG
  intWrkAr1 lcedg0_ = msh.get_iwork(10);
  intWrkAr1 lcfac0_ = msh.get_iwork(100);
  intWrkAr1 lctet0_ = msh.get_iwork(100);
  intAr1& lcedg0 = lcedg0_.get_array();
  intAr1& lcfac0 = lcfac0_.get_array();
  intAr1& lctet0 = lctet0_.get_array();
  cav.lcedg.copyTo(lcedg0);
  cav.lcfac.copyTo(lcfac0);
  cav.lctet.copyTo(lctet0);
  #endif



  //#ifdef NODELSURF
  //static int nwarn = 0;
  //// Disable surf
  //if(msh.get_tdim() < msh.idim && msh.param->iflag1 == 0 && idelaunay){
  //  if(nwarn++ < 10) MPRINTF("## WARNING DELAUNAY SURFACE DISABLED\n");
  //  idelaunay = false;
  //}
  //#endif

  msh.tag[ithrd1]++;
  if(idelaunay) msh.tag[ithrd2]++;

  // Tag entities and references
  for(int tdim = 1; tdim <= 3; tdim++){
    const intAr1& lcent = cav.lcent(tdim);
    for(int ientt : lcent){
      METRIS_ASSERT(ientt >= 0 && ientt < msh.nentt(tdim));
      METRIS_ASSERT(!isdeadent(ientt,msh.ent2poi(tdim)));
      msh.ent2tag(tdim)(ithrd1,ientt) = msh.tag[ithrd1];

      int iref = msh.ent2ref(tdim)[ientt];
      METRIS_ASSERT(iref >= 0);
      if(msh.ref2tag(tdim)(ithrd1,iref) < msh.tag[ithrd1]){
        CPRINTF1(" - ipins has edge ref {} \n",iref);
      }
      msh.ref2tag(tdim)(ithrd1,iref) = msh.tag[ithrd1];
    }
  }

  int pdim = msh.getpoitdim(cav.ipins);
  #ifndef NDEBUG
  {
  int cav_mindim = cav.lcedg.get_n() > 0 ? 1 :
                   cav.lcfac.get_n() > 0 ? 2 : 3;
  METRIS_ASSERT(pdim == cav_mindim);
  }
  #endif

  CPRINTF1("-- START increase_cavity ipins {} dim {} list initial cavity:\n", cav.ipins, pdim);
  cav.print(msh);


  // Get normal deviation of initial cavity.
  if(msh.nperiodic_face != 0){
    METRIS_THROW_MSG("TODO: ## CASE WITH PERIODIC FACES NOT HANDLED IN LOW_INCREASECAV")
    // I think the way to generalize this is not to go all in on generality as in reconnect_faccav,
    // but to keep this "happy path" centered approach and work around the exceptions locally.
    // It is rare in practice to have periodic faces and, even when some exist, most won't be, in real geoms.
    // Moreover, dealing with this is price paid for each surface insertion, not just on edges.
    // We left place in ref2nordev for a second entry, only for periodic refs.
  }

  dblWrkAr1 ref2nordev_ = msh.get_rwork(2*msh.CAD.ncadfa);
  dblAr2 ref2nordev(msh.CAD.ncadfa, 2, &ref2nordev_[0]);
  ref2nordev.fill(-1);

  if(msh.idim >= 3){
    for(int iface : cav.lcfac){
      INCVDEPTH(msh.param);
      int iref = msh.fac2ref[iface];
      double nordev;
      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
        nordev = getnordev<ideg>(msh,iface,NULL);
      }}CT_FOR1(ideg);
      ref2nordev(iref,0) = MAX(ref2nordev(iref,0) , nordev);
      CPRINTF1(" - iface {} nordev = {}\n",iface,nordev);
    }

    if(DOPRINTS1()){
      CPRINTF1(" - initial cavity nordev:\n");
      for(int iref = 0; iref < msh.CAD.ncadfa; iref++){
        double nordev = ref2nordev(iref,0);
        if(nordev < 0) continue;
        CPRINTF1(" - iref {} nordev = {}\n", iref, nordev);
      }
    }
  }



  int ent2pol[4];
  ent2pol[0] = cav.ipins;

  bool restart;
  int niter = 0;
  int ient0[2] = {0,0};


  int nnmet = (msh.idim * (msh.idim + 1)) / 2;
  double metl[6], lmet[6];
  double *metl_p;
  if(idelaunay){
    if(msh.met.getSpace() == MetSpace::Log){
      for(int ii = 0; ii < nnmet; ii++) lmet[ii] = msh.met(cav.ipins,ii);
      if(msh.idim == 2){
        getexpmet_cpy<2>(lmet, metl);
      }else{
        getexpmet_cpy<3>(lmet, metl);
      }
      metl_p = metl;
    }else{
      metl_p = msh.met[cav.ipins];
    }
  }

  do{

    restart = false;
    if(niter++ > 100){
      #ifndef NDEBUG
      MPRINTF("# Possibly infinite cavity ext iterations 100\n");
      printf("## WAIT\n");
      wait();
      #endif
      return 1;
    }

    CPRINTF1(" - inccav iter {} ifac0 {} itetr0 {} \n",niter,ient0[0],ient0[1]);


    int ient01[2] = {cav.lcfac.get_n(), cav.lctet.get_n()};

    for(int tdim = 2; tdim <= msh.get_tdim(); tdim++){

      intAr1 &lcent = cav.lcent(tdim);
      const intAr2 &ent2poi = msh.ent2poi(tdim);
      const intAr2 &ent2ent = msh.ent2ent(tdim);
      const intAr1 &ent2ref = msh.ent2ref(tdim);
      const intAr2 &ref2tag = msh.ref2tag(tdim);
      intAr2 &ent2tag = msh.ent2tag(tdim);
      const intAr2 &sub2tag = msh.ent2tag(tdim-1);

      CPRINTF1(" - inccav tdim {} ncent {}\n",tdim,lcent.get_n());

      // Note the bound is reeval'd, can't use range based
      for(int ientl = ient0[tdim-2]; ientl < lcent.get_n(); ientl++){
        INCVDEPTH(msh.param)
        int ientt = lcent[ientl];
        if(tdim == 2){
          CPRINTF1(" - inccav try {} / {} = {} ({},{},{}) \n",
                   ientl,lcent.get_n(),ientt,ent2poi(ientt,0),
                   ent2poi(ientt,1),ent2poi(ientt,2));
        }else{
          CPRINTF1(" - inccav try {} / {} = {} ({},{},{},{}) \n",
                   ientl,lcent.get_n(),ientt,ent2poi(ientt,0),
                   ent2poi(ientt,1),ent2poi(ientt,2),ent2poi(ientt,3));
        }

        double norCAD[3];
        if(msh.idim == 3 && tdim == 2){
          // If dimension 3 topo dim 2, get a normal for this face.
          if(msh.CAD()){
            getnorfacCAD(msh,ientt,norCAD);
          }else{
            getnorfacP1(ent2poi[ientt],msh.coord,norCAD);
          }
        }

        for(int inei = 0; inei < tdim + 1; inei++){

          bool iskip = false;
          for(int ii = 0; ii < tdim + 1; ii++){
            int ip = ent2poi(ientt, ii);
            if(ip != cav.ipins) continue;
            iskip = true;
            break;
          }
          if(iskip) continue;

          int ienei = ent2ent(ientt,inei);

          CPRINTF1("   - inei {} ienei = {}\n", inei, ienei);

          if(ienei >= 0){
            if(ent2tag(ithrd1,ienei) >= msh.tag[ithrd1]){
              CPRINTF1("   - ienei = {} is tagged {} >= {}\n",
                                 ienei,ent2tag(ithrd1,ienei),msh.tag[ithrd1]);
              continue;
            }
          }



          // tdim 2: if there's an edge here and it's in the cavity, then it will
          // be split and we'll get no face from it.
          // tdim 3: idem, faces.
          int isube = -1;
          if(tdim == 2){
            isube = msh.fac2edg(ientt,inei);
          }else{
            isube = msh.tet2fac(ientt,inei);
          }
          if(isube >= 0 && sub2tag(ithrd1,isube) >= msh.tag[ithrd1]){
            CPRINTF1("   - ientt {} -> isube {} is tagged, skip\n",ientt,isube);
            continue;
          }

          // New face is ipins, ip1, ip2
          if(tdim == 2){
            ent2pol[1] = ent2poi(ientt,lnoed2[inei][0]);
            ent2pol[2] = ent2poi(ientt,lnoed2[inei][1]);
          }else{
            ent2pol[lnofa3[0][0]] = ent2poi(ientt,lnofa3[inei][0]);
            ent2pol[lnofa3[0][1]] = ent2poi(ientt,lnofa3[inei][1]);
            ent2pol[lnofa3[0][2]] = ent2poi(ientt,lnofa3[inei][2]);
          }

          bool iflat;
          int nod2bpo[3];

          #ifndef NDEBUG
          try{
          #endif


            double nordev_tol = -1;
            if(msh.idim == 3 && tdim == 2){
              int iref = msh.fac2ref[ientt];
              nordev_tol = ref2nordev(iref,0);
              nod2bpo[0] = pdim == 2 ? msh.poi2ebp(cav.ipins, 2, ientt, iref) : -1;
              nod2bpo[1] = msh.poi2ebp(ent2pol[1], 2, ientt, iref);
              nod2bpo[2] = msh.poi2ebp(ent2pol[2], 2, ientt, iref);
              CPRINTF1(" - using nordevtol = {} for face ref {}\n", nordev_tol, iref);
              METRIS_ASSERT(nod2bpo[0] < 0 || msh.bpo2ibi(nod2bpo[0],1) == 2);
              METRIS_ASSERT(msh.bpo2ibi(nod2bpo[1],1) == 2);
              METRIS_ASSERT(msh.bpo2ibi(nod2bpo[2],1) == 2);
              METRIS_ASSERT(nod2bpo[0] < 0 || msh.fac2ref[msh.bpo2ibi(nod2bpo[0],2)] == iref);
              METRIS_ASSERT(msh.fac2ref[msh.bpo2ibi(nod2bpo[1],2)] == iref);
              METRIS_ASSERT(msh.fac2ref[msh.bpo2ibi(nod2bpo[2],2)] == iref);
            }


            // First, check if this is a sliver
            double meas0;
            bool ivalid = msh.idim == 2 ? isvalideltP1<2,2>(msh, ent2pol, NULL   , NULL, &meas0, nordev_tol)
                        :     tdim == 2 ? isvalideltP1<3,2>(msh, ent2pol, nod2bpo, NULL, &meas0, nordev_tol)
                                        : isvalideltP1<3,3>(msh, ent2pol, NULL   , NULL, &meas0, nordev_tol); // NORCAD
            iflat = !ivalid;
            CPRINTF1("   - inccav pdim {} tdim {} ent {} = {}\n",pdim,tdim,ientt,
                    intAr1(tdim+1,ent2pol));
            CPRINTF1("   - w/ vtol = {:e} got iflat = {} meas0 = {:15.7e} neighbour = {}\n",
                    msh.param->vtol,iflat,meas0,ienei);

          #ifndef NDEBUG
          }catch(const MetrisExcept& e){

            PRINTF("## isvalideltP1 threw for ientt {} tdim {}, nodes: {}\n",ientt,tdim,intAr1(tdim+1,ent2pol));
            if(msh.idim == 3 && tdim == 2){
              int iref = msh.fac2ref[ientt];
              int ibins = msh.poi2bpo[cav.ipins];
              PRINTF("## nod2bpo[0] using ipins {} poi2bpo = {}, bpo2ibi: {}\n",cav.ipins, ibins, intAr1(nibi,msh.bpo2ibi[ibins]));
              PRINTF("## List all ipins bpoi:\n");
              for(int ibpoi = ibins; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
                int ient0 = msh.bpo2ibi(ibpoi,2);
                int tdim0 = msh.bpo2ibi(ibpoi,1);
                int iref0 = -1;
                if(tdim0 > 0) iref0 = msh.ent2ref(tdim0)[ient0];
                PRINTF("##  {}: {}, entity ref {}\n", ibpoi, intAr1(nibi,msh.bpo2ibi[ibpoi]),iref0);
              }
              PRINTF("## USING pdim_ins = {}\n",pdim);
              PRINTF("## nod2bpo[1] using ipoin {} ientt {} iref {}, got ibpoi {}, bpo2ibi: {}\n",
                     ent2pol[1], ientt, iref, nod2bpo[1],
                     intAr1(nibi,msh.bpo2ibi[nod2bpo[1]]));
              PRINTF("## nod2bpo[2] using ipoin {} ientt {} iref {}, got ibpoi {}, bpo2ibi: {}\n",
                     ent2pol[2], ientt, iref, nod2bpo[2],
                     intAr1(nibi,msh.bpo2ibi[nod2bpo[2]]));
              PRINTF("Cavity:\n");
              cav.print(msh, 10);
              PRINTF("\nInitial:\n");
              MPRINTF(" - Edge cavity: \n");
              for(int iecav : lcedg0){
                MPRINTF("  {} : {}\n",iecav,intAr1(2,msh.edg2poi[iecav]));
              }
              MPRINTF(" - Face cavity: \n");
              for(int iecav : lcfac0){
                MPRINTF("  {} : {}\n",iecav,intAr1(2,msh.fac2poi[iecav]));
              }
              MPRINTF(" - Tetra cavity: \n");
              for(int iecav : lctet0){
                MPRINTF("  {} : {}\n",iecav,intAr1(2,msh.tet2poi[iecav]));
              }

            }
            throw(e);
          }
          #endif

          #if 0
          // Next check geodev
          // Actually not because adding more faces will only damage the cavity further
          // Do this in the future as pre reject, possibly.
          // Also depends on Pk etc. Probably best to leave in cav.
          double nrmal[3];
          if(msh.idim == 3 && pdim < 2){
            // Get the normal in the case we're on an edge in 3D, and get only
            // the correct side.
            int iref = msh.fac2ref[ientt];
            getnorballref<1>(msh,lcent,iref,nrmal);
          }
          #endif

          // if element created with this facet is negative, add the neighbour
          // to cavity.
          if(iflat){
            if(ienei >= 0){
              if(ref2tag(ithrd1,ent2ref[ienei]) < msh.tag[ithrd1]){
                CPRINTF1("   - ienei = {} is wrong ref {} -> cannot correct\n",
                         ienei,ent2ref[ienei]);
                return 1;
              }
              lcent.stack(ienei);
              ent2tag(ithrd1,ienei) = msh.tag[ithrd1];
              CPRINTF1("   - inccav added entt {} to stack \n", ienei);
              // If this is a face, we must also add the supported tets
              if(tdim == 2 && msh.nelem > 0){
                for(int ii = 0; ii < 2; ii++){
                  int ielem = msh.fac2tet(ienei, ii);
                  if(ielem < 0) continue;
                  if(msh.tet2tag(ithrd1,ielem) >= msh.tag[ithrd1]) continue;
                  int iref = msh.tet2ref[ielem];
                  if(msh.dom2tag(ithrd1,iref) < msh.tag[ithrd1]){
                    CPRINTF1("   - iface {} -> itetr {} is wrong ref {}\n",ienei,ielem,iref);
                    return 1;
                  }
                  cav.lctet.stack(ielem);
                  msh.tet2tag(ithrd1,ielem) = msh.tag[ithrd1];
                  CPRINTF1("   - inccav added tet {} to stack \n", ielem);
                }
              }
            }

            // There are two cases:
            // - ienei >= 0, then this entity is sandwiched and needs to be added
            // - ienei < 0, then the only hope of correction is adding this entity
            // Hence, in any case, if there is a subdim entity here, add it.
            if(isube >= 0){
              // If the point tdim is greater than this element's dim, it cannot
              // be added.
              if(pdim > tdim-1){
                CPRINTF1("   - ientt {} dim {} < ipins dim {}\n",isube,tdim-1,pdim);
                return 2;
              }
              // Add the boundary entity, but only if in allowed refs.
              int iref = msh.ent2ref(tdim-1)[isube];
              if(msh.ref2tag(tdim-1)(ithrd1,iref) < msh.tag[ithrd1]){
                CPRINTF1("   - ientt {} -> isube {} is wrong ref {}\n",ienei,isube,iref);
                return 1;
              }
              cav.lcent(tdim-1).stack(isube);
              msh.ent2tag(tdim-1)(ithrd1,isube) = msh.tag[ithrd1];
              CPRINTF1("   - inccav added dim {} ent {} to stack \n", tdim-1,
                       isube);
              // We added a lower dim entity, hence restart will be required.
              restart = true;
            }


            // If added due to validity, skip delaunay
            continue;
          }// if iflat

          // Only apply Delaunay on highest tdim elements.
          if(idelaunay && ienei >= 0 && tdim == msh.get_tdim()){
            if(ent2tag(ithrd2,ienei) >= msh.tag[ithrd2]){
              CPRINTF1("   - ienei = {} has already been checked for delaunay -> skip\n",
                ienei);
              continue;
            }
            ent2tag(ithrd2,ienei) = msh.tag[ithrd2];

            if(ref2tag(ithrd1,ent2ref[ienei]) < msh.tag[ithrd1]){
              CPRINTF1("   - ienei = {} is wrong ref {} -> skip Delaunay\n",
                       ienei,ent2ref[ienei]);
              continue;
            }


            // Check if Delaunay
            bool isinsph;
            try{
              if(tdim == 2){
                if(msh.idim == 2){
                  isinsph = indelsphere<2,2>(msh, msh.coord[cav.ipins], metl_p,
                                            ent2poi[ienei]);
                }else{
                  isinsph = indelsphere<3,2>(msh, msh.coord[cav.ipins], metl_p,
                                            ent2poi[ienei]);
                }
              }else{
                isinsph = indelsphere<3,3>(msh, msh.coord[cav.ipins], metl_p,
                                          ent2poi[ienei]);
              }
            }catch(const MetrisExcept& e){
              fmt::print("indelsphere call threw exception\n");
              fmt::print("with ienei = {} nodes {} ipins {}\n",
                         ienei,intAr1(tdim+1,ent2poi[ienei]),cav.ipins);
              //double meas0;
              //bool ivalid = isvalideltP1<3,3>(msh, ienei, NULL, &meas0);
              //fmt::print("elt measure {} valid {}\n",meas0,ivalid);
              throw(e);
            }

            if(isinsph){
              lcent.stack(ienei);
              ent2tag(ithrd1,ienei) = msh.tag[ithrd1];

              if(isube >= 0){
                // Add the boundary entity, but only if in allowed refs.
                int iref = msh.ent2ref(tdim-1)[isube];
                if(msh.ref2tag(tdim-1)(ithrd1,iref) < msh.tag[ithrd1]){
                  CPRINTF1("   - ientt {} -> isube {} is wrong ref {}\n",ienei,isube,iref);
                  return 1;
                }
                cav.lcent(tdim-1).stack(isube);
                msh.ent2tag(tdim-1)(ithrd1,isube) = msh.tag[ithrd1];
                CPRINTF1("   - inccav added dim {} ent {} to stack \n", tdim-1,
                         isube);
                // We added a lower dim entity, hence restart will be required.
                restart = true;
              }
            }

          }// if idelaunay

        } // for int inei

      } // for int ientl
    } // for int tdim

    ient0[0] = ient01[0];
    ient0[1] = ient01[1];

  }while(restart);

  CPRINTF1("-- END increase_cavity final cavity:\n");
  cav.print(msh);
  return 0;
}


template int increase_cavity(MeshMetric<MetricFieldAnalytical> &msh, MshCavity &cav,
                    bool idelaunay, int ithrd1, int ithrd2);
template int increase_cavity(MeshMetric<MetricFieldFE        > &msh, MshCavity &cav,
                    bool idelaunay, int ithrd1, int ithrd2);








// Increase for validity. Only allow same refs as ipins already has.
int increase_cavity_validity(MeshBase &msh, MshCavity &cav, int ithread){
  GETVDEPTH(msh.param);

  static int nwarnprt = 0;
  if(nwarnprt++ < 4){
    if(msh.param->iverb > 0) PRINTF("## Could move nordev checks to increase_cavity. For this, we need to precompute the ccos as in reconnect_faccav.\n");
  }
  METRIS_ASSERT(cav.ipins >= 0 && cav.ipins < msh.npoin);

  msh.tag[ithread]++;

  // Tag point's surface references if any. Filter entities
  aux_taginsrefs(msh,cav,ithread);

  for(int ielem : cav.lctet){
    METRIS_ASSERT(ielem >= 0 && ielem < msh.nelem);
    METRIS_ASSERT(!isdeadent(ielem,msh.tet2poi));
    msh.tet2tag(ithread,ielem) = msh.tag[ithread];
  }

  for(int iface : cav.lcfac){
    METRIS_ASSERT(iface >= 0 && iface < msh.nface);
    METRIS_ASSERT(!isdeadent(iface,msh.fac2poi));
    msh.fac2tag(ithread,iface) = msh.tag[ithread];
    if(!msh.isboundary_faces()) continue;

    int iref = msh.fac2ref[iface];
    METRIS_ASSERT(iref >= 0);
    METRIS_ASSERT(msh.cfa2tag(ithread,iref) <= msh.tag[ithread]);
    if(msh.cfa2tag(ithread,iref) < msh.tag[ithread]){
      CPRINTF1("## ERROR increase_cavity_validity: cavity face ref {} is not a ipins bdry ref\n",iref);
      return 2;
    }
  }

  for(int iedge : cav.lcedg){
    METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);
    METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));
    msh.edg2tag(ithread,iedge) = msh.tag[ithread];
    if(!msh.isboundary_edges()) continue;

    int iref = msh.edg2ref[iedge];
    METRIS_ASSERT(msh.ced2tag(ithread,iref) <= msh.tag[ithread]);
    if(msh.ced2tag(ithread,iref) < msh.tag[ithread]){
      CPRINTF1("## ERROR increase_cavity_validity: cavity edge is not a ipins bdry ref\n");
      return 2;
    }
  }

  CPRINTF1("-- START increase_cavity_validity ipins {} list initial cavity:\n", cav.ipins);
  if(DOPRINTS1()){
    if(cav.lcedg.get_n() > 0){
      CPRINTF1(" - Edge cavity: {}\n",cav.lcedg);
    }
    if(cav.lcfac.get_n() > 0){
      CPRINTF1(" - Face cavity: {}\n",cav.lcfac);
    }
    if(cav.lctet.get_n() > 0){
      CPRINTF1(" - Tetra cavity: {}\n",cav.lctet);
    }
  }
  if(DOPRINTS2()){
    for(int tdim = 1; tdim <= 3; tdim++){
      intAr1 &lcent = cav.lcent(tdim);
      int ncent = lcent.get_n();
      if(ncent <= 0) continue;
      intAr2 &ent2poi = msh.ent2poi(tdim);

      CPRINTF2(" - {} cavity:\n", tdim == 1 ? "Edge" : tdim == 2 ? "Face" : "Tetra");
      const int nnode = msh.nnode(tdim);
      for(int ientt : lcent)
        CPRINTF2("{} : {}\n",ientt,intAr1(nnode,ent2poi[ientt]));
    }
  }

  int ibins = msh.poi2bpo[cav.ipins];
  int pdim  = msh.get_tdim();
  if(ibins >= 0) pdim = msh.bpo2ibi(ibins,1);

  int ent2pol[4];
  ent2pol[0] = cav.ipins;

  bool restart;
  int niter = 0;
  int ient0[2] = {0,0};

  do{

    restart = false;
    if(niter++ > 100){
      #ifndef NDEBUG
      MPRINTF("# Possibly infinite cavity ext iterations 100\n");
      printf("## WAIT\n");
      wait();
      #endif
      return 1;
    }

    CPRINTF1(" - inccav iter {} ifac0 {} itetr0 {} \n",niter,
             ient0[0],ient0[1]);


    int ient01[2] = {cav.lcfac.get_n(), cav.lctet.get_n()};

    // Note the bound is reeval'd, can't use range based
    for(int tdim = 2; tdim <= 3; tdim++){

      intAr1 &lcent = cav.lcent(tdim);
      const intAr2 &ent2poi = msh.ent2poi(tdim);
      const intAr2 &ent2ent = msh.ent2ent(tdim);
      intAr2 &ent2tag = msh.ent2tag(tdim);


      for(int ientl = ient0[tdim-2]; ientl < lcent.get_n(); ientl++){
        INCVDEPTH(msh.param)
        int ientt = lcent[ientl];
        if(tdim == 2){
          CPRINTF1(" - inccav try {} / {} = {} ({},{},{}) \n",
                   ientl,lcent.get_n(),ientt,ent2poi(ientt,0),
                   ent2poi(ientt,1),ent2poi(ientt,2));
        }else{
          CPRINTF1(" - inccav try {} / {} = {} ({},{},{},{}) \n",
                   ientl,lcent.get_n(),ientt,ent2poi(ientt,0),
                   ent2poi(ientt,1),ent2poi(ientt,2),ent2poi(ientt,3));
        }

        double norCAD[3];
        if(msh.idim == 3 && tdim == 2){
          // If dimension 3 topo dim 2, get a normal for this face.
          if(msh.CAD()){
            getnorfacCAD(msh,ientt,norCAD);
          }else{
            getnorfacP1(ent2poi[ientt],msh.coord,norCAD);
          }
        }

        for(int inei = 0; inei < tdim + 1; inei++){

          bool iskip = false;
          for(int ii = 0; ii < tdim + 1; ii++){
            int ip = ent2poi(ientt, ii);
            if(ip != cav.ipins) continue;
            iskip = true;
            break;
          }
          if(iskip) continue;

          int ienei = ent2ent(ientt,inei);

          CPRINTF1("   - inei {} ienei = {}\n", inei, ienei);

          if(ienei >= 0){
            if(ent2tag(ithread,ienei) >= msh.tag[ithread]){
              CPRINTF1("   - ienei = {} is tagged {} >= {}\n",
                                 ienei,ent2tag(ithread,ienei),msh.tag[ithread]);
              continue;
            }
            if(tdim == 2){
              int iref = msh.fac2ref[ienei];
              if(msh.cfa2tag(ithread,iref) < msh.tag[ithread] && msh.isboundary_faces()){
                CPRINTF1("   - ienei = {} is wrong bdry ref {}\n",ienei,iref);
                continue;
              }
            }else{
              if(msh.tet2ref[ienei] != msh.tet2ref[ientt]){
                CPRINTF1("   - ienei {} ref = {} != ientt {} ref {} -> skip\n",
                ienei,msh.tet2ref[ienei],ientt,msh.tet2ref[ientt]);
                continue;
              }
            }
          }

          // tdim 2: if there's an edge here and it's in the cavity, then it will
          // be split and we'll get no face from it.
          // tdim 3: idem, faces.
          int iedge = -1, iface = -1;
          if(tdim == 2){
            iedge = msh.fac2edg(ientt,inei);
            if(iedge >= 0){
              if(msh.edg2tag(ithread,iedge) >= msh.tag[ithread]){
                CPRINTF1("   - iface {} -> iedge {} is tagged, skip\n",ientt,iedge);
                continue;
              }
              //int iref = msh.edg2ref[iedge];
              //if(msh.ced2tag(ithread,iref) < msh.tag[ithread] && msh.isboundary_edges()){
              //  CPRINTF1("   - iface {} -> iedge {} is wrong bdry ref {}\n",ienei,iedge,iref);
              //  continue;
              //}
            }
          }else{
            iface = msh.tet2fac(ientt,inei);
            if(iface >= 0){
              if(msh.fac2tag(ithread,iface) >= msh.tag[ithread]){
                CPRINTF1("   - itetr {} -> iface {} is tagged, skip\n",ientt,iface);
                continue;
              }
              //int iref = msh.fac2ref[iface];
              //if(msh.cfa2tag(ithread,iref) < msh.tag[ithread]){
              //  CPRINTF1("   - itetr {} -> iface {} is wrong bdry ref {}\n",ientt,iface,iref);
              //  continue;
              //}
            }
          }

          // New face is ipins, ip1, ip2
          if(tdim == 2){
            ent2pol[1] = ent2poi(ientt,lnoed2[inei][0]);
            ent2pol[2] = ent2poi(ientt,lnoed2[inei][1]);
          }else{
            ent2pol[lnofa3[0][0]] = ent2poi(ientt,lnofa3[inei][0]);
            ent2pol[lnofa3[0][1]] = ent2poi(ientt,lnofa3[inei][1]);
            ent2pol[lnofa3[0][2]] = ent2poi(ientt,lnofa3[inei][2]);
          }

          // First, check if this is a sliver
          int nod2bpo[3];
          if(msh.idim == 3 && tdim == 2){
            int iref = msh.fac2ref[ientt];
            nod2bpo[0] = ibins;
            nod2bpo[1] = msh.poi2ebp(ent2pol[1], 2, ientt, iref);
            nod2bpo[2] = msh.poi2ebp(ent2pol[2], 2, ientt, iref);
          }
          bool iflat;
          double meas0;
          meas0 = msh.idim == 2 ? getmeasentP1<2,2>(msh, ent2pol, nod2bpo, norCAD, &iflat, -1)
                :     tdim == 2 ? getmeasentP1<3,2>(msh, ent2pol, nod2bpo, norCAD, &iflat, -1)
                                : getmeasentP1<3,3>(msh, ent2pol, nod2bpo, norCAD, &iflat, -1);

          CPRINTF1("  - inccav pdim {} tdim {} ent {} = {}\n",
                   pdim,tdim,ientt,intAr1(tdim+1,ent2pol));
          CPRINTF1("  - w/ vtol = {} got iflat = {} meas0 = {:15.7e} neighbour = {}\n",
                   msh.param->vtol,iflat,meas0,ienei);

          #if 0
          // Next check geodev
          // Actually not because adding more faces will only damage the cavity further
          // Do this in the future as pre reject, possibly.
          // Also depends on Pk etc. Probably best to leave in cav.
          double nrmal[3];
          if(msh.idim == 3 && pdim < 2){
            // Get the normal in the case we're on an edge in 3D, and get only
            // the correct side.
            int iref = msh.fac2ref[ientt];
            getnorballref<1>(msh,lcent,iref,nrmal);
          }
          #endif
          // ignore ienei < 0 as it could be bdry -> edge remeshing
          if((iflat || meas0 < 0)){
            //if(ienei == -1) return 1;
            //// Cannot be corrected
            //if(ienei < 0){
            //  METRIS_ASSERT(iedge >= 0 && tdim == 2 || iface >= 0 && tdim ==3);
            //  CPRINTF1(" # abort flat no neighbour: meas {:23.15e}\n", meas0);
            //  return 1;
            //}

            if(ienei >= 0){
              lcent.stack(ienei);
              ent2tag(ithread,ienei) = msh.tag[ithread];
              CPRINTF1("   - inccav added entt {} to stack \n", ienei);
            }else{
              // Add the boundary entity, but only if in allowed refs.
              if(tdim == 2){
                int iref = msh.edg2ref[iedge];
                if(msh.ced2tag(ithread,iref) < msh.tag[ithread] && msh.isboundary_edges()){
                  CPRINTF1("   - iface {} -> iedge {} is wrong bdry ref {}\n",ienei,iedge,iref);
                  return 1;
                }
              }else{
                int iref = msh.fac2ref[iface];
                if(msh.cfa2tag(ithread,iref) < msh.tag[ithread]){
                  CPRINTF1("   - itetr {} -> iface {} is wrong bdry ref {}\n",ienei,iedge,iref);
                  return 1;
                }
              }
              restart = true;
            }

            // If a subdim entity was sandwiched here, we need to add it
            // Also true if no neighbour -> add bdry entity.
            if((tdim == 2 && iedge >= 0) || (tdim == 3 && iface >= 0)){
              if(tdim == 2){
                cav.lcedg.stack(iedge);
                msh.edg2tag(ithread,iedge) = msh.tag[ithread];
              }else{
                cav.lcfac.stack(iface);
                msh.fac2tag(ithread,iface) = msh.tag[ithread];
              }
              CPRINTF1("   - inccav added dim {} ent {} to stack \n", tdim - 1,
                       tdim == 2 ? iedge : iface);
            }

            // If this is a face, we must also add the supported tets
            if(tdim == 2 && msh.nelem > 0){
              for(int ii = 0; ii < 2; ii++){
                int ielem = msh.fac2tet(ientt, ii);
                if(ielem < 0) continue;
                if(msh.tet2tag(ithread,ielem) >= msh.tag[ithread]) continue;
                msh.tet2tag(ithread,ielem) = msh.tag[ithread];
                cav.lctet.stack(ielem);
                CPRINTF1("   - inccav added tet {} to stack \n", ielem);
              }
            }

          }

        } // for int inei

      } // for int ientl
    } // for int tdim

    ient0[0] = ient01[0];
    ient0[1] = ient01[1];

  }while(restart);

  return 0;
}




// Increase cavity for Delaunay criterion on ipoin
// Normal only needed in 3D case if cavity has faces
template<class MFT>
int increase_cavity_Delaunay(MeshMetric<MFT> &msh, MshCavity &cav, int tdim,
                             int ngrow, int ithread){

  if(tdim <= 1) return 0;

  GETVDEPTH(msh.param);
  METRIS_ASSERT(tdim <= cav.get_tdim());

  //#ifdef NODELSURF
  //static int nwarn = 0;

  //// Disable surf
  //if(msh.get_tdim() < msh.idim && msh.param->iflag1 == 0){
  //  if(nwarn++ < 10) MPRINTF("## WARNING DELAUNAY SURFACE DISABLED\n");
  //  return 0;
  //}
  //#endif


  //if(msh.get_tdim() == 3)
  //  METRIS_THROW_MSG("TODO: Unit test this for n = 3. Implement gettetfac instead of getfacedg");
  // Simply disable surface Delaunay for now

  int nnmet = (msh.idim * (msh.idim + 1)) / 2;

  msh.tag[ithread]++;

  // Tag point's surface references if any. Filter entities
  aux_taginsrefs(msh,cav,ithread);

  for(int ielem : cav.lctet){
    METRIS_ASSERT(ielem >= 0 && ielem < msh.nelem);
    METRIS_ASSERT(!isdeadent(ielem,msh.tet2poi));
    msh.tet2tag(ithread,ielem) = msh.tag[ithread];
  }

  for(int iface : cav.lcfac){
    METRIS_ASSERT(iface >= 0 && iface < msh.nface);
    METRIS_ASSERT(!isdeadent(iface,msh.fac2poi));
    msh.fac2tag(ithread,iface) = msh.tag[ithread];
  }

  for(int iedge : cav.lcedg){
    METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);
    METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));
    msh.edg2tag(ithread,iedge) = msh.tag[ithread];
  }

  // Actually do only the one dimension, keep this in for the future, maybe.
  //for(int tdim = msh.get_tdim(); tdim <= msh.get_tdim(); tdim++){
  intAr1 &lcent = cav.lcent(tdim);
  //if(lcent.get_n() == 0) continue;
  intAr1 &lcsub = cav.lcent(tdim-1);

  CPRINTF1("-- START increase_cavity_Delaunay {}\n",tdim);
  const intAr2&  ent2ent = msh.ent2ent(tdim);
  const intAr2&  ent2poi = msh.ent2poi(tdim);
        intAr2r& ent2tag = msh.ent2tag(tdim);
        intAr2r& sub2tag = msh.ent2tag(tdim-1);


  double metl[6], lmet[6];
  double *metl_p;
  if(msh.met.getSpace() == MetSpace::Log){
    for(int jj = 0; jj < nnmet; jj++) lmet[jj] = msh.met(cav.ipins,jj);
    if(msh.idim == 2){
      getexpmet_cpy<2>(lmet, metl);
    }else{
      getexpmet_cpy<3>(lmet, metl);
    }
    metl_p = metl;
  }else{
    metl_p = msh.met[cav.ipins];
  }


  int icen0 = 0, icen1 = lcent.get_n();
  for(int igrow = 0; igrow < ngrow || ngrow < 0; igrow++){

    for(int icent = icen0; icent < icen1; icent++){
      INCVDEPTH(msh.param);
      int ientt = lcent[icent];
      for(int jj = 0; jj < tdim + 1; jj++){
        int ienei = ent2ent(ientt,jj);
        if(ienei < 0) continue; // Non manifold skip
        CPRINTF1(" - check ienei = {} Delaunay\n",ienei);
        if(ent2tag(ithread,ienei) >= msh.tag[ithread]){
          CPRINTF1(" - ienei = {} is tagged {} >= {}\n",
                   ienei,ent2tag(ithread,ienei),msh.tag[ithread]);
          continue;
        }

        int isube = -1;
        if(tdim == 2){
          int iref2 = msh.fac2ref[ienei];
          if(msh.cfa2tag(ithread,iref2) < msh.tag[ithread] && msh.isboundary_faces()){
            CPRINTF1(" - ienei = {} is wrong bdry ref {}\n",ienei,iref2);
            continue;
          }
          int isube = msh.fac2edg(ientt,jj);
          if(isube >= 0){
            if(msh.edg2tag(ithread,isube) >= msh.tag[ithread]){
              CPRINTF1(" - iface {} -> iedge {} is tagged, skip\n",ientt,isube);
              continue;
            }
            int iref1 = msh.edg2ref[isube];
            if(msh.ced2tag(ithread,iref1) < msh.tag[ithread] && msh.isboundary_edges()){
              CPRINTF1(" - iface {} -> iedge {} is wrong bdry ref {}\n",ienei,isube,iref1);
              continue;
            }
          }
        }else{
          if(msh.tet2ref[ienei] != msh.tet2ref[ientt]){
            CPRINTF1(" - ienei {} ref = {} != ientt {} ref {} -> skip\n",
                     ienei,msh.tet2ref[ienei],ientt,msh.tet2ref[ientt]);
          }
          int isube = msh.tet2fac(ientt,jj);
          if(isube >= 0){
            if(msh.fac2tag(ithread,isube) >= msh.tag[ithread]){
              CPRINTF1(" - itetr {} -> iface {} is tagged, skip\n",ientt,isube);
              continue;
            }
            int iref1 = msh.fac2ref[isube];
            if(msh.ced2tag(ithread,iref1) < msh.tag[ithread]){
              CPRINTF1(" - itetr {} -> iface {} is wrong bdry ref {}\n",ienei,isube,iref1);
              continue;
            }
          }
        }


        ent2tag(ithread,ienei) = msh.tag[ithread];

        bool isinsph;
        try{
          if(tdim == 2){
            if(msh.idim == 2){
              isinsph = indelsphere<2,2>(msh, msh.coord[cav.ipins], metl_p,
                                        ent2poi[ienei]);
            }else{
              isinsph = indelsphere<3,2>(msh, msh.coord[cav.ipins], metl_p,
                                        ent2poi[ienei]);
            }
          }else{
            isinsph = indelsphere<3,3>(msh, msh.coord[cav.ipins], metl_p,
                                      ent2poi[ienei]);
          }
        }catch(const MetrisExcept& e){
          fmt::print("indelsphere call threw exception\n");
          fmt::print("with ienei = {} nodes {} ipins {}\n",
                     ienei,intAr1(tdim+1,ent2poi[ienei]),cav.ipins);
          //double meas0;
          //bool ivalid = isvalideltP1<3,3>(msh, ienei, NULL, &meas0);
          //fmt::print("elt measure {} valid {}\n",meas0,ivalid);
          throw(e);
        }
        if(isinsph){
          lcent.stack(ienei);
          CPRINTF1(" - stack dim {} ienei {}\n",tdim,ienei);
          if(isube >= 0){
            CPRINTF1(" - stack dim {} subent {}\n",tdim-1,isube);
            sub2tag(ithread,isube) = msh.tag[ithread];
            lcsub.stack(isube);
          }
          if(tdim == 2 && msh.get_tdim() >= 3){
            for(int ii = 0; ii < 2; ii++){
              int isupe = msh.fac2tet(ienei, ii);
              if(isupe < 0) continue;
              if(msh.tet2tag(ithread,isupe) >= msh.tag[ithread]) continue;
              CPRINTF1(" - stack dim {} supent {}\n",tdim+1,isupe);
              msh.tet2tag(ithread,isupe) = msh.tag[ithread];
              cav.lctet.stack(isupe);
            }
          }
        }

      }// for j = 0,tdim
    }// for icent


    icen0 = icen1;
    icen1 = lcent.get_n();
    CPRINTF1(" - del grow {} / {} + {} ent\n",igrow,ngrow,icen1-icen0);
    if(icen1 == icen0) break;
  }// for igrow
  //}// for tdim

  return 0;
}

template int increase_cavity_Delaunay(MeshMetric<MetricFieldAnalytical> &msh,
                                      MshCavity &cav, int tdim, int ngrow, int ithread);
template int increase_cavity_Delaunay(MeshMetric<MetricFieldFE        > &msh,
                                      MshCavity &cav, int tdim, int ngrow, int ithread);





template<class MFT>
int increase_cavity_lenedg(MeshMetric<MFT> &msh, MshCavity &cav,
                           const CavOprOpt &opts,
                           int ipins,int ithrd1, int ithrd2){
  int nprem = 0;
//  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
      nprem = increase_cavity_lenedg0<MFT,gdim>(msh,cav,opts,ipins,ithrd1,ithrd2);
    }}CT_FOR1(gdim);
//  }}CT_FOR1(ideg);
  return nprem;
}

template<class MFT, int gdim>
int increase_cavity_lenedg0(MeshMetric<MFT> &msh, MshCavity &cav,
                            const CavOprOpt &opts,
                            int ipins, int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);


  // Note ipins must be seeded with newbpotopo
  const int pdim_ipins = msh.getpoitdim(ipins);

  //const intAr2 &ent2ent = msh.ent2ent(tdim);
  msh.tag[ithrd1]++;
  for(int tdim = 1; tdim <= msh.get_tdim(); tdim++){
    for(int ientt : cav.lcent(tdim)){
      msh.ent2tag(tdim)(ithrd1,ientt) = msh.tag[ithrd1];
    }
  }

  // Tag point's surface references if any. Filter entities
  aux_taginsrefs(msh,cav,ithrd1);



  intAr1 lbtet(20), lbfac(20), lbedg(20);
  int iopen;

  int nprem = 0;

  int edg2pol[2];
  edg2pol[0] = ipins;
  double sz[2];

  //int ncomp = 0;
  //int ncav0 = lcent.get_n();

  // NB: loop bounds MUST be reevaluated ! don't range-for this
  int cdim = 0;
       if(cav.lctet.get_n() > 0) cdim = 3;
  else if(cav.lcfac.get_n() > 0) cdim = 2;
  else if(cav.lcedg.get_n() > 0) cdim = 1;
  const int nedgl = (cdim*(cdim+1))/2;
  const intAr2 lnoed(nedgl,2,cdim == 2 ? lnoed2[0] : lnoed3[0]);

  intAr1 &lcent = cav.lcent(cdim);
  for(int ii = 0; ii < lcent.get_n(); ii++){
    INCVDEPTH(msh.param);
    int ientt = lcent[ii];
    METRIS_ASSERT_MSG(!isdeadent(ientt, msh.ent2poi(cdim)),
      "entity {} tdim {} is dead", ientt, cdim);


    #if 0
    for(int ifa = 0; ifa < tdim + 1; ifa++){
      int ientn = ent2ent(ientt,ifa);
      if(ientn >= 0){
        if(ent2tag(ithrd1,ientn) >= msh.tag[ithrd1]) continue;
      }
      // Cavity boundary
      // Loop over face nodes
      int kk = -1;
      for(int ii = 0; ii < tdim; ii++){
        // Increment and skip when == to ifa (= not on facet)
        kk += 1 + ((kk + 1) == ifa);
        int ipoin = ent2poi(ientt,kk);
        if(ipoin == ipins) continue;
        if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];

        edg2pol[1] = ipoin;
        double len = getlenedg_geosz<MFT,gdim,1>(msh,edg2pol,sz);
        ncomp++;
        if(len > 1.0/sqrt(2)) continue;


        // Short edge

        if(!opts.allow_remove_points) return -1;
        if constexpr (tdim == 2){
          ball2(msh,ipoin,ientt,lbfac,dum,&iopen,&imani,ithrd2);
        }else{
          ball3(msh,ipoin,ientt,lbfac,&iopen,ithrd2);
        }
        nprem++;
        for(int ient2 : lbfac){
          if(ent2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          ent2tag(ithrd1,ient2) = msh.tag[ithrd1];
          lcent.stack(ient2);
        }
      }
    }
    #else
    for(int inode = 0; inode < cdim + 1; inode++){
      int ipoin = msh.ent2poi(cdim)(ientt,inode);
      if(ipoin == ipins) continue;
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
      msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];

      edg2pol[1] = ipoin;
      double len = getlenedg_geosz<MFT,gdim,1>(msh,edg2pol,sz);


      CPRINTF1(" - check len ipoin {} len = {} <? 1/sqrt(2) {}\n",
                ipoin,len,len <= 1.0/sqrt(2));

      if(len <= 1.0/sqrt(2)){
        int pdim = msh.getpoitdim(ipoin);

        if(pdim < pdim_ipins){
          CPRINTF1(" - short edge and other end has dim {} < {} = dim ipins -> reject\n",
            pdim, pdim_ipins);
          return -1;
        }

        if(pdim == pdim_ipins && !opts.allow_remove_points){
          CPRINTF1(" - short edge and other end has dim {} = {} = dim ipins "
                  "w/ opts.allow_remove_points == false -> reject\n",
                 pdim, pdim_ipins);
          return -1;
        }

        if(pdim > pdim_ipins && !opts.allow_remove_points_superdim){
          CPRINTF1(" - short edge and other end has dim {} > {} = dim ipins "
                  "w/ opts.allow_remove_points_superdim == false -> reject\n",
                 pdim, pdim_ipins);
          return -1;
        }

        lbedg.set_n(0);
        lbfac.set_n(0);
        lbtet.set_n(0);
        // ball can append while avoiding duplicates
        ball(msh, ipoin, lbedg, lbfac, lbtet, &iopen, true, ithrd2);
        //if(cdim == 2){
        //  ball2(msh,ipoin,ientt,lbfac,lbedg,&iopen,&imani,ithrd2);
        //}else{
        //  ball3(msh,ipoin,ientt,lbtet,&iopen,ithrd2);
        //  if(pdim <= 2){
        //    // Also get ball2 of point
        //    int iface = -1;
        //    if(pdim == 1){
        //      int iedge = msh.poi2ent(ipoin,0);
        //      iface = msh.edg2fac[iedge];
        //    }else{
        //      iface = msh.poi2ent(ipoin,0);
        //    }
        //    METRIS_ASSERT(iface >= 0 && iface < msh.nface);
        //    ball2(msh,ipoin,iface,lbfac,lbedg,&iopen,&imani,ithrd2);
        //  }
        //}
        int ncel0 = cav.lctet.get_n();
        int ncfa0 = cav.lcfac.get_n();
        int nced0 = cav.lcedg.get_n();

        bool ifail = false;
        for(int ient2 : lbedg){
          if(msh.edg2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          int iref = msh.edg2ref[ient2];
          if(msh.ced2tag(ithrd1,iref) < msh.tag[ithrd1]){
            ifail = true;
            goto failed;
          }
        }
        for(int ient2 : lbfac){
          if(msh.fac2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          int iref = msh.fac2ref[ient2];
          if(msh.cfa2tag(ithrd1,iref) < msh.tag[ithrd1]){
            ifail = true;
            goto failed;
          }
        }
        for(int ient2 : lbtet){
          if(msh.tet2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          int iref = msh.tet2ref[ient2];
          if(msh.dom2tag(ithrd1,iref) < msh.tag[ithrd1]){
            ifail = true;
            goto failed;
          }
        }

        for(int ient2 : lbedg){
          if(msh.edg2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          msh.edg2tag(ithrd1,ient2) = msh.tag[ithrd1];
          cav.lcedg.stack(ient2);
        }
        for(int ient2 : lbfac){
          if(msh.fac2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          msh.fac2tag(ithrd1,ient2) = msh.tag[ithrd1];
          cav.lcfac.stack(ient2);
        }
        for(int ient2 : lbtet){
          if(msh.tet2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          msh.tet2tag(ithrd1,ient2) = msh.tag[ithrd1];
          cav.lctet.stack(ient2);
        }

        nprem++;

        failed:
        if(ifail){
          CPRINTF1(" - Failed to add point {} to collapse\n",ipoin);
          cav.lcedg.set_n(nced0);
          cav.lcfac.set_n(ncfa0);
          cav.lctet.set_n(ncel0);
        }
      }
    }
    #endif

    //// Control height, only in dimension 2d.
    //if(tdim == 2){

    //}else{
    //  METRIS_THROW_MSG(
    //    "Implement height control in increase_cavity_lenedg 3D");
    //}
  }

  //printf("Debug ncavity init = {} final = {} ncomp = {} \n",ncav0,lcent.get_n(),ncomp);

  return nprem;
}

template int increase_cavity_lenedg(MeshMetric<MetricFieldAnalytical> &msh,
           MshCavity &cav, const CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);
template int increase_cavity_lenedg(MeshMetric<MetricFieldFE        > &msh,
           MshCavity &cav, const CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);


template int increase_cavity_lenedg0<MetricFieldAnalytical,2>(
                            MeshMetric<MetricFieldAnalytical> &msh,
           MshCavity &cav, const CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);
template int increase_cavity_lenedg0<MetricFieldFE        ,2>(
                            MeshMetric<MetricFieldFE        > &msh,
           MshCavity &cav, const CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);
template int increase_cavity_lenedg0<MetricFieldAnalytical,3>(
                            MeshMetric<MetricFieldAnalytical> &msh,
           MshCavity &cav, const CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);
template int increase_cavity_lenedg0<MetricFieldFE        ,3>(
                            MeshMetric<MetricFieldFE        > &msh,
           MshCavity &cav, const CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);


// Increase cavity based on quality
template<class MFT>
int increase_cavity_quality(Mesh<MFT> &msh, MshCavity &cav, int tdim,
                             int ngrow, int ithread){

  if(tdim <= 1) return 0;

  GETVDEPTH(msh.param);
  METRIS_ASSERT(tdim <= cav.get_tdim());

  int nnmet = (msh.idim * (msh.idim + 1)) / 2;

  msh.tag[ithread]++;

  // Tag point's surface references if any. Filter entities
  aux_taginsrefs(msh,cav,ithread);

  for(int ielem : cav.lctet){
    METRIS_ASSERT(ielem >= 0 && ielem < msh.nelem);
    METRIS_ASSERT(!isdeadent(ielem,msh.tet2poi));
    msh.tet2tag(ithread,ielem) = msh.tag[ithread];
  }

  for(int iface : cav.lcfac){
    METRIS_ASSERT(iface >= 0 && iface < msh.nface);
    METRIS_ASSERT(!isdeadent(iface,msh.fac2poi));
    msh.fac2tag(ithread,iface) = msh.tag[ithread];
  }

  for(int iedge : cav.lcedg){
    METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);
    METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));
    msh.edg2tag(ithread,iedge) = msh.tag[ithread];
  }

  // Actually do only the one dimension, keep this in for the future, maybe.
  //for(int tdim = msh.get_tdim(); tdim <= msh.get_tdim(); tdim++){
  intAr1 &lcent = cav.lcent(tdim);
  //if(lcent.get_n() == 0) continue;
  intAr1 &lcsub = cav.lcent(tdim-1);

  CPRINTF1("-- START increase_cavity_quality {}\n",tdim);
  const intAr2&  ent2ent = msh.ent2ent(tdim);
  const intAr2&  ent2poi = msh.ent2poi(tdim);
        intAr2r& ent2tag = msh.ent2tag(tdim);
        intAr2r& sub2tag = msh.ent2tag(tdim-1);

  int icen0 = 0, icen1 = lcent.get_n();
  for(int igrow = 0; igrow < ngrow || ngrow < 0; igrow++){

    // loop over current cavity entities
    for(int icent = icen0; icent < icen1; icent++){
      INCVDEPTH(msh.param);

      int ientt = lcent[icent]; // fetch entity ID

      // loop over neighbors of ientt (equivalently over boundary facets of ientt)
      for(int jj = 0; jj < tdim + 1; jj++){

        int ienei = ent2ent(ientt,jj); // fetch neighbor

        if(ienei < 0) continue; // Non manifold skip

        CPRINTF1(" - check ienei = {} quality\n",ienei);

        // if neighbor tagged means it belongs to cavity so skip
        if(ent2tag(ithread,ienei) >= msh.tag[ithread]){
          CPRINTF1(" - ienei = {} is tagged {} >= {}\n",
                   ienei,ent2tag(ithread,ienei),msh.tag[ithread]);
          continue;
        }


        // some checks that I don't fully understand, but I realize they must stay
        int isube = -1;
        if(tdim == 2){
          int iref2 = msh.fac2ref[ienei];
          if(msh.cfa2tag(ithread,iref2) < msh.tag[ithread] && msh.isboundary_faces()){
            CPRINTF1(" - ienei = {} is wrong bdry ref {}\n",ienei,iref2);
            continue;
          }
          isube = msh.fac2edg(ientt,jj);
          if(isube >= 0){
            if(msh.edg2tag(ithread,isube) >= msh.tag[ithread]){
              CPRINTF1(" - iface {} -> iedge {} is tagged, skip\n",ientt,isube);
              continue;
            }
            int iref1 = msh.edg2ref[isube];
            if(msh.ced2tag(ithread,iref1) < msh.tag[ithread] && msh.isboundary_edges()){
              CPRINTF1(" - iface {} -> iedge {} is wrong bdry ref {}\n",ienei,isube,iref1);
              continue;
            }
          }
        }else{
          if(msh.tet2ref[ienei] != msh.tet2ref[ientt]){
            CPRINTF1(" - ienei {} ref = {} != ientt {} ref {} -> skip\n",
                     ienei,msh.tet2ref[ienei],ientt,msh.tet2ref[ientt]);
                     continue;
          }
          isube = msh.tet2fac(ientt,jj);
          if(isube >= 0){
            if(msh.fac2tag(ithread,isube) >= msh.tag[ithread]){
              CPRINTF1(" - itetr {} -> iface {} is tagged, skip\n",ientt,isube);
              continue;
            }
            int iref1 = msh.fac2ref[isube];
            if(msh.ced2tag(ithread,iref1) < msh.tag[ithread]){
              CPRINTF1(" - itetr {} -> iface {} is wrong bdry ref {}\n",ienei,isube,iref1);
              continue;
            }
          }
        }

        // at this point, we have that ienei is an OUTSIDE element to ientt across local facet jj

        // next thing is to identify all the cavity elements that ienei is neighbor of
        // those plus ienei itself form the current configuration of the patch
        intAr1 entInCurrentConfig;
        for(int kk = 0; kk < tdim + 1; kk++){

          int ieneinei = ent2ent(ienei,kk); // fetch neighbor of ienei

          if(ieneinei < 0) continue; // Non manifold skip

          // if neighbor tagged means it belongs to cavity so stack it
          if(ent2tag(ithread,ieneinei) >= msh.tag[ithread]) entInCurrentConfig.stack(ieneinei);
        }
        entInCurrentConfig.stack(ienei);

        // compute aggregate quality of the current configuration
        dblAr1 quaOfEach;
        double quaCurrentConfig = -1.;

        constexpr AsDeg asdmsh = AsDeg::P1;
        constexpr AsDeg asdmet = AsDeg::P1;
        if (tdim == 2){
          if (msh.idim == 2)
            for (const int iele : entInCurrentConfig) quaOfEach.stack(metqua<MFT,2,2,QuaFun::SizeShape>(msh,asdmsh,asdmet,iele,1.));
          else
            for (const int iele : entInCurrentConfig) quaOfEach.stack(metqua<MFT,3,2,QuaFun::SizeShape>(msh,asdmsh,asdmet,iele,1.));
        }
        else
          for (const int iele : entInCurrentConfig) quaOfEach.stack(metqua<MFT,3,3,QuaFun::SizeShape>(msh,asdmsh,asdmet,iele,1.));

        // set quality of the current configuration as maximum among the elements in the configuration
        // (recall here quality is actually distortion, so max is the worst case)
        for (const double qual : quaOfEach) if (qual > quaCurrentConfig) quaCurrentConfig = qual;

        // now we need to obtain the quality of the "would-be" elements if we add ienei to the cavity
        // for this, we create a new local cavity with all the elements in the current configuration

        MshCavity cav_loc(4,4,2);
        cav_loc.reset();
        cav_loc.ipins = cav.ipins;
        cav_loc.inewp = 1;

        if (tdim == 2){
          // add the faces (elements) in current configuration
          for (const int iface : entInCurrentConfig) cav_loc.lcfac.stack(iface);

          // if mesh is surface in 3D, add the tets incident to the all faces in the current config
          if (msh.idim == 3){
            msh.tag[ithread]++;
            for (const int iface : entInCurrentConfig){
              for (int it = 0; it < 2; it++) {
                int itet = msh.fac2tet(iface, it);
                if (itet < 0) continue;
                if (msh.tet2tag(ithread, itet) == msh.tag[ithread]) continue;
                msh.tet2tag(ithread, itet) = msh.tag[ithread];
                cav_loc.lctet.stack(itet);
              }
            }
          }
        }
        else {
          // add the tets in current configuration
          for (const int iele : entInCurrentConfig) cav_loc.lctet.stack(iele);
        }

        // dry-run cavity operator on local cavity
        // this reconnects local cavity and compute the max quality
        // if we have an invalid element in the resulting reconnection it returns non-zero
        CavOprOpt opts_loc;
        opts_loc.allow_topological_correction = true;
        opts_loc.skip_topo_checks = false;
        opts_loc.dryrun = true;
        opts_loc.qmax_nec = -1;
        opts_loc.qmax_suf = -1;
        opts_loc.qmax_iff = -1;

        CavOprInfo info_loc;
        CavWrkArrs work_loc;

        int ierr = 0;
        CT_FOR0_INC(1, METRIS_MAX_DEG, ideg){ if (msh.curdeg == ideg) {
          ierr = cavity_operator<MFT, ideg>(msh, cav_loc, opts_loc, work_loc, info_loc, ithread);
        }} CT_FOR1(ideg);
        if (ierr != 0){
          // this will catch, among other things, if one of the would-be elements is invalid
          CPRINTF1(" - local cavity reconnection failed for ientt {}, jj {}, ienei {}\n", ientt,jj,ienei);
          continue;
        }

        double quanew = info_loc.qmax_end;

        bool improvesQua = quaCurrentConfig > quanew + 1e-12;

        // add ienei to cavity if that improves quality
        if(improvesQua){
          lcent.stack(ienei);
          ent2tag(ithread,ienei) = msh.tag[ithread];
          CPRINTF1(" - stack dim {} ienei {}\n",tdim,ienei);
          if(isube >= 0){
            CPRINTF1(" - stack dim {} subent {}\n",tdim-1,isube);
            sub2tag(ithread,isube) = msh.tag[ithread];
            lcsub.stack(isube);
          }
          if(tdim == 2 && msh.get_tdim() >= 3){
            for(int ii = 0; ii < 2; ii++){
              int isupe = msh.fac2tet(ienei, ii);
              if(isupe < 0) continue;
              if(msh.tet2tag(ithread,isupe) >= msh.tag[ithread]) continue;
              CPRINTF1(" - stack dim {} supent {}\n",tdim+1,isupe);
              msh.tet2tag(ithread,isupe) = msh.tag[ithread];
              cav.lctet.stack(isupe);
            }
          }
        }

      }// for jj = 0,tdim
    }// for icent

    icen0 = icen1;
    icen1 = lcent.get_n();
    CPRINTF1(" - del grow {} / {} + {} ent\n",igrow,ngrow,icen1-icen0);
    if(icen1 == icen0) break;
  }// for igrow
  //}// for tdim

  return 0;
}

template int increase_cavity_quality(Mesh<MetricFieldAnalytical> &msh,
                                      MshCavity &cav, int tdim, int ngrow, int ithread);
template int increase_cavity_quality(Mesh<MetricFieldFE        > &msh,
                                      MshCavity &cav, int tdim, int ngrow, int ithread);


void aux_taginsrefs(MeshBase &msh, MshCavity &cav, int ithread){
  GETVDEPTH(msh.param);
  METRIS_ASSERT_MSG(ithread >= 0, "ithread = {} < 0", ithread);
  CPRINTF2("-- START aux_taginsrefs\n");
  for(int iedge : cav.lcedg){
    int iref = msh.edg2ref[iedge];
    METRIS_ASSERT(iref >= 0);
    if(msh.ced2tag(ithread,iref) < msh.tag[ithread]){
      CPRINTF3(" - ipins has edge ref {} \n",iref);
    }
    msh.ced2tag(ithread,iref) = msh.tag[ithread];
  }
  for(int iface : cav.lcfac){
    int iref = msh.fac2ref[iface];
    METRIS_ASSERT(iref >= 0);
    if(msh.ced2tag(ithread,iref) < msh.tag[ithread]){
      CPRINTF3(" - ipins has face ref {} \n",iref);
    }
    msh.cfa2tag(ithread,iref) = msh.tag[ithread];
  }
  for(int ielem : cav.lctet){
    int iref = msh.tet2ref[ielem];
    METRIS_ASSERT_MSG(iref >= 0, "ielem = {} invalid iref = {}", ielem, iref);
    if(msh.dom2tag(ithread,iref) < msh.tag[ithread]){
      CPRINTF3(" - ipins has tetra ref {} \n",iref);
    }
    msh.dom2tag(ithread,iref) = msh.tag[ithread];
  }
  #if 0
  for(int ibpoi = msh.poi2bpo[cav.ipins]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
    int bdim = msh.bpo2ibi(ibpoi,1);
    if(bdim == 0) continue;
    int ientt = msh.bpo2ibi(ibpoi,2);
    if(bdim == 1){
      int iref = msh.edg2ref[ientt];
      METRIS_ASSERT(iref >= 0);
      if(msh.ced2tag(ithread,iref) < msh.tag[ithread]){
        CPRINTF1(" - ipins has edge ref {} \n",iref);
      }
      msh.ced2tag(ithread,iref) = msh.tag[ithread];
    }else{
      int iref = msh.fac2ref[ientt];
      METRIS_ASSERT(iref >= 0);
      if(msh.cfa2tag(ithread,iref) < msh.tag[ithread]){
        CPRINTF1(" - ipins has face ref {} \n",iref);
      }
      msh.cfa2tag(ithread,iref) = msh.tag[ithread];
    }
  }
  #endif
}



} // end namespace
