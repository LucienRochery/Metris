//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_insert.hxx"
#include "aux_insert.hxx"
#include "EdgeSeed.hxx"

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
#include "../../quality/low_metqua.hxx"

#include "../../msh_checktopo.hxx"

#include <cmath>

namespace Metris{

#ifdef TESTQUALITYALGO
template<class MFT, QuaFun iquaf>
double completed_p2_element_objective(Mesh<MFT>& msh,
                                      int tdim,
                                      int element){
  METRIS_ASSERT(tdim == 2 || tdim == 3);
  return tdim == 2
      ? metqua<MFT,2,2,iquaf>(msh,AsDeg::Pk,AsDeg::P1,element,1.)
      : metqua<MFT,3,3,iquaf>(msh,AsDeg::Pk,AsDeg::P1,element,1.);
}

template<class MFT, QuaFun iquaf>
void configure_completed_p2_insertion_acceptance(
    Mesh<MFT>& msh,
    const MshCavity& cav,
    int tdim,
    const BadEntHandler& handler,
    CavOprOpt& opts){
  double old_numerator = 0.;
  int old_element_count = 0;
  const intAr2& element_to_point = msh.ent2poi(tdim);
  for(int element : cav.lcent(tdim)){
    if(isdeadent(element,element_to_point)) continue;
    old_numerator += completed_p2_element_objective<MFT,iquaf>(
        msh,tdim,element);
    old_element_count++;
  }

  METRIS_ENFORCE(old_element_count > 0);
  if constexpr(iquaf == QuaFun::StepDistance){
    if(msh.param->step_distance_cavity_target_average){
      const StepDistanceObjectiveState global_objective = tdim == 2
          ? step_distance_global_objective_state<MFT,2,2>(
                msh,AsDeg::Pk,AsDeg::P1)
          : step_distance_global_objective_state<MFT,3,3>(
                msh,AsDeg::Pk,AsDeg::P1);
      opts.accept_completed_elements =
          [&msh,tdim,old_numerator,old_element_count,global_objective]
          (int candidate_tdim, int first_new, int end_new){
            if(candidate_tdim != tdim) return false;
            double new_numerator = 0.;
            int new_element_count = 0;
            const intAr2& candidate_to_point = msh.ent2poi(tdim);
            for(int element = first_new; element < end_new; element++){
              if(isdeadent(element,candidate_to_point)) continue;
              const double value
                  = completed_p2_element_objective<MFT,iquaf>(
                        msh,tdim,element);
              if(!std::isfinite(value)) return false;
              new_numerator += value;
              new_element_count++;
            }
            return new_element_count > 0
                && global_objective.accepts_replacement(
                       old_numerator,old_element_count,old_element_count,
                       new_numerator,new_element_count,new_element_count);
          };
      return;
    }
  }

  opts.accept_completed_elements =
      [&msh,&handler,tdim,old_numerator]
      (int candidate_tdim, int first_new, int end_new){
        if(candidate_tdim != tdim) return false;
        double new_numerator = 0.;
        int new_element_count = 0;
        const intAr2& candidate_to_point = msh.ent2poi(tdim);
        for(int element = first_new; element < end_new; element++){
          if(isdeadent(element,candidate_to_point)) continue;
          const double value = completed_p2_element_objective<MFT,iquaf>(
              msh,tdim,element);
          if(!std::isfinite(value)) return false;
          new_numerator += value;
          new_element_count++;
        }
        return new_element_count > 0
            && handler.checkSuccess(new_numerator,old_numerator);
      };
}
#endif

// Return 0 if done nothing, 1 if error, -1 if done swap
template<class MFT, QuaFun iquaf>
int insertEdge(Mesh<MFT>& msh,
               const EdgeSeed &insertionSeed,
               double lenqua_short_max, // maximum quality (error) a new short edge can have
               bool icollapse,
               MshCavity &cav, CavWrkArrs &work,
               intAr1 &lerro,
               #ifdef TESTQUALITYALGO
               BadEntHandler& handler,
               const bool lengthBased,
               const double worsenPctg,
               #endif
               int ithrd1, int ithrd2){

  int iverb0   = msh.param->iverb;
  int ivdepth0 = msh.param->ivdepth;
  int ierro_cavity;

  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1   != ithrd2);

  CavOprOpt opts;
  CavOprInfo info;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = false;
  opts.dryrun = false;
  opts.allow_remove_points = true;
  opts.allow_remove_points_superdim = true; // For boundary
  opts.qmax_nec = -1;
  opts.qmax_suf = -1;
  opts.qmax_iff = -1;
  bool idbg = false;
  cav.inewp = 1;

  int mgrow = 5;

  CPRINTF1("-- START insertEdge tdimp = {} iseed = {}\n",
           insertionSeed.tdimp,insertionSeed.iseed);

  // cav.print(msh);
  CPRINTF1(" - cavity seed nedge {} nface {} ntetr {}\n",
           cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());
  METRIS_ASSERT(cav.lcedg.get_n() > 0 || cav.lcfac.get_n() > 0 || cav.lctet.get_n() > 0);

  int ierro = 0;
  bool irestart_cav;

  int nced0 = cav.lcedg.get_n();
  int ncfa0 = cav.lcfac.get_n();
  int ncte0 = cav.lctet.get_n();

  // save number of elements in mesh originally, before cavity operator modifies it
  const int nentt0 = msh.nentt(insertionSeed.tdim_adp);

  double algnd[3];

  // Create the point, set info for localization
  //cav.ipins = msh.newpoitopo(insertionSeed.tdimp, insertionSeed.iseed);
  //int ibins = -1;
  //if(msh.isboundary_tdim(insertionSeed.tdimp))
  //  ibins = msh.newbpotopo(cav.ipins,insertionSeed.tdimp,insertionSeed.iseed);

  // Proper surface seeding
  cav.ipins = msh.newpoint(PointType::Vertex, insertionSeed.tdimp, insertionSeed.iseed);

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
  if(nwarnprt++ < 10) PRINTF("## WARNING REMOVE STATI FROM NOCOMP\n");
  #endif
  static std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp;


  CPRINTF1(" - create ipins {} tdim = {} seed {} ref {} icollapse {}\n",cav.ipins,insertionSeed.tdimp,insertionSeed.iseed, insertionSeed.iref,icollapse);

  bool imoved_point = false;

  ierro = aux_bisecPointLen(msh, insertionSeed, msh.poi2bpo[cav.ipins], icollapse, cav);
  if(ierro != 0){
    CPRINTF1(" # Failed aux_bisecPointLen ierro = {}\n",ierro);
    goto cleanup;
  }

  // Seed the cavity properly
  #ifndef NDEBUG
  try{
  #endif
    ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
  #ifndef NDEBUG
  }catch(const MetrisExcept& exc){
    fmt::print("## increase_cavity failed, tdim_adp {} tdimp {} iseed {} iref {}\n",
              insertionSeed.tdim_adp,insertionSeed.tdimp,
              insertionSeed.iseed,insertionSeed.iref);
    throw(exc);
  }
  #endif

  if (ierro > 0) goto cleanup;

  // if boundary insertion, check normal deviation: if fails abort operation
  if (msh.get_tdim() == 3 && insertionSeed.tdimp < 3 && !icollapse){

    constexpr bool usemax = false;
    bool nordevOK = seedCav_nordevOK<usemax>(msh,cav,ithrd1);

    if (!nordevOK){
      ierro = 1;
      goto cleanup;
    }
  }

  if(icollapse){
    #ifdef TESTQUALITYALGO
    ierro = increase_cavity_quality<MFT,iquaf>(msh,cav,insertionSeed.tdim_adp,5,handler,ithrd1);
    if(ierro != 0){
      CPRINTF1(" # increase_cavity_quality could not improve patch\n");
      ierro = INS2D_ERR_NOQUALIMPROV;
      goto cleanup;
    }
    #else
    ierro = collrejcav_lenqua(msh, cav, false, false, false, -1, nocomp, ithrd2);
    if(ierro > 0){
      ierro = INS2D_ERR_SHORTEDG6;
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
    #endif
    goto call_cavity;
  }

  // -- This section only if !icollapse

  #ifdef TESTQUALITYALGO
  if (lengthBased){
    ierro = setCavityInsertion3(msh,cav,opts,insertionSeed,mgrow,lenqua_short_max,nocomp,ithrd1,ithrd2);
    if (ierro != 0) goto cleanup;
    ierro = checkCavityQuality(msh,cav,insertionSeed.tdim_adp,5,handler,worsenPctg,ithrd1);
    if (ierro != 0) ierro = INS2D_ERR_NOQUALIMPROV;
  }
  else{
    ierro = setCavityInsertionQuality<MFT,iquaf>(msh,cav,opts,insertionSeed,mgrow,handler,lenqua_short_max,nocomp,ithrd1,ithrd2);
  }
  #else
  ierro = setCavityInsertion3(msh,cav,opts,insertionSeed,mgrow,lenqua_short_max,nocomp,ithrd1,ithrd2);
  #endif

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

  #ifdef TESTQUALITYALGO
  if(!icollapse && msh.curdeg == 2){
    configure_completed_p2_insertion_acceptance<MFT,iquaf>(
        msh,cav,insertionSeed.tdim_adp,handler,opts);
  }else{
    opts.accept_completed_elements = {};
  }
  #endif

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
    #ifdef TESTQUALITYALGO
    ierro = increase_cavity_quality<MFT,iquaf>(msh,cav,insertionSeed.tdim_adp,5,handler,ithrd1);
    if(ierro != 0){
      CPRINTF1(" - +cav error {}\n",ierro);
      ierro = INS2D_ERR_NOQUALIMPROV;
      goto cleanup;
    }
    #else
    ierro = increase_cavity_Delaunay(msh, cav, insertionSeed.tdim_adp, -1, ithrd1);
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
    #endif

    goto restart_cavity;
  }

  if(ierro > 0) lerro[ierro-1]++;


  if(ierro != 0) ierro = INS2D_ERR_CAVITYOPERATOR;

  if(info.done){

    #ifdef TESTQUALITYALGO
    // tell handler about killed and new entities
    const int tdim_adp = insertionSeed.tdim_adp;

    // killed entities
    const intAr1& deadEntts = cav.lcent(tdim_adp);
    for (const auto& ideadEntt : deadEntts) handler.deadEntts.push_back(ideadEntt);

    // new entities with their qualities
    const int nenttNew = msh.nentt(tdim_adp);
    double difto = 1.;
    for (int ienttNew = nentt0; ienttNew < nenttNew; ienttNew++) {

      double quael;
      if (tdim_adp == 2){
        quael = metqua<MFT,2,2,iquaf>(msh,AsDeg::Pk,AsDeg::P1,ienttNew,difto);
      }
      else {
        quael = metqua<MFT,3,3,iquaf>(msh,AsDeg::Pk,AsDeg::P1,ienttNew,difto);
      }
      handler.affectedEnttsAlive[ienttNew] = quael;
    }

    #endif
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



#ifdef TESTQUALITYALGO
#define INSTANTIATE_INSERT_EDGE(MFT_VAL, QUAFUN_VAL) \
template int insertEdge<MFT_VAL,QUAFUN_VAL>(Mesh<MFT_VAL>& msh, \
                         const EdgeSeed &insertionSeed, \
                         double lenqua_short_max, bool icollapse, \
                         MshCavity &cav, CavWrkArrs &work, \
                         intAr1 &lerro, BadEntHandler& handler, \
                         const bool lengthBased, const double worsenPctg, \
                         int ithrd1, int ithrd2);
#else
#define INSTANTIATE_INSERT_EDGE(MFT_VAL, QUAFUN_VAL) \
template int insertEdge<MFT_VAL,QUAFUN_VAL>(Mesh<MFT_VAL>& msh, \
                         const EdgeSeed &insertionSeed, \
                         double lenqua_short_max, bool icollapse, \
                         MshCavity &cav, CavWrkArrs &work, intAr1 &lerro, \
                         int ithrd1, int ithrd2);
#endif

INSTANTIATE_INSERT_EDGE(MetricFieldAnalytical, QuaFun::SizeShape)
INSTANTIATE_INSERT_EDGE(MetricFieldAnalytical, QuaFun::StepDistance)
INSTANTIATE_INSERT_EDGE(MetricFieldFE, QuaFun::SizeShape)
INSTANTIATE_INSERT_EDGE(MetricFieldFE, QuaFun::StepDistance)

#undef INSTANTIATE_INSERT_EDGE




} // end namespace
