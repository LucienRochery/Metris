//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_swap.hxx"
#include "msh_swap.hxx" // for swapOptions

#include "../Mesh/Mesh.hxx"

#include "../linalg/det.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../utils/aux_misc.hxx"
#include "../low_geo/lenedg.hxx"
#include "../utils/mprintf.hxx"
#include "../low_geo/normal.hxx"
#include "../io_libmeshb.hxx"
#include "../quality/low_metqua.hxx"

#include <limits>

namespace Metris{

template<class MFT>
void accumulate_tet_objective(Mesh<MFT>& msh,
                              const intAr1& region,
                              double& numerator,
                              int& elementCount,
                              double& targetWeight){
  #ifdef STEPDISTANCE
  constexpr QuaFun iquaf = QuaFun::StepDistance;
  #else
  constexpr QuaFun iquaf = QuaFun::SizeShape;
  #endif
  numerator = 0.;
  elementCount = 0;
  targetWeight = 0.;
  for(const int itet : region){
    if(isdeadent(itet,msh.tet2poi)) continue;
    const double quality = metqua<MFT,3,3,iquaf>(
        msh,AsDeg::P1,AsDeg::P1,itet,1.0);
    #ifdef STEPDISTANCE
    if(msh.param->step_distance_cavity_target_average){
      const double weight = step_distance_element_target_weight<MFT,3,3>(
          msh,AsDeg::P1,itet);
      numerator += weight*quality;
      targetWeight += weight;
    }else{
      numerator += quality;
    }
    #else
    numerator += quality;
    #endif
    elementCount++;
  }
}

template<class MFT, int ideg>
int attempt_tet_objective_reconnection(
    Mesh<MFT>& msh,
    MshCavity& cav,
    CavOprOpt& opts,
    CavWrkArrs& work,
    double oldNumerator,
    int oldCount,
    double oldTargetWeight,
    StepDistanceObjectiveState *globalObjective,
    double& qnrm1,
    int ithread){
  opts.dryrun = true;
  opts.qmax_nec = std::numeric_limits<double>::max();
  CavOprInfo probeInfo;
  int ierro = cavity_operator<MFT,ideg>(
      msh,cav,opts,work,probeInfo,ithread);
  qnrm1 = probeInfo.qmax_end;

  bool accepted = probeInfo.objective_element_count_end > 0;
  if(accepted && globalObjective != nullptr){
    accepted = globalObjective->accepts_replacement(
        oldNumerator,oldCount,oldTargetWeight,
        probeInfo.objective_numerator_end,
        probeInfo.objective_element_count_end,
        probeInfo.objective_target_weight_end);
  }else if(accepted){
    accepted = objective_strictly_improves(
            probeInfo.objective_numerator_end,oldNumerator);
  }
  if(!accepted) return ierro == CAV_ERR_DRYFAIL1 ? 0 : ierro;

  opts.dryrun = false;
  CavOprInfo commitInfo;
  ierro = cavity_operator<MFT,ideg>(
      msh,cav,opts,work,commitInfo,ithread);
  opts.dryrun = true;
  qnrm1 = commitInfo.qmax_end;
  if(!commitInfo.done) return ierro;

  if(globalObjective != nullptr){
    globalObjective->replace(
        oldNumerator,oldCount,oldTargetWeight,
        commitInfo.objective_numerator_end,
        commitInfo.objective_element_count_end,
        commitInfo.objective_target_weight_end);
  }
  return -1;
}

// Return > 0 if error
//          0 if nothing done
//         -1 if face -> edge swap done
//         -2 if edge -> face(s) swap done
// tet swaps: try face and edge swaps.
// if lazy, take first quality improving operation.
// Tet swaps can only be quality, not length, based (for now and perhaps ever).
template<class MFT, int ideg>
int swaptetra(Mesh<MFT>& msh, int itetr, swapOptions opt,
             MshCavity &cav, CavWrkArrs &work,
             double *qnrm0_, double *qnrm1_,
             StepDistanceObjectiveState *globalObjective,
             int ithrd1, int ithrd2){
  INCVDEPTH(msh.param);
  constexpr int tdim = 3;
  constexpr int gdim = 3;
  constexpr AsDeg asdmet = AsDeg::P1;

  const bool ilazy = true;
  METRIS_ENFORCE_MSG(ilazy == true, "Non-lazy swaptetra not implemented");


  static int nwarnprt = 0;
  if(nwarnprt++ < 10 && opt.swap_norm != msh.param->opt_power)
    CPRINTF1("## WARNING forced qpnorm = {}, provided {}\n",msh.param->opt_pnorm,opt.swap_norm);

  CavOprOpt opts;
  opts.allow_topological_correction = false;
  opts.skip_topo_checks = false;
  opts.allow_remove_points = false;
  opts.allow_remove_points_superdim = false;
  opts.cache_tetra_quality = true;
  // No cavity extension means no issues should arise.
  opts.skip_topo_checks = true;

  // Precompute initial tetra quality and store in cavity hash table.
  cav.reset(); // this resets the quality hash table
  #ifdef TESTQUALITYALGO
  #ifdef STEPDISTANCE
  double quael = metqua<MFT,gdim,tdim,QuaFun::StepDistance>(msh,AsDeg::P1,asdmet,itetr,1.0);
  #else
  double quael = metqua<MFT,gdim,tdim,QuaFun::SizeShape>(msh,AsDeg::P1,asdmet,itetr,1.0);
  #endif
  #else
  double quael = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,itetr,1.0);
  #endif

  auto key = stup4(msh.tet2poi(itetr,0),msh.tet2poi(itetr,1),
                   msh.tet2poi(itetr,2),msh.tet2poi(itetr,3));
  cav.qtetr[key] = quael;

  CPRINTF1("-- START swaptetra itetr {} quael {} ilazy {}\n", itetr, quael,ilazy);

  opts.dryrun = true;
  opts.qmax_nec = std::numeric_limits<double>::max();


  for(int ifa = 0; ifa < 4; ifa++){
    int ierro = aux_swaptetface<MFT,ideg>(msh, opt, itetr, ifa, quael, cav, opts, work,
                                          qnrm0_, qnrm1_,globalObjective,ithrd1);

    CPRINTF1(" - tried face {} ierro {} got qual {} -> {}\n",ifa, ierro, *qnrm0_, *qnrm1_);
    if(ierro < 0){
      CPRINTF1(" - Accepted tet {} swap face {} quality {} -> {}\n",
               itetr, ifa, *qnrm0_, *qnrm1_);
      return -1;
    }
  }

  if(!msh.param->opt_swap_tet_expensive) return 0;

  for(int ied = 0; ied < 6; ied++){
    int ierro = aux_swaptetedge<MFT,ideg>(msh, opt, itetr, ied, quael, cav, opts, work,
                                          qnrm0_, qnrm1_,globalObjective,
                                          ithrd1, ithrd2);
    CPRINTF1(" - tried edge {} ierro {} got qual {} -> {}\n",ied, ierro, *qnrm0_, *qnrm1_);
    if(ierro < 0){
      CPRINTF1(" - Accepted tet {} swap edge {} quality {} -> {}\n",
               itetr, ied, *qnrm0_, *qnrm1_);
      return -2;
    }
  }

  return 0;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template int swaptetra<MetricFieldAnalytical,n>(Mesh<MetricFieldAnalytical>& msh, \
                                    int itetr, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1,\
                                    StepDistanceObjectiveState*, int ithrd1, int ithrd2);\
template int swaptetra<MetricFieldFE        ,n>(Mesh<MetricFieldFE        >& msh, \
                                    int itetr, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1,\
                                    StepDistanceObjectiveState*, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


// Returns > 0 if error, 0 if no error and no operation, -1 if operation done.
template<class MFT, int ideg>
int aux_swaptetface(Mesh<MFT>& msh, swapOptions opt, int itetr, int ifacl, double quae1,
                    MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,
                    double *qnrm0_, double *qnrm1_,
                    StepDistanceObjectiveState *globalObjective,
                    int ithread){
  GETVDEPTH(msh.param);
  constexpr int tdim = 3;
  constexpr int gdim = 3;
  constexpr AsDeg asdmet = AsDeg::P1;

  static int nwarnprt = 0;
  if(nwarnprt++ < 10 && opt.swap_norm != msh.param->opt_pnorm)
    CPRINTF1("## WARNING forced qpnorm = {}, provided {}\n",msh.param->opt_pnorm,opt.swap_norm);


  double &qnrm0 = *qnrm0_;
  double &qnrm1 = *qnrm1_;

  int itet2 = msh.tet2tet(itetr, ifacl);
  if(itet2 < 0){
    CPRINTF1("# END aux_swaptetface -> no neighbour\n");
    return 1;
  }

  int idom1 = msh.tet2ref[itetr];
  if(idom1 != msh.tet2ref[itet2]){
    CPRINTF1("# END aux_swaptetface: crosses refs {} != {}\n",idom1, msh.tet2ref[itet2]);
    return 1;
  }


  #ifndef NDEBUG
  int iface = getfacglo(msh, msh.tet2poi(itetr,lnofa3[ifacl][0]),
                             msh.tet2poi(itetr,lnofa3[ifacl][1]),
                             msh.tet2poi(itetr,lnofa3[ifacl][2]));
  METRIS_ASSERT_MSG(iface < 0, "# END aux_swaptetface: found face between two same-domn elements\n"
                    " itet1 = {} itet2 = {}, doms {}, {}\n"
                    " iface = {}\n",itetr, itet2, idom1, msh.tet2ref[itet2], iface);
  #endif



  // Compute second tetrahedron quality for computing qnrm0.
  // We cache this quality for the following cases:
  // - The caller had done dryruns and is now effecting the operation -> skip computation
  // - The caller will eventually call an edge-based swap, its qnrm0 will involve
  // this element if the edge is on the face.
  auto key = stup4(msh.tet2poi(itet2,0),msh.tet2poi(itet2,1),
                   msh.tet2poi(itet2,2),msh.tet2poi(itet2,3));
  auto tt = cav.qtetr.find(key);
  double quae2;
  if(tt != cav.qtetr.end()){
    quae2 = tt->second;
    CPRINTF2(" - found cached quality for neighbour tet {}: {}\n",itet2,quae2);
  }else{
    #ifdef TESTQUALITYALGO
    #ifdef STEPDISTANCE
    quae2 = metqua<MFT,gdim,tdim,QuaFun::StepDistance>(msh,AsDeg::P1,asdmet,itet2,1.0);
    #else
    quae2 = metqua<MFT,gdim,tdim,QuaFun::SizeShape>(msh,AsDeg::P1,asdmet,itet2,1.0);
    #endif
    #else
    quae2 = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,itet2,1.0);
    #endif
    cav.qtetr[key] = quae2;
  }

  qnrm0 = MAX(quae1,quae2);
  int ifa2 = -1;
  for(int itmp = 0; itmp < 4; itmp++){
    if(msh.tet2tet(itet2,itmp) != itetr) continue;
    ifa2 = itmp;
    break;
  }
  METRIS_ASSERT(ifa2 >= 0);

  // Get point opposite face in second tetra.
  int ipopp = msh.tet2poi(itet2,ifa2);

  // Do not call reset()! We need the quality hash table to persist between calls.
  cav.lcedg.set_n(0);
  cav.lcfac.set_n(0);
  cav.lctet.set_n(0);

  cav.ipins = ipopp;
  cav.inewp = 0;

  cav.lctet.stack(itetr);
  cav.lctet.stack(itet2);

  double oldNumerator,oldTargetWeight;
  int oldCount;
  accumulate_tet_objective(
      msh,cav.lctet,oldNumerator,oldCount,oldTargetWeight);
  const int ierro = attempt_tet_objective_reconnection<MFT,ideg>(
      msh,cav,opts,work,
      oldNumerator,oldCount,oldTargetWeight,
      globalObjective,qnrm1,ithread);
  CPRINTF1("- aux_swaptetface objective cavity result {} qnrm1 = {}\n",
           ierro,qnrm1);
  return ierro;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template int aux_swaptetface<MetricFieldAnalytical,n>(Mesh<MetricFieldAnalytical>& msh, \
                        swapOptions opt, int itetr, int ifacl, double quae1,\
                        MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,\
                        double *qnrm0_, double *qnrm1_,\
                        StepDistanceObjectiveState*, int ithread);\
template int aux_swaptetface<MetricFieldFE,n>(Mesh<MetricFieldFE>& msh, \
                        swapOptions opt, int itetr, int ifacl, double quae1,\
                        MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,\
                        double *qnrm0_, double *qnrm1_,\
                        StepDistanceObjectiveState*, int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()







// Returns > 0 if error, 0 if no error and no operation, -1 if operation done.
template<class MFT, int ideg>
int aux_swaptetedge(Mesh<MFT>& msh, swapOptions opt, int itetr, int iedgl, double quae1,
                    MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,
                    double *qnrm0_, double *qnrm1_,
                    StepDistanceObjectiveState *globalObjective,
                    int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);
  constexpr int tdim = 3;
  constexpr int gdim = 3;
  constexpr AsDeg asdmet = AsDeg::P1;


  static int nwarnprt = 0;
  if(nwarnprt++ < 10 && opt.swap_norm != msh.param->opt_pnorm)
    CPRINTF1("## WARNING forced qpnorm = {}, provided {}\n",msh.param->opt_pnorm,opt.swap_norm);


  double &qnrm0 = *qnrm0_;
  double &qnrm1 = *qnrm1_;

  int ipoi1 = msh.tet2poi(itetr, lnoed3[iedgl][0]);
  int ipoi2 = msh.tet2poi(itetr, lnoed3[iedgl][1]);

  auto skey = stup2(ipoi1, ipoi2);
  auto tt = msh.edgHshTab.find(skey);
  if(tt != msh.edgHshTab.end()){
    CPRINTF1("# END aux_swaptetedge: edge {} {} is geometric\n",ipoi1,ipoi2);
    return 1;
  }

  // Do not call reset()! We need the quality hash table to persist between calls.
  cav.lcedg.set_n(0);
  cav.lcfac.set_n(0);
  cav.lctet.set_n(0);
  cav.lcedg.allocate(1);
  cav.lcfac.allocate(10);
  cav.lctet.allocate(10);
  cav.inewp = 0;

  int iopen;
  shell(msh, ipoi1, ipoi2, 3, itetr, cav.lcedg, cav.lcfac, cav.lctet, &iopen);
  if(iopen != 0){
    CPRINTF1("# END aux_swaptetedge: open shell\n");
    return 2;
  }
  if(cav.lcedg.get_n() > 0 || cav.lcfac.get_n() > 0){
    CPRINTF1("# END aux_swaptetedge: {} edges and {} faces in shell\n",
      cav.lcedg.get_n(), cav.lcfac.get_n());
    return 4;
  }

  int idom1 = msh.tet2ref[itetr];
  for(int ielem : cav.lctet){
    int idom2 = msh.tet2ref[ielem];
    if(idom2 != idom1){
      CPRINTF1("# END aux_swaptetedge: multi-ref shell {} != {}\n",idom1,idom2);
      return 3;
    }
  }

  // Compute quality of initial cavity. Qualities have been cached by prior calls
  // including with other edges, and calls to aux_swaptetface
  qnrm0 = quae1;
  for(int ielem : cav.lctet){
    if(ielem == itetr) continue;
    double quael;
    auto key = stup4(msh.tet2poi(ielem,0),msh.tet2poi(ielem,1),
                     msh.tet2poi(ielem,2),msh.tet2poi(ielem,3));
    auto tt = cav.qtetr.find(key);
    if(tt != cav.qtetr.end()){
      quael = tt->second;
      CPRINTF2(" - found cached quality for shell tet {}: {}\n",ielem,quael);
    }else{
      #ifdef TESTQUALITYALGO
      #ifdef STEPDISTANCE
      quael = metqua<MFT,gdim,tdim,QuaFun::StepDistance>(msh,AsDeg::P1,asdmet,ielem,1.0);
      #else
      quael = metqua<MFT,gdim,tdim,QuaFun::SizeShape>(msh,AsDeg::P1,asdmet,ielem,1.0);
      #endif
      #else
      quael = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,ielem,1.0);
      #endif
      CPRINTF2(" - computed quality for shell tet {}: {}\n",ielem,quael);
      cav.qtetr[key] = quael;
    }
    qnrm0 = MAX(qnrm0, quael);
  }
  double oldNumerator,oldTargetWeight;
  int oldCount;
  accumulate_tet_objective(
      msh,cav.lctet,oldNumerator,oldCount,oldTargetWeight);

  // A shell to face(s) (= n -> 2(n-2) ) swap can be made in as many ways as
  // there exist triangulations of the n vertices in the shell but not on the edge
  // We can't exhaust them all using our cavity operator, only those that are
  // starrings of one of these vertices. We try these combinations.
  // For the case n = 3 (single face) and n = 4, there are resp. 1/2 unique
  // combinations. Otherwise, loop over all i = 1,n.

  int ierro = 0;
  int nshell = cav.lctet.get_n();
  if(nshell == 3){ // 3 -> 2 swap
    CPRINTF1(" - aux_swaptetedge special case 3 -> 2\n");
    cav.ipins = -1;
    // Grab any point not on the edge.
    for(int iver = 0; iver < 4; iver++){
      if(iver == lnoed3[iedgl][0]) continue;
      if(iver == lnoed3[iedgl][1]) continue;
      cav.ipins = msh.tet2poi(itetr, iver);
      break;
    }
    ierro = attempt_tet_objective_reconnection<MFT,ideg>(
        msh,cav,opts,work,
        oldNumerator,oldCount,oldTargetWeight,
        globalObjective,qnrm1,ithrd2);
    CPRINTF1("- aux_swaptetedge objective cavity result {} qnrm1 = {}\n",
             ierro,qnrm1);
    if(ierro < 0) return -1;
  }else if(nshell == 4){ // 4 -> 4 swap
    // Get the edge opposite iedgl.
    int iedgo = -1;
    for(iedgo = 0; iedgo < 6; iedgo++){
      if(lnoed3[iedgo][0] == lnoed3[iedgl][0]) continue;
      if(lnoed3[iedgo][0] == lnoed3[iedgl][1]) continue;
      if(lnoed3[iedgo][1] == lnoed3[iedgl][0]) continue;
      if(lnoed3[iedgo][1] == lnoed3[iedgl][1]) continue;
      break;
    }
    METRIS_ASSERT(iedgo >= 0 && iedgo < 6);

    CPRINTF1(" - aux_swaptetedge special case 4 -> 4, candidates {} {}\n",
      msh.tet2poi(itetr, lnoed3[iedgo][0]),msh.tet2poi(itetr, lnoed3[iedgo][1]));
    // Two possible configs, using either point.
    for(int icfg = 0; icfg < 2; icfg++){
      cav.ipins = msh.tet2poi(itetr, lnoed3[iedgo][icfg]);
      ierro = attempt_tet_objective_reconnection<MFT,ideg>(
          msh,cav,opts,work,
          oldNumerator,oldCount,oldTargetWeight,
          globalObjective,qnrm1,ithrd2);
      CPRINTF1("- aux_swaptetedge objective cavity result {} qnrm1 = {}\n",
               ierro,qnrm1);
      if(ierro < 0) return -1;
    }
  }else{ // general n -> 2(n-2) swap
    CPRINTF1(" - aux_swaptetedge general {} -> {}\n",nshell, 2*(nshell-2));
    cav.iwrk1.allocate(nshell);
    cav.iwrk1.set_n(0);
    msh.tag[ithrd1]++;
    for(int ielem : cav.lctet){
      for(int iver = 0; iver < 4; iver++){
        int ipoin = msh.tet2poi(ielem,iver);
        if(ipoin == ipoi1) continue;
        if(ipoin == ipoi2) continue;
        if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];
        CPRINTF1(" - {} -> {} swap try point {} from elt {}\n",nshell, 2*(nshell-2), ipoin,ielem);
        cav.ipins = ipoin;
        ierro = attempt_tet_objective_reconnection<MFT,ideg>(
            msh,cav,opts,work,
            oldNumerator,oldCount,oldTargetWeight,
            globalObjective,qnrm1,ithrd2);
        CPRINTF1("- aux_swaptetedge objective cavity result {} qnrm1 = {}\n",
                 ierro,qnrm1);
        if(ierro < 0) return -1;
      }
    }
    return 0;
  }


  return 0; // 0 = nothing done
}

#define BOOST_PP_LOCAL_MACRO(n)\
template int aux_swaptetedge<MetricFieldAnalytical,n>(Mesh<MetricFieldAnalytical>& msh, \
                        swapOptions opt, int itetr, int ifacl, double quae1,\
                        MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,\
                        double *qnrm0_, double *qnrm1_,\
                        StepDistanceObjectiveState*, int ithrd1, int ithrd2);\
template int aux_swaptetedge<MetricFieldFE,n>(Mesh<MetricFieldFE>& msh, \
                        swapOptions opt, int itetr, int ifacl, double quae1,\
                        MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,\
                        double *qnrm0_, double *qnrm1_,\
                        StepDistanceObjectiveState*, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


} // end namespace
