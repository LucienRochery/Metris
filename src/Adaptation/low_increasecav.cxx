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
#include "../cavity/reconnect_geometry.hxx"
#include "../low_geo/validity.hxx"
#include "../smoothing/msh_smooball.hxx"
#include "../Mesh/Mesh.hxx"
#include "../io_libmeshb.hxx"
#include "../smoothing/low_smoolen.hxx"

#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"

#include <cmath>
#include <limits>

//#define NODELSURF

namespace Metris{

namespace {

template<int tdim>
constexpr int simplex_edge_count(){
  static_assert(tdim == 2 || tdim == 3);
  return tdim == 2 ? 3 : 6;
}

template<int tdim>
int simplex_edge_vertex(int edge, int endpoint){
  static_assert(tdim == 2 || tdim == 3);
  if constexpr(tdim == 2) return lnoed2[edge][endpoint];
  return lnoed3[edge][endpoint];
}

template<int tdim>
int quadratic_edge_node(int edge){
  int index[tdim + 1] = {};
  index[simplex_edge_vertex<tdim>(edge,0)] = 1;
  index[simplex_edge_vertex<tdim>(edge,1)] = 1;
  return mul2nod(tdim,index);
}

template<class MFT>
bool evaluate_common_cad_midpoint(const MeshMetric<MFT>& msh,
                                  int point0,
                                  int point1,
                                  int minimum_owner_dimension,
                                  double* coefficient){
  if(!msh.CAD() || msh.getBasis() != FEBasis::Lagrange) return false;

  for(int owner_dimension = minimum_owner_dimension;
      owner_dimension <= 2 && owner_dimension < msh.idim;
      owner_dimension++){
    for(int boundary_point = msh.poi2bpo[point0]; boundary_point >= 0;
        boundary_point = msh.bpo2ibi(boundary_point,3)){
      if(msh.bpo2ibi(boundary_point,1) != owner_dimension) continue;
      const int owner_entity = msh.bpo2ibi(boundary_point,2);
      if(owner_entity < 0 || owner_entity >= msh.nentt(owner_dimension)){
        continue;
      }
      const int owner_reference = msh.ent2ref(owner_dimension)[owner_entity];
      const int other_boundary_point = msh.poi2ebp(
          point1,owner_dimension,-1,owner_reference);
      if(other_boundary_point < 0) continue;

      double parameter[2] = {0.,0.};
      for(int component = 0; component < owner_dimension; component++){
        parameter[component]
            = 0.5*(msh.bpo2rbi(boundary_point,component)
                  + msh.bpo2rbi(other_boundary_point,component));
      }
      const ego owner = owner_dimension == 1
          ? msh.CAD.cad2edg[owner_reference]
          : msh.CAD.cad2fac[owner_reference];
      double result[18];
      const int error = EG_evaluate(owner,parameter,result);
      METRIS_ENFORCE_MSG(error == 0,
          "EG_evaluate failed with error {} for local P2 cavity probe",
          error);
      for(int component = 0; component < msh.idim; component++){
        coefficient[component] = result[component];
      }
      return true;
    }
  }
  return false;
}

template<class MFT, int gdim, int tdim>
void complete_quadratic_cone_element(MeshMetric<MFT>& msh,
                                     const MshCavity& cav,
                                     int source_element,
                                     int candidate_element){
  static_assert(gdim == tdim);
  METRIS_ASSERT(msh.curdeg == 2);
  intAr2& element_to_point = msh.ent2poi(tdim);

  for(int candidate_edge = 0;
      candidate_edge < simplex_edge_count<tdim>(); candidate_edge++){
    const int local0 = simplex_edge_vertex<tdim>(candidate_edge,0);
    const int local1 = simplex_edge_vertex<tdim>(candidate_edge,1);
    const int point0 = element_to_point(candidate_element,local0);
    const int point1 = element_to_point(candidate_element,local1);
    const int coefficient_point
        = element_to_point(candidate_element,
                           quadratic_edge_node<tdim>(candidate_edge));

    bool copied_existing_edge = false;
    if(point0 != cav.ipins && point1 != cav.ipins){
      for(int source_edge = 0;
          source_edge < simplex_edge_count<tdim>(); source_edge++){
        const int source0 = element_to_point(
            source_element,simplex_edge_vertex<tdim>(source_edge,0));
        const int source1 = element_to_point(
            source_element,simplex_edge_vertex<tdim>(source_edge,1));
        if(!((source0 == point0 && source1 == point1)
          || (source0 == point1 && source1 == point0))) continue;
        const int source_coefficient = element_to_point(
            source_element,quadratic_edge_node<tdim>(source_edge));
        for(int component = 0; component < gdim; component++){
          msh.coord(coefficient_point,component)
              = msh.coord(source_coefficient,component);
        }
        copied_existing_edge = true;
        break;
      }
      METRIS_ENFORCE_MSG(copied_existing_edge,
          "Local P2 cavity probe could not find inherited edge {} {} "
          "in source element {}",point0,point1,source_element);
      continue;
    }

    const int other_point = point0 == cav.ipins ? point1 : point0;
    const bool is_parent_endpoint
        = cav.split_edge_points.get_n() >= 2
       && (other_point == cav.split_edge_points[0]
        || other_point == cav.split_edge_points[1]);

    // A CAD curve is authoritative only for a true child of the split edge.
    // A common CAD surface also owns newly created boundary spokes.
    const int minimum_cad_dimension = is_parent_endpoint ? 1 : 2;
    bool initialized = evaluate_common_cad_midpoint(
        msh,cav.ipins,other_point,minimum_cad_dimension,
        msh.coord[coefficient_point]);
    if(!initialized){
      initialized = initialize_quadratic_split_child_coefficient<MFT,gdim>(
          msh,cav,point0,point1,msh.coord[coefficient_point]);
    }
    if(!initialized){
      for(int component = 0; component < gdim; component++){
        msh.coord(coefficient_point,component)
            = 0.5*(msh.coord(point0,component)
                  + msh.coord(point1,component));
      }
    }
  }
}

template<class MFT, int gdim, int tdim>
class TemporaryCompletedP2Cone {
public:
  TemporaryCompletedP2Cone(MeshMetric<MFT>& mesh,
                           const MshCavity& cavity,
                           int source_element,
                           const int* vertices)
    : msh(mesh), original_element_count(mesh.nentt(tdim)),
      original_point_count(mesh.npoin), element(original_element_count){
    static_assert(gdim == tdim);
    METRIS_ASSERT(msh.curdeg == 2);
    msh.set_nentt(tdim,original_element_count + 1);
    for(int local = 0; local < tdim + 1; local++){
      msh.ent2poi(tdim)(element,local) = vertices[local];
    }
    for(int edge = 0; edge < simplex_edge_count<tdim>(); edge++){
      const int point = msh.newpoitopo(PointType::CtrlPt,tdim,element);
      msh.ent2poi(tdim)(element,quadratic_edge_node<tdim>(edge)) = point;
    }
    complete_quadratic_cone_element<MFT,gdim,tdim>(
        msh,cavity,source_element,element);
  }

  TemporaryCompletedP2Cone(const TemporaryCompletedP2Cone&) = delete;
  TemporaryCompletedP2Cone& operator=(const TemporaryCompletedP2Cone&) = delete;

  ~TemporaryCompletedP2Cone(){
    msh.set_nentt(tdim,original_element_count);
    msh.set_npoin(original_point_count);
  }

  int index() const noexcept { return element; }

private:
  MeshMetric<MFT>& msh;
  int original_element_count;
  int original_point_count;
  int element;
};

template<class MFT, int gdim, int tdim>
bool completed_p2_cone_is_valid(MeshMetric<MFT>& msh,
                                const MshCavity& cav,
                                int source_element,
                                const int* vertices){
  TemporaryCompletedP2Cone<MFT,gdim,tdim> candidate(
      msh,cav,source_element,vertices);
  return classify_element_validity<gdim,2>(msh,candidate.index())
      .accepted_conservatively();
}

} // namespace

template<class MFT, int gdim, QuaFun iquaf>
double evaluate_completed_p2_cavity_cone(
    Mesh<MFT>& msh,
    const MshCavity& cav,
    int source_element,
    const int* vertices,
    bool* valid){
  static_assert(gdim == 2 || gdim == 3);
  TemporaryCompletedP2Cone<MFT,gdim,gdim> candidate(
      msh,cav,source_element,vertices);
  const bool is_valid = classify_element_validity<gdim,2>(
      msh,candidate.index()).accepted_conservatively();
  if(valid != nullptr) *valid = is_valid;
  if(!is_valid) return std::numeric_limits<double>::infinity();
  return metqua<MFT,gdim,gdim,iquaf>(
      msh,AsDeg::Pk,AsDeg::P1,candidate.index(),1.0);
}

template<QuaFun iquaf, class MFT, int gdim, int tdim>
double cavity_element_contribution(Mesh<MFT>& msh,
                                   AsDeg asdmet,
                                   int ientt,
                                   double element_value,
                                   double& target_weight_sum){
  if constexpr(iquaf == QuaFun::StepDistance){
    if(msh.param->step_distance_cavity_target_average){
      (void)asdmet;
      (void)ientt;
      target_weight_sum += 1.0;
      return step_distance_region_contribution(
          element_value,1.0,true);
    }
  }
  return element_value;
}

template<QuaFun iquaf, class MFT>
double cavity_region_objective(const Mesh<MFT>& msh,
                               double elemental_sum,
                               int element_count,
                               double target_weight_sum = 0.0){
  (void)target_weight_sum;
  if constexpr(iquaf == QuaFun::StepDistance){
    return step_distance_region_objective(
        elemental_sum,element_count,
        msh.param->step_distance_cavity_target_average);
  }
  return elemental_sum;
}

template<QuaFun iquaf, class MFT>
void cavity_replacement_objectives(
    const Mesh<MFT>& msh,
    const BadEntHandler& handler,
    double old_elemental_sum,
    int old_element_count,
    double new_elemental_sum,
    int new_element_count,
    double old_target_weight_sum,
    double new_target_weight_sum,
    double& old_local_objective,
    double& new_local_objective,
    double& old_global_objective,
    double& new_global_objective){
  old_local_objective = cavity_region_objective<iquaf>(
      msh,old_elemental_sum,old_element_count,old_target_weight_sum);
  new_local_objective = cavity_region_objective<iquaf>(
      msh,new_elemental_sum,new_element_count,new_target_weight_sum);

  if constexpr(iquaf == QuaFun::StepDistance){
    if(msh.param->step_distance_cavity_target_average){
      old_global_objective = step_distance_region_objective(
          handler.getQualitySum(),handler.getQualityCount(),true);
      new_global_objective = step_distance_replaced_region_objective(
          handler.getQualitySum(),old_elemental_sum,new_elemental_sum,
          handler.getQualityCount(),old_element_count,new_element_count);
      return;
    }
  }

  old_global_objective = old_local_objective;
  new_global_objective = new_local_objective;
}

template<QuaFun iquaf, class MFT>
bool cavity_replacement_accepts(
    const Mesh<MFT>& msh,
    const BadEntHandler& handler,
    double old_global_objective,
    double new_global_objective){
  if constexpr(iquaf == QuaFun::StepDistance){
    if(msh.param->step_distance_cavity_target_average){
      return objective_strictly_improves(
          new_global_objective,old_global_objective);
    }
  }
  return handler.checkSuccess(new_global_objective,old_global_objective);
}

template<class MFT, QuaFun iquaf>
bool edge_cavity_length_objective_nonworsening_2d(
    Mesh<MFT>& msh, MshCavity& cav, int ithread,
    double& objective_old, double& objective_new){
  objective_old = 0.0;
  objective_new = 0.0;

  const bool protect_edges = msh.idim == 2
                          && msh.get_tdim() == 2
                          && msh.CAD()
                          && msh.param->adp_line_adapt
                          && cav.lcedg.get_n() > 0;
  if(!protect_edges) return true;

  // Match reconnect_lincav's definition of the old 1D cavity and of its
  // reconnected boundary.  In particular, a non-manifold incidence is kept
  // as a distinct new edge.
  for(const int iedge : cav.lcedg){
    msh.edg2tag(ithread,iedge) = msh.tag[ithread];
    objective_old += metqua1_length<MFT,2,iquaf>(
        msh,msh.edg2poi[iedge]);
  }

  int number_new_edges = 0;
  for(const int iedge : cav.lcedg){
    for(int inei = 0; inei < 2; inei++){
      const int iedge_neighbor = msh.edg2edg(iedge,inei);
      if(iedge_neighbor >= 0
         && msh.edg2tag(ithread,iedge_neighbor) >= msh.tag[ithread]){
        continue;
      }

      const int ipseed = msh.edg2poi(iedge,1-inei);
      if(ipseed == cav.ipins) continue;

      const int edge_points[2] = {cav.ipins,ipseed};
      objective_new += metqua1_length<MFT,2,iquaf>(msh,edge_points);
      number_new_edges++;
    }
  }

  if(number_new_edges == 0) return false;
  const double tolerance = 64.0*std::numeric_limits<double>::epsilon()
                         * MAX(1.0,std::abs(objective_old));
  return objective_new <= objective_old + tolerance;
}


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
template<class MFT, QuaFun iquaf>
int setCavityInsertionQuality(Mesh<MFT>& msh, MshCavity &cav, const CavOprOpt &opts,
                       const EdgeSeed &insertionSeed, int mgrow, BadEntHandler& handler, double lenqua_short_max,
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

  ierro = increase_cavity_quality<MFT,iquaf>(msh,cav,tdim,5,handler,ithrd1);
  if(DOPRINTS2()){
    // writeMeshCavity("insert_cavity1",msh,cav);
  }
  // if(ierro > 0){
  //   CPRINTF1(" # increase_cavity_quality error {}\n",ierro);
  //   return INS2D_ERR_INCCAVDEL;
  // }
  if(ierro < 0){
    CPRINTF1(" # insertion: increase_cavity_quality could not improve patch\n");
    return INS2D_ERR_NOQUALIMPROV;
  }

  return 0;
}

#define INSTANTIATE_SET_CAVITY_QUALITY(MFT_VAL, QUAFUN_VAL) \
template int setCavityInsertionQuality<MFT_VAL,QUAFUN_VAL>( \
    Mesh<MFT_VAL>& msh, MshCavity &cav, const CavOprOpt &opts, \
    const EdgeSeed &insertionSeed, int mgrow, BadEntHandler& handler, \
    double lenqua_short_max, \
    std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp, \
    int ithrd1, int ithrd2);

INSTANTIATE_SET_CAVITY_QUALITY(MetricFieldAnalytical, QuaFun::SizeShape)
INSTANTIATE_SET_CAVITY_QUALITY(MetricFieldAnalytical, QuaFun::StepDistance)
INSTANTIATE_SET_CAVITY_QUALITY(MetricFieldFE, QuaFun::SizeShape)
INSTANTIATE_SET_CAVITY_QUALITY(MetricFieldFE, QuaFun::StepDistance)

#undef INSTANTIATE_SET_CAVITY_QUALITY



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
  // cav.print(msh);


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
              // nordev_tol = ref2nordev(iref,0);
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


            // The full-dimensional reconnection must already satisfy the
            // completed-P2 validity contract before objective-driven growth
            // begins. Embedded boundary faces retain their separate P1
            // compatibility test until the deferred surface certificate is
            // implemented.
            double meas0 = std::numeric_limits<double>::quiet_NaN();
            bool ivalid;
            if(msh.curdeg == 2 && tdim == msh.idim){
              ivalid = msh.idim == 2
                  ? completed_p2_cone_is_valid<MFT,2,2>(
                        msh,cav,ientt,ent2pol)
                  : completed_p2_cone_is_valid<MFT,3,3>(
                        msh,cav,ientt,ent2pol);
            }else{
              ivalid = msh.idim == 2
                  ? isvalideltP1<2,2>(msh,ent2pol,NULL,NULL,&meas0,nordev_tol)
                  : tdim == 2
                  ? isvalideltP1<3,2>(msh,ent2pol,nod2bpo,NULL,&meas0,nordev_tol)
                  : isvalideltP1<3,3>(msh,ent2pol,NULL,NULL,&meas0,nordev_tol);
            }
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
              // cav.print(msh, 10);
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
  // cav.print(msh);

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
template<class MFT, QuaFun iquaf>
int increase_cavity_quality(Mesh<MFT> &msh, MshCavity &cav, int tdim,
                             int ngrow, BadEntHandler& handler, int ithread){

  if(tdim <= 1) return 0;

  GETVDEPTH(msh.param);
  METRIS_ASSERT(tdim <= cav.get_tdim());
  METRIS_ASSERT(tdim == 2 || tdim == 3);
  METRIS_ASSERT_MSG(msh.idim + tdim != 5, "Qual-based algo not yet supported for 2D mesh in 3D space. I guess..., it might work");

  int nnmet = (msh.idim * (msh.idim + 1)) / 2;

  msh.tag[ithread]++;

  #ifdef DEBUGCAV
  std::cout << "INITIAL CAVITY" << std::endl;
  for (int iedge : cav.lcedg){
    std::cout << "iedge = " << iedge << std::endl;
  }
  for (int ifac : cav.lcfac){
    std::cout << "ifac = " << ifac << std::endl;
  }
  for (int itet : cav.lctet){
    std::cout << "itet = " << itet << std::endl;
  }

  int ipdbg =  msh.newpoitopo(PointType::Vertex,-1,-1);
  int ibdbg =  msh.newbpotopo(Vertex{ipdbg},0,ipdbg);
  for(int ii = 0; ii < msh.idim; ii++)
    msh.coord(ipdbg,ii) = msh.coord(cav.ipins,ii);

  writeMesh("mshCavSingDet.meshb",msh);
  writeMeshCavity("cavSingDet_ini.meshb",msh,cav);
  #endif

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

  intAr1 &lcent = cav.lcent(tdim);
  intAr1 &lcsub = cav.lcent(tdim-1);

  const bool ipinsOnBnd = lcsub.get_n() > 0;

  int iedins = -1;
  const bool ipinsOnEdge = cav.lcedg.get_n() > 0;
  if (ipinsOnEdge){
    METRIS_ASSERT(cav.lcedg.get_n() <= 3); // 1 if insertion, and could be 2 or 3 if collapse
    iedins = cav.lcedg[0];
  }

  CPRINTF1("-- START increase_cavity_quality {}\n",tdim);
  const intAr2&  ent2ent = msh.ent2ent(tdim);
        intAr2&  ent2poi = msh.ent2poi(tdim);
        intAr2r& ent2tag = msh.ent2tag(tdim);
        intAr2r& sub2tag = msh.ent2tag(tdim-1);

  // start by computing the quality of the initial cav
  // for both configurations: current and reconnected, and for both entities and subentities

  const int nentt0 = msh.nentt(tdim);
  msh.set_nentt(tdim,nentt0+1);
  const int tmpEntt = nentt0; // index for temporary entity to construct would-be elements
  const int npoiProbe0 = msh.npoin;
  if(msh.curdeg == 2){
    const int number_of_edges = tdim == 2
        ? simplex_edge_count<2>() : simplex_edge_count<3>();
    for(int edge = 0; edge < number_of_edges; edge++){
      const int point = msh.newpoitopo(PointType::CtrlPt,tdim,tmpEntt);
      const int node = tdim == 2
          ? quadratic_edge_node<2>(edge)
          : quadratic_edge_node<3>(edge);
      msh.ent2poi(tdim)(tmpEntt,node) = point;
    }
  }

  const int ipins = cav.ipins;

  double quaCav0 = 0.; // for current config (sum of all qual)
  double quaMax0 = -1.; // worst qual for current config
  double quaCav1 = 0.; // for reconnected config (sum of all qual)
  double quaMax1 = 1.; // worst qual for reconnected config
  int nQuaCav0 = 0;
  int nQuaCav1 = 0;
  double targetWeightCav0 = 0.;
  double targetWeightCav1 = 0.;

  // same as above but for subentities (faces when tdim == 3)
  double quaSub0 = 0.;
  double quaMaxSub0 = -1.;
  double quaSub1 = 0.;
  double quaMaxSub1 = -1.;
  int nQuaSub0 = 0;
  int nQuaSub1 = 0;
  double targetWeightSub0 = 0.;
  double targetWeightSub1 = 0.;
  const int nsube0 = msh.nentt(tdim-1);
  msh.set_nentt(tdim-1,nsube0+1); // in case we need to create temporary subentities
  const int tmpSubEntt = nsube0;

  double difto = 1.;

  const AsDeg growth_geometry_degree
      = msh.curdeg == 1 ? AsDeg::P1 : AsDeg::Pk;

  auto existingObjective2D = [&](int element){
    return metqua<MFT,2,2,iquaf>(
        msh,growth_geometry_degree,AsDeg::P1,element,difto);
  };
  auto existingObjective3D = [&](int element){
    return metqua<MFT,3,3,iquaf>(
        msh,growth_geometry_degree,AsDeg::P1,element,difto);
  };
  auto completedCandidate2D = [&](int source_element,
                                  double& objective) -> bool{
    if(msh.curdeg == 1){
      double measure;
      if(!isvalideltP1<2,2>(msh,tmpEntt,NULL,&measure)) return false;
    }else{
      complete_quadratic_cone_element<MFT,2,2>(
          msh,cav,source_element,tmpEntt);
      if(!classify_element_validity<2,2>(msh,tmpEntt)
              .accepted_conservatively()) return false;
    }
    objective = metqua<MFT,2,2,iquaf>(
        msh,growth_geometry_degree,AsDeg::P1,tmpEntt,difto);
    return std::isfinite(objective);
  };
  auto completedCandidate3D = [&](int source_element,
                                  double& objective) -> bool{
    if(msh.curdeg == 1){
      double measure;
      if(!isvalideltP1<3,3>(msh,tmpEntt,NULL,&measure)) return false;
    }else{
      complete_quadratic_cone_element<MFT,3,3>(
          msh,cav,source_element,tmpEntt);
      if(!classify_element_validity<3,2>(msh,tmpEntt)
              .accepted_conservatively()) return false;
    }
    objective = metqua<MFT,3,3,iquaf>(
        msh,growth_geometry_degree,AsDeg::P1,tmpEntt,difto);
    return std::isfinite(objective);
  };

  auto restoreIpins = [&](const double* coor0, const double* met0,
                          const intAr1& lbpoi0, const dblAr2& rbpoi0){
    for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipins,ii) = coor0[ii];
    for(int ii = 0; ii < nnmet; ii++) msh.met(ipins,ii) = met0[ii];

    for(int ib = 0; ib < lbpoi0.get_n(); ib++){
      int ibpoi = lbpoi0[ib];
      for(int jj = 0; jj < nrbi; jj++) msh.bpo2rbi(ibpoi,jj) = rbpoi0(ib,jj);
    }
  };

  auto retagCavity = [&](){
    aux_taginsrefs(msh,cav,ithread);

    for(int ielem : cav.lctet) msh.tet2tag(ithread,ielem) = msh.tag[ithread];
    for(int iface : cav.lcfac) msh.fac2tag(ithread,iface) = msh.tag[ithread];
    for(int iedge : cav.lcedg) msh.edg2tag(ithread,iedge) = msh.tag[ithread];
  };

  auto stackCavityBoundaryEdges = [&]([[maybe_unused]] int iface){
    if(tdim != 2) return;
    if(!ipinsOnEdge) return;

    for(int kk = 0; kk < 3; kk++){

      int iedge = msh.fac2edg(iface,kk);
      if(iedge < 0) continue;

      int iref = msh.edg2ref[iedge];
      if(msh.ced2tag(ithread,iref) < msh.tag[ithread]) continue;

      if(msh.edg2tag(ithread,iedge) >= msh.tag[ithread]) continue;
      lcsub.stack(iedge);
      sub2tag(ithread,iedge) = msh.tag[ithread];
      CPRINTF1(" - stack dim {} subent {}\n",tdim-1,iedge);
    }
  };

  auto getCavityQuality2D = [&](double& qua0, double& qua1,
                                double& qmax0, double& qmax1,
                                int& nqua0, int& nqua1,
                                double& targetWeight0,
                                double& targetWeight1) -> bool{
    qua0 = 0.;
    qua1 = 0.;
    qmax0 = -1.;
    qmax1 = -1.;
    nqua0 = 0;
    nqua1 = 0;
    targetWeight0 = 0.;
    targetWeight1 = 0.;

    for(const int ienttCav : lcent){

      double qua = existingObjective2D(ienttCav);
      qua0 += cavity_element_contribution<iquaf,MFT,2,2>(
          msh,AsDeg::P1,ienttCav,qua,targetWeight0);
      nqua0++;
      if(qua > qmax0) qmax0 = qua;

      int ent2pol[3];
      for(int jj = 0; jj < 3; jj++){

        const int ienei = ent2ent(ienttCav,jj);
        if(ienei >= 0 && ent2tag(ithread,ienei) >= msh.tag[ithread]) continue;

        if(ipinsOnEdge){
          int iedgeGlobal = msh.fac2edg(ienttCav,jj);
          if(iedgeGlobal >= 0){
            if(msh.edg2ref[iedgeGlobal] == msh.edg2ref[iedins]) continue;
          }
        }

        ent2pol[0] = ipins;
        ent2pol[lnoed2[0][0]] = ent2poi(ienttCav,lnoed2[jj][0]);
        ent2pol[lnoed2[0][1]] = ent2poi(ienttCav,lnoed2[jj][1]);

        if(ent2pol[1] == ipins || ent2pol[2] == ipins) return false;

        ent2poi(tmpEntt,0) = ent2pol[0];
        ent2poi(tmpEntt,1) = ent2pol[1];
        ent2poi(tmpEntt,2) = ent2pol[2];

        if(!completedCandidate2D(ienttCav,qua)) return false;
        qua1 += cavity_element_contribution<iquaf,MFT,2,2>(
            msh,AsDeg::P1,tmpEntt,qua,targetWeight1);
        nqua1++;
        if(qua > qmax1) qmax1 = qua;
      }
    }

    if constexpr(iquaf == QuaFun::StepDistance){
      if(msh.param->step_distance_cavity_target_average){
        return nqua0 > 0 && nqua1 > 0;
      }
    }
    return true;
  };

  #ifdef DEBUGCAV
  std::cout << "GROWING CAVITY" << std::endl;
  #endif

  #ifdef CAVGROWTH

  bool useGlobalCavityGrowth = false;
  if constexpr(iquaf == QuaFun::StepDistance){
    useGlobalCavityGrowth =
        msh.param->step_distance_cavity_target_average;
  }
  StepDistanceObjectiveState cavityGrowthObjective;

  // The growth probe normally compares only the patch changed by adding one
  // outside neighbor. CavityTargetAverage is different because its arithmetic
  // mean changes denominator when the reconnection changes the element count.
  // Build the current tentative reconnection once, then maintain its global
  // numerator/count with each accepted local replacement.
  double quaCav1BeforeGrowth = 0.;
  double quaMax1BeforeGrowth = -1.;
  double targetWeightBeforeGrowth = 0.;
  int nQuaBeforeGrowth = 0;

  for (const int ienttCav : lcent){

    int ent2pol[4];
    for(int jj = 0; jj < tdim + 1; jj++){

      const int ienei = ent2ent(ienttCav,jj);

      if (ienei >= 0 && ent2tag(ithread,ienei) >= msh.tag[ithread]) continue;

      if (tdim == 2){

        if (ipinsOnEdge){
          int iedgeGlobal = msh.fac2edg(ienttCav,jj);
          if (iedgeGlobal >= 0){
            if (msh.edg2ref[iedgeGlobal] == msh.edg2ref[iedins]) continue;
          }
        }

        ent2pol[0] = ipins;
        ent2pol[lnoed2[0][0]] = ent2poi(ienttCav,lnoed2[jj][0]);
        ent2pol[lnoed2[0][1]] = ent2poi(ienttCav,lnoed2[jj][1]);

        if (ent2pol[1] == ipins || ent2pol[2] == ipins) continue;

        ent2poi(tmpEntt,0) = ent2pol[0];
        ent2poi(tmpEntt,1) = ent2pol[1];
        ent2poi(tmpEntt,2) = ent2pol[2];

        double qua;
        const bool isValid = completedCandidate2D(ienttCav,qua);
        METRIS_ENFORCE_MSG(isValid,
            "Initial cavity has invalid completed element when reconnected");
        quaCav1BeforeGrowth +=
            cavity_element_contribution<iquaf,MFT,2,2>(
                msh,AsDeg::P1,tmpEntt,qua,targetWeightBeforeGrowth);
        nQuaBeforeGrowth++;
        if (qua > quaMax1BeforeGrowth) quaMax1BeforeGrowth = qua;
      }
      else{

        int ifaceGlobal = msh.tet2fac(ienttCav, jj);
        if(ifaceGlobal >= 0 && msh.fac2tag(ithread, ifaceGlobal) >= msh.tag[ithread]) continue;

        ent2pol[0] = ipins;
        ent2pol[lnofa3[0][0]] = ent2poi(ienttCav, lnofa3[jj][0]);
        ent2pol[lnofa3[0][1]] = ent2poi(ienttCav, lnofa3[jj][1]);
        ent2pol[lnofa3[0][2]] = ent2poi(ienttCav, lnofa3[jj][2]);

        if(ent2pol[1]==ipins || ent2pol[2]==ipins || ent2pol[3]==ipins) continue;

        ent2poi(tmpEntt,0)=ent2pol[0];
        ent2poi(tmpEntt,1)=ent2pol[1];
        ent2poi(tmpEntt,2)=ent2pol[2];
        ent2poi(tmpEntt,3)=ent2pol[3];

        double qua;
        const bool isValid = completedCandidate3D(ienttCav,qua);
        METRIS_ENFORCE_MSG(isValid,
            "Initial cavity has invalid completed element when reconnected");
        quaCav1BeforeGrowth +=
            cavity_element_contribution<iquaf,MFT,3,3>(
                msh,AsDeg::P1,tmpEntt,qua,targetWeightBeforeGrowth);
        nQuaBeforeGrowth++;
        if(qua > quaMax1BeforeGrowth) quaMax1BeforeGrowth = qua;
      }
    }
  }

  if(useGlobalCavityGrowth){
    double quaCav0BeforeGrowth = 0.;
    int nQuaCav0BeforeGrowth = 0;
    for(const int ienttCav : lcent){
      double qua;
      if(tdim == 2){
        qua = existingObjective2D(ienttCav);
      }else{
        qua = existingObjective3D(ienttCav);
      }
      quaCav0BeforeGrowth += qua;
      nQuaCav0BeforeGrowth++;
    }

    cavityGrowthObjective.numerator =
        handler.getQualitySum()
        - quaCav0BeforeGrowth + quaCav1BeforeGrowth;
    cavityGrowthObjective.element_count =
        handler.getQualityCount()
        - nQuaCav0BeforeGrowth + nQuaBeforeGrowth;
    METRIS_ENFORCE(cavityGrowthObjective.element_count > 0);
    cavityGrowthObjective.target_weight =
        cavityGrowthObjective.element_count;
  }

  #ifdef CAVSMOOTHING
  // High-order cavity-targeted smoothing is Phase 7 step 5. Keep the P1
  // behavior unchanged, but do not enter its P1-only variable path on P2.
  if(msh.curdeg == 1){
    double quaCav1AfterInitialSmoo;
    double quaMax1AfterInitialSmoo;
    double targetWeightAfterInitialSmoo;
    const double objCav1BeforeInitialSmoo =
        cavity_region_objective<iquaf>(
            msh,quaCav1BeforeGrowth,nQuaBeforeGrowth,
            targetWeightBeforeGrowth);
    [[maybe_unused]] double statInitialSmooCav = smoothCavity(
        msh,cav,handler,iquaf,
        quaCav1BeforeGrowth,quaMax1BeforeGrowth,targetWeightBeforeGrowth,
        quaCav1AfterInitialSmoo,quaMax1AfterInitialSmoo,
        targetWeightAfterInitialSmoo,ithread,ithread);
    retagCavity();

    const double objCav1AfterInitialSmoo =
        cavity_region_objective<iquaf>(
            msh,quaCav1AfterInitialSmoo,nQuaBeforeGrowth,
            targetWeightAfterInitialSmoo);
    METRIS_ENFORCE_MSG(
        objCav1AfterInitialSmoo <= objCav1BeforeInitialSmoo,
        "Initial cavity smoothing worsen quality!");

    if(useGlobalCavityGrowth){
      cavityGrowthObjective.numerator +=
          quaCav1AfterInitialSmoo - quaCav1BeforeGrowth;
    }
    quaCav1BeforeGrowth = quaCav1AfterInitialSmoo;
    quaMax1BeforeGrowth = quaMax1AfterInitialSmoo;
    targetWeightBeforeGrowth = targetWeightAfterInitialSmoo;
  }
  #endif

  int icen0 = 0, icen1 = lcent.get_n();
  for(int igrow = 0; igrow < ngrow || ngrow < 0; igrow++){

    #ifdef DEBUGCAV
    std::cout << "  igrow = " << igrow << std::endl;
    #endif

    // loop over current cavity entities
    for(int icent = icen0; icent < icen1; icent++){
      INCVDEPTH(msh.param);

      int ientt = lcent[icent]; // fetch entity ID

      #ifdef DEBUGCAV
      std::cout << "    ientt = " << ientt << std::endl;
      #endif

      // loop over boundary facets of ientt to fetch neighbors
      for(int jj = 0; jj < tdim + 1; jj++){

        int ieneijj = ent2ent(ientt,jj); // fetch neighbor

        #ifdef DEBUGCAV
        std::cout << "      ieneijj = " << ieneijj << std::endl;
        #endif

        if(ieneijj < 0) continue; // cannot grow in this direction, no entt there!

        // if neighbor tagged means it belongs to cavity so skip
        if(ent2tag(ithread,ieneijj) >= msh.tag[ithread]) continue;

        #ifdef DEBUGCAV
        std::cout << "      outside nei identified" << std::endl;
        #endif

        // at this point, we have that ieneijj is an OUTSIDE entt neighbor to ientt across local facet jj

        // next thing is to identify all the cavity elements that ieneijj is neighbor of

        // those plus ieneijj itself form the current configuration of the local patch (think of it as a mini cavity)

        int facetsOfNeiTouchingCav[4];
        for (int kk = 0; kk < 4; kk++) facetsOfNeiTouchingCav[kk] = -1; // initially mark all facets as not touching the cavity and interior to the domain
        // the encoding then will be as follow:
        // facetsOfNeiTouchingCav[kk] == -1, facet kk does not touch the cavity and is interior face in the domain (i.e. not a mesh face on a CAD face)
        // facetsOfNeiTouchingCav[kk] == -2, facet kk does not touch the cavity and is a boundary face of the domain
        // facetsOfNeiTouchingCav[kk] ==  1, facet kk touches the cavity

        for(int kk = 0; kk < tdim + 1; kk++){

          int ieneinei = ent2ent(ieneijj,kk); // fetch neighbor of ienei

          if(ieneinei < 0){
            facetsOfNeiTouchingCav[kk] = -2;
            continue;
          }

          // if neighbor tagged means it belongs to cavity
          if(ent2tag(ithread,ieneinei) >= msh.tag[ithread]) facetsOfNeiTouchingCav[kk] = 1;

        }

        #ifdef DEBUGCAV
        for(int kk = 0; kk < tdim + 1; kk++){
          std::cout << "      face identifier kk = " << facetsOfNeiTouchingCav[kk] << std::endl;
        }
        #endif

        // now we need to compare current local patch configuration to the reconnected configuration

        if (tdim == 2){

          double quaLocalReconnect = 0;
          double quaMaxLocalReconnect = -1;
          double quaLocal = 0;
          double quaMaxLocal = -1;
          int nQuaLocalReconnect = 0;
          int nQuaLocal = 0;
          double targetWeightLocalReconnect = 0.;
          double targetWeightLocal = 0.;

          // loop over faces of ieneijj
          bool invalid = false;
          METRIS_ASSERT(facetsOfNeiTouchingCav[3] == -1);
          for (int iedge = 0; iedge < 3; iedge++){

            if (facetsOfNeiTouchingCav[iedge] > 0){ // this edge is touching cavity, so it adds to the current config

              // fetch the neighbor inside the cavity
              int ienttcav = ent2ent(ieneijj,iedge);
              METRIS_ASSERT(msh.fac2tag(ithread,ienttcav) >= msh.tag[ithread]);

              int iedgeFromInside = -1;
              for (int kk = 0; kk < 3; kk++){
                if (ent2ent(ienttcav,kk) == ieneijj){
                  iedgeFromInside = kk;
                  break;
                }
              }
              METRIS_ASSERT(iedgeFromInside >= 0);

              int ent2pol[3];
              ent2pol[0] = ipins;
              ent2pol[lnoed2[0][0]] = ent2poi(ienttcav,lnoed2[iedgeFromInside][0]);
              ent2pol[lnoed2[0][1]] = ent2poi(ienttcav,lnoed2[iedgeFromInside][1]);
              if (ent2pol[1] == ipins || ent2pol[2] == ipins){

                invalid = true;
                break;
              }

              ent2poi(tmpEntt,0) = ent2pol[0];
              ent2poi(tmpEntt,1) = ent2pol[1];
              ent2poi(tmpEntt,2) = ent2pol[2];

              double qua;
              if(!completedCandidate2D(ienttcav,qua)){
                METRIS_THROW_MSG(
                    "Current completed-P2 cavity configuration is invalid")
              }

              quaLocal += cavity_element_contribution<iquaf,MFT,2,2>(
                  msh,AsDeg::P1,tmpEntt,qua,targetWeightLocal);
              nQuaLocal++;
              if (qua > quaMaxLocal) quaMaxLocal = qua;

            }else{ // this edge is NOT touching cavity, so it might add to the reconnected config if not on same boundary as insertion edge

              // if the insertion point is on a boundary edge, we need to check if iedge is on same boundary, to not attemp to create a triangle in that case
              if (ipinsOnEdge){

                if (facetsOfNeiTouchingCav[iedge] == -2){ // iedge on boundary as well
                  int iedgeGlobal = msh.fac2edg(ieneijj,iedge);
                  if (msh.edg2ref[iedgeGlobal] == msh.edg2ref[iedins]) continue; // skip if on same boundary as insertion edge
                }
              }

              int ent2pol[3];
              ent2pol[0] = ipins;
              ent2pol[lnoed2[0][0]] = ent2poi(ieneijj,lnoed2[iedge][0]);
              ent2pol[lnoed2[0][1]] = ent2poi(ieneijj,lnoed2[iedge][1]);
              if (ent2pol[1] == ipins || ent2pol[2] == ipins){

                invalid = true;
                break;
              }

              ent2poi(tmpEntt,0) = ent2pol[0];
              ent2poi(tmpEntt,1) = ent2pol[1];
              ent2poi(tmpEntt,2) = ent2pol[2];

              double qua;
              if(!completedCandidate2D(ieneijj,qua)){
                invalid = true;
                break;
              }

              quaLocalReconnect +=
                  cavity_element_contribution<iquaf,MFT,2,2>(
                      msh,AsDeg::P1,tmpEntt,qua,
                      targetWeightLocalReconnect);
              nQuaLocalReconnect++;
              if (qua > quaMaxLocalReconnect) quaMaxLocalReconnect = qua;

            }
          }

          if (invalid) continue;
          if constexpr(iquaf == QuaFun::StepDistance){
            if(msh.param->step_distance_cavity_target_average
               && nQuaLocalReconnect == 0) continue;
          }

          double quaLocalInside = quaLocal;

          // outside element
          double quaOutsideEntt = existingObjective2D(ieneijj);

          quaLocal = quaLocalInside
              + cavity_element_contribution<iquaf,MFT,2,2>(
                  msh,AsDeg::P1,ieneijj,quaOutsideEntt,
                  targetWeightLocal);
          nQuaLocal++;
          quaMaxLocal = MAX(quaMaxLocal,quaOutsideEntt);

          const double objLocalReconnect =
              cavity_region_objective<iquaf>(
                  msh,quaLocalReconnect,nQuaLocalReconnect,
                  targetWeightLocalReconnect);
          const double objLocal =
              cavity_region_objective<iquaf>(
                  msh,quaLocal,nQuaLocal,targetWeightLocal);

          double candidateGlobalNumerator = cavityGrowthObjective.numerator;
          int candidateGlobalElementCount =
              cavityGrowthObjective.element_count;
          double objGrowthCurrent = objLocal;
          double objGrowthCandidate = objLocalReconnect;
          bool improveGrowthSum = cavity_replacement_accepts<iquaf>(
              msh,handler,objLocal,objLocalReconnect);
          if(useGlobalCavityGrowth){
            candidateGlobalNumerator +=
                quaLocalReconnect - quaLocal;
            candidateGlobalElementCount +=
                nQuaLocalReconnect - nQuaLocal;
            METRIS_ENFORCE(candidateGlobalElementCount > 0);
            objGrowthCurrent = cavityGrowthObjective.value();
            objGrowthCandidate = step_distance_region_objective(
                candidateGlobalNumerator,candidateGlobalElementCount,true);
            improveGrowthSum = objective_strictly_improves(
                objGrowthCandidate,objGrowthCurrent);
          }
          int incidentCavityElements = 0;
          for(int edge = 0; edge < 3; edge++){
            if(facetsOfNeiTouchingCav[edge] > 0){
              incidentCavityElements++;
            }
          }
          METRIS_ENFORCE(nQuaLocal == incidentCavityElements + 1);
          if(cav.inspect_growth_probe){
            CavityGrowthProbeInfo probe;
            probe.topological_dimension = 2;
            probe.geometry_degree = msh.curdeg;
            probe.outside_element = ieneijj;
            probe.current_cavity_element_count = lcent.get_n();
            probe.incident_cavity_element_count = incidentCavityElements;
            probe.current_local_element_count = nQuaLocal;
            probe.enlarged_local_element_count = nQuaLocalReconnect;
            probe.current_configuration_valid = true;
            probe.enlarged_configuration_valid = true;
            probe.current_objective = objGrowthCurrent;
            probe.enlarged_objective = objGrowthCandidate;
            cav.inspect_growth_probe(probe);
          }
          bool improveLocalMax = true;
          #ifdef IMPROVEMAXQUAL
          improveLocalMax = quaMaxLocalReconnect <= quaMaxLocal;
          #endif

          const double nearMissRel = 5.0e-1;
          bool nearMissGrowthSum =
              objGrowthCandidate <= (1.0 + nearMissRel)*objGrowthCurrent;
          bool acceptCandidate = improveGrowthSum && improveLocalMax;

          #ifdef CAVSMOOTHING
          if (msh.curdeg == 1
          && !acceptCandidate && nearMissGrowthSum && improveLocalMax){

            double coor0[3] = {};
            double met0[6] = {};
            for(int ii = 0; ii < msh.idim; ii++) coor0[ii] = msh.coord(ipins,ii);
            for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipins,ii);

            int nbpoi0 = 0;
            for(int ibpoi = msh.poi2bpo[ipins]; ibpoi >= 0;
                ibpoi = msh.bpo2ibi(ibpoi,3)) nbpoi0++;

            intAr1 lbpoi0(nbpoi0);
            dblAr2 rbpoi0(nbpoi0,nrbi);
            int ib = 0;
            for(int ibpoi = msh.poi2bpo[ipins]; ibpoi >= 0;
                ibpoi = msh.bpo2ibi(ibpoi,3)){
              lbpoi0[ib] = ibpoi;
              for(int jj = 0; jj < nrbi; jj++) rbpoi0(ib,jj) = msh.bpo2rbi(ibpoi,jj);
              ib++;
            }

            int ncent0 = lcent.get_n();
            int nsub0  = lcsub.get_n();

            lcent.stack(ieneijj);
            ent2tag(ithread,ieneijj) = msh.tag[ithread];
            stackCavityBoundaryEdges(ieneijj);

            double quaNear0, quaNear1, quaMaxNear0, quaMaxNear1;
            int nQuaNear0, nQuaNear1;
            double targetWeightNear0, targetWeightNear1;
            bool validNear = getCavityQuality2D(quaNear0,quaNear1,
                                                quaMaxNear0,quaMaxNear1,
                                                nQuaNear0,nQuaNear1,
                                                targetWeightNear0,
                                                targetWeightNear1);
            if(validNear){
              const double quaNear1BeforeSmoo = quaNear1;
              const int nQuaNear1BeforeSmoo = nQuaNear1;
              double quaNearAfterSmoo;
              double quaMaxNearAfterSmoo;
              double targetWeightNearAfterSmoo;
              double statNearSmooCav = smoothCavity(
                  msh,cav,handler,iquaf,
                  quaNear1,quaMaxNear1,targetWeightNear1,
                  quaNearAfterSmoo,quaMaxNearAfterSmoo,
                  targetWeightNearAfterSmoo,ithread,ithread);
              retagCavity();
              targetWeightNear1 = targetWeightNearAfterSmoo;

              if(useGlobalCavityGrowth){
                validNear = getCavityQuality2D(
                    quaNear0,quaNearAfterSmoo,
                    quaMaxNear0,quaMaxNearAfterSmoo,
                    nQuaNear0,nQuaNear1,
                    targetWeightNear0,targetWeightNear1);
              }

              if(validNear){
                const double objNearAfterSmoo =
                    cavity_region_objective<iquaf>(
                        msh,quaNearAfterSmoo,nQuaNear1,
                        targetWeightNear1);
                const double objNear0 =
                    cavity_region_objective<iquaf>(
                        msh,quaNear0,nQuaNear0,targetWeightNear0);
                bool improveNearSum;
                if(useGlobalCavityGrowth){
                  candidateGlobalNumerator +=
                      quaNearAfterSmoo - quaNear1BeforeSmoo;
                  candidateGlobalElementCount +=
                      nQuaNear1 - nQuaNear1BeforeSmoo;
                  METRIS_ENFORCE(candidateGlobalElementCount > 0);
                  objGrowthCurrent = cavityGrowthObjective.value();
                  objGrowthCandidate = step_distance_region_objective(
                      candidateGlobalNumerator,
                      candidateGlobalElementCount,true);
                  improveNearSum = objective_strictly_improves(
                      objGrowthCandidate,objGrowthCurrent);
                }else{
                  objGrowthCurrent = objNear0;
                  objGrowthCandidate = objNearAfterSmoo;
                  improveNearSum =
                      handler.checkSuccess(objNearAfterSmoo,objNear0);
                }
                bool improveNearMax = true;
                #ifdef IMPROVEMAXQUAL
                improveNearMax = quaMaxNearAfterSmoo <= quaMaxNear0;
                #endif

                acceptCandidate = improveNearSum && improveNearMax;
                if(acceptCandidate){
                  double edgeObjectiveNear0 = 0.0;
                  double edgeObjectiveNear1 = 0.0;
                  acceptCandidate =
                      edge_cavity_length_objective_nonworsening_2d<MFT,iquaf>(
                          msh,cav,ithread,
                          edgeObjectiveNear0,edgeObjectiveNear1);
                  if(cav.lcedg.get_n() > 0 && msh.idim == 2
                     && msh.CAD() && msh.param->adp_line_adapt){
                    CPRINTF2(" - growth 1D objective {} -> {}, accepted {}\n",
                             edgeObjectiveNear0,edgeObjectiveNear1,
                             acceptCandidate);
                  }
                }
                CPRINTF1(" - near miss smoothing dim {} ieneijj {} objective {} -> {} accepted {}\n",
                         tdim,ieneijj,objGrowthCurrent,
                         objGrowthCandidate,acceptCandidate);
              }else{
                acceptCandidate = false;
              }
            }

            if(!acceptCandidate){
              restoreIpins(coor0,met0,lbpoi0,rbpoi0);

              for(int ii = nsub0; ii < lcsub.get_n(); ii++){
                int isub = lcsub[ii];
                sub2tag(ithread,isub) = 0;
              }
              lcsub.set_n(nsub0);

              ent2tag(ithread,ieneijj) = 0;
              lcent.set_n(ncent0);
            }
          }
          #endif

          // A directly accepted 2D growth candidate has not yet been stacked.
          // Probe the complete would-be CAD-edge subcavity before committing
          // the growth.  A near-miss candidate was stacked and checked above.
          if(acceptCandidate
             && ent2tag(ithread,ieneijj) < msh.tag[ithread]){
            const int nsub0 = lcsub.get_n();
            stackCavityBoundaryEdges(ieneijj);

            double edgeObjectiveCandidate0 = 0.0;
            double edgeObjectiveCandidate1 = 0.0;
            const bool edgeCandidateOK =
                edge_cavity_length_objective_nonworsening_2d<MFT,iquaf>(
                    msh,cav,ithread,
                    edgeObjectiveCandidate0,edgeObjectiveCandidate1);
            if(cav.lcedg.get_n() > 0 && msh.idim == 2
               && msh.CAD() && msh.param->adp_line_adapt){
              CPRINTF2(" - growth 1D objective {} -> {}, accepted {}\n",
                       edgeObjectiveCandidate0,edgeObjectiveCandidate1,
                       edgeCandidateOK);
            }
            if(!edgeCandidateOK){
              for(int ii = nsub0; ii < lcsub.get_n(); ii++){
                sub2tag(ithread,lcsub[ii]) = 0;
              }
              lcsub.set_n(nsub0);
              acceptCandidate = false;
            }
          }

          if (acceptCandidate){

            if(useGlobalCavityGrowth){
              cavityGrowthObjective.numerator =
                  candidateGlobalNumerator;
              cavityGrowthObjective.element_count =
                  candidateGlobalElementCount;
              cavityGrowthObjective.target_weight =
                  candidateGlobalElementCount;
            }

            // first thing: add outside element to cavity
            if(ent2tag(ithread,ieneijj) < msh.tag[ithread]){
              lcent.stack(ieneijj);
              ent2tag(ithread,ieneijj) = msh.tag[ithread];
              CPRINTF1(" - stack dim {} ieneijj {}\n",tdim,ieneijj);
            }

            // stack also boundary edges that need to be split
            stackCavityBoundaryEdges(ieneijj);

            // update quality of cavity for both configs

            // quaCav0 += quaOutsideEntt;
            // quaCav1 += quaLocalReconnect - quaLocalInside;

            // if (quaOutsideEntt > quaMax0) quaMax0 = quaOutsideEntt;
            // if (quaMaxLocalReconnect > quaMax1) quaMax1 = quaMaxLocalReconnect;

          }

        }
        else{ // tdim == 3

          double quaLocalReconnect = 0;
          double quaMaxLocalReconnect = -1;
          double quaLocal = 0;
          double quaMaxLocal = -1;
          int nQuaLocalReconnect = 0;
          int nQuaLocal = 0;
          double targetWeightLocalReconnect = 0.;
          double targetWeightLocal = 0.;

          bool worsenFaces = false;
          double quaSub0Backup = quaSub0;
          double quaSub1Backup = quaSub1;
          double quaMaxSub0Backup = quaMaxSub0;
          double quaMaxSub1Backup = quaMaxSub1;

          #ifdef DEBUGCAV
          std::cout << "      LOOP OVER FACES" << std::endl;
          #endif

          // loop over faces of ieneijj
          bool invalid = false;
          for (int iface = 0; iface < 4; iface++){

            #ifdef DEBUGCAV
            std::cout << "        iface = " << iface << std::endl;
            #endif

            if (facetsOfNeiTouchingCav[iface] > 0){ // this face is touching cavity, so it adds to the current config

              #ifdef DEBUGCAV
              std::cout << "        touching cav" << std::endl;
              #endif

              // fetch the neighbor inside the cavity
              int ienttcav = ent2ent(ieneijj,iface);
              METRIS_ASSERT(msh.tet2tag(ithread,ienttcav) >= msh.tag[ithread]);

              int ifaceFromInside = -1;
              for (int kk = 0; kk < 4; kk++){
                if (ent2ent(ienttcav,kk) == ieneijj){
                  ifaceFromInside = kk;
                  break;
                }
              }
              METRIS_ASSERT(ifaceFromInside >= 0);

              int ent2pol[4];
              ent2pol[0] = ipins;
              ent2pol[lnofa3[0][0]] = ent2poi(ienttcav,lnofa3[ifaceFromInside][0]);
              ent2pol[lnofa3[0][1]] = ent2poi(ienttcav,lnofa3[ifaceFromInside][1]);
              ent2pol[lnofa3[0][2]] = ent2poi(ienttcav,lnofa3[ifaceFromInside][2]);
              if (ent2pol[1] == ipins || ent2pol[2] == ipins || ent2pol[3]== ipins) continue;

              ent2poi(tmpEntt,0) = ent2pol[0];
              ent2poi(tmpEntt,1) = ent2pol[1];
              ent2poi(tmpEntt,2) = ent2pol[2];
              ent2poi(tmpEntt,3) = ent2pol[3];

              double qua;
              const bool isValid = completedCandidate3D(ienttcav,qua);
              METRIS_ENFORCE_MSG(isValid,
                  "Current completed-P2 cavity configuration is invalid");

              quaLocal += cavity_element_contribution<iquaf,MFT,3,3>(
                  msh,AsDeg::P1,tmpEntt,qua,targetWeightLocal);
              nQuaLocal++;
              if (qua > quaMaxLocal) quaMaxLocal = qua;

            }else{ // this face is NOT touching cavity, so it might add to the reconnected config or else we need to see if it affects faces cavity

              #ifdef DEBUGCAV
              std::cout << "        not touching cav" << std::endl;
              #endif

              bool coneFace = true;
              if (facetsOfNeiTouchingCav[iface] == -2){ // facet is on boundary

                // we check first if the insertion point is on a boundary
                if (ipinsOnBnd){ // insertion point on boundary

                  int ifaceOutside = msh.tet2fac(ieneijj,iface); // notice ifaceOutside is the face OUTISDE the faces cavity
                  METRIS_ASSERT(ifaceOutside >= 0);

                  // next we check if our boundary face reference is the same as one of the existing faces in the cavity
                  for (int ifaceCav : lcsub){
                    if (msh.fac2ref[ifaceOutside] == msh.fac2ref[ifaceCav]){
                      coneFace = false; // this means the facet shares boundary with the insertion point -> do not cone it
                      break;
                    }
                  }

                  #ifdef CHECKSUBENTTQUAL
                  // if we found out the facet shares boundary with the insertion point, we do not cone the facet
                  // instead, we need to compare qualities of current faces vs the would-be faces if we add this facet to faces cavity (lcsub)
                  if (!coneFace){

                    double quaFaceLocalReconnect = 0.;
                    double quaMaxFaceLocalReconnect = 0.;
                    double quaFaceLocal = 0.;
                    double quaMaxFaceLocal = 0.;
                    for (int iedge = 0; iedge < 3; iedge++){

                      // skip if this edge is in the same CAD edge as insertion edge
                      int iedgeGlobal = msh.fac2edg(ifaceOutside,iedge);
                      if (iedgeGlobal >= 0 && msh.edg2ref[iedgeGlobal] == msh.edg2ref[iedins]) continue;

                      int ifacenei = msh.fac2fac(ifaceOutside,iedge);

                      #ifndef NDEBUG
                      bool touchesCav = false;
                      #endif

                      if (ifacenei >=0 && msh.fac2tag(ithread,ifacenei) >= msh.tag[ithread]){ // this edge touches the faces cavity -> contributes to current config

                        #ifndef NDEBUG
                        touchesCav = true;
                        #endif

                        // ifacenei is an INSIDE face of faces cavity
                        int ifaceInside = ifacenei; // just for more clarity

                        METRIS_ASSERT(msh.fac2ref[ifaceInside] == msh.fac2ref[ifaceOutside]); // at this point both faces must be in SAME CAD face

                        // need to find the inside edge that leads to ifaceOutside
                        int iedgeFromInside = -1;
                        for (int ll = 0; ll < 3; ll++){
                          if (fac2fac(ifaceInside,ll) == ifaceOutside){
                            iedgeFromInside = ll;
                            break;
                          }
                        }
                        METRIS_ASSERT(iedgeFromInside >= 0);

                        int ent2pol[3];
                        ent2pol[0] = ipins;
                        ent2pol[lnoed2[0][0]] = fac2poi(ifaceInside,lnoed2[iedgeFromInside][0]);
                        ent2pol[lnoed2[0][1]] = fac2poi(ifaceInside,lnoed2[iedgeFromInside][1]);

                        fac2poi(tmpSubEntt,0) = ent2pol[0];
                        fac2poi(tmpSubEntt,1) = ent2pol[1];
                        fac2poi(tmpSubEntt,2) = ent2pol[2];

                        msh.fac2ref[tmpSubEntt] = msh.fac2ref[ifaceInside];

                        // no need to check validity, current cavity should be valid at all times by construction

                        #ifndef NDEBUG
                        // put this for debug build anyways
                        double meas;
                        bool isValid = isvalideltP1<3,2>(msh, tmpSubEntt, NULL, &meas);
                        METRIS_ASSERT_MSG(isValid, "Current cavity config invalid. Shouldn't ever happen");
                        #endif

                        double qua = metqua<MFT,3,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,tmpSubEntt,difto);
                        quaFaceLocal += qua;
                        if (qua > quaMaxFaceLocal) quaMaxFaceLocal = qua;

                      }
                      else{ // this edge does not touch the faces cavity -> contributes to reconnection

                        int ent2pol[3];
                        ent2pol[0] = ipins;
                        ent2pol[lnoed2[0][0]] = fac2poi(ifaceOutside,lnoed2[iedge][0]);
                        ent2pol[lnoed2[0][1]] = fac2poi(ifaceOutside,lnoed2[iedge][1]);

                        fac2poi(tmpSubEntt,0) = ent2pol[0];
                        fac2poi(tmpSubEntt,1) = ent2pol[1];
                        fac2poi(tmpSubEntt,2) = ent2pol[2];

                        msh.fac2ref[tmpSubEntt] = msh.fac2ref[ifaceOutside];

                        double meas;
                        if(!isvalideltP1<3,2>(msh, tmpSubEntt, NULL, &meas)){
                          invalid = true;
                          break; // iedge loop (edges of the outside face)
                        }

                        double qua = metqua<MFT,3,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,tmpSubEntt,difto);
                        quaFaceLocalReconnect += qua;
                        if (qua > quaMaxFaceLocalReconnect) quaMaxFaceLocalReconnect = qua;

                      } // if edge touches cavity else not
                    } // for iedge
                    METRIS_ASSERT_MSG(touchesCav,"Boundary face of outside cavity element is not neighbor of any face in cavity!");

                    if (invalid) break; // iface loop (faces of outside tet ieneijj)

                    double quaFaceLocalInside = quaFaceLocal;
                    double quaOutsideFace = metqua<MFT,3,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,ifaceOutside,difto);
                    quaFaceLocal = quaFaceLocalInside + quaOutsideFace;
                    quaMaxFaceLocal = MAX(quaMaxFaceLocal,quaOutsideFace);

                    // here we already check if this operation improves the faces quality. if not, we skip the tet ieneijj

                    if (quaFaceLocalReconnect < quaFaceLocal){
                      // NOTE: might seem risky to modify these without knowing if the volume will also imrove, but we have the backup values in case the volume gets worse

                      quaSub0 += quaOutsideFace;
                      quaSub1 += quaFaceLocalReconnect - quaFaceLocalInside;

                      if (quaMaxFaceLocal > quaMaxSub0) quaMaxSub0 = quaMaxFaceLocal;
                      if (quaMaxFaceLocalReconnect > quaMaxSub1) quaMaxSub1 = quaMaxFaceLocalReconnect;
                    }
                    else{
                      worsenFaces = true;
                    }
                  } // if !coneFace

                  #endif
                } // if insertion point in boundary

                // if insertion point not in boundary, then it is safe to cone the facet, nothing else to do here

              } // if facet is on boundary

              #ifdef DEBUGCAV
              std::cout << "        coneFace = " << coneFace << std::endl;
              std::cout << "        facetsOfNeiTouchingCav[iface] = " << facetsOfNeiTouchingCav[iface] << std::endl;
              #endif

              if (facetsOfNeiTouchingCav[iface] == -1 || coneFace){ // interior facet or not sharing boundary with ipins

                int ent2pol[4];
                ent2pol[0] = ipins;
                ent2pol[lnofa3[0][0]] = ent2poi(ieneijj,lnofa3[iface][0]);
                ent2pol[lnofa3[0][1]] = ent2poi(ieneijj,lnofa3[iface][1]);
                ent2pol[lnofa3[0][2]] = ent2poi(ieneijj,lnofa3[iface][2]);
                if (ent2pol[1] == ipins || ent2pol[2] == ipins || ent2pol[3]== ipins) continue;

                ent2poi(tmpEntt,0) = ent2pol[0];
                ent2poi(tmpEntt,1) = ent2pol[1];
                ent2poi(tmpEntt,2) = ent2pol[2];
                ent2poi(tmpEntt,3) = ent2pol[3];

                #ifdef DEBUGCAV
                std::cout << "          probing new element with outside face " << iface << std::endl;
                for (int oo = 0; oo < tdim+1; oo++){
                  std::cout << "          point " << oo << ": " << ent2pol[oo] << std::endl;
                }
                #endif

                double qua;
                if(!completedCandidate3D(ieneijj,qua)){
                  invalid = true;
                  break; // iface loop (faces of outside tet ieneijj)
                }


                quaLocalReconnect +=
                    cavity_element_contribution<iquaf,MFT,3,3>(
                        msh,AsDeg::P1,tmpEntt,qua,
                        targetWeightLocalReconnect);
                nQuaLocalReconnect++;
                if (qua > quaMaxLocalReconnect) quaMaxLocalReconnect = qua;
              }
            }
          } // for iface (faces of outside tet ieneijj)

          if (invalid){
            quaSub0 = quaSub0Backup;
            quaSub1 = quaSub1Backup;
            quaMaxSub0 = quaMaxSub0Backup;
            quaMaxSub1 = quaMaxSub1Backup;
            continue; // jj loop, neighbors of ientt (current entity in cavity lcent[icent])
          }

          #ifdef CHECKSUBENTTQUAL
          if (worsenFaces){
            quaSub0 = quaSub0Backup;
            quaSub1 = quaSub1Backup;
            quaMaxSub0 = quaMaxSub0Backup;
            quaMaxSub1 = quaMaxSub1Backup;
            continue; // jj loop, neighbors of ientt (current entity in cavity lcent[icent])
          }
          #endif

          if constexpr(iquaf == QuaFun::StepDistance){
            if(msh.param->step_distance_cavity_target_average
               && nQuaLocalReconnect == 0) continue;
          }

          double quaLocalInside = quaLocal;

          // outside tet
          double quaOutsideEntt = existingObjective3D(ieneijj);

          quaLocal = quaLocalInside
              + cavity_element_contribution<iquaf,MFT,3,3>(
                  msh,AsDeg::P1,ieneijj,quaOutsideEntt,
                  targetWeightLocal);
          nQuaLocal++;
          quaMaxLocal = MAX(quaMaxLocal,quaOutsideEntt);

          const double objLocalReconnect =
              cavity_region_objective<iquaf>(
                  msh,quaLocalReconnect,nQuaLocalReconnect,
                  targetWeightLocalReconnect);
          const double objLocal =
              cavity_region_objective<iquaf>(
                  msh,quaLocal,nQuaLocal,targetWeightLocal);
          double candidateGlobalNumerator = cavityGrowthObjective.numerator;
          int candidateGlobalElementCount =
              cavityGrowthObjective.element_count;
          double objGrowthCurrent = objLocal;
          double objGrowthCandidate = objLocalReconnect;
          bool improveGrowthSum = cavity_replacement_accepts<iquaf>(
              msh,handler,objLocal,objLocalReconnect);
          if(useGlobalCavityGrowth){
            candidateGlobalNumerator +=
                quaLocalReconnect - quaLocal;
            candidateGlobalElementCount +=
                nQuaLocalReconnect - nQuaLocal;
            METRIS_ENFORCE(candidateGlobalElementCount > 0);
            objGrowthCurrent = cavityGrowthObjective.value();
            objGrowthCandidate = step_distance_region_objective(
                candidateGlobalNumerator,
                candidateGlobalElementCount,true);
            improveGrowthSum = objective_strictly_improves(
                objGrowthCandidate,objGrowthCurrent);
          }
          int incidentCavityElements = 0;
          for(int face = 0; face < 4; face++){
            if(facetsOfNeiTouchingCav[face] > 0){
              incidentCavityElements++;
            }
          }
          METRIS_ENFORCE(nQuaLocal == incidentCavityElements + 1);
          if(cav.inspect_growth_probe){
            CavityGrowthProbeInfo probe;
            probe.topological_dimension = 3;
            probe.geometry_degree = msh.curdeg;
            probe.outside_element = ieneijj;
            probe.current_cavity_element_count = lcent.get_n();
            probe.incident_cavity_element_count = incidentCavityElements;
            probe.current_local_element_count = nQuaLocal;
            probe.enlarged_local_element_count = nQuaLocalReconnect;
            probe.current_configuration_valid = true;
            probe.enlarged_configuration_valid = true;
            probe.current_objective = objGrowthCurrent;
            probe.enlarged_objective = objGrowthCandidate;
            cav.inspect_growth_probe(probe);
          }
          bool improveLocalMax = true;
          #ifdef IMPROVEMAXQUAL
          improveLocalMax = quaMaxLocalReconnect <= quaMaxLocal;
          #endif

          if (improveGrowthSum && improveLocalMax){

            if(useGlobalCavityGrowth){
              cavityGrowthObjective.numerator =
                  candidateGlobalNumerator;
              cavityGrowthObjective.element_count =
                  candidateGlobalElementCount;
              cavityGrowthObjective.target_weight =
                  candidateGlobalElementCount;
            }

            #ifdef DEBUGCAV
            std::cout << "      ADDING NEI = " << ieneijj << std::endl;
            #endif

            // first thing: add outside tet to cavity
            lcent.stack(ieneijj);
            ent2tag(ithread,ieneijj) = msh.tag[ithread];

            #ifdef DEBUGCAV
            std::string cavName = "cavSingDet_ntet_" + std::to_string(cav.lctet.get_n()) + ".meshb";
            writeMeshCavity(cavName,msh,cav);
            #endif

            CPRINTF1(" - stack dim {} ieneijj {}\n",tdim,ieneijj);
            // if(isube >= 0){
            //   CPRINTF1(" - stack dim {} subent {}\n",tdim-1,isube);
            //   sub2tag(ithread,isube) = msh.tag[ithread];
            //   lcsub.stack(isube);
            // }

            // stack also boundary faces that need to be split
            if (ipinsOnBnd){
              for (int kk = 0; kk < 4; kk++){

                int iface = msh.tet2fac(ieneijj,kk);
                if (iface >= 0) {
                  for (int ifaceCav : lcsub){
                    if (msh.fac2ref[ifaceCav] == msh.fac2ref[iface] && msh.fac2tag(ithread,iface) < msh.tag[ithread]){

                      lcsub.stack(iface);
                      msh.fac2tag(ithread,iface) = msh.tag[ithread];
                      break;
                    }
                  }
                }

              }
            }

            // update quality of cavity for both configs

            // quaCav0 += quaOutsideEntt;
            // quaCav1 += quaLocalReconnect - quaLocalInside;

            // if (quaOutsideEntt > quaMax0) quaMax0 = quaOutsideEntt;
            // if (quaMaxLocalReconnect > quaMax1) quaMax1 = quaMaxLocalReconnect;

          } // if improves quality

        } // if tdim == 2 else 3
      }// for jj neighbor of icent
    } // for icent
    icen0 = icen1;
    icen1 = lcent.get_n();
    CPRINTF1(" - del grow {} / {} + {} ent\n",igrow,ngrow,icen1-icen0);
    if(icen1 == icen0) break; // igrow loop
  } // for igrow
  #endif

  // now compute final state of cavity: compare current config vs reconnected config
  // first the entities
  for (const int ienttCav : lcent){

    // first add entt qual for current config
    double qua;
    if (tdim == 2) qua = existingObjective2D(ienttCav);
    else           qua = existingObjective3D(ienttCav);
    if(tdim == 2){
      quaCav0 += cavity_element_contribution<iquaf,MFT,2,2>(
          msh,AsDeg::P1,ienttCav,qua,targetWeightCav0);
    }else{
      quaCav0 += cavity_element_contribution<iquaf,MFT,3,3>(
          msh,AsDeg::P1,ienttCav,qua,targetWeightCav0);
    }
    nQuaCav0++;
    if (qua > quaMax0) quaMax0 = qua;

    // now identy boundary facets in this element to construct
    // and compute quality of reconnected config

    int ent2pol[4];
    for(int jj = 0; jj < tdim + 1; jj++){

      const int ienei = ent2ent(ienttCav,jj);

      // if neighbor tagged it is in cavity -> skip
      if (ienei >= 0 && ent2tag(ithread,ienei) >= msh.tag[ithread]) continue;

      // at this point, facet jj is at boundary of the cavity
      // need to probe quality of the reconnected element

      if (tdim == 2){

        if (ipinsOnEdge){

          // check that facet jj, if also on boundary, is not in same boundary as the insertion edge. otherwise the new triangle would be flat
          int iedgeGlobal = msh.fac2edg(ienttCav,jj);
          if (iedgeGlobal >= 0){
            if (msh.edg2ref[iedgeGlobal] == msh.edg2ref[iedins]) continue;
          }
        }

        // put together new triangle
        ent2pol[0] = ipins;
        ent2pol[lnoed2[0][0]] = ent2poi(ienttCav,lnoed2[jj][0]);
        ent2pol[lnoed2[0][1]] = ent2poi(ienttCav,lnoed2[jj][1]);

        if (ent2pol[1] == ipins || ent2pol[2] == ipins) continue;

        ent2poi(tmpEntt,0) = ent2pol[0];
        ent2poi(tmpEntt,1) = ent2pol[1];
        ent2poi(tmpEntt,2) = ent2pol[2];

        const bool isValid = completedCandidate2D(ienttCav,qua);
        METRIS_ENFORCE_MSG(isValid,
            "Final tentative cavity has an invalid completed-P2 element");
        quaCav1 += cavity_element_contribution<iquaf,MFT,2,2>(
            msh,AsDeg::P1,tmpEntt,qua,targetWeightCav1);
        nQuaCav1++;
        if (qua > quaMax1) quaMax1 = qua;
      }
      else{ // tdim == 3

        // if boundary face itself is in cavity (tagged), it will be split => no single cone tet
        int ifaceGlobal = msh.tet2fac(ienttCav, jj);
        if(ifaceGlobal >= 0 && msh.fac2tag(ithread, ifaceGlobal) >= msh.tag[ithread]) continue;

        int ent2pol[4];
        ent2pol[0] = ipins;
        ent2pol[lnofa3[0][0]] = ent2poi(ienttCav, lnofa3[jj][0]);
        ent2pol[lnofa3[0][1]] = ent2poi(ienttCav, lnofa3[jj][1]);
        ent2pol[lnofa3[0][2]] = ent2poi(ienttCav, lnofa3[jj][2]);

        if(ent2pol[1]==ipins || ent2pol[2]==ipins || ent2pol[3]==ipins) continue;

        // copy into tmpEntt for metqua
        ent2poi(tmpEntt,0)=ent2pol[0];
        ent2poi(tmpEntt,1)=ent2pol[1];
        ent2poi(tmpEntt,2)=ent2pol[2];
        ent2poi(tmpEntt,3)=ent2pol[3];

        double qua;
        const bool isValid = completedCandidate3D(ienttCav,qua);
        METRIS_ENFORCE_MSG(isValid,
            "Final tentative cavity has an invalid completed-P2 element");
        quaCav1 += cavity_element_contribution<iquaf,MFT,3,3>(
            msh,AsDeg::P1,tmpEntt,qua,targetWeightCav1);
        nQuaCav1++;
        if(qua > quaMax1) quaMax1 = qua;
      } // if tdim == 2 else 3
    } // for jj (bnd facets of ienttCav)
  } // for ienttCav

  #ifdef CHECKSUBENTTQUAL
  // second the subentities
  const intAr2& fac2fac = msh.ent2ent(2);
        intAr2& fac2poi = msh.ent2poi(2);
  for (int isubentt : lcsub){
    if (tdim == 2) break; // skip for edges in 2D (2D already working good without this, might implement it in the future)

    double qua;
    qua = metqua<MFT,3,2,iquaf>(msh, AsDeg::P1, AsDeg::P1, isubentt, difto);
    quaSub0 += cavity_element_contribution<iquaf,MFT,3,2>(
        msh,AsDeg::P1,isubentt,qua,targetWeightSub0);
    nQuaSub0++;
    if (qua > quaMaxSub0) quaMaxSub0 = qua;

    // loop over edges of this face
    int ent2pol[3];
    for (int jj = 0; jj < 3; jj++){

      int ifnei = fac2fac(isubentt,jj);
      if (ifnei >= 0 && msh.fac2tag(ithread,ifnei) >= msh.tag[ithread]) continue; // neighbor face in cavity

      ent2pol[0] = ipins;
      ent2pol[lnoed2[0][0]] = fac2poi(isubentt,lnoed2[jj][0]);
      ent2pol[lnoed2[0][1]] = fac2poi(isubentt,lnoed2[jj][1]);

      fac2poi(tmpSubEntt,0) = ent2pol[0];
      fac2poi(tmpSubEntt,1) = ent2pol[1];
      fac2poi(tmpSubEntt,2) = ent2pol[2];

      msh.fac2ref[tmpSubEntt] = msh.fac2ref[isubentt];

      // no need to check validity, reconnection of seed cavity should yield valid config
      // by construction of seed cavity

      #ifndef NDEBUG
      // put this for debug build anyways
      double meas;
      bool isValid = isvalideltP1<3,2>(msh, tmpSubEntt, NULL, &meas);
      METRIS_ASSERT_MSG(isValid, "Initial cavity has invalid element when reconnected. Shouldn't ever happen");
      #endif

      qua = metqua<MFT,3,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,tmpSubEntt,difto);
      quaSub1 += cavity_element_contribution<iquaf,MFT,3,2>(
          msh,AsDeg::P1,tmpSubEntt,qua,targetWeightSub1);
      nQuaSub1++;
      if (qua > quaMaxSub1) quaMaxSub1 = qua;

    }
  }
  #endif

  #if defined(CAVGROWTH) && !defined(NDEBUG)
  if constexpr(iquaf == QuaFun::StepDistance){
    if(msh.param->step_distance_cavity_target_average){
      const double expectedGrowthNumerator =
          handler.getQualitySum() - quaCav0 + quaCav1;
      const int expectedGrowthElementCount =
          handler.getQualityCount() - nQuaCav0 + nQuaCav1;
      const double growthScale =
          MAX(1.,std::abs(expectedGrowthNumerator));
      METRIS_ASSERT(
          std::abs(cavityGrowthObjective.numerator
                   - expectedGrowthNumerator)
          <= 1.e-11*growthScale);
      METRIS_ASSERT(cavityGrowthObjective.element_count
                    == expectedGrowthElementCount);
    }
  }
  #endif

  #ifdef CAVSMOOTHING
  if(msh.curdeg == 1){
    double quaCav1AfterSmoo;
    double quaMax1AfterSmoo;
    double targetWeightCav1AfterSmoo;
    const double objCav1BeforeSmoo =
        cavity_region_objective<iquaf>(
            msh,quaCav1,nQuaCav1,targetWeightCav1);
    [[maybe_unused]] double statSmooCav = smoothCavity(
        msh,cav,handler,iquaf,
        quaCav1,quaMax1,targetWeightCav1,
        quaCav1AfterSmoo,quaMax1AfterSmoo,targetWeightCav1AfterSmoo,
        ithread,ithread);
    retagCavity();

    quaCav1 = quaCav1AfterSmoo;
    quaMax1 = quaMax1AfterSmoo;
    targetWeightCav1 = targetWeightCav1AfterSmoo;
    if constexpr(iquaf == QuaFun::StepDistance){
      if(tdim == 2 && msh.param->step_distance_cavity_target_average){
        const bool validAfterSmoo = getCavityQuality2D(
            quaCav0,quaCav1,quaMax0,quaMax1,
            nQuaCav0,nQuaCav1,targetWeightCav0,targetWeightCav1);
        METRIS_ENFORCE(validAfterSmoo);
      }
    }
    const double objCav1AfterSmoo =
        cavity_region_objective<iquaf>(
            msh,quaCav1,nQuaCav1,targetWeightCav1);
    METRIS_ENFORCE_MSG(objCav1AfterSmoo <= objCav1BeforeSmoo,
                       "Cavity smoothing worsen quality!");
  }
  #endif

  double edgeObjective0 = 0.0;
  double edgeObjective1 = 0.0;
  const bool improveEdgeCavity =
      edge_cavity_length_objective_nonworsening_2d<MFT,iquaf>(
          msh,cav,ithread,edgeObjective0,edgeObjective1);
  if(cav.lcedg.get_n() > 0 && msh.idim == 2
     && msh.CAD() && msh.param->adp_line_adapt){
    CPRINTF2(" - 1D length objective {} -> {}, accepted {}\n",
             edgeObjective0,edgeObjective1,improveEdgeCavity);
  }

  // restore to original number of entities in mesh
  msh.set_nentt(tdim,nentt0);
  msh.set_nentt(tdim-1,nsube0);
  msh.set_npoin(npoiProbe0);

  double objCav0;
  double objCav1;
  double objGlobal0;
  double objGlobal1;
  cavity_replacement_objectives<iquaf>(
      msh,handler,
      quaCav0,nQuaCav0,
      quaCav1,nQuaCav1,
      targetWeightCav0,targetWeightCav1,
      objCav0,objCav1,objGlobal0,objGlobal1);
  CPRINTF2(" - replacement objective local {} -> {}, global {} -> {}; "
           "numerator {} -> {}, unit weight {} -> {}, elements {} -> {}\n",
           objCav0,objCav1,objGlobal0,objGlobal1,
           quaCav0,quaCav1,targetWeightCav0,targetWeightCav1,
           nQuaCav0,nQuaCav1);
  const bool improveEnttsSum =
      cavity_replacement_accepts<iquaf>(
          msh,handler,objGlobal0,objGlobal1);
  bool improveEnttsMax = true;
  #ifdef IMPROVEMAXQUAL
  improveEnttsMax = quaMax1 <= quaMax0;
  #endif

  bool improveSubEnttsSum = true;
  bool improveSubEnttsMax = true;
  #ifdef CHECKSUBENTTQUAL
  if (ipinsOnBnd){
    const double objSub0 =
        cavity_region_objective<iquaf>(
            msh,quaSub0,nQuaSub0,targetWeightSub0);
    const double objSub1 =
        cavity_region_objective<iquaf>(
            msh,quaSub1,nQuaSub1,targetWeightSub1);
    improveSubEnttsSum = handler.checkSuccess(objSub1,objSub0);
    #ifdef IMPROVEMAXQUAL
    improveSubEnttsMax = quaMaxSub1 <= quaMaxSub0;
    #endif
  }
  #endif

  if (improveEnttsSum && improveEnttsMax && improveSubEnttsSum
      && improveSubEnttsMax && improveEdgeCavity) return 0; // reconnected cavity has better quality than original config
  else return -1;                                       // original config is better
}

#define INSTANTIATE_INCREASE_CAVITY_QUALITY(MFT_VAL, QUAFUN_VAL) \
template int increase_cavity_quality<MFT_VAL,QUAFUN_VAL>( \
    Mesh<MFT_VAL> &msh, MshCavity &cav, int tdim, int ngrow, \
    BadEntHandler& handler, int ithread);

INSTANTIATE_INCREASE_CAVITY_QUALITY(MetricFieldAnalytical, QuaFun::SizeShape)
INSTANTIATE_INCREASE_CAVITY_QUALITY(MetricFieldAnalytical, QuaFun::StepDistance)
INSTANTIATE_INCREASE_CAVITY_QUALITY(MetricFieldFE, QuaFun::SizeShape)
INSTANTIATE_INCREASE_CAVITY_QUALITY(MetricFieldFE, QuaFun::StepDistance)

#undef INSTANTIATE_INCREASE_CAVITY_QUALITY


// Check cavity quality
template<class MFT>
int checkCavityQuality(Mesh<MFT> &msh, MshCavity &cav, int tdim,
                       int ngrow, BadEntHandler& handler, const double worsenPctg, int ithread){

  if(tdim <= 1) return 0;

  #ifdef STEPDISTANCE
  constexpr QuaFun iquaf = QuaFun::StepDistance;
  #else
  constexpr QuaFun iquaf = QuaFun::SizeShape;
  #endif

  GETVDEPTH(msh.param);
  METRIS_ASSERT(tdim <= cav.get_tdim());
  METRIS_ASSERT(tdim == 2 || tdim == 3);
  METRIS_ASSERT_MSG(msh.idim + tdim != 5, "Qual-based algo not yet supported for 2D mesh in 3D space. I guess..., it might work");

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

  intAr1 &lcent = cav.lcent(tdim);
  intAr1 &lcsub = cav.lcent(tdim-1);

  const bool ipinsOnBnd = lcsub.get_n() > 0;

  int iedins = -1;
  const bool ipinsOnEdge = cav.lcedg.get_n() > 0;
  if (ipinsOnEdge){
    METRIS_ASSERT(cav.lcedg.get_n() <= 3); // 1 if insertion, and could be 2 or 3 if collapse
    iedins = cav.lcedg[0];
  }

  CPRINTF1("-- START checkCavityQuality {}\n",tdim);
  const intAr2&  ent2ent = msh.ent2ent(tdim);
        intAr2&  ent2poi = msh.ent2poi(tdim);
        intAr2r& ent2tag = msh.ent2tag(tdim);
        intAr2r& sub2tag = msh.ent2tag(tdim-1);

  // start by computing the quality of the initial cav
  // for both configurations: current and reconnected, and for both entities and subentities

  const int nentt0 = msh.nentt(tdim);
  msh.set_nentt(tdim,nentt0+1);
  const int tmpEntt = nentt0; // index for temporary entity to construct would-be elements

  const int ipins = cav.ipins;

  double quaCav0 = 0.; // for current config (sum of all qual)
  double quaMax0 = -1.; // worst qual for current config
  double quaCav1 = 0.; // for reconnected config (sum of all qual)
  double quaMax1 = 1.; // worst qual for reconnected config
  int nQuaCav0 = 0;
  int nQuaCav1 = 0;
  double targetWeightCav0 = 0.;
  double targetWeightCav1 = 0.;

  // same as above but for subentities (faces when tdim == 3)
  double quaSub0 = 0.;
  double quaMaxSub0 = -1.;
  double quaSub1 = 0.;
  double quaMaxSub1 = -1.;
  int nQuaSub0 = 0;
  int nQuaSub1 = 0;
  double targetWeightSub0 = 0.;
  double targetWeightSub1 = 0.;
  const int nsube0 = msh.nentt(tdim-1);
  msh.set_nentt(tdim-1,nsube0+1); // in case we need to create temporary subentities
  const int tmpSubEntt = nsube0;

  double difto = 1.;

  // compare current config vs reconnected config
  // first the entities
  for (const int ienttCav : lcent){

    // first add entt qual for current config
    double qua;
    if (tdim == 2) qua = metqua<MFT,2,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,ienttCav,difto);
    else           qua = metqua<MFT,3,3,iquaf>(msh,AsDeg::P1,AsDeg::P1,ienttCav,difto);
    if(tdim == 2){
      quaCav0 += cavity_element_contribution<iquaf,MFT,2,2>(
          msh,AsDeg::P1,ienttCav,qua,targetWeightCav0);
    }else{
      quaCav0 += cavity_element_contribution<iquaf,MFT,3,3>(
          msh,AsDeg::P1,ienttCav,qua,targetWeightCav0);
    }
    nQuaCav0++;
    if (qua > quaMax0) quaMax0 = qua;

    // now identy boundary facets in this element to construct
    // and compute quality of reconnected config

    int ent2pol[4];
    for(int jj = 0; jj < tdim + 1; jj++){

      const int ienei = ent2ent(ienttCav,jj);

      // if neighbor tagged it is in cavity -> skip
      if (ienei >= 0 && ent2tag(ithread,ienei) >= msh.tag[ithread]) continue;

      // at this point, facet jj is at boundary of the cavity
      // need to probe quality of the reconnected element

      if (tdim == 2){

        if (ipinsOnEdge){

          // check that facet jj, if also on boundary, is not in same boundary as the insertion edge. otherwise the new triangle would be flat
          int iedgeGlobal = msh.fac2edg(ienttCav,jj);
          if (iedgeGlobal >= 0){
            if (msh.edg2ref[iedgeGlobal] == msh.edg2ref[iedins]) continue;
          }
        }

        // put together new triangle
        ent2pol[0] = ipins;
        ent2pol[lnoed2[0][0]] = ent2poi(ienttCav,lnoed2[jj][0]);
        ent2pol[lnoed2[0][1]] = ent2poi(ienttCav,lnoed2[jj][1]);

        if (ent2pol[1] == ipins || ent2pol[2] == ipins) continue;

        ent2poi(tmpEntt,0) = ent2pol[0];
        ent2poi(tmpEntt,1) = ent2pol[1];
        ent2poi(tmpEntt,2) = ent2pol[2];

        // no need to check validity, must be valid by construction

        #ifndef NDEBUG
        // put this for debug build anyways
        double meas;
        bool isValid = isvalideltP1<2,2>(msh, tmpEntt, NULL, &meas);
        METRIS_ASSERT_MSG(isValid, "Final cavity has invalid element when reconnected. Shouldn't ever happen");
        #endif

        qua = metqua<MFT,2,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,tmpEntt,difto);
        quaCav1 += cavity_element_contribution<iquaf,MFT,2,2>(
            msh,AsDeg::P1,tmpEntt,qua,targetWeightCav1);
        nQuaCav1++;
        if (qua > quaMax1) quaMax1 = qua;
      }
      else{ // tdim == 3

        // if boundary face itself is in cavity (tagged), it will be split => no single cone tet
        int ifaceGlobal = msh.tet2fac(ienttCav, jj);
        if(ifaceGlobal >= 0 && msh.fac2tag(ithread, ifaceGlobal) >= msh.tag[ithread]) continue;

        int ent2pol[4];
        ent2pol[0] = ipins;
        ent2pol[lnofa3[0][0]] = ent2poi(ienttCav, lnofa3[jj][0]);
        ent2pol[lnofa3[0][1]] = ent2poi(ienttCav, lnofa3[jj][1]);
        ent2pol[lnofa3[0][2]] = ent2poi(ienttCav, lnofa3[jj][2]);

        if(ent2pol[1]==ipins || ent2pol[2]==ipins || ent2pol[3]==ipins) continue;

        // copy into tmpEntt for metqua
        ent2poi(tmpEntt,0)=ent2pol[0];
        ent2poi(tmpEntt,1)=ent2pol[1];
        ent2poi(tmpEntt,2)=ent2pol[2];
        ent2poi(tmpEntt,3)=ent2pol[3];

        // no need to check validity, must be valid by construction

        #ifndef NDEBUG
        // put this for debug build anyways
        double meas;
        bool isValid = isvalideltP1<3,3>(msh, tmpEntt, NULL, &meas);
        METRIS_ASSERT_MSG(isValid, "Final cavity has invalid element when reconnected. Shouldn't ever happen");
        #endif

        double qua = metqua<MFT,3,3,iquaf>(msh, AsDeg::P1, AsDeg::P1, tmpEntt, difto);
        quaCav1 += cavity_element_contribution<iquaf,MFT,3,3>(
            msh,AsDeg::P1,tmpEntt,qua,targetWeightCav1);
        nQuaCav1++;
        if(qua > quaMax1) quaMax1 = qua;
      } // if tdim == 2 else 3
    } // for jj (bnd facets of ienttCav)
  } // for ienttCav

  #ifdef CHECKSUBENTTQUAL
  // second the subentities
  const intAr2& fac2fac = msh.ent2ent(2);
        intAr2& fac2poi = msh.ent2poi(2);
  for (int isubentt : lcsub){
    if (tdim == 2) break; // skip for edges in 2D (2D already working good without this, might implement it in the future)

    double qua;
    qua = metqua<MFT,3,2,iquaf>(msh, AsDeg::P1, AsDeg::P1, isubentt, difto);
    quaSub0 += cavity_element_contribution<iquaf,MFT,3,2>(
        msh,AsDeg::P1,isubentt,qua,targetWeightSub0);
    nQuaSub0++;
    if (qua > quaMaxSub0) quaMaxSub0 = qua;

    // loop over edges of this face
    int ent2pol[3];
    for (int jj = 0; jj < 3; jj++){

      int ifnei = fac2fac(isubentt,jj);
      if (ifnei >= 0 && msh.fac2tag(ithread,ifnei) >= msh.tag[ithread]) continue; // neighbor face in cavity

      ent2pol[0] = ipins;
      ent2pol[lnoed2[0][0]] = fac2poi(isubentt,lnoed2[jj][0]);
      ent2pol[lnoed2[0][1]] = fac2poi(isubentt,lnoed2[jj][1]);

      fac2poi(tmpSubEntt,0) = ent2pol[0];
      fac2poi(tmpSubEntt,1) = ent2pol[1];
      fac2poi(tmpSubEntt,2) = ent2pol[2];

      msh.fac2ref[tmpSubEntt] = msh.fac2ref[isubentt];

      // no need to check validity, reconnection of seed cavity should yield valid config
      // by construction of seed cavity

      #ifndef NDEBUG
      // put this for debug build anyways
      double meas;
      bool isValid = isvalideltP1<3,2>(msh, tmpSubEntt, NULL, &meas);
      METRIS_ASSERT_MSG(isValid, "Initial cavity has invalid element when reconnected. Shouldn't ever happen");
      #endif

      qua = metqua<MFT,3,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,tmpSubEntt,difto);
      quaSub1 += cavity_element_contribution<iquaf,MFT,3,2>(
          msh,AsDeg::P1,tmpSubEntt,qua,targetWeightSub1);
      nQuaSub1++;
      if (qua > quaMaxSub1) quaMaxSub1 = qua;

    }
  }
  #endif

  double edgeObjective0 = 0.0;
  double edgeObjective1 = 0.0;
  const bool improveEdgeCavity =
      edge_cavity_length_objective_nonworsening_2d<MFT,iquaf>(
          msh,cav,ithread,edgeObjective0,edgeObjective1);
  if(cav.lcedg.get_n() > 0 && msh.idim == 2
     && msh.CAD() && msh.param->adp_line_adapt){
    CPRINTF2(" - 1D length objective {} -> {}, accepted {}\n",
             edgeObjective0,edgeObjective1,improveEdgeCavity);
  }

  // restore to original number of entities in mesh
  msh.set_nentt(tdim,nentt0);
  msh.set_nentt(tdim-1,nsube0);

  double objCav0;
  double objCav1;
  double objGlobal0;
  double objGlobal1;
  cavity_replacement_objectives<iquaf>(
      msh,handler,
      quaCav0,nQuaCav0,
      quaCav1,nQuaCav1,
      targetWeightCav0,targetWeightCav1,
      objCav0,objCav1,objGlobal0,objGlobal1);
  CPRINTF2(" - replacement objective local {} -> {}, global {} -> {}; "
           "numerator {} -> {}, unit weight {} -> {}, elements {} -> {}\n",
           objCav0,objCav1,objGlobal0,objGlobal1,
           quaCav0,quaCav1,targetWeightCav0,targetWeightCav1,
           nQuaCav0,nQuaCav1);
  // A committed topology change is always a descent step for its configured
  // acceptance objective. worsenPctg may be used while probing/growing a
  // candidate cavity, but it must not relax the final test. For
  // CavityTargetAverage, the configured objective is the mesh-wide arithmetic
  // mean; no separate local descent condition is imposed.
  (void)worsenPctg;
  const bool cavAccepted =
      cavity_replacement_accepts<iquaf>(
          msh,handler,objGlobal0,objGlobal1);

  bool improveEnttsMax = true;
  #ifdef IMPROVEMAXQUAL
  improveEnttsMax = quaMax1 <= quaMax0;
  #endif

  bool improveSubEnttsSum = true;
  bool improveSubEnttsMax = true;
  #ifdef CHECKSUBENTTQUAL
  if (ipinsOnBnd){
    const double objSub0 =
        cavity_region_objective<iquaf>(
            msh,quaSub0,nQuaSub0,targetWeightSub0);
    const double objSub1 =
        cavity_region_objective<iquaf>(
            msh,quaSub1,nQuaSub1,targetWeightSub1);
    improveSubEnttsSum = handler.checkSuccess(objSub1,objSub0);
    #ifdef IMPROVEMAXQUAL
    improveSubEnttsMax = quaMaxSub1 <= quaMaxSub0;
    #endif
  }
  #endif

  if (cavAccepted && improveEnttsMax && improveSubEnttsSum
      && improveSubEnttsMax && improveEdgeCavity) return 0; // reconnected cavity has an acceptable quality
  else return -1;                                       // original config must be kept
}

template int checkCavityQuality(Mesh<MetricFieldAnalytical> &msh,
                                      MshCavity &cav, int tdim, int ngrow, BadEntHandler& handler, const double worsenPctg, int ithread);
template int checkCavityQuality(Mesh<MetricFieldFE        > &msh,
                                      MshCavity &cav, int tdim, int ngrow, BadEntHandler& handler, const double worsenPctg, int ithread);

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

#define INSTANTIATE_COMPLETED_P2_CAVITY_CONE(MFT_VALUE,DIMENSION,OBJECTIVE) \
template double evaluate_completed_p2_cavity_cone< \
    MFT_VALUE,DIMENSION,OBJECTIVE>( \
        Mesh<MFT_VALUE>&,const MshCavity&,int,const int*,bool*);

INSTANTIATE_COMPLETED_P2_CAVITY_CONE(
    MetricFieldAnalytical,2,QuaFun::SizeShape)
INSTANTIATE_COMPLETED_P2_CAVITY_CONE(
    MetricFieldAnalytical,2,QuaFun::StepDistance)
INSTANTIATE_COMPLETED_P2_CAVITY_CONE(
    MetricFieldAnalytical,3,QuaFun::SizeShape)
INSTANTIATE_COMPLETED_P2_CAVITY_CONE(
    MetricFieldAnalytical,3,QuaFun::StepDistance)
INSTANTIATE_COMPLETED_P2_CAVITY_CONE(
    MetricFieldFE,2,QuaFun::SizeShape)
INSTANTIATE_COMPLETED_P2_CAVITY_CONE(
    MetricFieldFE,2,QuaFun::StepDistance)
INSTANTIATE_COMPLETED_P2_CAVITY_CONE(
    MetricFieldFE,3,QuaFun::SizeShape)
INSTANTIATE_COMPLETED_P2_CAVITY_CONE(
    MetricFieldFE,3,QuaFun::StepDistance)

#undef INSTANTIATE_COMPLETED_P2_CAVITY_CONE



} // end namespace
