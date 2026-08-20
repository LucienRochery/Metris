// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_stepdistance_smoothing_scan

#include "common_setup.hxx"

#include "MetrisRunner/MetrisRunner.hxx"
#include "aux_topo.hxx"
#include "ho_constants.hxx"
#include "io_libmeshb.hxx"
#include "low_eval.hxx"
#include "low_geo/measure.hxx"
#include "quality/aux_volumeMeasure.hxx"
#include "quality/low_metqua.hxx"
#include "smoothing/low_smooballdiff.hxx"
#include "smoothing/msh_smooball.hxx"
#include "smoothing/smoothing_progress.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>

using namespace Metris;

#if defined(STEPDISTANCE)

namespace{

constexpr int gdim = 2;
constexpr int tdim = 2;
constexpr int ideg = 1;
constexpr int nnmet = 3;

enum class PointKind{Interior, Edge, Face};

const char* point_kind_name(PointKind kind){
  switch(kind){
  case PointKind::Interior: return "interior";
  case PointKind::Edge: return "edge";
  case PointKind::Face: return "face";
  }
  return "unknown";
}

struct MeshQualityStats{
  int elements = 0;
  int invalid = 0;
  double sum = 0.;
  double average = 0.;
  double minimum = std::numeric_limits<double>::infinity();
  double maximum = -std::numeric_limits<double>::infinity();
  double minimum_measure = std::numeric_limits<double>::infinity();
  double maximum_measure = 0.;
  double minimum_metric_volume = std::numeric_limits<double>::infinity();
  double maximum_metric_volume = 0.;
};

struct PointAttempt{
  int ipoin = -1;
  int seed = -1;
  int ball_size = 0;
  PointKind kind = PointKind::Interior;
  int error = 0;
  bool solver_success = false;
  bool accepted = false;
  bool substantive = false;
  double qsum0 = 0.;
  double qsum1 = 0.;
  double qmax0 = 0.;
  double qmax1 = 0.;
  double displacement = 0.;
  std::array<double,gdim> initial{{0.,0.}};
  std::array<double,gdim> final{{0.,0.}};
};

struct ScanStats{
  std::string phase;
  bool commit = false;
  int active_points = 0;
  int skipped_dead = 0;
  int skipped_corner = 0;
  int setup_failures = 0;
  int attempted = 0;
  int attempted_interior = 0;
  int attempted_edge = 0;
  int attempted_face = 0;
  int solver_success = 0;
  int accepted = 0;
  int substantive = 0;
  int accepted_local_max_worse = 0;
  int rejected_local = 0;
  int rejected_global_objective = 0;
  int rejected_global_max = 0;
  int exceptions = 0;
  std::map<int,int> errors;
  MeshQualityStats before;
  MeshQualityStats after;
  std::vector<PointAttempt> rows;
};

struct BoundarySnapshot{
  int ibpoi = -1;
  double r0 = 0.;
  double r1 = 0.;
};

std::string case_dir(){
  if(const char* value = std::getenv("METRIS_SMOOTHING_CASE_DIR")){
    return value;
  }
  return "/Users/renat/MIT/HOMeshing/"
         "L2Project_MOESS_ConfBasedXClassicInsertions/"
         "L2Project_MOESS_CavSmooth_0Pctg/CornSing_CG/"
         "L2Project_2D_StepDistance/Corner/Metris_CG_000500_P1";
}

int case_iteration(){
  if(const char* value = std::getenv("METRIS_SMOOTHING_ITERATION")){
    std::string iteration(value);
    if(!iteration.empty() && iteration.front() == 'a') iteration.erase(0,1);
    return std::atoi(iteration.c_str());
  }
  return 6;
}

int max_attempts(){
  if(const char* value = std::getenv("METRIS_SMOOTHING_MAX_ATTEMPTS")){
    return std::max(0,std::atoi(value));
  }
  return 0;
}

std::filesystem::path output_dir(){
  if(const char* value = std::getenv("METRIS_SMOOTHING_OUTPUT_DIR")){
    return value;
  }
  return std::filesystem::path(METRIS_ROOT_DIR) / "build" / "codex_diagnostic"
       / "stepdistance_smoothing_scan";
}

std::string input_command(){
  const std::string dir = case_dir();
  const std::string suffix = "a" + std::to_string(case_iteration());
  const std::string mesh = dir + "/mesh_MOESS_initial_" + suffix + ".meshb";
  const std::string metric = dir + "/met_MOESS_initial_" + suffix + ".solb";
  std::string cad;
  if (const char* value = std::getenv("METRIS_SMOOTHING_CAD"))
  {
    cad = value;
  }
  else
  {
    cad = dir + "/CAD_MOESS.egads";
  }
  if(!std::filesystem::exists(cad)){
    cad = "/Users/renat/MIT/HOMeshing/"
          "L2Project_MOESS_ConfBasedXClassicInsertions/CAD_MOESS.egads";
  }

  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(mesh),"Missing mesh: " + mesh);
  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(metric),"Missing metric: " + metric);
  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(cad),"Missing CAD: " + cad);

  std::string command = "-in " + mesh + " -met " + metric + " -cad " + cad
      + " -verb 0 -vdepth 0 -adapt 0 -opt-niter 0 -adp-opt-niter 0"
      + " -opt-pnorm 1 -opt-power 1 -opt-smoo-niter 1"
      ;

  if(const char* step_args =
      std::getenv("METRIS_SMOOTHING_STEP_DISTANCE_ARGS")){
    command += " ";
    command += step_args;
  }else{
    command += " --step-distance-p 1"
               " --step-distance-regularization 1e-8"
               " --step-distance-barrier-rho0 0.7"
               " --step-distance-barrier-beta 2";
  }

  if(const char* extra = std::getenv("METRIS_SMOOTHING_EXTRA_ARGS")){
    command += " ";
    command += extra;
  }
  return command;
}

template<class MFT>
MeshQualityStats mesh_quality(Mesh<MFT>& msh){
  const auto quafun = get_quafun<MFT,gdim,tdim>(QuaFun::StepDistance);
  const intAr2& ent2poi = msh.ent2poi(tdim);
  MeshQualityStats stats;

  for(int iface = 0; iface < msh.nface; iface++){
    if(isdeadent(iface,ent2poi)) continue;
    if(!isvalideltP1<gdim,tdim>(msh,iface)) stats.invalid++;
    const double quality = quafun(msh,AsDeg::P1,AsDeg::P1,iface,1.);
    const double measure =
        getmeasentP1<gdim>(ent2poi[iface],msh.coord);
    stats.elements++;
    stats.sum += quality;
    stats.minimum = std::min(stats.minimum,quality);
    stats.maximum = std::max(stats.maximum,quality);
    stats.minimum_measure = std::min(stats.minimum_measure,measure);
    stats.maximum_measure = std::max(stats.maximum_measure,measure);

    constexpr int nquad = tdim+2;
    for(int iquad = 0; iquad < nquad; iquad++){
      double bary[tdim+1] = {0.,0.,0.};
      if(iquad < tdim+1){
        bary[iquad] = 1.;
      }else{
        for(double& value : bary) value = 1./(tdim+1);
      }

      double coordinate[gdim];
      double jmat[tdim*gdim];
      eval2<gdim,ideg>(msh.coord,ent2poi[iface],msh.getBasis(),
                       DifVar::Bary,DifVar::None,
                       bary,coordinate,jmat,NULL);

      double metric[nnmet];
      if(iquad < tdim+1){
        const int ipoin = ent2poi(iface,iquad);
        for(int imet = 0; imet < nnmet; imet++)
          metric[imet] = msh.met(ipoin,imet);
      }else{
        msh.met.getMetBary(AsDeg::P1,DifVar::None,msh.met.getSpace(),
                           ent2poi[iface],tdim,bary,metric,NULL);
      }

      double jreg_t[tdim*gdim];
      for(int i = 0; i < tdim; i++){
        for(int a = 0; a < gdim; a++){
          jreg_t[i*gdim+a] = 0.;
          for(int k = 0; k < tdim; k++){
            jreg_t[i*gdim+a] +=
                Constants::invtJ_0[hana::type_c<double>][tdim][i*tdim+k]
               *jmat[k*gdim+a];
          }
        }
      }

      double rho = 0.;
      double barrier = 0.;
      VolumeMeasureHelpers::eval_metric_volume_barrier_fixed_metric_grad<
          gdim,tdim,double>(jreg_t,metric,NULL,0.,0.,&rho,&barrier,NULL);
      stats.minimum_metric_volume =
          std::min(stats.minimum_metric_volume,rho);
      stats.maximum_metric_volume =
          std::max(stats.maximum_metric_volume,rho);
    }
  }
  if(stats.elements > 0) stats.average = stats.sum/stats.elements;
  return stats;
}

template<class MFT>
std::vector<BoundarySnapshot> snapshot_boundary(Mesh<MFT>& msh, int ipoin){
  std::vector<BoundarySnapshot> snapshot;
  int ibpoi = msh.poi2bpo[ipoin];
  int guard = 0;
  while(ibpoi >= 0 && msh.bpo2ibi(ibpoi,0) == ipoin){
    snapshot.push_back({ibpoi,msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1)});
    ibpoi = msh.bpo2ibi(ibpoi,3);
    METRIS_ENFORCE_MSG(++guard <= msh.nbpoi,
                       "Cycle in boundary records for point {}",ipoin);
  }
  return snapshot;
}

template<class MFT>
void restore_boundary(Mesh<MFT>& msh,
                      const std::vector<BoundarySnapshot>& snapshot){
  for(const BoundarySnapshot& entry : snapshot){
    msh.bpo2rbi(entry.ibpoi,0) = entry.r0;
    msh.bpo2rbi(entry.ibpoi,1) = entry.r1;
  }
}

template<class MFT>
void update_edge_face_parameters(Mesh<MFT>& msh, int ipoin, int ibpoin){
  const double topt = msh.bpo2rbi(ibpoin,0);
  const int iedge = msh.bpo2ibi(ibpoin,2);
  METRIS_ENFORCE(iedge >= 0 && iedge < msh.nedge);
  const int iref_edge = msh.edg2ref[iedge];
  ego cad_edge = msh.CAD.cad2edg[iref_edge];

  int record = ibpoin;
  int guard = 0;
  while(record >= 0 && msh.bpo2ibi(record,0) == ipoin){
    if(msh.bpo2ibi(record,1) == 2){
      const int iface = msh.bpo2ibi(record,2);
      METRIS_ENFORCE(iface >= 0 && iface < msh.nface);
      ego cad_face = msh.CAD.cad2fac[msh.fac2ref[iface]];
      double uv[2];
      const int error = EG_getEdgeUV(cad_face,cad_edge,0,topt,uv);
      METRIS_ENFORCE_MSG(error == EGADS_SUCCESS,
                         "EG_getEdgeUV failed: {}",error);
      msh.bpo2rbi(record,0) = uv[0];
      msh.bpo2rbi(record,1) = uv[1];
    }
    record = msh.bpo2ibi(record,3);
    METRIS_ENFORCE_MSG(++guard <= msh.nbpoi,
                       "Cycle in boundary records for point {}",ipoin);
  }
}

double median(std::vector<double> values){
  if(values.empty()) return 0.;
  const size_t middle = values.size()/2;
  std::nth_element(values.begin(),values.begin()+middle,values.end());
  if(values.size()%2 == 1) return values[middle];
  const double upper = values[middle];
  std::nth_element(values.begin(),values.begin()+middle-1,values.begin()+middle);
  return 0.5*(values[middle-1]+upper);
}

void write_rows(const ScanStats& stats, const std::filesystem::path& file){
  std::ofstream out(file);
  BOOST_REQUIRE_MESSAGE(out.good(),"Could not open " + file.string());
  out << std::setprecision(17);
  out << "phase,commit,ipoin,kind,seed,ball_size,error,solver_success,"
         "accepted,substantive,qsum0,qsum1,qsum_decrease,qmax0,qmax1,"
         "qmax_decrease,x0,y0,x1,y1,displacement\n";
  for(const PointAttempt& row : stats.rows){
    out << stats.phase << ',' << stats.commit << ',' << row.ipoin << ','
        << point_kind_name(row.kind) << ',' << row.seed << ',' << row.ball_size
        << ',' << row.error << ',' << row.solver_success << ',' << row.accepted
        << ',' << row.substantive << ',' << row.qsum0 << ',' << row.qsum1
        << ',' << row.qsum0-row.qsum1 << ',' << row.qmax0 << ',' << row.qmax1
        << ',' << row.qmax0-row.qmax1 << ',' << row.initial[0] << ','
        << row.initial[1] << ',' << row.final[0] << ',' << row.final[1] << ','
        << row.displacement << '\n';
  }
}

void write_summary(const ScanStats& stats, const std::filesystem::path& file){
  std::vector<double> displacements;
  std::vector<double> decreases;
  for(const PointAttempt& row : stats.rows){
    if(!row.accepted) continue;
    displacements.push_back(row.displacement);
    decreases.push_back(row.qsum0-row.qsum1);
  }

  std::ofstream out(file);
  BOOST_REQUIRE_MESSAGE(out.good(),"Could not open " + file.string());
  out << std::setprecision(17);
  out << "phase " << stats.phase << '\n'
      << "commit " << stats.commit << '\n'
      << "active_points " << stats.active_points << '\n'
      << "skipped_dead " << stats.skipped_dead << '\n'
      << "skipped_corner " << stats.skipped_corner << '\n'
      << "setup_failures " << stats.setup_failures << '\n'
      << "attempted " << stats.attempted << '\n'
      << "attempted_interior " << stats.attempted_interior << '\n'
      << "attempted_edge " << stats.attempted_edge << '\n'
      << "attempted_face " << stats.attempted_face << '\n'
      << "solver_success " << stats.solver_success << '\n'
      << "accepted " << stats.accepted << '\n'
      << "substantive " << stats.substantive << '\n'
      << "accepted_local_max_worse " << stats.accepted_local_max_worse << '\n'
      << "rejected_local " << stats.rejected_local << '\n'
      << "rejected_global_objective "
      << stats.rejected_global_objective << '\n'
      << "rejected_global_max " << stats.rejected_global_max << '\n'
      << "exceptions " << stats.exceptions << '\n'
      << "success_fraction "
      << (stats.attempted > 0 ? double(stats.accepted)/stats.attempted : 0.)
      << '\n'
      << "substantive_fraction "
      << (stats.attempted > 0 ? double(stats.substantive)/stats.attempted : 0.)
      << '\n'
      << "quality_before_elements " << stats.before.elements << '\n'
      << "quality_before_invalid " << stats.before.invalid << '\n'
      << "quality_before_sum " << stats.before.sum << '\n'
      << "quality_before_average " << stats.before.average << '\n'
      << "quality_before_minimum " << stats.before.minimum << '\n'
      << "quality_before_maximum " << stats.before.maximum << '\n'
      << "measure_before_minimum " << stats.before.minimum_measure << '\n'
      << "measure_before_maximum " << stats.before.maximum_measure << '\n'
      << "metric_volume_before_minimum "
      << stats.before.minimum_metric_volume << '\n'
      << "metric_volume_before_maximum "
      << stats.before.maximum_metric_volume << '\n'
      << "quality_after_invalid " << stats.after.invalid << '\n'
      << "quality_after_sum " << stats.after.sum << '\n'
      << "quality_after_average " << stats.after.average << '\n'
      << "quality_after_minimum " << stats.after.minimum << '\n'
      << "quality_after_maximum " << stats.after.maximum << '\n'
      << "measure_after_minimum " << stats.after.minimum_measure << '\n'
      << "measure_after_maximum " << stats.after.maximum_measure << '\n'
      << "metric_volume_after_minimum "
      << stats.after.minimum_metric_volume << '\n'
      << "metric_volume_after_maximum "
      << stats.after.maximum_metric_volume << '\n'
      << "minimum_measure_ratio "
      << stats.after.minimum_measure/stats.before.minimum_measure << '\n'
      << "quality_sum_decrease " << stats.before.sum-stats.after.sum << '\n'
      << "accepted_displacement_median " << median(displacements) << '\n'
      << "accepted_displacement_maximum "
      << (displacements.empty() ? 0. : *std::max_element(displacements.begin(),
                                                         displacements.end()))
      << '\n'
      << "accepted_local_sum_decrease_median " << median(decreases) << '\n';
  for(const auto& error : stats.errors){
    out << "error_" << error.first << ' ' << error.second << '\n';
  }
}

void print_summary(const ScanStats& stats){
  std::cout << std::setprecision(8)
            << "\nStepDistance smoothing scan: " << stats.phase << '\n'
            << "  attempted: " << stats.attempted
            << " (interior " << stats.attempted_interior
            << ", edge " << stats.attempted_edge
            << ", face " << stats.attempted_face << ")\n"
            << "  solver success: " << stats.solver_success << '\n'
            << "  accepted: " << stats.accepted
            << " (" << (stats.attempted > 0
                          ? 100.*stats.accepted/stats.attempted : 0.) << "%)\n"
            << "  substantive: " << stats.substantive
            << " (" << (stats.attempted > 0
                          ? 100.*stats.substantive/stats.attempted : 0.) << "%)\n"
            << "  accepted with worse local maximum: "
            << stats.accepted_local_max_worse << '\n'
            << "  rejected locally: " << stats.rejected_local
            << ", rejected by global objective: "
            << stats.rejected_global_objective
            << ", rejected by global max: " << stats.rejected_global_max
            << ", setup failures: " << stats.setup_failures
            << ", exceptions: " << stats.exceptions << '\n'
            << "  whole-mesh sum: " << stats.before.sum << " -> "
            << stats.after.sum << '\n'
            << "  whole-mesh max: " << stats.before.maximum << " -> "
            << stats.after.maximum << '\n';
}

template<class MFT>
struct RegionTotals{
  double numerator = 0.;
  double target_weight = 0.;
  int count = 0;
};

template<class MFT>
RegionTotals<MFT> region_totals(Mesh<MFT>& msh, const intAr1& region){
  const auto quafun = get_quafun<MFT,gdim,tdim>(QuaFun::StepDistance);
  RegionTotals<MFT> totals;
  for(const int iface : region){
    if(isdeadent(iface,msh.fac2poi)) continue;
    const double quality =
        quafun(msh,AsDeg::P1,AsDeg::P1,iface,1.);
    if(msh.param->step_distance_cavity_target_average){
      const double weight =
          step_distance_element_target_weight<MFT,gdim,tdim>(
              msh,AsDeg::P1,iface);
      totals.numerator += weight*quality;
      totals.target_weight += weight;
    }else{
      totals.numerator += quality;
    }
    totals.count++;
  }
  return totals;
}

template<class MFT>
ScanStats run_scan(const std::string& command,
                   const std::string& phase,
                   bool commit,
                   bool require_global_improvement = false){
  cargHandler arg(command);
  MetrisRunner run(arg.c,arg.v);
  BOOST_REQUIRE(run.metricFE);
  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
  BOOST_REQUIRE_EQUAL(msh.idim,gdim);
  BOOST_REQUIRE_EQUAL(msh.get_tdim(),tdim);
  BOOST_REQUIRE_EQUAL(msh.curdeg,ideg);
  msh.met.setSpace(MetSpace::Exp);
  msh.setBasis(FEBasis::Lagrange);

  ScanStats stats;
  stats.phase = phase;
  stats.commit = commit;
  stats.before = mesh_quality(msh);
  const double global_qmax = stats.before.maximum;
  StepDistanceObjectiveState global_objective;
  if(require_global_improvement){
    global_objective =
        step_distance_global_objective_state<MFT,gdim,tdim>(
            msh,AsDeg::P1,AsDeg::P1);
  }

  std::vector<double> initial_coordinates(msh.npoin*gdim);
  std::vector<double> initial_metrics(msh.npoin*nnmet);
  std::vector<double> initial_boundary(msh.nbpoi*2);
  for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
    for(int i = 0; i < gdim; i++)
      initial_coordinates[ipoin*gdim+i] = msh.coord(ipoin,i);
    for(int i = 0; i < nnmet; i++)
      initial_metrics[ipoin*nnmet+i] = msh.met(ipoin,i);
  }
  for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
    initial_boundary[2*ibpoi] = msh.bpo2rbi(ibpoi,0);
    initial_boundary[2*ibpoi+1] = msh.bpo2rbi(ibpoi,1);
  }

  const int ithrd_ball = 2;
  intAr1 lball(100);
  const int limit = max_attempts();

  for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
    if(msh.isdeadpoint(ipoin)){
      stats.skipped_dead++;
      continue;
    }
    stats.active_points++;

    const int ibpoin = msh.poi2bpo[ipoin];
    PointKind kind = PointKind::Interior;
    int cad_dimension = 0;
    if(ibpoin >= 0){
      BOOST_REQUIRE_EQUAL(msh.bpo2ibi(ibpoin,0),ipoin);
      cad_dimension = msh.bpo2ibi(ibpoin,1);
      if(cad_dimension == 0){
        stats.skipped_corner++;
        continue;
      }
      if(cad_dimension == 1) kind = PointKind::Edge;
      else if(cad_dimension == 2) kind = PointKind::Face;
      else{
        stats.setup_failures++;
        continue;
      }
    }

    if(limit > 0 && stats.attempted >= limit) break;

    const int seed = getpoient(msh,ipoin,tdim);
    if(seed < 0 || seed >= msh.nface){
      stats.setup_failures++;
      continue;
    }
    const int local_vertex = msh.template getverfac<ideg>(seed,ipoin);
    if(local_vertex < 0 || local_vertex >= tdim+1){
      stats.setup_failures++;
      continue;
    }

    int iopen = 0;
    bool manifold = false;
    intAr1 dummy;
    const int ball_error = ball2(msh,ipoin,seed,lball,dummy,
                                 &iopen,&manifold,ithrd_ball);
    if(ball_error != 0 || !manifold || lball.get_n() <= 0){
      stats.setup_failures++;
      stats.errors[1000+ball_error]++;
      continue;
    }

    PointAttempt row;
    row.ipoin = ipoin;
    row.seed = seed;
    row.ball_size = lball.get_n();
    row.kind = kind;
    for(int i = 0; i < gdim; i++) row.initial[i] = msh.coord(ipoin,i);

    double coordinate0[gdim];
    double metric0[nnmet];
    for(int i = 0; i < gdim; i++) coordinate0[i] = msh.coord(ipoin,i);
    for(int i = 0; i < nnmet; i++) metric0[i] = msh.met(ipoin,i);
    const std::vector<BoundarySnapshot> boundary0 = snapshot_boundary(msh,ipoin);
    const RegionTotals<MFT> old_region = region_totals(msh,lball);

    stats.attempted++;
    if(kind == PointKind::Interior) stats.attempted_interior++;
    if(kind == PointKind::Edge) stats.attempted_edge++;
    if(kind == PointKind::Face) stats.attempted_face++;

    try{
      if(kind == PointKind::Interior){
        row.error = smooballdiff<MFT,gdim,ideg>(
            msh,ipoin,lball,&row.qsum0,&row.qmax0,&row.qsum1,&row.qmax1,
            QuaFun::StepDistance);
      }else{
        row.error = smooballdiff_boundary<MFT,gdim,ideg>(
            msh,ipoin,cad_dimension,lball,
            &row.qsum0,&row.qmax0,&row.qsum1,&row.qmax1,
            QuaFun::StepDistance);
      }
    }catch(...){
      row.error = -999;
      stats.exceptions++;
    }

    row.solver_success = row.error == 0;
    if(row.solver_success) stats.solver_success++;

    RegionTotals<MFT> new_region;
    bool global_improves = true;
    if(row.solver_success && require_global_improvement){
      new_region = region_totals(msh,lball);
      global_improves = global_objective.accepts_replacement(
          old_region.numerator,old_region.count,old_region.target_weight,
          new_region.numerator,new_region.count,new_region.target_weight);
    }

    if(row.solver_success && !global_improves){
      row.error = -3;
      stats.rejected_global_objective++;
    }else if(row.solver_success && row.qmax1 > global_qmax){
      row.error = -2;
      stats.rejected_global_max++;
    }else if(row.solver_success){
      row.accepted = true;
      stats.accepted++;
      row.substantive = smoothing_neighborhood_should_be_reactivated(
          row.qsum0,row.qsum1,row.qmax0,row.qmax1,
          msh.param->opt_smoo_tol,msh.param->opt_smoo_tol);
      if(row.substantive) stats.substantive++;
      if(row.qmax1 > row.qmax0) stats.accepted_local_max_worse++;
      if(require_global_improvement){
        global_objective.replace(
            old_region.numerator,old_region.count,old_region.target_weight,
            new_region.numerator,new_region.count,new_region.target_weight);
      }
    }else{
      stats.rejected_local++;
    }

    for(int i = 0; i < gdim; i++) row.final[i] = msh.coord(ipoin,i);
    for(int i = 0; i < gdim; i++){
      const double delta = row.final[i]-row.initial[i];
      row.displacement += delta*delta;
    }
    row.displacement = std::sqrt(row.displacement);

    if(!commit || !row.accepted){
      for(int i = 0; i < gdim; i++) msh.coord(ipoin,i) = coordinate0[i];
      for(int i = 0; i < nnmet; i++) msh.met(ipoin,i) = metric0[i];
      restore_boundary(msh,boundary0);
    }else if(kind == PointKind::Edge){
      update_edge_face_parameters(msh,ipoin,ibpoin);
    }

    if(!row.accepted) stats.errors[row.error]++;
    stats.rows.push_back(row);
  }

  stats.after = mesh_quality(msh);

  if(!commit){
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      for(int i = 0; i < gdim; i++)
        BOOST_CHECK_EQUAL(msh.coord(ipoin,i),initial_coordinates[ipoin*gdim+i]);
      for(int i = 0; i < nnmet; i++)
        BOOST_CHECK_EQUAL(msh.met(ipoin,i),initial_metrics[ipoin*nnmet+i]);
    }
    for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
      BOOST_CHECK_EQUAL(msh.bpo2rbi(ibpoi,0),initial_boundary[2*ibpoi]);
      BOOST_CHECK_EQUAL(msh.bpo2rbi(ibpoi,1),initial_boundary[2*ibpoi+1]);
    }
  }else{
    const std::filesystem::path dir = output_dir();
    std::filesystem::create_directories(dir);
    const std::string suffix = "a" + std::to_string(case_iteration());
    writeMesh((dir/(phase+"_"+suffix+".meshb")).string(),msh,false);
    msh.met.writeMetricFile((dir/(phase+"_"+suffix+".solb")).string());
  }

  return stats;
}

struct ProductionPassStats{
  int substantive = 0;
  MeshQualityStats before;
  MeshQualityStats after;
};

template<class MFT>
ProductionPassStats run_production_pass(const std::string& command){
  cargHandler arg(command);
  MetrisRunner run(arg.c,arg.v);
  BOOST_REQUIRE(run.metricFE);
  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
  BOOST_REQUIRE_EQUAL(msh.idim,gdim);
  BOOST_REQUIRE_EQUAL(msh.get_tdim(),tdim);
  BOOST_REQUIRE_EQUAL(msh.curdeg,ideg);
  msh.met.setSpace(MetSpace::Exp);
  msh.setBasis(FEBasis::Lagrange);

  ProductionPassStats stats;
  stats.before = mesh_quality(msh);
  const double normalized_operations =
      smoothInterior_Ball<MFT>(msh,QuaFun::StepDistance,1,2);
  stats.substantive = int(std::llround(normalized_operations*msh.nface));
  stats.after = mesh_quality(msh);
  return stats;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_stepdistance_smoothing_scan)
{
  const std::string command = input_command();
  const std::string cavity_command =
      command + " --step-distance-cavity-target-average";
  const std::filesystem::path dir = output_dir();
  std::filesystem::create_directories(dir);

  const ScanStats isolated =
      run_scan<MetricFieldFE>(command,"isolated",false);
  const ScanStats local_committed =
      run_scan<MetricFieldFE>(command,"local_committed",true);
  const ScanStats global_committed =
      run_scan<MetricFieldFE>(cavity_command,"global_committed",true,true);
  const ProductionPassStats production =
      run_production_pass<MetricFieldFE>(cavity_command);

  write_rows(isolated,dir/"isolated_attempts.csv");
  write_rows(local_committed,dir/"local_committed_attempts.csv");
  write_rows(global_committed,dir/"global_committed_attempts.csv");
  write_summary(isolated,dir/"isolated_summary.txt");
  write_summary(local_committed,dir/"local_committed_summary.txt");
  write_summary(global_committed,dir/"global_committed_summary.txt");
  print_summary(isolated);
  print_summary(local_committed);
  print_summary(global_committed);
  std::cout << "  production-pass substantive moves: "
            << production.substantive << '\n'
            << "  production-pass whole-mesh sum: " << production.before.sum
            << " -> " << production.after.sum << '\n';

  BOOST_REQUIRE_GT(isolated.attempted,0);
  BOOST_REQUIRE_GT(local_committed.attempted,0);
  BOOST_REQUIRE_GT(global_committed.attempted,0);
  BOOST_CHECK_EQUAL(isolated.before.invalid,0);
  BOOST_CHECK_EQUAL(isolated.after.invalid,0);
  BOOST_CHECK_EQUAL(local_committed.after.invalid,0);
  BOOST_CHECK_EQUAL(global_committed.after.invalid,0);
  BOOST_CHECK_SMALL(isolated.after.sum-isolated.before.sum,
                    1.e-12*std::max(1.,isolated.before.sum));
  BOOST_CHECK_LE(local_committed.after.sum,
                 local_committed.before.sum
                   + 1.e-10*std::max(1.,local_committed.before.sum));
  if(max_attempts() == 0 || global_committed.substantive == 0){
    BOOST_CHECK_EQUAL(production.substantive,global_committed.substantive);
    BOOST_CHECK_CLOSE(production.after.sum,global_committed.after.sum,1.e-10);
    BOOST_CHECK_CLOSE(production.after.maximum,
                      global_committed.after.maximum,1.e-10);
  }
}

#else

BOOST_AUTO_TEST_CASE(test_stepdistance_smoothing_scan_requires_stepdistance_build)
{
  BOOST_TEST_MESSAGE(
      "Skipping diagnostic: build with STEPDISTANCE.");
}

#endif
