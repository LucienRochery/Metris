//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_quality_insert_divergence

#include "common_setup.hxx"

#include "Adaptation/Insertion/EdgeSeed.hxx"
#include "Adaptation/Insertion/aux_insert.hxx"
#include "Adaptation/Insertion/insert_errors.hxx"
#include "Adaptation/Insertion/low_insert.hxx"
#include "Adaptation/low_collapse.hxx"
#include "Adaptation/low_cavqual.hxx"
#include "Adaptation/low_increasecav.hxx"
#include "Adaptation/msh_lineadapt.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "aux_badEntHandler.hxx"
#include "cavity/msh_cavity.hxx"
#include "io_libmeshb.hxx"
#include "low_geo/lenedg.hxx"
#include "msh_lenedg.hxx"
#include "quality/msh_metqua.hxx"
#include "quality/aux_volumeMeasure.hxx"
#include "quality/quafun_sizeshape.hxx"
#include "ho_constants.hxx"
#include "linalg/eigen.hxx"
#include "low_eval.hxx"
#include "smoothing/low_smoolen.hxx"
#include "smoothing/msh_smooball.hxx"
#include "utils/CT_loop.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

using namespace Metris;

namespace{

struct EdgeOpCandidate{
  int ientt;
  int ied;
  double qentt;
  double len;
  double dev;
};

struct InsertAttemptResult{
  int ierro = 0;
  int ipins = -1;
  int ncedg = 0;
  int ncfac = 0;
  int nctet = 0;
};

struct ObjectiveInsertCandidate{
  int ientt = -1;
  int ied = -1;
  int ip1 = -1;
  int ip2 = -1;
  double len = 0.;
  double size_shape_quality = 0.;
  double step_distance_quality = 0.;
};

struct ObjectiveInsertComparison{
  ObjectiveInsertCandidate candidate;
  InsertAttemptResult size_shape;
  InsertAttemptResult step_distance;
  InsertAttemptResult step_distance_length;
  std::string category;
};

struct ObjectiveInsertSummary{
  int n_candidates = 0;
  int n_both = 0;
  int n_size_shape_only = 0;
  int n_step_distance_only = 0;
  int n_neither = 0;
  int n_size_shape_noop = 0;
  int n_step_distance_noop = 0;
  int n_step_distance_length_success = 0;
  std::map<int,int> size_shape_errors;
  std::map<int,int> step_distance_errors;
  std::map<int,int> step_distance_length_errors;
  std::vector<ObjectiveInsertComparison> comparisons;
};

enum class CommonScanObjective{
  SizeShapeP1,
  StepDistanceP1,
  StepDistanceP2
};

struct CommonQuadratureValues{
  double size_shape_p1 = 0.;
  double step_distance_p1 = 0.;
  double step_distance_p2 = 0.;
  double step_distance_barrier = 0.;
  double min_metric_volume = std::numeric_limits<double>::infinity();
};

struct CommonCavityValues{
  CommonQuadratureValues original;
  CommonQuadratureValues reconnected;
  bool valid = true;
};

struct CommonObjectiveOptimization{
  double original = std::numeric_limits<double>::infinity();
  double initial = std::numeric_limits<double>::infinity();
  double optimized = std::numeric_limits<double>::infinity();
  double x = 0.;
  double y = 0.;
  double parameter = 0.;
  double min_vertex_distance_ratio = 0.;
  double min_child_area_fraction = 0.;
  double min_metric_volume = std::numeric_limits<double>::infinity();
  double original_barrier_unscaled = 0.;
  double optimized_barrier_unscaled = 0.;
  int probes = 0;
};

struct CommonPScanRow{
  ObjectiveInsertCandidate candidate;
  int setup_error = 0;
  CommonObjectiveOptimization size_shape;
  CommonObjectiveOptimization step_p1;
  CommonObjectiveOptimization step_p2;
};

struct CommonPScanSummary{
  int n_candidates = 0;
  int n_setup_failures = 0;
  int n_valid = 0;
  int n_size_before = 0;
  int n_step_p1_before = 0;
  int n_step_p2_before = 0;
  int n_size_after = 0;
  int n_step_p1_after = 0;
  int n_step_p2_after = 0;
  int n_size_step_p1_both = 0;
  int n_size_only_vs_step_p1 = 0;
  int n_step_p1_only_vs_size = 0;
  int n_size_step_p1_neither = 0;
  int n_step_p1_rescues_step_p2 = 0;
  int n_step_p2_only_vs_step_p1 = 0;
  std::vector<CommonPScanRow> rows;
};

struct DivergenceCase{
  bool found = false;
  int tdim = -1;
  int ientt = -1;
  int ied = -1;
  int ip1 = -1;
  int ip2 = -1;
  double qentt = -1;
  double len = -1;
  int ierro_quality = 0;
  int ierro_length = 0;
  InsertAttemptResult quality;
  InsertAttemptResult length;
};

struct StatefulAdaptResult{
  bool found = false;
  int iter = 0;
  int ntrySmoothing = 0;
  int nSuccessSmoothing = 0;
  int ntryInsert = 0;
  int nSuccessInsert = 0;
  int ntryInsertLength = 0;
  int nSuccessInsertLength = 0;
  int nSuccessInsertLengthInterior = 0;
  int nSuccessInsertLengthBoundary = 0;
  int ntryCollapse = 0;
  int nSuccessCollapse = 0;
  int ientt = -1;
  int ied = -1;
  int ip1 = -1;
  int ip2 = -1;
  double qentt = -1;
  double len = -1;
  int ierro_quality = 0;
  int ierro_length = 0;
  int ipins = -1;
  int ncedg = 0;
  int ncfac = 0;
  int nctet = 0;
  std::string pre_mesh;
  std::string pre_met;
};

struct CavityQualityStats{
  int nqua0 = 0;
  int nqua1 = 0;
  double qmin0 = 0;
  double qmin1 = 0;
  double qmax0 = 0;
  double qmax1 = 0;
  double qsum0 = 0;
  double qsum1 = 0;
};

const char* insertion_error_name(int error){
  switch(error){
  case -1: return "success";
  case INS2D_NOERR: return "no operation";
  case INS2D_ERR_CAVITYOPERATOR: return "cavity operator";
  case INS2D_ERR_LENQUA: return "length quality";
  case INS2D_ERR_NOQUALIMPROV: return "no objective improvement";
  default: return "other";
  }
}

CavOprOpt diagnostic_cavity_options();

std::string case_dir(){
  if(const char* env_dir = std::getenv("METRIS_QUALITY_INSERT_CASE_DIR")){
    return std::string(env_dir);
  }
  return "/Users/renat/MIT/HOMeshing/L2Project_MOESS_ConfBasedXClassicInsertions/"
         "L2Project_MOESS_0Pctg/CornSing_CG/L2Project_2D_QualBased/Corner/"
         "Metris_CG_000500_P1";
}

std::string cad_file(const std::string& dir){
  if(const char* env_cad = std::getenv("METRIS_QUALITY_INSERT_CAD")){
    return std::string(env_cad);
  }
  return dir + "/CAD_MOESS.egads";
}

int adaptive_iteration(){
  if(const char* env_iter = std::getenv("METRIS_QUALITY_INSERT_ADAPT_ITERATION")){
    std::string iter(env_iter);
    if(!iter.empty() && iter[0] == 'a') iter.erase(iter.begin());
    return std::atoi(iter.c_str());
  }
  return 1;
}

int diagnostic_verbosity(){
  if(const char* env_verb =
      std::getenv("METRIS_QUALITY_INSERT_VERBOSITY")){
    return std::atoi(env_verb);
  }
  return 0;
}

int diagnostic_vdepth(){
  if(const char* env_vdepth =
      std::getenv("METRIS_QUALITY_INSERT_VDEPTH")){
    return std::atoi(env_vdepth);
  }
  return 20;
}

int trace_ientt(){
  if(const char* env_ientt = std::getenv("METRIS_QUALITY_INSERT_TRACE_IENTT")){
    return std::atoi(env_ientt);
  }
  return 71;
}

int trace_ied(){
  if(const char* env_ied = std::getenv("METRIS_QUALITY_INSERT_TRACE_IED")){
    return std::atoi(env_ied);
  }
  return 0;
}

int trace_skip_face(){
  if(const char* env_face = std::getenv("METRIS_QUALITY_INSERT_TRACE_SKIP_FACE")){
    return std::atoi(env_face);
  }
  return -1;
}

bool scan_smooth_initial_cavity(){
  if(const char* env_smoo = std::getenv("METRIS_QUALITY_INSERT_SMOOTH_INITIAL")){
    return std::atoi(env_smoo) != 0;
  }
  return false;
}

int stateful_adapt_max_iter(){
  if(const char* env_iter = std::getenv("METRIS_QUALITY_ADAPT_MAX_ITER")){
    return std::atoi(env_iter);
  }
  return 10000;
}

bool stateful_adapt_require_boundary_divergence(){
  if(const char* env_boundary =
      std::getenv("METRIS_QUALITY_ADAPT_REQUIRE_BOUNDARY")){
    return std::atoi(env_boundary) != 0;
  }
  return false;
}

int stateful_adapt_trace_interval(){
  if(const char* env_iter = std::getenv("METRIS_QUALITY_ADAPT_TRACE_INTERVAL")){
    return std::atoi(env_iter);
  }
  return 0;
}

int stateful_adapt_trace_start(){
  if(const char* env_iter = std::getenv("METRIS_QUALITY_ADAPT_TRACE_START")){
    return std::atoi(env_iter);
  }
  return 0;
}

int stateful_adapt_save_iteration(){
  if(const char* env_iter = std::getenv("METRIS_QUALITY_ADAPT_SAVE_ITER")){
    return std::atoi(env_iter);
  }
  return 0;
}

std::string stateful_adapt_save_dir(){
  if(const char* env_dir = std::getenv("METRIS_QUALITY_ADAPT_SAVE_DIR")){
    return std::string(env_dir);
  }
  return {};
}

std::string input_command(){
  const std::string dir = case_dir();
  const std::string adapt_suffix = "a" + std::to_string(adaptive_iteration());
  const std::string mesh = dir + "/mesh_MOESS_initial_" + adapt_suffix + ".meshb";
  const std::string met  = dir + "/met_MOESS_initial_" + adapt_suffix + ".solb";
  const std::string cad  = cad_file(dir);

  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(mesh),
                        "Missing diagnostic mesh: " + mesh);
  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(met),
                        "Missing diagnostic metric: " + met);
  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(cad),
                        "Missing diagnostic CAD: " + cad);

  std::string command = "-in " + mesh + " -met " + met + " -cad " + cad
      + " -verb " + std::to_string(diagnostic_verbosity())
      + " -vdepth " + std::to_string(diagnostic_vdepth())
      + " -adapt 0 -opt-niter 0 -adp-opt-niter 0";
  if(const char* extra = std::getenv("METRIS_QUALITY_INSERT_EXTRA_ARGS")){
    command += " ";
    command += extra;
  }
  return command;
}

std::string input_command_from_files(const std::string& mesh,
                                     const std::string& met,
                                     const std::string& cad){
  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(mesh),
                        "Missing diagnostic mesh: " + mesh);
  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(met),
                        "Missing diagnostic metric: " + met);
  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(cad),
                        "Missing diagnostic CAD: " + cad);

  std::string command = "-in " + mesh + " -met " + met + " -cad " + cad
      + " -verb 0 -vdepth 0 -adapt 0 -opt-niter 0 -adp-opt-niter 0";
  if(const char* extra = std::getenv("METRIS_OBJECTIVE_COMPARE_EXTRA_ARGS")){
    command += " ";
    command += extra;
  }
  return command;
}

std::string objective_compare_case_dir(){
  if(const char* env_dir = std::getenv("METRIS_OBJECTIVE_COMPARE_CASE_DIR")){
    return std::string(env_dir);
  }
  return "/Users/renat/MIT/HOMeshing/"
         "L2Project_MOESS_ConfBasedXClassicInsertions/"
         "L2Project_MOESS_CavSmooth_0Pctg/CornSing_CG/"
         "L2Project_2D_QualBased/Corner/Metris_CG_000500_P1";
}

int objective_compare_iteration(){
  if(const char* env_iter = std::getenv("METRIS_OBJECTIVE_COMPARE_ITERATION")){
    std::string iter(env_iter);
    if(!iter.empty() && iter[0] == 'a') iter.erase(iter.begin());
    return std::atoi(iter.c_str());
  }
  return 6;
}

double objective_compare_length_threshold(){
  if(const char* env_len = std::getenv("METRIS_OBJECTIVE_COMPARE_LENGTH_THRESHOLD")){
    return std::atof(env_len);
  }
  return std::sqrt(2.);
}

int objective_compare_max_candidates(){
  if(const char* env_max = std::getenv("METRIS_OBJECTIVE_COMPARE_MAX_CANDIDATES")){
    return std::atoi(env_max);
  }
  return 0;
}

double objective_compare_min_child_area_fraction(){
  if(const char* env_min =
      std::getenv("METRIS_OBJECTIVE_COMPARE_MIN_CHILD_AREA_FRACTION")){
    return std::atof(env_min);
  }
  return 0.;
}

double objective_compare_step_barrier_beta(){
  if(const char* env_beta =
      std::getenv("METRIS_OBJECTIVE_COMPARE_STEP_BARRIER_BETA")){
    return std::max(0.,std::atof(env_beta));
  }
  return 0.;
}

double objective_compare_step_barrier_rho0(){
  if(const char* env_rho0 =
      std::getenv("METRIS_OBJECTIVE_COMPARE_STEP_BARRIER_RHO0")){
    return std::max(0.,std::atof(env_rho0));
  }
  return 0.1;
}

std::string objective_compare_cad_file(const std::string& dir){
  if(const char* env_cad = std::getenv("METRIS_OBJECTIVE_COMPARE_CAD")){
    return std::string(env_cad);
  }
  const std::string local_cad = dir + "/CAD_MOESS.egads";
  if(std::filesystem::exists(local_cad)) return local_cad;
  return "/Users/renat/MIT/HOMeshing/"
         "L2Project_MOESS_ConfBasedXClassicInsertions/CAD_MOESS.egads";
}

std::string objective_compare_input_command(){
  const std::string dir = objective_compare_case_dir();
  const std::string suffix = "a" + std::to_string(objective_compare_iteration());
  return input_command_from_files(dir + "/mesh_MOESS_initial_" + suffix + ".meshb",
                                  dir + "/met_MOESS_initial_" + suffix + ".solb",
                                  objective_compare_cad_file(dir));
}

std::filesystem::path objective_compare_output_dir(){
  if(const char* output =
      std::getenv("METRIS_OBJECTIVE_COMPARE_OUTPUT_DIR")){
    return output;
  }
  return std::filesystem::current_path() / "build" / "codex_diagnostic"
                                      / "objective_insert_comparison";
}

std::filesystem::path trace_dir(){
  return std::filesystem::current_path() / "build" / "codex_diagnostic" / "quality_insert_divergence_trace";
}

std::filesystem::path vizir_dir(){
  return std::filesystem::current_path() / "build" / "codex_diagnostic" / "quality_insert_divergence_vizir";
}

std::filesystem::path runmetris_dir(){
  return std::filesystem::current_path() / "build" / "codex_diagnostic"
                                      / "quality_insert_runmetris";
}

std::string production_run_command(){
  const std::string dir = case_dir();
  const std::string adapt_suffix = "a" + std::to_string(adaptive_iteration());
  const std::string mesh = dir + "/mesh_MOESS_initial_" + adapt_suffix + ".meshb";
  const std::string met  = dir + "/met_MOESS_initial_" + adapt_suffix + ".solb";
  const std::string cad  = cad_file(dir);
  const std::filesystem::path outdir = runmetris_dir();

  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(mesh),
                        "Missing diagnostic mesh: " + mesh);
  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(met),
                        "Missing diagnostic metric: " + met);
  BOOST_REQUIRE_MESSAGE(std::filesystem::exists(cad),
                        "Missing diagnostic CAD: " + cad);

  std::filesystem::create_directories(outdir);
  std::filesystem::remove(outdir / ("outputAdaptStats_MOESS_a"
                                  + std::to_string(adaptive_iteration())
                                  + ".txt"));

  return "-in " + mesh + " -met " + met + " -cad " + cad
       + " -prefix " + outdir.string() + "/"
       + " -out debug_Metris_L2proj"
       + " -verb 0 -vdepth 0"
       + " -opt-niter 5 -opt-pnorm 1 -opt-power 1"
       + " -adp-opt-niter 1"
       + " -MOESS_adapt_it " + std::to_string(adaptive_iteration());
}

template<class MFT, int gdim, int ideg>
void setup_quality_mesh(Mesh<MFT>& msh){
  BOOST_REQUIRE_EQUAL(msh.idim, gdim);
  BOOST_REQUIRE_EQUAL(msh.curdeg, ideg);
  msh.met.setSpace(MetSpace::Exp);
  msh.setBasis(FEBasis::Lagrange);
}

template<class MFT>
void log_point_position(const std::string& label,
                        const Mesh<MFT>& msh,
                        int ipoin,
                        std::ofstream& log){
  log << label << " ipoin " << ipoin << " coord [";
  for(int ii = 0; ii < msh.idim; ii++){
    if(ii > 0) log << " ";
    log << msh.coord(ipoin,ii);
  }
  log << "]\n";
}

template<class MFT>
void write_cavity_trace(const std::filesystem::path& dir,
                        const std::string& stage,
                        Mesh<MFT>& msh,
                        const MshCavity& cav,
                        std::ofstream& log){
  const std::filesystem::path mesh_path = dir / (stage + ".meshb");
  writeMeshCavity(mesh_path.string(), msh, cav);

  log << stage
      << " ipins " << cav.ipins
      << " nedge " << cav.lcedg.get_n()
      << " nface " << cav.lcfac.get_n()
      << " ntet " << cav.lctet.get_n()
      << "\n";
  log << "  lcedg " << cav.lcedg << "\n";
  log << "  lcfac " << cav.lcfac << "\n";
  log << "  lctet " << cav.lctet << "\n";
  log_point_position("  ipins_position", msh, cav.ipins, log);
  log << "  mesh " << mesh_path << "\n";
}

template<class MFT>
void write_painted_point_meshes(const std::filesystem::path& dir,
                                const std::string& stage,
                                Mesh<MFT>& msh,
                                const MshCavity& cav,
                                std::ofstream& log){
  std::filesystem::create_directories(dir);

  log_point_position(stage + " final_ipins_position", msh, cav.ipins, log);

  const int ipdbg = msh.newpoitopo(PointType::Vertex,-1,-1);
  msh.newbpotopo(Vertex{ipdbg},0,ipdbg);
  for(int ii = 0; ii < msh.idim; ii++){
    msh.coord(ipdbg,ii) = msh.coord(cav.ipins,ii);
  }

  const std::filesystem::path mesh_path = dir / (stage + "_painted_point.meshb");
  const std::filesystem::path cavity_path = dir / (stage + "_painted_point_cavity.meshb");
  writeMesh(mesh_path.string(), msh, false);
  writeMeshCavity(cavity_path.string(), msh, cav);

  log << stage
      << " painted_point ipdbg " << ipdbg
      << " from_ipins " << cav.ipins
      << "\n";
  log_point_position("  painted_point_position", msh, ipdbg, log);
  log << "  painted_mesh " << mesh_path << "\n";
  log << "  painted_cavity_mesh " << cavity_path << "\n";

  msh.killpoint(ipdbg);
}

template<class MFT>
void write_painted_cavity_context(const std::filesystem::path& dir,
                                  const std::string& stage,
                                  Mesh<MFT>& msh,
                                  const MshCavity& cav,
                                  std::ofstream& log){
  std::filesystem::create_directories(dir);

  intAr1 old_edg_refs(cav.lcedg.get_n());
  intAr1 old_fac_refs(cav.lcfac.get_n());
  old_edg_refs.set_n(cav.lcedg.get_n());
  old_fac_refs.set_n(cav.lcfac.get_n());

  for(int ii = 0; ii < cav.lcedg.get_n(); ii++){
    const int iedge = cav.lcedg[ii];
    old_edg_refs[ii] = msh.edg2ref[iedge];
    msh.edg2ref[iedge] = 9001;
  }
  for(int ii = 0; ii < cav.lcfac.get_n(); ii++){
    const int iface = cav.lcfac[ii];
    old_fac_refs[ii] = msh.fac2ref[iface];
    msh.fac2ref[iface] = 9002 + ii;
  }

  const std::filesystem::path mesh_path = dir / (stage + "_cavity_context.meshb");
  writeMesh(mesh_path.string(), msh, false);

  for(int ii = 0; ii < cav.lcedg.get_n(); ii++){
    msh.edg2ref[cav.lcedg[ii]] = old_edg_refs[ii];
  }
  for(int ii = 0; ii < cav.lcfac.get_n(); ii++){
    msh.fac2ref[cav.lcfac[ii]] = old_fac_refs[ii];
  }

  log << stage
      << " painted_cavity_context " << mesh_path
      << " edge_ref 9001 face_refs_start 9002"
      << "\n";
}

const char* cavity_error_name(int ierro){
  switch(ierro){
  case CAV_NOERR: return "CAV_NOERR";
  case CAV_ERR_NOBPO: return "CAV_ERR_NOBPO";
  case CAV_ERR_TDIMN: return "CAV_ERR_TDIMN";
  case CAV_ERR_DUPEDG: return "CAV_ERR_DUPEDG";
  case CAV_ERR_INTEDG: return "CAV_ERR_INTEDG";
  case CAV_ERR_FLATFAC: return "CAV_ERR_FLATFAC";
  case CAV_ERR_NEGFAC: return "CAV_ERR_NEGFAC";
  case CAV_ERR_DUPEDG2: return "CAV_ERR_DUPEDG2";
  case CAV_ERR_QMAXNEC: return "CAV_ERR_QMAXNEC";
  case CAV_ERR_QMAXIFF: return "CAV_ERR_QMAXIFF";
  case CAV_ERR_QFACNEG: return "CAV_ERR_QFACNEG";
  case CAV_ERR_CADFAR: return "CAV_ERR_CADFAR";
  case CAV_ERR_INCORRECTIBLE: return "CAV_ERR_INCORRECTIBLE";
  case CAV_ERR_DRYFAIL1: return "CAV_ERR_DRYFAIL1";
  case CAV_ERR_DRYFAIL2: return "CAV_ERR_DRYFAIL2";
  case CAV_ERR_LINETOPO: return "CAV_ERR_LINETOPO";
  case CAV_ERR_GEODEVLIN: return "CAV_ERR_GEODEVLIN";
  case CAV_ERR_FLATEDG: return "CAV_ERR_FLATEDG";
  case CAV_ERR_DUPFAC: return "CAV_ERR_DUPFAC";
  case CAV_ERR_LINFAC: return "CAV_ERR_LINFAC";
  case CAV_ERR_FLATTET: return "CAV_ERR_FLATTET";
  case CAV_ERR_MULTIREFTET: return "CAV_ERR_MULTIREFTET";
  case CAV_ERR_NOFACTET: return "CAV_ERR_NOFACTET";
  case CAV_ERR_INTERPMETBACK: return "CAV_ERR_INTERPMETBACK";
  case CAV_ERR_REMPT: return "CAV_ERR_REMPT";
  case CAV_ERR_NOEDGFAC: return "CAV_ERR_NOEDGFAC";
  case CAV_ERR_INTFAC: return "CAV_ERR_INTFAC";
  case CAV_ERR_BDRYTET: return "CAV_ERR_BDRYTET";
  case CAV_ERR_BDRYTET2: return "CAV_ERR_BDRYTET2";
  case CAV_ERR_CORRECTCAV: return "CAV_ERR_CORRECTCAV";
  case CAV_ERR_OBJECTIVE: return "CAV_ERR_OBJECTIVE";
  default: return "CAV_ERR_UNKNOWN";
  }
}

template<class MFT, int ideg>
void trace_raw_cavity_operator(Mesh<MFT>& msh,
                               MshCavity& cav,
                               const std::string& label,
                               std::ofstream& log,
                               int ithrd){
  CavOprOpt opts = diagnostic_cavity_options();
  CavWrkArrs work;
  CavOprInfo info;
  cav.inewp = 1;

  const int ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd);
  log << label
      << " raw_cavity_operator ierro " << ierro
      << " " << cavity_error_name(ierro)
      << " info.done " << info.done
      << " qmax_end " << info.qmax_end
      << " lbad_n " << work.lbad.get_n()
      << "\n";
  if(work.lbad.get_n() > 0){
    log << label << " raw_cavity_operator lbad";
    for(int ii = 0; ii < work.lbad.get_n(); ii++){
      log << " [" << work.lbad(ii,0) << " " << work.lbad(ii,1) << "]";
    }
    log << "\n";
  }
}

template<class MFT>
bool point_has_edge_ref(Mesh<MFT>& msh, int ipoin, int iref){
  for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0;
      ibpoi = msh.bpo2ibi(ibpoi,3)){
    if(msh.bpo2ibi(ibpoi,1) != 1) continue;
    const int iedge = msh.bpo2ibi(ibpoi,2);
    if(iedge >= 0 && msh.edg2ref[iedge] == iref) return true;
  }
  return false;
}

template<class MFT>
void trace_linfac_candidates(Mesh<MFT>& msh,
                             MshCavity& cav,
                             const std::string& label,
                             std::ofstream& log,
                             int ithrd){
  const int tag0 = ++msh.tag[ithrd];
  for(int iedge : cav.lcedg) msh.edg2tag(ithrd,iedge) = tag0;
  for(int iface : cav.lcfac) msh.fac2tag(ithrd,iface) = tag0;

  const int ibins = msh.poi2bpo[cav.ipins];
  const int iedins = ibins >= 0 ? msh.bpo2ibi(ibins,2) : -1;
  const int irefins = iedins >= 0 ? msh.edg2ref[iedins] : -1;

  log << label << " linfac_probe ipins " << cav.ipins
      << " edge_ref " << irefins << "\n";

  for(int iface : cav.lcfac){
    for(int iedge = 0; iedge < 3; iedge++){
      const int ifac2 = msh.fac2fac(iface,iedge);
      if(ifac2 >= 0 && msh.fac2tag(ithrd,ifac2) >= tag0) continue;
      const int iedgeGlobal = msh.facedg2glo(iface,iedge);
      if(iedgeGlobal >= 0 && msh.edg2tag(ithrd,iedgeGlobal) >= tag0) continue;

      const int ip1 = msh.fac2poi(iface,lnoed2[iedge][0]);
      const int ip2 = msh.fac2poi(iface,lnoed2[iedge][1]);
      if(ip1 == cav.ipins || ip2 == cav.ipins) continue;

      const int iref = iedgeGlobal >= 0 ? msh.edg2ref[iedgeGlobal] : irefins;
      const bool same_ref =
        iref >= 0
        && point_has_edge_ref(msh,cav.ipins,iref)
        && point_has_edge_ref(msh,ip1,iref)
        && point_has_edge_ref(msh,ip2,iref);

      log << label
          << " boundary_new_face from_face " << iface
          << " local_edge " << iedge
          << " global_edge " << iedgeGlobal
          << " edge_ref " << iref
          << " vertices [" << cav.ipins << " " << ip1 << " " << ip2 << "]"
          << " all_on_same_edge_ref " << same_ref
          << "\n";
    }
  }
}

template<class MFT, QuaFun iquaf>
CavityQualityStats compute_cavity_quality_2d(Mesh<MFT>& msh,
                                             MshCavity& cav,
                                             int ithrd,
                                             const std::string& label,
                                             std::ofstream* log){
  auto update_stats = [](double qua, int& nqua, double& qmin,
                         double& qmax, double& qsum){
    qmin = MIN(qmin, qua);
    qmax = MAX(qmax, qua);
    qsum += qua;
    nqua++;
  };

  const int tdim = 2;
  intAr1& lcent = cav.lcent(tdim);
  intAr2& ent2poi = msh.ent2poi(tdim);
  const intAr2& ent2ent = msh.ent2ent(tdim);
  intAr2r& ent2tag = msh.ent2tag(tdim);

  CavityQualityStats stats;
  stats.qmin0 =  std::numeric_limits<double>::max();
  stats.qmin1 =  std::numeric_limits<double>::max();
  stats.qmax0 = -std::numeric_limits<double>::max();
  stats.qmax1 = -std::numeric_limits<double>::max();
  double difto = 1.;

  for(int iface : lcent){
    double qua = metqua<MFT,2,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,iface,difto);
    if(log != nullptr){
      *log << label
           << " original_face " << iface
           << " nodes [" << ent2poi(iface,0) << " "
           << ent2poi(iface,1) << " " << ent2poi(iface,2) << "]"
           << " quality " << qua << "\n";
    }
    update_stats(qua, stats.nqua0, stats.qmin0, stats.qmax0, stats.qsum0);
  }

  msh.tag[ithrd]++;
  const int tag = msh.tag[ithrd];
  for(int iface : lcent){
    ent2tag(ithrd, iface) = tag;
  }

  const int nface0 = msh.nface;
  msh.set_nface(nface0 + 1);
  const int tmpFace = nface0;
  ent2poi(tmpFace,0) = cav.ipins;

  for(int iface : lcent){
    for(int iedge = 0; iedge < 3; iedge++){
      int ineighbor = ent2ent(iface, iedge);
      if(ineighbor >= 0 && ent2tag(ithrd, ineighbor) >= tag) continue;

      ent2poi(tmpFace,1) = ent2poi(iface,lnoed2[iedge][0]);
      ent2poi(tmpFace,2) = ent2poi(iface,lnoed2[iedge][1]);
      if(ent2poi(tmpFace,1) == cav.ipins || ent2poi(tmpFace,2) == cav.ipins){
        if(log != nullptr){
          *log << label
               << " skipped_reconnected_face boundary_of_face " << iface
               << " edge " << iedge
               << " because it already contains ipins\n";
        }
        continue;
      }

      double meas = 0;
      if(!isvalideltP1<2,2>(msh, tmpFace, NULL, &meas)){
        if(log != nullptr){
          *log << label
               << " invalid_reconnected_face boundary_of_face " << iface
               << " edge " << iedge
               << " meas " << meas << "\n";
        }
        continue;
      }

      double qua = metqua<MFT,2,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,tmpFace,difto);
      if(log != nullptr){
        *log << label
             << " reconnected_face from_face " << iface
             << " boundary_edge " << iedge
             << " nodes [" << ent2poi(tmpFace,0) << " "
             << ent2poi(tmpFace,1) << " " << ent2poi(tmpFace,2) << "]"
             << " measure " << meas
             << " quality " << qua << "\n";
      }
      update_stats(qua, stats.nqua1, stats.qmin1, stats.qmax1, stats.qsum1);
    }
  }

  ent2poi(tmpFace,0) = -1;
  msh.set_nface(nface0);

  if(stats.nqua0 == 0){
    stats.qmin0 = stats.qmax0 = stats.qsum0 = 0;
  }
  if(stats.nqua1 == 0){
    stats.qmin1 = stats.qmax1 = stats.qsum1 = 0;
  }
  return stats;
}

inline constexpr std::array<double,12> step_distance_p_values = {
  1.,1.1,1.2,1.25,1.3,1.35,1.4,1.5,1.75,2.,3.,4.
};

struct TriangleIntegrationValues{
  double size_shape_one_point = 0.;
  double size_shape_vertex_barycenter = 0.;
  double step_distance_one_point = 0.;
  double step_distance_vertex_barycenter = 0.;
  std::array<double,step_distance_p_values.size()> step_p_one_point = {};
  std::array<double,step_distance_p_values.size()> step_p_vertex_barycenter = {};
};

struct CavityIntegrationValues{
  TriangleIntegrationValues original;
  TriangleIntegrationValues reconnected;
};

void add_triangle_integration_values(TriangleIntegrationValues& sum,
                                     const TriangleIntegrationValues& value){
  sum.size_shape_one_point += value.size_shape_one_point;
  sum.size_shape_vertex_barycenter += value.size_shape_vertex_barycenter;
  sum.step_distance_one_point += value.step_distance_one_point;
  sum.step_distance_vertex_barycenter += value.step_distance_vertex_barycenter;
  for(std::size_t ip = 0; ip < step_distance_p_values.size(); ip++){
    sum.step_p_one_point[ip] += value.step_p_one_point[ip];
    sum.step_p_vertex_barycenter[ip] += value.step_p_vertex_barycenter[ip];
  }
}

template<class MFT>
TriangleIntegrationValues write_triangle_objective_decomposition_2d(
    Mesh<MFT>& msh,
    int iface,
    const std::string& label,
    std::ofstream& log){
  const int* nodes = msh.fac2poi[iface];
  double meas = 0.;
  BOOST_REQUIRE((isvalideltP1<2,2>(msh,iface,nullptr,&meas)));

  double bary_center[3] = {1./3.,1./3.,1./3.};
  double metric_center[3] = {};
  for(int inode = 0; inode < 3; inode++){
    for(int imet = 0; imet < 3; imet++){
      metric_center[imet] += msh.met(nodes[inode],imet) / 3.;
    }
  }
  const double size_shape_raw =
      quafun_sizeshape<MFT,2,2,double>(msh,AsDeg::P1,AsDeg::P1,
                                      nodes,bary_center,metric_center);
  const double sqrt_det_metric_center = std::sqrt(
      metric_center[0]*metric_center[2]
      - metric_center[1]*metric_center[1]);
  const double size_shape_integrand = std::abs(size_shape_raw - 1.);
  const double size_shape_contribution =
      size_shape_integrand * meas * sqrt_det_metric_center;

  log << label << " nodes [" << nodes[0] << " " << nodes[1] << " "
      << nodes[2] << "] measure " << meas
      << " sizeshape_raw " << size_shape_raw
      << " sizeshape_abs_error " << size_shape_integrand
      << " sizeshape_sqrt_det_metric " << sqrt_det_metric_center
      << " sizeshape_contribution " << size_shape_contribution << "\n";

  TriangleIntegrationValues values;
  for(int iquad = 0; iquad < 4; iquad++){
    double bary[3] = {};
    if(iquad < 3){
      bary[iquad] = 1.;
    }else{
      bary[0] = bary[1] = bary[2] = 1./3.;
    }

    double coopr[2];
    double jmat[4];
    eval2<2,1>(msh.coord,nodes,msh.getBasis(),
               DifVar::Bary,DifVar::None,bary,coopr,jmat,nullptr);

    double metric[3];
    if(iquad < 3){
      for(int imet = 0; imet < 3; imet++){
        metric[imet] = msh.met(nodes[iquad],imet);
      }
    }else if constexpr(std::is_same<MFT,MetricFieldAnalytical>::value){
      msh.met.getMetPhys(DifVar::None,msh.met.getSpace(),
                         coopr,metric,nullptr);
    }else{
      msh.met.getMetBary(AsDeg::P1,DifVar::None,msh.met.getSpace(),
                         nodes,2,bary,metric,nullptr);
    }

    double jreg_t[4] = {};
    for(int ii = 0; ii < 2; ii++){
      for(int aa = 0; aa < 2; aa++){
        for(int kk = 0; kk < 2; kk++){
          jreg_t[2*ii+aa] +=
              Constants::invtJ_0[hana::type_c<double>][2][2*ii+kk]
              * jmat[2*kk+aa];
        }
      }
    }

    double afull[4] = {};
    for(int ii = 0; ii < 2; ii++){
      for(int jj = 0; jj < 2; jj++){
        afull[2*ii+jj] =
            jreg_t[2*ii] * metric[0] * jreg_t[2*jj]
          + jreg_t[2*ii] * metric[1] * jreg_t[2*jj+1]
          + jreg_t[2*ii+1] * metric[1] * jreg_t[2*jj]
          + jreg_t[2*ii+1] * metric[2] * jreg_t[2*jj+1];
      }
    }
    double apacked[3] = {afull[0],afull[1],afull[3]};
    double eigenvalues[2];
    double eigenvectors[4];
    geteigsym<2,double>(apacked,eigenvalues,eigenvectors);
    const double log0 = std::log(eigenvalues[0]);
    const double log1 = std::log(eigenvalues[1]);
    const double phi = log0*log0 + log1*log1;
    const double distance = std::sqrt(phi);
    const double size_shape_sample =
        quafun_sizeshape<MFT,2,2,double>(msh,AsDeg::P1,AsDeg::P1,
                                        nodes,bary,metric);
    const double size_shape_error = std::abs(size_shape_sample - 1.);
    double theta = 0.;
    VolumeMeasureHelpers::eval_theta_fixed_metric_grad<2,2,double>(
        jreg_t,metric,nullptr,&theta,nullptr);
    const double contribution = phi*theta/4.;
    const double size_shape_contribution_vb = size_shape_error*theta/4.;
    values.step_distance_vertex_barycenter += contribution;
    values.size_shape_vertex_barycenter += size_shape_contribution_vb;
    for(std::size_t ip = 0; ip < step_distance_p_values.size(); ip++){
      const double distance_power = std::pow(distance,step_distance_p_values[ip]);
      values.step_p_vertex_barycenter[ip] += distance_power*theta/4.;
      if(iquad == 3){
        values.step_p_one_point[ip] = distance_power*theta;
      }
    }
    if(iquad == 3){
      values.step_distance_one_point = phi*theta;
      values.size_shape_one_point = size_shape_error*theta;
    }

    log << label << " step_sample " << iquad
        << " bary [" << bary[0] << " " << bary[1] << " " << bary[2] << "]"
        << " eigenvalues [" << eigenvalues[0] << " " << eigenvalues[1] << "]"
        << " phi " << phi
        << " sizeshape_raw " << size_shape_sample
        << " sizeshape_abs_error " << size_shape_error
        << " theta " << theta
        << " step_weighted_contribution " << contribution
        << " sizeshape_weighted_contribution "
        << size_shape_contribution_vb << "\n";
  }
  log << label
      << " integration_values"
      << " sizeshape_one_point " << values.size_shape_one_point
      << " sizeshape_vertex_barycenter "
      << values.size_shape_vertex_barycenter
      << " stepdistance_one_point " << values.step_distance_one_point
      << " stepdistance_vertex_barycenter "
      << values.step_distance_vertex_barycenter << "\n";
  return values;
}

template<class MFT>
CavityIntegrationValues write_cavity_objective_decomposition_2d(
    Mesh<MFT>& msh,
    MshCavity& cav,
    int ithrd,
    const std::string& label,
    std::ofstream& log){
  intAr1& cavity_faces = cav.lcfac;
  intAr2& fac2poi = msh.fac2poi;
  const intAr2& fac2fac = msh.fac2fac;
  intAr2r& fac2tag = msh.fac2tag;

  const int tag = ++msh.tag[ithrd];
  for(int iface : cavity_faces) fac2tag(ithrd,iface) = tag;

  CavityIntegrationValues values;
  for(int iface : cavity_faces){
    const TriangleIntegrationValues triangle =
        write_triangle_objective_decomposition_2d(
            msh,iface,label + " original_face_" + std::to_string(iface),log);
    add_triangle_integration_values(values.original,triangle);
  }

  const int nface0 = msh.nface;
  msh.set_nface(nface0 + 1);
  const int tmp_face = nface0;
  fac2poi(tmp_face,0) = cav.ipins;
  for(int iface : cavity_faces){
    for(int iedge = 0; iedge < 3; iedge++){
      const int neighbor = fac2fac(iface,iedge);
      if(neighbor >= 0 && fac2tag(ithrd,neighbor) >= tag) continue;
      fac2poi(tmp_face,1) = fac2poi(iface,lnoed2[iedge][0]);
      fac2poi(tmp_face,2) = fac2poi(iface,lnoed2[iedge][1]);
      if(fac2poi(tmp_face,1) == cav.ipins
          || fac2poi(tmp_face,2) == cav.ipins) continue;
      const TriangleIntegrationValues triangle =
          write_triangle_objective_decomposition_2d(
              msh,tmp_face,
              label + " reconnected_from_" + std::to_string(iface)
              + "_edge_" + std::to_string(iedge),log);
      add_triangle_integration_values(values.reconnected,triangle);
    }
  }
  fac2poi(tmp_face,0) = -1;
  msh.set_nface(nface0);
  log << label << " cavity_integration_summary"
      << " sizeshape_one_point "
      << values.original.size_shape_one_point << " -> "
      << values.reconnected.size_shape_one_point
      << " sizeshape_vertex_barycenter "
      << values.original.size_shape_vertex_barycenter << " -> "
      << values.reconnected.size_shape_vertex_barycenter
      << " stepdistance_one_point "
      << values.original.step_distance_one_point << " -> "
      << values.reconnected.step_distance_one_point
      << " stepdistance_vertex_barycenter "
      << values.original.step_distance_vertex_barycenter << " -> "
      << values.reconnected.step_distance_vertex_barycenter << "\n";
  for(std::size_t ip = 0; ip < step_distance_p_values.size(); ip++){
    log << label << " stepdistance_p_summary p "
        << step_distance_p_values[ip]
        << " one_point " << values.original.step_p_one_point[ip]
        << " -> " << values.reconnected.step_p_one_point[ip]
        << " vertex_barycenter "
        << values.original.step_p_vertex_barycenter[ip]
        << " -> " << values.reconnected.step_p_vertex_barycenter[ip]
        << "\n";
  }
  return values;
}

template<class MFT>
bool cavity_reconnection_is_valid_2d(Mesh<MFT>& msh,
                                      MshCavity& cav,
                                      int ithrd){
  const int tag = ++msh.tag[ithrd];
  for(int iface : cav.lcfac) msh.fac2tag(ithrd,iface) = tag;

  const int nface0 = msh.nface;
  msh.set_nface(nface0 + 1);
  const int tmp_face = nface0;
  msh.fac2poi(tmp_face,0) = cav.ipins;
  bool valid = true;
  for(int iface : cav.lcfac){
    for(int iedge = 0; iedge < 3; iedge++){
      const int neighbor = msh.fac2fac(iface,iedge);
      if(neighbor >= 0 && msh.fac2tag(ithrd,neighbor) >= tag) continue;
      msh.fac2poi(tmp_face,1) = msh.fac2poi(iface,lnoed2[iedge][0]);
      msh.fac2poi(tmp_face,2) = msh.fac2poi(iface,lnoed2[iedge][1]);
      if(msh.fac2poi(tmp_face,1) == cav.ipins
          || msh.fac2poi(tmp_face,2) == cav.ipins) continue;
      double measure = 0.;
      if(!isvalideltP1<2,2>(msh,tmp_face,nullptr,&measure)){
        valid = false;
        break;
      }
    }
    if(!valid) break;
  }
  msh.fac2poi(tmp_face,0) = -1;
  msh.set_nface(nface0);
  return valid;
}

struct StepPOptimizationResult{
  double original = 0.;
  double initial = 0.;
  double optimized = 0.;
  double x = 0.;
  double y = 0.;
};

template<class MFT>
StepPOptimizationResult optimize_step_p_point_2d(Mesh<MFT>& msh,
                                                  MshCavity& cav,
                                                  int ithrd,
                                                  std::size_t ip_value,
                                                  bool vertex_barycenter){
  static std::ofstream null_log("/dev/null");
  const double x0 = msh.coord(cav.ipins,0);
  const double y0 = msh.coord(cav.ipins,1);

  auto evaluate = [&](double x, double y,
                      double* original_value = nullptr) -> double{
    msh.coord(cav.ipins,0) = x;
    msh.coord(cav.ipins,1) = y;
    if(!cavity_reconnection_is_valid_2d(msh,cav,ithrd)){
      return std::numeric_limits<double>::infinity();
    }
    const CavityIntegrationValues values =
        write_cavity_objective_decomposition_2d(
            msh,cav,ithrd,"optimization_probe",null_log);
    if(original_value != nullptr){
      *original_value = vertex_barycenter
          ? values.original.step_p_vertex_barycenter[ip_value]
          : values.original.step_p_one_point[ip_value];
    }
    return vertex_barycenter
        ? values.reconnected.step_p_vertex_barycenter[ip_value]
        : values.reconnected.step_p_one_point[ip_value];
  };

  double xmin = std::numeric_limits<double>::max();
  double xmax = -std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max();
  double ymax = -std::numeric_limits<double>::max();
  std::set<int> cavity_points;
  for(int iface : cav.lcfac){
    for(int inode = 0; inode < 3; inode++){
      cavity_points.insert(msh.fac2poi(iface,inode));
    }
  }
  for(int ipoin : cavity_points){
    xmin = std::min(xmin,msh.coord(ipoin,0));
    xmax = std::max(xmax,msh.coord(ipoin,0));
    ymin = std::min(ymin,msh.coord(ipoin,1));
    ymax = std::max(ymax,msh.coord(ipoin,1));
  }

  StepPOptimizationResult result;
  result.x = x0;
  result.y = y0;
  result.initial = evaluate(x0,y0,&result.original);
  result.optimized = result.initial;

  double step = 0.1*std::max(xmax-xmin,ymax-ymin);
  constexpr double directions[8][2] = {
    {1.,0.},{-1.,0.},{0.,1.},{0.,-1.},
    {1.,1.},{1.,-1.},{-1.,1.},{-1.,-1.}
  };
  for(int iteration = 0; iteration < 200 && step > 1.e-12; iteration++){
    bool improved = false;
    for(const auto& direction : directions){
      const double x = result.x + step*direction[0];
      const double y = result.y + step*direction[1];
      const double value = evaluate(x,y);
      if(value < result.optimized){
        result.optimized = value;
        result.x = x;
        result.y = y;
        improved = true;
      }
    }
    if(!improved) step *= 0.5;
  }

  msh.coord(cav.ipins,0) = x0;
  msh.coord(cav.ipins,1) = y0;
  return result;
}

double select_common_objective(const CommonQuadratureValues& values,
                               CommonScanObjective objective){
  switch(objective){
  case CommonScanObjective::SizeShapeP1:
    return values.size_shape_p1;
  case CommonScanObjective::StepDistanceP1:
    return values.step_distance_p1
         + objective_compare_step_barrier_beta()
             *values.step_distance_barrier;
  case CommonScanObjective::StepDistanceP2:
    return values.step_distance_p2;
  }
  return std::numeric_limits<double>::infinity();
}

void add_common_quadrature_values(CommonQuadratureValues& sum,
                                  const CommonQuadratureValues& value){
  sum.size_shape_p1 += value.size_shape_p1;
  sum.step_distance_p1 += value.step_distance_p1;
  sum.step_distance_p2 += value.step_distance_p2;
  sum.step_distance_barrier += value.step_distance_barrier;
  sum.min_metric_volume =
      std::min(sum.min_metric_volume,value.min_metric_volume);
}

template<class MFT>
CommonQuadratureValues evaluate_triangle_common_vb_2d(Mesh<MFT>& msh,
                                                       int iface){
  const int* nodes = msh.fac2poi[iface];
  CommonQuadratureValues values;

  for(int iquad = 0; iquad < 4; iquad++){
    double bary[3] = {};
    if(iquad < 3){
      bary[iquad] = 1.;
    }else{
      bary[0] = bary[1] = bary[2] = 1./3.;
    }

    double coopr[2];
    double jmat[4];
    eval2<2,1>(msh.coord,nodes,msh.getBasis(),
               DifVar::Bary,DifVar::None,bary,coopr,jmat,nullptr);

    double metric[3];
    if(iquad < 3){
      for(int imet = 0; imet < 3; imet++){
        metric[imet] = msh.met(nodes[iquad],imet);
      }
    }else if constexpr(std::is_same<MFT,MetricFieldAnalytical>::value){
      msh.met.getMetPhys(DifVar::None,msh.met.getSpace(),
                         coopr,metric,nullptr);
    }else{
      msh.met.getMetBary(AsDeg::P1,DifVar::None,msh.met.getSpace(),
                         nodes,2,bary,metric,nullptr);
    }

    double jreg_t[4] = {};
    for(int ii = 0; ii < 2; ii++){
      for(int aa = 0; aa < 2; aa++){
        for(int kk = 0; kk < 2; kk++){
          jreg_t[2*ii+aa] +=
              Constants::invtJ_0[hana::type_c<double>][2][2*ii+kk]
              * jmat[2*kk+aa];
        }
      }
    }

    const double det_jreg = jreg_t[0]*jreg_t[3] - jreg_t[1]*jreg_t[2];
    const double det_metric = metric[0]*metric[2] - metric[1]*metric[1];
    if(!std::isfinite(det_jreg) || !std::isfinite(det_metric)
        || std::abs(det_jreg) <= 1.e-15 || det_metric <= 0.){
      values.size_shape_p1 = std::numeric_limits<double>::infinity();
      values.step_distance_p1 = std::numeric_limits<double>::infinity();
      values.step_distance_p2 = std::numeric_limits<double>::infinity();
      return values;
    }

    double afull[4] = {};
    for(int ii = 0; ii < 2; ii++){
      for(int jj = 0; jj < 2; jj++){
        afull[2*ii+jj] =
            jreg_t[2*ii] * metric[0] * jreg_t[2*jj]
          + jreg_t[2*ii] * metric[1] * jreg_t[2*jj+1]
          + jreg_t[2*ii+1] * metric[1] * jreg_t[2*jj]
          + jreg_t[2*ii+1] * metric[2] * jreg_t[2*jj+1];
      }
    }
    double apacked[3] = {afull[0],afull[1],afull[3]};
    double eigenvalues[2];
    double eigenvectors[4];
    geteigsym<2,double>(apacked,eigenvalues,eigenvectors);
    if(eigenvalues[0] <= 0. || eigenvalues[1] <= 0.){
      values.size_shape_p1 = std::numeric_limits<double>::infinity();
      values.step_distance_p1 = std::numeric_limits<double>::infinity();
      values.step_distance_p2 = std::numeric_limits<double>::infinity();
      return values;
    }

    const double log0 = std::log(eigenvalues[0]);
    const double log1 = std::log(eigenvalues[1]);
    const double distance2 = log0*log0 + log1*log1;
    const double distance1 = std::sqrt(distance2);
    // This is quafun_sizeshape for tdim=2 and opt_power=1, expressed
    // directly through the already available eigenvalues. Besides avoiding
    // duplicate work, it lets an almost singular trial return a very large
    // finite/infinite cost instead of throwing inside quafun_tradet.
    const double trace_a = eigenvalues[0] + eigenvalues[1];
    const double det_a = eigenvalues[0]*eigenvalues[1];
    const double size_shape =
        trace_a*trace_a*(1. + 1./(det_a*det_a))/8.;
    const double size_shape_error = std::abs(size_shape - 1.);

    // In the full-dimensional 2D case this is exactly
    // sqrt(det(J^T J))*sqrt(det(M)), evaluated without forming the Gram
    // matrix. The determinant form avoids cancellation for nearly flat
    // infeasible optimization probes.
    const double theta = std::abs(det_jreg)*std::sqrt(det_metric);
    const double weight = theta/4.;
    const double rho0 = objective_compare_step_barrier_rho0();
    if(rho0 > 0. && theta < rho0){
      const double log_ratio = std::log(rho0/theta);
      const double log_ratio2 = log_ratio*log_ratio;
      values.step_distance_barrier += log_ratio2*log_ratio2/4.;
    }
    values.min_metric_volume = std::min(values.min_metric_volume,theta);
    values.size_shape_p1 += weight*size_shape_error;
    values.step_distance_p1 += weight*distance1;
    values.step_distance_p2 += weight*distance2;
  }
  return values;
}

template<class MFT>
CommonCavityValues evaluate_cavity_common_vb_2d(Mesh<MFT>& msh,
                                                 MshCavity& cav,
                                                 int ithrd){
  CommonCavityValues values;
  const int tag = ++msh.tag[ithrd];
  for(int iface : cav.lcfac) msh.fac2tag(ithrd,iface) = tag;

  for(int iface : cav.lcfac){
    add_common_quadrature_values(values.original,
                                 evaluate_triangle_common_vb_2d(msh,iface));
  }

  const bool point_on_edge = cav.lcedg.get_n() > 0;
  int insertion_ref = -1;
  if(point_on_edge){
    const int insertion_edge = cav.lcedg[0];
    if(insertion_edge >= 0) insertion_ref = msh.edg2ref[insertion_edge];
  }

  const int nface0 = msh.nface;
  msh.set_nface(nface0 + 1);
  const int tmp_face = nface0;
  msh.fac2poi(tmp_face,0) = cav.ipins;

  for(int iface : cav.lcfac){
    for(int iedge = 0; iedge < 3; iedge++){
      const int neighbor = msh.fac2fac(iface,iedge);
      if(neighbor >= 0 && msh.fac2tag(ithrd,neighbor) >= tag) continue;

      if(point_on_edge){
        const int global_edge = msh.fac2edg(iface,iedge);
        if(global_edge >= 0 && msh.edg2ref[global_edge] == insertion_ref){
          continue;
        }
      }

      msh.fac2poi(tmp_face,1) = msh.fac2poi(iface,lnoed2[iedge][0]);
      msh.fac2poi(tmp_face,2) = msh.fac2poi(iface,lnoed2[iedge][1]);
      if(msh.fac2poi(tmp_face,1) == cav.ipins
          || msh.fac2poi(tmp_face,2) == cav.ipins){
        continue;
      }

      double measure = 0.;
      if(!isvalideltP1<2,2>(msh,tmp_face,nullptr,&measure)){
        values.valid = false;
        break;
      }
      add_common_quadrature_values(
          values.reconnected,evaluate_triangle_common_vb_2d(msh,tmp_face));
    }
    if(!values.valid) break;
  }

  msh.fac2poi(tmp_face,0) = -1;
  msh.set_nface(nface0);
  return values;
}

template<class MFT>
double common_cavity_min_child_area_fraction_2d(Mesh<MFT>& msh,
                                                 const MshCavity& cav,
                                                 double x,
                                                 double y){
  std::set<int> cavity_faces;
  double cavity_area = 0.;
  for(int iface : cav.lcfac){
    cavity_faces.insert(iface);
    const int ip0 = msh.fac2poi(iface,0);
    const int ip1 = msh.fac2poi(iface,1);
    const int ip2 = msh.fac2poi(iface,2);
    const double twice_area =
        (msh.coord(ip1,0) - msh.coord(ip0,0))
          *(msh.coord(ip2,1) - msh.coord(ip0,1))
      - (msh.coord(ip1,1) - msh.coord(ip0,1))
          *(msh.coord(ip2,0) - msh.coord(ip0,0));
    cavity_area += 0.5*std::abs(twice_area);
  }

  const bool inserted_on_edge = cav.lcedg.get_n() > 0;
  int insertion_ref = -1;
  if(inserted_on_edge && cav.lcedg[0] >= 0){
    insertion_ref = msh.edg2ref[cav.lcedg[0]];
  }
  double min_child_area = std::numeric_limits<double>::infinity();
  for(int iface : cav.lcfac){
    for(int iedge = 0; iedge < 3; iedge++){
      const int neighbor = msh.fac2fac(iface,iedge);
      if(neighbor >= 0 && cavity_faces.count(neighbor) > 0) continue;
      if(inserted_on_edge){
        const int global_edge = msh.fac2edg(iface,iedge);
        if(global_edge >= 0 && msh.edg2ref[global_edge] == insertion_ref){
          continue;
        }
      }
      const int ip1 = msh.fac2poi(iface,lnoed2[iedge][0]);
      const int ip2 = msh.fac2poi(iface,lnoed2[iedge][1]);
      const double twice_area =
          (msh.coord(ip1,0) - x)*(msh.coord(ip2,1) - y)
        - (msh.coord(ip1,1) - y)*(msh.coord(ip2,0) - x);
      min_child_area = std::min(min_child_area,0.5*std::abs(twice_area));
    }
  }
  return cavity_area > 0. && std::isfinite(min_child_area)
      ? min_child_area/cavity_area : 0.;
}

template<class MFT>
CommonObjectiveOptimization optimize_common_objective_point_2d(
    Mesh<MFT>& msh,
    MshCavity& cav,
    int ithrd,
    CommonScanObjective objective){
  constexpr int nnmet = 3;
  const int ipins = cav.ipins;
  const double x0 = msh.coord(ipins,0);
  const double y0 = msh.coord(ipins,1);
  double metric0[nnmet];
  for(int ii = 0; ii < nnmet; ii++) metric0[ii] = msh.met(ipins,ii);

  const int ibpoin = msh.poi2bpo[ipins];
  const bool point_on_edge = ibpoin >= 0 && msh.bpo2ibi(ibpoin,1) == 1;
  const double parameter0 = point_on_edge ? msh.bpo2rbi(ibpoin,0) : 0.;
  const double min_child_area_fraction =
      objective_compare_min_child_area_fraction();

  CommonObjectiveOptimization result;
  result.x = x0;
  result.y = y0;
  result.parameter = parameter0;

  auto restore = [&](){
    msh.coord(ipins,0) = x0;
    msh.coord(ipins,1) = y0;
    for(int ii = 0; ii < nnmet; ii++) msh.met(ipins,ii) = metric0[ii];
    if(point_on_edge) msh.bpo2rbi(ibpoin,0) = parameter0;
  };

  auto evaluate_current = [&]() -> double{
    result.probes++;
    if(min_child_area_fraction > 0.
        && common_cavity_min_child_area_fraction_2d(
             msh,cav,msh.coord(ipins,0),msh.coord(ipins,1))
             < min_child_area_fraction){
      return std::numeric_limits<double>::infinity();
    }
    const CommonCavityValues values =
        evaluate_cavity_common_vb_2d(msh,cav,ithrd);
    if(!values.valid) return std::numeric_limits<double>::infinity();
    result.original = select_common_objective(values.original,objective);
    return select_common_objective(values.reconnected,objective);
  };

  auto refresh_metric = [&]() -> bool{
    return msh.interpMetBack(ipins) <= 0;
  };

  result.initial = evaluate_current();
  result.optimized = result.initial;

  if(point_on_edge){
    const int insertion_edge = msh.bpo2ibi(ibpoin,2);
    const int insertion_ref = msh.edg2ref[insertion_edge];
    ego cad_edge = msh.CAD.cad2edg[insertion_ref];
    double range[2];
    int periodic = 0;
    if(EG_getRange(cad_edge,range,&periodic) == EGADS_SUCCESS){
      if(range[0] > range[1]) std::swap(range[0],range[1]);
      double step = 0.1*(range[1] - range[0]);
      for(int iteration = 0; iteration < 160 && step > 1.e-12; iteration++){
        bool improved = false;
        for(double direction : {-1.,1.}){
          const double parameter = result.parameter + direction*step;
          if(parameter < range[0] || parameter > range[1]) continue;
          double eg_parameter[2] = {parameter,0.};
          double evaluation[18];
          if(EG_evaluate(cad_edge,eg_parameter,evaluation) != EGADS_SUCCESS){
            continue;
          }
          msh.bpo2rbi(ibpoin,0) = parameter;
          msh.coord(ipins,0) = evaluation[0];
          msh.coord(ipins,1) = evaluation[1];
          if(!refresh_metric()) continue;
          const double value = evaluate_current();
          if(value < result.optimized){
            result.optimized = value;
            result.parameter = parameter;
            result.x = msh.coord(ipins,0);
            result.y = msh.coord(ipins,1);
            improved = true;
          }
        }
        if(!improved) step *= 0.5;
      }
    }
  }else{
    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();
    std::set<int> cavity_points;
    for(int iface : cav.lcfac){
      for(int inode = 0; inode < 3; inode++){
        cavity_points.insert(msh.fac2poi(iface,inode));
      }
    }
    for(int ipoin : cavity_points){
      xmin = std::min(xmin,msh.coord(ipoin,0));
      xmax = std::max(xmax,msh.coord(ipoin,0));
      ymin = std::min(ymin,msh.coord(ipoin,1));
      ymax = std::max(ymax,msh.coord(ipoin,1));
    }
    double step = 0.1*std::max(xmax-xmin,ymax-ymin);
    constexpr double directions[8][2] = {
      {1.,0.},{-1.,0.},{0.,1.},{0.,-1.},
      {1.,1.},{1.,-1.},{-1.,1.},{-1.,-1.}
    };
    for(int iteration = 0; iteration < 160 && step > 1.e-12; iteration++){
      bool improved = false;
      for(const auto& direction : directions){
        msh.coord(ipins,0) = result.x + step*direction[0];
        msh.coord(ipins,1) = result.y + step*direction[1];
        if(!refresh_metric()) continue;
        const double value = evaluate_current();
        if(value < result.optimized){
          result.optimized = value;
          result.x = msh.coord(ipins,0);
          result.y = msh.coord(ipins,1);
          improved = true;
        }
      }
      if(!improved) step *= 0.5;
    }
  }

  // Geometry diagnostics for detecting optimizers that obtain a low integral
  // by approaching an existing vertex or a cavity-boundary line.
  std::set<int> cavity_points;
  for(int iface : cav.lcfac){
    for(int inode = 0; inode < 3; inode++){
      cavity_points.insert(msh.fac2poi(iface,inode));
    }
  }
  double diameter = 0.;
  double min_distance = std::numeric_limits<double>::infinity();
  for(int ipoin : cavity_points){
    const double dx = result.x - msh.coord(ipoin,0);
    const double dy = result.y - msh.coord(ipoin,1);
    min_distance = std::min(min_distance,std::sqrt(dx*dx + dy*dy));
    for(int jpoin : cavity_points){
      const double ex = msh.coord(ipoin,0) - msh.coord(jpoin,0);
      const double ey = msh.coord(ipoin,1) - msh.coord(jpoin,1);
      diameter = std::max(diameter,std::sqrt(ex*ex + ey*ey));
    }
  }
  result.min_vertex_distance_ratio = diameter > 0.
      ? min_distance/diameter : 0.;
  result.min_child_area_fraction =
      common_cavity_min_child_area_fraction_2d(msh,cav,result.x,result.y);

  msh.coord(ipins,0) = result.x;
  msh.coord(ipins,1) = result.y;
  if(point_on_edge) msh.bpo2rbi(ibpoin,0) = result.parameter;
  if(refresh_metric()){
    const CommonCavityValues values =
        evaluate_cavity_common_vb_2d(msh,cav,ithrd);
    if(values.valid){
      result.min_metric_volume = values.reconnected.min_metric_volume;
      result.original_barrier_unscaled =
          values.original.step_distance_barrier;
      result.optimized_barrier_unscaled =
          values.reconnected.step_distance_barrier;
    }
  }

  restore();
  return result;
}

template<class MFT, QuaFun iquaf>
void write_cavity_quality_2d(Mesh<MFT>& msh,
                             MshCavity& cav,
                             int ithrd,
                             const std::string& label,
                             std::ofstream& log){
  CavityQualityStats stats =
    compute_cavity_quality_2d<MFT,iquaf>(msh, cav, ithrd, label, &log);
  double qavg0 = stats.nqua0 > 0 ? stats.qsum0 / stats.nqua0 : 0;
  double qavg1 = stats.nqua1 > 0 ? stats.qsum1 / stats.nqua1 : 0;

  log << label
      << " nqua " << stats.nqua0 << " -> " << stats.nqua1
      << " qsum " << stats.qsum0 << " -> " << stats.qsum1
      << " qmin " << stats.qmin0 << " -> " << stats.qmin1
      << " qmax " << stats.qmax0 << " -> " << stats.qmax1
      << " qavg " << qavg0 << " -> " << qavg1
      << "\n";
}

template<class MFT, QuaFun iquaf>
void sample_boundary_cavity_quality_2d(Mesh<MFT>& msh,
                                       MshCavity& cav,
                                       int ithrd,
                                       const std::string& label,
                                       std::ofstream& log){
  const int ipins = cav.ipins;
  const int ibpoin = msh.poi2bpo[ipins];
  if(ibpoin < 0 || msh.bpo2ibi(ibpoin,1) != 1) return;

  constexpr int nnmet = (2*(2+1))/2;
  double coord0[2];
  double met0[nnmet];
  for(int ii = 0; ii < 2; ii++) coord0[ii] = msh.coord(ipins,ii);
  for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipins,ii);

  const double t0 = msh.bpo2rbi(ibpoin,0);
  const int iedge = msh.bpo2ibi(ibpoin,2);
  const int iref = msh.edg2ref[iedge];
  ego cadEdge = msh.CAD.cad2edg[iref];

  const double dt[9] = {-0.1,-0.05,-0.025,-0.01,0.,0.01,0.025,0.05,0.1};
  log << label << " samples t0 " << t0 << "\n";
  for(double dti : dt){
    double egParam[2] = {t0 + dti, 0.};
    double evalResult[18];
    int estat = EG_evaluate(cadEdge, egParam, evalResult);
    if(estat != EGADS_SUCCESS){
      log << "  dt " << dti << " t " << egParam[0]
          << " egads_error " << estat << "\n";
      continue;
    }

    msh.bpo2rbi(ibpoin,0) = egParam[0];
    for(int ii = 0; ii < 2; ii++) msh.coord(ipins,ii) = evalResult[ii];
    int ierro = msh.interpMetBack(ipins);
    if(ierro > 0){
      log << "  dt " << dti << " t " << egParam[0]
          << " interp_error " << ierro << "\n";
      continue;
    }

    CavityQualityStats stats =
      compute_cavity_quality_2d<MFT,iquaf>(msh, cav, ithrd,
                                           label + " sample", nullptr);
    log << "  dt " << dti
        << " t " << egParam[0]
        << " coord [" << msh.coord(ipins,0) << " " << msh.coord(ipins,1) << "]"
        << " current_qsum " << stats.qsum0
        << " reconnect_qsum " << stats.qsum1
        << " reconnect_qmax " << stats.qmax1
        << "\n";
  }

  msh.bpo2rbi(ibpoin,0) = t0;
  for(int ii = 0; ii < 2; ii++) msh.coord(ipins,ii) = coord0[ii];
  for(int ii = 0; ii < nnmet; ii++) msh.met(ipins,ii) = met0[ii];
}

template<class MFT>
BadEntHandler* build_handler_for_mesh(Mesh<MFT>& msh, int tdim, dblAr1*& lquae_store){
  bool iinva = false;
  double qmin = 0, qmax = 0, qavg = 0;
  lquae_store = new dblAr1(msh.nentt(tdim));
  dblAr1& lquae = *lquae_store;

  #ifdef STEPDISTANCE
  getmetquamesh<MFT, QuaFun::StepDistance>(msh,tdim,AsDeg::P1,AsDeg::P1,
                                           &iinva,&qmin,&qmax,&qavg,&lquae);
  #else
  getmetquamesh<MFT, QuaFun::SizeShape>(msh,tdim,AsDeg::P1,AsDeg::P1,
                                        &iinva,&qmin,&qmax,&qavg,&lquae);
  #endif

  BadEntHandler* handler = new BadEntHandler(tdim, 100., 0.00001);
  const intAr2& ent2poi = msh.ent2poi(tdim);
  handler->setCallbacks([&lquae](int ientt){ return lquae[ientt]; },
                        [&](int ientt){ return isdeadent(ientt,ent2poi); });

  std::vector<int> sorted_ids(msh.nentt(tdim));
  std::iota(sorted_ids.begin(), sorted_ids.end(), 0);
  std::sort(sorted_ids.begin(), sorted_ids.end(),
            [&](int a, int b){ return lquae[a] > lquae[b]; });
  handler->seedFromSortedIDs(sorted_ids);
  return handler;
}

template<class MFT, QuaFun iquaf>
void trace_quality_growth_2d(Mesh<MFT>& msh,
                             MshCavity& cav,
                             BadEntHandler& handler,
                             int ithrd,
                             const std::string& label,
                             std::ofstream& log){
  const int tdim = 2;
  aux_taginsrefs(msh, cav, ithrd);
  msh.tag[ithrd]++;
  const int tag = msh.tag[ithrd];

  for(int iface : cav.lcfac){
    msh.fac2tag(ithrd, iface) = tag;
  }
  for(int iedge : cav.lcedg){
    msh.edg2tag(ithrd, iedge) = tag;
  }

  intAr1& lcent = cav.lcent(tdim);
  const intAr2& ent2ent = msh.ent2ent(tdim);
  intAr2& ent2poi = msh.ent2poi(tdim);
  intAr2r& ent2tag = msh.ent2tag(tdim);

  const int nface0 = msh.nface;
  msh.set_nface(nface0 + 1);
  const int tmpFace = nface0;
  const int ipins = cav.ipins;
  const bool ipinsOnEdge = cav.lcedg.get_n() > 0;
  const int iedins = ipinsOnEdge ? cav.lcedg[0] : -1;
  const int skip_face = trace_skip_face();
  double difto = 1.;

  log << label << " local growth probes\n";
  int icen0 = 0;
  int icen1 = lcent.get_n();
  for(int igrow = 0; igrow < 5; igrow++){
    bool grew = false;
    log << " grow_step " << igrow
        << " cavity_faces_before " << intAr1(lcent.get_n(), &lcent[0])
        << "\n";

    for(int icent = icen0; icent < icen1; icent++){
      int ientt = lcent[icent];
      for(int jj = 0; jj < 3; jj++){
        int ienei = ent2ent(ientt,jj);
        if(ienei < 0) continue;
        if(ent2tag(ithrd, ienei) >= tag) continue;
        if(ienei == skip_face){
          log << "  candidate_neighbor " << ienei
              << " from_face " << ientt
              << " across_local_edge " << jj
              << " skipped_by_diagnostic 1\n";
          continue;
        }

        int facetsTouching[4] = {-1,-1,-1,-1};
        for(int kk = 0; kk < 3; kk++){
          int ieneinei = ent2ent(ienei,kk);
          if(ieneinei < 0){
            facetsTouching[kk] = -2;
          }else if(ent2tag(ithrd, ieneinei) >= tag){
            facetsTouching[kk] = 1;
          }
        }

        double quaLocalReconnect = 0;
        double quaMaxLocalReconnect = -1;
        double quaLocal = 0;
        double quaMaxLocal = -1;
        bool invalid = false;

        for(int iedge = 0; iedge < 3; iedge++){
          if(facetsTouching[iedge] > 0){
            int ienttcav = ent2ent(ienei,iedge);
            int iedgeFromInside = -1;
            for(int kk = 0; kk < 3; kk++){
              if(ent2ent(ienttcav,kk) == ienei){
                iedgeFromInside = kk;
                break;
              }
            }

            int ent2pol[3];
            ent2pol[0] = ipins;
            ent2pol[lnoed2[0][0]] = ent2poi(ienttcav,lnoed2[iedgeFromInside][0]);
            ent2pol[lnoed2[0][1]] = ent2poi(ienttcav,lnoed2[iedgeFromInside][1]);
            if(ent2pol[1] == ipins || ent2pol[2] == ipins){
              invalid = true;
              break;
            }

            ent2poi(tmpFace,0) = ent2pol[0];
            ent2poi(tmpFace,1) = ent2pol[1];
            ent2poi(tmpFace,2) = ent2pol[2];

            double meas = 0;
            if(!isvalideltP1<2,2>(msh, tmpFace, NULL, &meas)){
              invalid = true;
              break;
            }
            double qua = metqua<MFT,2,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,tmpFace,difto);
            quaLocal += qua;
            quaMaxLocal = MAX(quaMaxLocal, qua);
          }else{
            if(ipinsOnEdge && facetsTouching[iedge] == -2){
              int iedgeGlobal = msh.fac2edg(ienei,iedge);
              if(iedgeGlobal >= 0 && msh.edg2ref[iedgeGlobal] == msh.edg2ref[iedins]){
                log << "  candidate_neighbor " << ienei
                    << " skip_edge " << iedge
                    << " same_boundary_ref " << msh.edg2ref[iedgeGlobal]
                    << "\n";
                continue;
              }
            }

            int ent2pol[3];
            ent2pol[0] = ipins;
            ent2pol[lnoed2[0][0]] = ent2poi(ienei,lnoed2[iedge][0]);
            ent2pol[lnoed2[0][1]] = ent2poi(ienei,lnoed2[iedge][1]);
            if(ent2pol[1] == ipins || ent2pol[2] == ipins){
              invalid = true;
              break;
            }

            ent2poi(tmpFace,0) = ent2pol[0];
            ent2poi(tmpFace,1) = ent2pol[1];
            ent2poi(tmpFace,2) = ent2pol[2];

            double meas = 0;
            if(!isvalideltP1<2,2>(msh, tmpFace, NULL, &meas)){
              invalid = true;
              break;
            }
            double qua = metqua<MFT,2,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,tmpFace,difto);
            quaLocalReconnect += qua;
            quaMaxLocalReconnect = MAX(quaMaxLocalReconnect, qua);
          }
        }

        if(invalid){
          log << "  candidate_neighbor " << ienei
              << " from_face " << ientt
              << " across_local_edge " << jj
              << " invalid_reconnect 1\n";
          continue;
        }

        double quaLocalInside = quaLocal;
        double quaOutside = metqua<MFT,2,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,ienei,difto);
        quaLocal += quaOutside;
        quaMaxLocal = MAX(quaMaxLocal, quaOutside);
        bool improveLocalSum = quaLocalReconnect <= quaLocal;
        bool improveLocalMax = true;
        #ifdef IMPROVEMAXQUAL
        improveLocalMax = quaMaxLocalReconnect <= quaMaxLocal;
        #endif

        log << "  candidate_neighbor " << ienei
            << " from_face " << ientt
            << " across_local_edge " << jj
            << " current_inside_sum " << quaLocalInside
            << " outside_quality " << quaOutside
            << " current_sum " << quaLocal
            << " reconnect_sum " << quaLocalReconnect
            << " current_max " << quaMaxLocal
            << " reconnect_max " << quaMaxLocalReconnect
            << " improve_sum " << improveLocalSum
            << " improve_max " << improveLocalMax
            << "\n";

        if(improveLocalSum && improveLocalMax){
          lcent.stack(ienei);
          ent2tag(ithrd, ienei) = tag;
          grew = true;
          log << "   -> quality path would stack " << ienei << "\n";
        }
      }
    }

    icen0 = icen1;
    icen1 = lcent.get_n();
    if(!grew || icen0 == icen1) break;
  }

  ent2poi(tmpFace,0) = -1;
  msh.set_nface(nface0);
}

CavOprOpt diagnostic_cavity_options(){
  CavOprOpt opts;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = false;
  opts.dryrun = false;
  opts.allow_remove_points = true;
  opts.allow_remove_points_superdim = true;
  opts.qmax_nec = -1;
  opts.qmax_suf = -1;
  opts.qmax_iff = -1;
  return opts;
}

template<class MFT, int gdim, int ideg>
void initialize_insertion_cavity(const std::string& cmd,
                                 int tdim,
                                 int ientt,
                                 int ied,
                                 MetrisRunner*& run,
                                 Mesh<MFT>*& msh,
                                 MshCavity*& cav,
                                 EdgeSeed*& insertion_seed,
                                 BadEntHandler*& handler,
                                 dblAr1*& lquae_store,
                                 std::ofstream& log,
                                 const std::filesystem::path& dir,
                                 const std::string& prefix){
  cargHandler arg(cmd);
  run = new MetrisRunner(arg.c, arg.v);
  msh = &static_cast<Mesh<MFT>&>(*run->msh_g);
  setup_quality_mesh<MFT,gdim,ideg>(*msh);

  handler = build_handler_for_mesh(*msh, tdim, lquae_store);
  cav = new MshCavity(100,100,1);
  insertion_seed = new EdgeSeed(*msh, *cav, tdim, tdim, ientt, ied);

  cav->ipins = msh->newpoint(PointType::Vertex,
                             insertion_seed->tdimp,
                             insertion_seed->iseed);
  log << prefix << " newpoint ipins " << cav->ipins
      << " tdimp " << insertion_seed->tdimp
      << " iseed " << insertion_seed->iseed
      << " iref " << insertion_seed->iref
      << " endpoints " << insertion_seed->ipedg[0] << " "
      << insertion_seed->ipedg[1] << "\n";
  write_cavity_trace(dir, prefix + "_00_seed_shell", *msh, *cav, log);

  int ierro = aux_bisecPointLen(*msh, *insertion_seed,
                                msh->poi2bpo[cav->ipins], false, *cav);
  log << prefix << " aux_bisecPointLen ierro " << ierro << "\n";
  write_cavity_trace(dir, prefix + "_01_after_bisec_point", *msh, *cav, log);
  BOOST_REQUIRE_EQUAL(ierro, 0);

  ierro = increase_cavity(*msh, *cav, false, 1, 2);
  log << prefix << " initial increase_cavity ierro " << ierro << "\n";
  write_cavity_trace(dir, prefix + "_02_after_initial_increase", *msh, *cav, log);
  BOOST_REQUIRE(ierro <= 0);
}

template<class MFT>
void cleanup_trace_objects(MetrisRunner* run,
                           MshCavity* cav,
                           EdgeSeed* insertion_seed,
                           BadEntHandler* handler,
                           dblAr1* lquae_store){
  delete insertion_seed;
  delete cav;
  delete handler;
  delete lquae_store;
  delete run;
}

template<class MFT, int gdim, int ideg>
void trace_exact_case_impl(const std::string& cmd,
                           int ientt_in = -1,
                           int ied_in = -1){
  const int tdim = 2;
  const int ientt = ientt_in >= 0 ? ientt_in : trace_ientt();
  const int ied = ied_in >= 0 ? ied_in : trace_ied();
  const int ithrd1 = 1;
  const int ithrd2 = 2;

  std::filesystem::create_directories(trace_dir());
  std::filesystem::create_directories(vizir_dir());
  std::ofstream log(trace_dir() / "trace.log");
  BOOST_REQUIRE(log.good());

  {
    cargHandler arg(cmd);
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
    setup_quality_mesh<MFT,gdim,ideg>(msh);
    double sz[2];
    double elen = getlenedg_geosz<MFT,gdim,ideg>(msh, ientt, tdim, ied, sz);
    log << "case adaptive_iteration " << adaptive_iteration()
        << " tdim " << tdim
        << " ientt " << ientt
        << " ied " << ied
        << " edge_length " << elen
        << " sz0 " << sz[0]
        << " sz1 " << sz[1]
        << "\n";
    log << "face_nodes " << intAr1(msh.nnode(tdim), msh.fac2poi[ientt]) << "\n";
  }

  lenStat lenstat0;
  {
    cargHandler arg(cmd);
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
    setup_quality_mesh<MFT,gdim,ideg>(msh);
    intAr2 ilned;
    dblAr1 rlned;
    getLengthEdges<MFT>(msh,tdim,-1,ilned,rlned,lenstat0);
  }
  const double lenqua_short_max = (lenstat0.qua_short + lenstat0.qua_long) / 2;
  log << "lenqua_short_max " << lenqua_short_max
      << " qua_short " << lenstat0.qua_short
      << " qua_long " << lenstat0.qua_long << "\n";

  {
    cargHandler arg(cmd);
    MetrisRunner run_insert(arg.c, arg.v);
    Mesh<MFT>& msh_insert = static_cast<Mesh<MFT>&>(*run_insert.msh_g);
    setup_quality_mesh<MFT,gdim,ideg>(msh_insert);

    bool iinva = false;
    double qmin = 0, qmax = 0, qavg = 0;
    dblAr1 lquae_insert(msh_insert.nentt(tdim));
    #ifdef STEPDISTANCE
    getmetquamesh<MFT, QuaFun::StepDistance>(msh_insert,tdim,AsDeg::P1,AsDeg::P1,
                                             &iinva,&qmin,&qmax,&qavg,&lquae_insert);
    #else
    getmetquamesh<MFT, QuaFun::SizeShape>(msh_insert,tdim,AsDeg::P1,AsDeg::P1,
                                          &iinva,&qmin,&qmax,&qavg,&lquae_insert);
    #endif

    BadEntHandler handler_insert(tdim, 100., 0.00001);
    const intAr2& ent2poi_insert = msh_insert.ent2poi(tdim);
    handler_insert.setCallbacks([&](int ientt_){ return lquae_insert[ientt_]; },
                                [&](int ientt_){ return isdeadent(ientt_,ent2poi_insert); });
    std::vector<int> sorted_ids(msh_insert.nentt(tdim));
    std::iota(sorted_ids.begin(), sorted_ids.end(), 0);
    std::sort(sorted_ids.begin(), sorted_ids.end(),
              [&](int a, int b){ return lquae_insert[a] > lquae_insert[b]; });
    handler_insert.seedFromSortedIDs(sorted_ids);

    MshCavity cav_insert(100,100,1);
    EdgeSeed seed_insert(msh_insert, cav_insert, tdim, tdim, ientt, ied);

    CavWrkArrs work_insert;
    intAr1 lcaverr(CAV_ERR_NERROR);
    lcaverr.set_n(CAV_ERR_NERROR);
    for(int ii = 0; ii < CAV_ERR_NERROR; ii++) lcaverr[ii] = 0;

    int ierro_insert = insertEdge(msh_insert, seed_insert, lenqua_short_max,
                                  false, cav_insert, work_insert, lcaverr,
                                  handler_insert, false, 0., ithrd1, ithrd2);

    log << "quality_insert_actual insertEdge ierro " << ierro_insert
        << " ipins " << cav_insert.ipins
        << " nedge " << cav_insert.lcedg.get_n()
        << " nface " << cav_insert.lcfac.get_n()
        << " ntet " << cav_insert.lctet.get_n()
        << "\n";
    log << "quality_insert_actual lcedg " << cav_insert.lcedg << "\n";
    log << "quality_insert_actual lcfac " << cav_insert.lcfac << "\n";
    log << "quality_insert_actual lcaverr";
    for(int ii = 0; ii < lcaverr.get_n(); ii++){
      if(lcaverr[ii] > 0){
        log << " [" << ii + 1 << " " << cavity_error_name(ii + 1)
            << " " << lcaverr[ii] << "]";
      }
    }
    log << "\n";
  }

  {
    MetrisRunner* run_qop = nullptr;
    Mesh<MFT>* msh_qop = nullptr;
    MshCavity* cav_qop = nullptr;
    EdgeSeed* seed_qop = nullptr;
    BadEntHandler* handler_qop = nullptr;
    dblAr1* lquae_qop = nullptr;
    initialize_insertion_cavity<MFT,gdim,ideg>(cmd, tdim, ientt, ied,
                                               run_qop, msh_qop,
                                               cav_qop, seed_qop,
                                               handler_qop, lquae_qop,
                                               log, trace_dir(),
                                               "quality_operator_input");

    CavOprOpt opts_qop = diagnostic_cavity_options();
    std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp_qop;
    int ierro_qop = setCavityInsertionQuality(*msh_qop, *cav_qop,
                                              opts_qop, *seed_qop, 5,
                                              *handler_qop, lenqua_short_max,
                                              nocomp_qop, ithrd1, ithrd2);

    log << "quality_operator_input setCavityInsertionQuality ierro "
        << ierro_qop
        << " ipins " << cav_qop->ipins
        << " nedge " << cav_qop->lcedg.get_n()
        << " nface " << cav_qop->lcfac.get_n()
        << " ntet " << cav_qop->lctet.get_n()
        << "\n";
    log << "quality_operator_input lcedg " << cav_qop->lcedg << "\n";
    log << "quality_operator_input lcfac " << cav_qop->lcfac << "\n";
    log_point_position("quality_operator_input position",
                       *msh_qop, cav_qop->ipins, log);
    write_cavity_trace(trace_dir(), "quality_operator_input_final",
                       *msh_qop, *cav_qop, log);
    write_painted_point_meshes(vizir_dir(), "quality_operator_input_final",
                               *msh_qop, *cav_qop, log);
    write_painted_cavity_context(vizir_dir(), "quality_operator_input_final",
                                 *msh_qop, *cav_qop, log);
    trace_linfac_candidates(*msh_qop, *cav_qop,
                            "quality_operator_input_final",
                            log, ithrd1);

    cleanup_trace_objects<MFT>(run_qop, cav_qop, seed_qop,
                               handler_qop, lquae_qop);
  }

  MetrisRunner* run_quality = nullptr;
  Mesh<MFT>* msh_quality = nullptr;
  MshCavity* cav_quality = nullptr;
  EdgeSeed* seed_quality = nullptr;
  BadEntHandler* handler_quality = nullptr;
  dblAr1* lquae_quality = nullptr;
  initialize_insertion_cavity<MFT,gdim,ideg>(cmd, tdim, ientt, ied,
                                             run_quality, msh_quality,
                                             cav_quality, seed_quality,
                                             handler_quality, lquae_quality,
                                             log, trace_dir(),
                                             "quality");

  write_cavity_objective_decomposition_2d(
      *msh_quality,*cav_quality,ithrd1,
      "common_midpoint integration_comparison",log);

  #ifdef STEPDISTANCE
  CavityQualityStats quality_initial_unsmoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::StepDistance>(*msh_quality, *cav_quality,
                                                        ithrd1,
                                                        "quality initial smoothing before",
                                                        &log);
  double quality_initial_smooth_qsum = quality_initial_unsmoothed_stats.qsum1;
  double quality_initial_smooth_qmax = quality_initial_unsmoothed_stats.qmax1;
  constexpr QuaFun quality_smooth_fun = QuaFun::StepDistance;
  #else
  CavityQualityStats quality_initial_unsmoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality, *cav_quality,
                                                     ithrd1,
                                                     "quality initial smoothing before",
                                                     &log);
  double quality_initial_smooth_qsum = quality_initial_unsmoothed_stats.qsum1;
  double quality_initial_smooth_qmax = quality_initial_unsmoothed_stats.qmax1;
  constexpr QuaFun quality_smooth_fun = QuaFun::SizeShape;
  #endif

  log << "quality initial smoothing diagnostic before"
      << " current_qsum " << quality_initial_unsmoothed_stats.qsum0
      << " reconnect_qsum " << quality_initial_unsmoothed_stats.qsum1
      << " reconnect_qmax " << quality_initial_unsmoothed_stats.qmax1
      << "\n";
  log_point_position("quality initial smoothing before position",
                     *msh_quality, cav_quality->ipins, log);
  #ifdef STEPDISTANCE
  sample_boundary_cavity_quality_2d<MFT,QuaFun::StepDistance>(*msh_quality,
                                                              *cav_quality,
                                                              ithrd1,
                                                              "quality initial smoothing",
                                                              log);
  #else
  sample_boundary_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality,
                                                           *cav_quality,
                                                           ithrd1,
                                                           "quality initial smoothing",
                                                           log);
  #endif

  double initial_smooth_noper = 0;
  double quality_initial_smooth_target_weight = 0.;
  const int nface_before_initial_smoothing = msh_quality->nface;
  if(!msh_quality->param->step_distance_cavity_target_average){
    msh_quality->set_nface(nface_before_initial_smoothing + 1);
    initial_smooth_noper =
      smoothCavity(*msh_quality, *cav_quality, *handler_quality,
                   quality_smooth_fun,
                   quality_initial_unsmoothed_stats.qsum1,
                   quality_initial_unsmoothed_stats.qmax1,
                   0.,
                   quality_initial_smooth_qsum,
                   quality_initial_smooth_qmax,
                   quality_initial_smooth_target_weight,
                   ithrd1, ithrd2);
    msh_quality->fac2poi(nface_before_initial_smoothing,0) = -1;
    msh_quality->set_nface(nface_before_initial_smoothing);
  }else{
    log << "quality initial smoothing diagnostic skipped: raw-quality "
           "decomposition has no cavity target-weight accumulator\n";
  }

  #ifdef STEPDISTANCE
  CavityQualityStats quality_initial_smoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::StepDistance>(*msh_quality, *cav_quality,
                                                        ithrd1,
                                                        "quality initial smoothing after",
                                                        &log);
  #else
  CavityQualityStats quality_initial_smoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality, *cav_quality,
                                                     ithrd1,
                                                     "quality initial smoothing after",
                                                     &log);
  #endif

  log << "quality initial smoothing diagnostic after"
      << " noper " << initial_smooth_noper
      << " smooth_return_qsum " << quality_initial_smooth_qsum
      << " recomputed_current_qsum " << quality_initial_smoothed_stats.qsum0
      << " recomputed_reconnect_qsum " << quality_initial_smoothed_stats.qsum1
      << " smooth_return_qmax " << quality_initial_smooth_qmax
      << "\n";
  log_point_position("quality initial smoothing after position",
                     *msh_quality, cav_quality->ipins, log);
  write_cavity_trace(trace_dir(), "quality_03_after_initial_cavity_smoothing",
                     *msh_quality, *cav_quality, log);

  #ifdef STEPDISTANCE
  write_cavity_quality_2d<MFT,QuaFun::StepDistance>(*msh_quality, *cav_quality,
                                                    ithrd1,
                                                    "quality before growth",
                                                    log);
  #else
  write_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality, *cav_quality,
                                                 ithrd1,
                                                 "quality before growth",
                                                 log);
  #endif
  #ifdef STEPDISTANCE
  trace_quality_growth_2d<MFT,QuaFun::StepDistance>(*msh_quality, *cav_quality,
                                                    *handler_quality,
                                                    ithrd1, "quality", log);
  #else
  trace_quality_growth_2d<MFT,QuaFun::SizeShape>(*msh_quality, *cav_quality,
                                                 *handler_quality,
                                                 ithrd1, "quality", log);
  #endif
  write_cavity_trace(trace_dir(), "quality_03_after_quality_growth_trace",
                     *msh_quality, *cav_quality, log);
  #ifdef STEPDISTANCE
  write_cavity_quality_2d<MFT,QuaFun::StepDistance>(*msh_quality, *cav_quality,
                                                    ithrd1,
                                                    "quality after growth trace",
                                                    log);
  #else
  write_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality, *cav_quality,
                                                 ithrd1,
                                                 "quality after growth trace",
                                                 log);
  #endif
  write_painted_point_meshes(vizir_dir(), "quality_final",
                             *msh_quality, *cav_quality, log);

  #ifdef STEPDISTANCE
  CavityQualityStats quality_unsmoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::StepDistance>(*msh_quality, *cav_quality,
                                                        ithrd1,
                                                        "quality smoothing before",
                                                        &log);
  double quality_smooth_qsum = quality_unsmoothed_stats.qsum1;
  double quality_smooth_qmax = quality_unsmoothed_stats.qmax1;
  #else
  CavityQualityStats quality_unsmoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality, *cav_quality,
                                                     ithrd1,
                                                     "quality smoothing before",
                                                     &log);
  double quality_smooth_qsum = quality_unsmoothed_stats.qsum1;
  double quality_smooth_qmax = quality_unsmoothed_stats.qmax1;
  #endif

  log << "quality smoothing diagnostic before"
      << " current_qsum " << quality_unsmoothed_stats.qsum0
      << " reconnect_qsum " << quality_unsmoothed_stats.qsum1
      << " reconnect_qmax " << quality_unsmoothed_stats.qmax1
      << "\n";
  log_point_position("quality smoothing before position",
                     *msh_quality, cav_quality->ipins, log);

  double smooth_noper = 0;
  double quality_smooth_target_weight = 0.;
  const int nface_before_smoothing = msh_quality->nface;
  if(!msh_quality->param->step_distance_cavity_target_average){
    msh_quality->set_nface(nface_before_smoothing + 1);
    smooth_noper =
      smoothCavity(*msh_quality, *cav_quality, *handler_quality,
                   quality_smooth_fun,
                   quality_unsmoothed_stats.qsum1,
                   quality_unsmoothed_stats.qmax1,
                   0.,
                   quality_smooth_qsum,
                   quality_smooth_qmax,
                   quality_smooth_target_weight,
                   ithrd1, ithrd2);
    msh_quality->fac2poi(nface_before_smoothing,0) = -1;
    msh_quality->set_nface(nface_before_smoothing);
  }else{
    log << "quality smoothing diagnostic skipped: raw-quality decomposition "
           "has no cavity target-weight accumulator\n";
  }

  #ifdef STEPDISTANCE
  CavityQualityStats quality_smoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::StepDistance>(*msh_quality, *cav_quality,
                                                        ithrd1,
                                                        "quality smoothing after",
                                                        &log);
  #else
  CavityQualityStats quality_smoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality, *cav_quality,
                                                     ithrd1,
                                                     "quality smoothing after",
                                                     &log);
  #endif

  const bool quality_smoothing_passes_original =
    handler_quality->checkSuccess(quality_smooth_qsum,
                                  quality_unsmoothed_stats.qsum0);
  log << "quality smoothing diagnostic after"
      << " noper " << smooth_noper
      << " smooth_return_qsum " << quality_smooth_qsum
      << " recomputed_current_qsum " << quality_smoothed_stats.qsum0
      << " recomputed_reconnect_qsum " << quality_smoothed_stats.qsum1
      << " smooth_return_qmax " << quality_smooth_qmax
      << " passes_original_qsum " << quality_smoothing_passes_original
      << "\n";
  log_point_position("quality smoothing after position",
                     *msh_quality, cav_quality->ipins, log);
  write_cavity_objective_decomposition_2d(
      *msh_quality,*cav_quality,ithrd1,
      "quality final_decomposition",log);
  write_cavity_trace(trace_dir(), "quality_04_after_cavity_smoothing",
                     *msh_quality, *cav_quality, log);
  write_painted_point_meshes(vizir_dir(), "quality_smoothed_final",
                             *msh_quality, *cav_quality, log);

  cleanup_trace_objects<MFT>(run_quality, cav_quality, seed_quality,
                             handler_quality, lquae_quality);

  MetrisRunner* run_quality_sizeshape = nullptr;
  Mesh<MFT>* msh_quality_sizeshape = nullptr;
  MshCavity* cav_quality_sizeshape = nullptr;
  EdgeSeed* seed_quality_sizeshape = nullptr;
  BadEntHandler* handler_quality_sizeshape = nullptr;
  dblAr1* lquae_quality_sizeshape = nullptr;
  initialize_insertion_cavity<MFT,gdim,ideg>(cmd, tdim, ientt, ied,
                                             run_quality_sizeshape,
                                             msh_quality_sizeshape,
                                             cav_quality_sizeshape,
                                             seed_quality_sizeshape,
                                             handler_quality_sizeshape,
                                             lquae_quality_sizeshape,
                                             log, trace_dir(),
                                             "quality_sizeshape");

  CavityQualityStats quality_sizeshape_initial_unsmoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality_sizeshape,
                                                     *cav_quality_sizeshape,
                                                     ithrd1,
                                                     "quality_sizeshape initial smoothing before",
                                                     &log);
  double quality_sizeshape_initial_smooth_qsum =
    quality_sizeshape_initial_unsmoothed_stats.qsum1;
  double quality_sizeshape_initial_smooth_qmax =
    quality_sizeshape_initial_unsmoothed_stats.qmax1;

  log << "quality_sizeshape initial smoothing diagnostic before"
      << " current_qsum " << quality_sizeshape_initial_unsmoothed_stats.qsum0
      << " reconnect_qsum " << quality_sizeshape_initial_unsmoothed_stats.qsum1
      << " reconnect_qmax " << quality_sizeshape_initial_unsmoothed_stats.qmax1
      << "\n";
  log_point_position("quality_sizeshape initial smoothing before position",
                     *msh_quality_sizeshape, cav_quality_sizeshape->ipins, log);
  sample_boundary_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality_sizeshape,
                                                           *cav_quality_sizeshape,
                                                           ithrd1,
                                                           "quality_sizeshape initial smoothing",
                                                           log);

  double sizeshape_initial_smooth_noper = 0;
  double sizeshape_initial_smooth_target_weight = 0.;
  const int nface_before_sizeshape_initial_smoothing = msh_quality_sizeshape->nface;
  msh_quality_sizeshape->set_nface(nface_before_sizeshape_initial_smoothing + 1);
  sizeshape_initial_smooth_noper =
    smoothCavity(*msh_quality_sizeshape, *cav_quality_sizeshape,
                 *handler_quality_sizeshape,
                 QuaFun::SizeShape,
                 quality_sizeshape_initial_unsmoothed_stats.qsum1,
                 quality_sizeshape_initial_unsmoothed_stats.qmax1,
                 0.,
                 quality_sizeshape_initial_smooth_qsum,
                 quality_sizeshape_initial_smooth_qmax,
                 sizeshape_initial_smooth_target_weight,
                 ithrd1, ithrd2);
  msh_quality_sizeshape->fac2poi(nface_before_sizeshape_initial_smoothing,0) = -1;
  msh_quality_sizeshape->set_nface(nface_before_sizeshape_initial_smoothing);

  CavityQualityStats quality_sizeshape_initial_smoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality_sizeshape,
                                                     *cav_quality_sizeshape,
                                                     ithrd1,
                                                     "quality_sizeshape initial smoothing after",
                                                     &log);

  log << "quality_sizeshape initial smoothing diagnostic after"
      << " noper " << sizeshape_initial_smooth_noper
      << " smooth_return_qsum " << quality_sizeshape_initial_smooth_qsum
      << " recomputed_current_qsum " << quality_sizeshape_initial_smoothed_stats.qsum0
      << " recomputed_reconnect_qsum " << quality_sizeshape_initial_smoothed_stats.qsum1
      << " smooth_return_qmax " << quality_sizeshape_initial_smooth_qmax
      << "\n";
  log_point_position("quality_sizeshape initial smoothing after position",
                     *msh_quality_sizeshape, cav_quality_sizeshape->ipins, log);
  write_cavity_trace(trace_dir(), "quality_sizeshape_03_after_initial_cavity_smoothing",
                     *msh_quality_sizeshape, *cav_quality_sizeshape, log);

  write_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality_sizeshape,
                                                 *cav_quality_sizeshape,
                                                 ithrd1,
                                                 "quality_sizeshape before growth",
                                                 log);
  trace_quality_growth_2d<MFT,QuaFun::SizeShape>(*msh_quality_sizeshape,
                                                 *cav_quality_sizeshape,
                                                 *handler_quality_sizeshape,
                                                 ithrd1,
                                                 "quality_sizeshape",
                                                 log);
  write_cavity_trace(trace_dir(), "quality_sizeshape_03_after_quality_growth_trace",
                     *msh_quality_sizeshape, *cav_quality_sizeshape, log);
  write_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality_sizeshape,
                                                 *cav_quality_sizeshape,
                                                 ithrd1,
                                                 "quality_sizeshape after growth trace",
                                                 log);
  write_painted_point_meshes(vizir_dir(), "quality_sizeshape_final",
                             *msh_quality_sizeshape, *cav_quality_sizeshape, log);
  write_painted_cavity_context(vizir_dir(), "quality_sizeshape_final",
                               *msh_quality_sizeshape, *cav_quality_sizeshape, log);

  CavityQualityStats quality_sizeshape_unsmoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality_sizeshape,
                                                     *cav_quality_sizeshape,
                                                     ithrd1,
                                                     "quality_sizeshape smoothing before",
                                                     &log);
  double quality_sizeshape_smooth_qsum = quality_sizeshape_unsmoothed_stats.qsum1;
  double quality_sizeshape_smooth_qmax = quality_sizeshape_unsmoothed_stats.qmax1;

  log << "quality_sizeshape smoothing diagnostic before"
      << " current_qsum " << quality_sizeshape_unsmoothed_stats.qsum0
      << " reconnect_qsum " << quality_sizeshape_unsmoothed_stats.qsum1
      << " reconnect_qmax " << quality_sizeshape_unsmoothed_stats.qmax1
      << "\n";
  log_point_position("quality_sizeshape smoothing before position",
                     *msh_quality_sizeshape, cav_quality_sizeshape->ipins, log);

  double sizeshape_smooth_noper = 0;
  double sizeshape_smooth_target_weight = 0.;
  const int nface_before_sizeshape_smoothing = msh_quality_sizeshape->nface;
  msh_quality_sizeshape->set_nface(nface_before_sizeshape_smoothing + 1);
  sizeshape_smooth_noper =
    smoothCavity(*msh_quality_sizeshape, *cav_quality_sizeshape,
                 *handler_quality_sizeshape,
                 QuaFun::SizeShape,
                 quality_sizeshape_unsmoothed_stats.qsum1,
                 quality_sizeshape_unsmoothed_stats.qmax1,
                 0.,
                 quality_sizeshape_smooth_qsum,
                 quality_sizeshape_smooth_qmax,
                 sizeshape_smooth_target_weight,
                 ithrd1, ithrd2);
  msh_quality_sizeshape->fac2poi(nface_before_sizeshape_smoothing,0) = -1;
  msh_quality_sizeshape->set_nface(nface_before_sizeshape_smoothing);

  CavityQualityStats quality_sizeshape_smoothed_stats =
    compute_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_quality_sizeshape,
                                                     *cav_quality_sizeshape,
                                                     ithrd1,
                                                     "quality_sizeshape smoothing after",
                                                     &log);

  const bool quality_sizeshape_smoothing_passes_original =
    handler_quality_sizeshape->checkSuccess(quality_sizeshape_smooth_qsum,
                                            quality_sizeshape_unsmoothed_stats.qsum0);
  log << "quality_sizeshape smoothing diagnostic after"
      << " noper " << sizeshape_smooth_noper
      << " smooth_return_qsum " << quality_sizeshape_smooth_qsum
      << " recomputed_current_qsum " << quality_sizeshape_smoothed_stats.qsum0
      << " recomputed_reconnect_qsum " << quality_sizeshape_smoothed_stats.qsum1
      << " smooth_return_qmax " << quality_sizeshape_smooth_qmax
      << " passes_original_qsum " << quality_sizeshape_smoothing_passes_original
      << "\n";
  log_point_position("quality_sizeshape smoothing after position",
                     *msh_quality_sizeshape, cav_quality_sizeshape->ipins, log);
  write_cavity_objective_decomposition_2d(
      *msh_quality_sizeshape,*cav_quality_sizeshape,ithrd1,
      "quality_sizeshape final_decomposition",log);
  write_cavity_trace(trace_dir(), "quality_sizeshape_04_after_cavity_smoothing",
                     *msh_quality_sizeshape, *cav_quality_sizeshape, log);
  write_painted_point_meshes(vizir_dir(), "quality_sizeshape_smoothed_final",
                             *msh_quality_sizeshape, *cav_quality_sizeshape, log);
  write_painted_cavity_context(vizir_dir(), "quality_sizeshape_smoothed_final",
                               *msh_quality_sizeshape, *cav_quality_sizeshape, log);
  trace_linfac_candidates(*msh_quality_sizeshape, *cav_quality_sizeshape,
                          "quality_sizeshape_smoothed_final",
                          log, ithrd1);
  trace_raw_cavity_operator<MFT,ideg>(*msh_quality_sizeshape,
                                      *cav_quality_sizeshape,
                                      "quality_sizeshape_smoothed_final",
                                      log, ithrd1);

  cleanup_trace_objects<MFT>(run_quality_sizeshape, cav_quality_sizeshape,
                             seed_quality_sizeshape, handler_quality_sizeshape,
                             lquae_quality_sizeshape);

  MetrisRunner* run_length = nullptr;
  Mesh<MFT>* msh_length = nullptr;
  MshCavity* cav_length = nullptr;
  EdgeSeed* seed_length = nullptr;
  BadEntHandler* handler_length = nullptr;
  dblAr1* lquae_length = nullptr;
  initialize_insertion_cavity<MFT,gdim,ideg>(cmd, tdim, ientt, ied,
                                             run_length, msh_length,
                                             cav_length, seed_length,
                                             handler_length, lquae_length,
                                             log, trace_dir(),
                                             "length");

  CavOprOpt opts = diagnostic_cavity_options();
  std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp;

  int ierro = movePointCavLen<MFT>(*msh_length, *cav_length, 5, ithrd1);
  log << "length movePointCavLen ierro " << ierro << "\n";
  write_cavity_trace(trace_dir(), "length_03_after_movePointCavLen",
                     *msh_length, *cav_length, log);
  BOOST_REQUIRE(ierro <= 0);

  ierro = increase_cavity_Delaunay(*msh_length, *cav_length, tdim, 5, ithrd1);
  log << "length increase_cavity_Delaunay ierro " << ierro << "\n";
  write_cavity_trace(trace_dir(), "length_04_after_delaunay",
                     *msh_length, *cav_length, log);
  BOOST_REQUIRE(ierro <= 0);

  ierro = increase_cavity_validity(*msh_length, *cav_length, ithrd1);
  log << "length increase_cavity_validity ierro " << ierro << "\n";
  write_cavity_trace(trace_dir(), "length_05_after_validity",
                     *msh_length, *cav_length, log);
  BOOST_REQUIRE(ierro <= 0);

  intWrkAr1 lrempoi = msh_length->get_iwork(10);
  check_cavity_rempoint(*msh_length, *cav_length, opts,
                        lrempoi.get_array(), true, ithrd1);
  log << "length check_cavity_rempoint nrempoi " << lrempoi.get_n()
      << " lrempoi " << lrempoi.get_array() << "\n";

  ierro = collrejcav_lenqua(*msh_length, *cav_length,
                            true, false, true, lenqua_short_max,
                            nocomp, ithrd2);
  log << "length collrejcav_lenqua ierro " << ierro << "\n";
  write_cavity_trace(trace_dir(), "length_06_after_lenqua_check",
                     *msh_length, *cav_length, log);
  BOOST_REQUIRE_EQUAL(ierro, 0);

  #ifdef STEPDISTANCE
  write_cavity_quality_2d<MFT,QuaFun::StepDistance>(*msh_length, *cav_length,
                                                    ithrd1,
                                                    "length final cavity quality before checkCavityQuality",
                                                    log);
  #else
  write_cavity_quality_2d<MFT,QuaFun::SizeShape>(*msh_length, *cav_length,
                                                 ithrd1,
                                                 "length final cavity quality before checkCavityQuality",
                                                 log);
  #endif

  ierro = checkCavityQuality(*msh_length, *cav_length, tdim, 5,
                             *handler_length, 0., ithrd1);
  log << "length checkCavityQuality worsenPctg=0 ierro " << ierro << "\n";
  write_cavity_trace(trace_dir(), "length_07_after_quality_check",
                     *msh_length, *cav_length, log);
  write_painted_point_meshes(vizir_dir(), "length_final",
                             *msh_length, *cav_length, log);
  write_painted_cavity_context(vizir_dir(), "length_final",
                               *msh_length, *cav_length, log);
  trace_linfac_candidates(*msh_length, *cav_length,
                          "length_final",
                          log, ithrd1);
  trace_raw_cavity_operator<MFT,ideg>(*msh_length, *cav_length,
                                      "length_final",
                                      log, ithrd1);

  log << "trace_dir " << trace_dir() << "\n";
  log << "vizir_dir " << vizir_dir() << "\n";
  log.flush();

  // The trace has already written every artifact by this point. Keeping the
  // length-side runner alive avoids teardown noise from the diagnostic path.
  (void)run_length;
  (void)msh_length;
  (void)cav_length;
  (void)seed_length;
  (void)handler_length;
  (void)lquae_length;
}

template<class MFT>
void dispatch_trace_exact_case(const std::string& cmd){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  MeshBase& msh = *run.msh_g;

  bool dispatched = false;
  CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      trace_exact_case_impl<MFT,gdim,ideg>(cmd);
      dispatched = true;
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
  BOOST_REQUIRE(dispatched);
}

template<class MFT>
void dispatch_trace_exact_case(const std::string& cmd, int ientt, int ied){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  MeshBase& msh = *run.msh_g;

  bool dispatched = false;
  CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      trace_exact_case_impl<MFT,gdim,ideg>(cmd, ientt, ied);
      dispatched = true;
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
  BOOST_REQUIRE(dispatched);
}

template<class MFT, int gdim, int ideg>
CavityIntegrationValues compare_exact_case_integration_schemes_impl(
    const std::string& cmd){
  constexpr int tdim = 2;
  constexpr int ientt = 214;
  constexpr int ied = 0;
  constexpr int ithrd1 = 1;

  const std::filesystem::path dir =
      std::filesystem::current_path() / "build" / "codex_diagnostic"
                                      / "integration_method_comparison";
  std::filesystem::create_directories(dir);
  std::ofstream log(dir / "a6_face214_edge0.log");
  BOOST_REQUIRE(log.good());

  MetrisRunner* run = nullptr;
  Mesh<MFT>* msh = nullptr;
  MshCavity* cav = nullptr;
  EdgeSeed* seed = nullptr;
  BadEntHandler* handler = nullptr;
  dblAr1* lquae = nullptr;
  initialize_insertion_cavity<MFT,gdim,ideg>(
      cmd,tdim,ientt,ied,run,msh,cav,seed,handler,lquae,
      log,dir,"common_midpoint");

  const CavityIntegrationValues values =
      write_cavity_objective_decomposition_2d(
          *msh,*cav,ithrd1,"common_midpoint",log);

  auto print_comparison = [](const char* objective,
                             const char* integration,
                             double old_value,
                             double new_value){
    const double pct = 100.*(new_value/old_value - 1.);
    fmt::print("   {:12s} {:19s}: {:.8e} -> {:.8e} ({:+.3f}%) {}\n",
               objective,integration,old_value,new_value,pct,
               new_value < old_value ? "ACCEPT" : "REJECT");
  };

  fmt::print("\n-- Exact a6 integration-scheme comparison\n");
  fmt::print("   face 214 local edge 0 endpoints 84-148; common midpoint\n");
  print_comparison("SizeShape","one point",
                   values.original.size_shape_one_point,
                   values.reconnected.size_shape_one_point);
  print_comparison("StepDistance","one point",
                   values.original.step_distance_one_point,
                   values.reconnected.step_distance_one_point);
  print_comparison("SizeShape","vertices+barycenter",
                   values.original.size_shape_vertex_barycenter,
                   values.reconnected.size_shape_vertex_barycenter);
  print_comparison("StepDistance","vertices+barycenter",
                   values.original.step_distance_vertex_barycenter,
                   values.reconnected.step_distance_vertex_barycenter);
  fmt::print("\n   StepDistance p scan on the same geometry:\n");
  fmt::print("   {:>6s} {:>13s} {:>10s} {:>13s} {:>10s}\n",
             "p","one-point %","decision","vert+bary %","decision");
  for(std::size_t ip = 0; ip < step_distance_p_values.size(); ip++){
    const double pct_one = 100.*(
        values.reconnected.step_p_one_point[ip]
        / values.original.step_p_one_point[ip] - 1.);
    const double pct_vb = 100.*(
        values.reconnected.step_p_vertex_barycenter[ip]
        / values.original.step_p_vertex_barycenter[ip] - 1.);
    fmt::print("   {:6.2f} {:+12.3f}% {:>10s} {:+12.3f}% {:>10s}\n",
               step_distance_p_values[ip],pct_one,
               pct_one < 0. ? "ACCEPT" : "REJECT",
               pct_vb,pct_vb < 0. ? "ACCEPT" : "REJECT");
  }

  fmt::print("\n   StepDistance after p-specific derivative-free point optimization:\n");
  fmt::print("   {:>6s} {:>13s} {:>10s} {:>13s} {:>10s}\n",
             "p","one-point %","decision","vert+bary %","decision");
  for(std::size_t ip = 0; ip < step_distance_p_values.size(); ip++){
    const StepPOptimizationResult optimized_one =
        optimize_step_p_point_2d(*msh,*cav,ithrd1,ip,false);
    const StepPOptimizationResult optimized_vb =
        optimize_step_p_point_2d(*msh,*cav,ithrd1,ip,true);
    const double pct_one =
        100.*(optimized_one.optimized/optimized_one.original - 1.);
    const double pct_vb =
        100.*(optimized_vb.optimized/optimized_vb.original - 1.);
    if(step_distance_p_values[ip] == 1.){
      BOOST_TEST(optimized_one.optimized < optimized_one.original);
      BOOST_TEST(optimized_vb.optimized < optimized_vb.original);
    }
    if(step_distance_p_values[ip] == 2.){
      BOOST_TEST(optimized_one.optimized > optimized_one.original);
      BOOST_TEST(optimized_vb.optimized > optimized_vb.original);
    }
    fmt::print("   {:6.2f} {:+12.3f}% {:>10s} {:+12.3f}% {:>10s}\n",
               step_distance_p_values[ip],pct_one,
               pct_one < 0. ? "ACCEPT" : "REJECT",
               pct_vb,pct_vb < 0. ? "ACCEPT" : "REJECT");
    log << "optimized_stepdistance p " << step_distance_p_values[ip]
        << " one_point " << optimized_one.original << " -> "
        << optimized_one.optimized
        << " point [" << optimized_one.x << " " << optimized_one.y << "]"
        << " vertex_barycenter " << optimized_vb.original << " -> "
        << optimized_vb.optimized
        << " point [" << optimized_vb.x << " " << optimized_vb.y << "]\n";
  }
  fmt::print("   detailed log: {}\n",(dir / "a6_face214_edge0.log").string());

  cleanup_trace_objects<MFT>(run,cav,seed,handler,lquae);
  return values;
}

template<class MFT>
CavityIntegrationValues dispatch_compare_exact_case_integration_schemes(
    const std::string& cmd){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c,arg.v);
  MeshBase& msh = *run.msh_g;

  CavityIntegrationValues values;
  bool dispatched = false;
  CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      values = compare_exact_case_integration_schemes_impl<MFT,gdim,ideg>(cmd);
      dispatched = true;
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
  BOOST_REQUIRE(dispatched);
  return values;
}

template<class MFT, int gdim, int ideg,
         QuaFun iquaf = DefaultQualityFunction>
InsertAttemptResult attempt_insert(const std::string& cmd,
                                   int tdim,
                                   int ientt,
                                   int ied,
                                   bool length_based,
                                   double worsen_pctg,
                                   double lenqua_short_max){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
  setup_quality_mesh<MFT,gdim,ideg>(msh);

  const int ithrd1 = 1;
  const int ithrd2 = 2;

  bool iinva = false;
  double qmin = 0, qmax = 0, qavg = 0;
  dblAr1 lquae(msh.nentt(tdim));

  getmetquamesh<MFT,iquaf>(msh,tdim,AsDeg::P1,AsDeg::P1,
                            &iinva,&qmin,&qmax,&qavg,&lquae);

  BadEntHandler handler(tdim, 100., 0.00001);
  const intAr2& ent2poi = msh.ent2poi(tdim);
  handler.setCallbacks([&](int ientt_){ return lquae[ientt_]; },
                       [&](int ientt_){ return isdeadent(ientt_,ent2poi); });
  std::vector<int> sorted_ids(msh.nentt(tdim));
  std::iota(sorted_ids.begin(), sorted_ids.end(), 0);
  std::sort(sorted_ids.begin(), sorted_ids.end(),
            [&](int a, int b){ return lquae[a] > lquae[b]; });
  handler.seedFromSortedIDs(sorted_ids);

  MshCavity cav(100,100,1);
  CavWrkArrs work;
  intAr1 lcaverr(CAV_ERR_NERROR);
  EdgeSeed insertion_seed(msh, cav, tdim, tdim, ientt, ied);

  InsertAttemptResult result;
  result.ierro = insertEdge<MFT,iquaf>(msh, insertion_seed,
                                      lenqua_short_max, false,
                                      cav, work, lcaverr, handler,
                                      length_based, worsen_pctg,
                                      ithrd1, ithrd2);
  result.ipins = cav.ipins;
  result.ncedg = cav.lcedg.get_n();
  result.ncfac = cav.lcfac.get_n();
  result.nctet = cav.lctet.get_n();
  return result;
}

template<class MFT, int gdim, int ideg>
ObjectiveInsertSummary compare_objective_insertions(const std::string& cmd){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
  setup_quality_mesh<MFT,gdim,ideg>(msh);

  const int tdim = msh.get_tdim();
  BOOST_REQUIRE_EQUAL(tdim,2);

  bool iinva_size = false;
  bool iinva_step = false;
  double qmin_size = 0., qmax_size = 0., qavg_size = 0.;
  double qmin_step = 0., qmax_step = 0., qavg_step = 0.;
  dblAr1 lquae_size(msh.nentt(tdim));
  dblAr1 lquae_step(msh.nentt(tdim));
  getmetquamesh<MFT,QuaFun::SizeShape>(
      msh,tdim,AsDeg::P1,AsDeg::P1,
      &iinva_size,&qmin_size,&qmax_size,&qavg_size,&lquae_size);
  getmetquamesh<MFT,QuaFun::StepDistance>(
      msh,tdim,AsDeg::P1,AsDeg::P1,
      &iinva_step,&qmin_step,&qmax_step,&qavg_step,&lquae_step);
  BOOST_REQUIRE(!iinva_size);
  BOOST_REQUIRE(!iinva_step);

  constexpr double legacy_smoothing_threshold = 0.01;
  const intAr2& ent2poi = msh.ent2poi(tdim);
  int alive_elements = 0;
  int size_shape_below_threshold = 0;
  int step_distance_below_threshold = 0;
  for (int ientt = 0; ientt < msh.nentt(tdim); ientt++)
  {
    if (isdeadent(ientt,ent2poi)) continue;
    alive_elements++;
    if (lquae_size[ientt] < legacy_smoothing_threshold)
    {
      size_shape_below_threshold++;
    }
    if (lquae_step[ientt] < legacy_smoothing_threshold)
    {
      step_distance_below_threshold++;
    }
  }

  lenStat lenstat0;
  intAr2 ilned;
  dblAr1 rlned;
  getLengthEdges<MFT>(msh,tdim,-1,ilned,rlned,lenstat0);
  const double lenqua_short_max =
      (lenstat0.qua_short + lenstat0.qua_long) / 2.;
  const double length_threshold = objective_compare_length_threshold();

  const int nedgl = tdim * (tdim + 1) / 2;
  const intAr2 lnoed(nedgl,2,lnoed2[0]);
  std::set<std::pair<int,int>> seen_edges;
  std::vector<ObjectiveInsertCandidate> candidates;

  for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
    if(isdeadent(ientt,ent2poi)) continue;
    for(int ied = 0; ied < nedgl; ied++){
      int ip1 = ent2poi(ientt,lnoed(ied,0));
      int ip2 = ent2poi(ientt,lnoed(ied,1));
      if(ip2 < ip1) std::swap(ip1,ip2);
      if(!seen_edges.emplace(ip1,ip2).second) continue;

      double sz[2];
      const double len = getlenedg_geosz<MFT,gdim,ideg>(
          msh,ientt,tdim,ied,sz);
      if(len <= length_threshold) continue;

      candidates.push_back({ientt,ied,ip1,ip2,len,
                            lquae_size[ientt],lquae_step[ientt]});
    }
  }

  std::sort(candidates.begin(),candidates.end(),
            [](const ObjectiveInsertCandidate& left,
               const ObjectiveInsertCandidate& right){
              if(left.len != right.len) return left.len > right.len;
              if(left.ip1 != right.ip1) return left.ip1 < right.ip1;
              return left.ip2 < right.ip2;
            });
  const int max_candidates = objective_compare_max_candidates();
  if(max_candidates > 0 && static_cast<int>(candidates.size()) > max_candidates){
    candidates.resize(max_candidates);
  }

  ObjectiveInsertSummary summary;
  summary.n_candidates = candidates.size();
  summary.comparisons.reserve(candidates.size());

  fmt::print("\n-- SizeShape / StepDistance initial insertion scan\n");
  fmt::print("   unique long-edge candidates : {}\n",candidates.size());
  fmt::print("   metric length threshold     : {:.16e}\n",length_threshold);
  fmt::print("   lenqua_short_max            : {:.16e}\n",lenqua_short_max);
  fmt::print("   SizeShape quality min/avg/max    : {:.8e} {:.8e} {:.8e}\n",
             qmin_size,qavg_size,qmax_size);
  fmt::print("   StepDistance quality min/avg/max : {:.8e} {:.8e} {:.8e}\n",
             qmin_step,qavg_step,qmax_step);
  fmt::print("   Elements below legacy 0.01 gate : "
             "SizeShape {} ({:.2f}%), StepDistance {} ({:.2f}%)\n",
             size_shape_below_threshold,
             100.0*size_shape_below_threshold/alive_elements,
             step_distance_below_threshold,
             100.0*step_distance_below_threshold/alive_elements);

  for(std::size_t ii = 0; ii < candidates.size(); ii++){
    const ObjectiveInsertCandidate& candidate = candidates[ii];
    ObjectiveInsertComparison comparison;
    comparison.candidate = candidate;
    comparison.size_shape =
      attempt_insert<MFT,gdim,ideg,QuaFun::SizeShape>(
          cmd,tdim,candidate.ientt,candidate.ied,
          false,0.,lenqua_short_max);
    comparison.step_distance =
      attempt_insert<MFT,gdim,ideg,QuaFun::StepDistance>(
          cmd,tdim,candidate.ientt,candidate.ied,
          false,0.,lenqua_short_max);
    comparison.step_distance_length =
      attempt_insert<MFT,gdim,ideg,QuaFun::StepDistance>(
          cmd,tdim,candidate.ientt,candidate.ied,
          true,0.,lenqua_short_max);

    const bool size_success = comparison.size_shape.ierro < 0;
    const bool step_success = comparison.step_distance.ierro < 0;
    if(size_success && step_success){
      comparison.category = "both";
      summary.n_both++;
    }else if(size_success){
      comparison.category = "size_shape_only";
      summary.n_size_shape_only++;
    }else if(step_success){
      comparison.category = "step_distance_only";
      summary.n_step_distance_only++;
    }else{
      comparison.category = "neither";
      summary.n_neither++;
    }
    if(comparison.size_shape.ierro == 0) summary.n_size_shape_noop++;
    if(comparison.step_distance.ierro == 0) summary.n_step_distance_noop++;
    if(comparison.step_distance_length.ierro < 0){
      summary.n_step_distance_length_success++;
    }else{
      summary.step_distance_length_errors[
          comparison.step_distance_length.ierro]++;
    }
    if(!size_success) summary.size_shape_errors[comparison.size_shape.ierro]++;
    if(!step_success) summary.step_distance_errors[comparison.step_distance.ierro]++;
    summary.comparisons.push_back(comparison);

    if((ii + 1) % 25 == 0 || ii + 1 == candidates.size()){
      fmt::print("   scanned {}/{} candidates\n",ii + 1,candidates.size());
    }
  }

  const std::filesystem::path outdir = objective_compare_output_dir();
  std::filesystem::create_directories(outdir);
  const std::filesystem::path csv =
      outdir / ("a" + std::to_string(objective_compare_iteration()) + ".csv");
  std::ofstream fout(csv);
  BOOST_REQUIRE_MESSAGE(fout.good(),"Could not open comparison CSV: " + csv.string());
  fout << "ientt,ied,ip1,ip2,length,size_shape_quality,step_distance_quality,"
          "size_shape_ierro,step_distance_ierro,step_distance_length_ierro,"
          "category,size_shape_ipins,step_distance_ipins,"
          "step_distance_length_ipins,size_shape_nface,step_distance_nface,"
          "step_distance_length_nface\n";
  fout.precision(17);
  for(const ObjectiveInsertComparison& comparison : summary.comparisons){
    const auto& c = comparison.candidate;
    fout << c.ientt << ',' << c.ied << ',' << c.ip1 << ',' << c.ip2 << ','
         << c.len << ',' << c.size_shape_quality << ','
         << c.step_distance_quality << ','
         << comparison.size_shape.ierro << ','
         << comparison.step_distance.ierro << ','
         << comparison.step_distance_length.ierro << ','
         << comparison.category << ','
         << comparison.size_shape.ipins << ','
         << comparison.step_distance.ipins << ','
         << comparison.step_distance_length.ipins << ','
         << comparison.size_shape.ncfac << ','
         << comparison.step_distance.ncfac << ','
         << comparison.step_distance_length.ncfac << '\n';
  }

  fmt::print("\n-- Initial insertion comparison summary\n");
  fmt::print("   both succeed       : {}\n",summary.n_both);
  fmt::print("   SizeShape only     : {}\n",summary.n_size_shape_only);
  fmt::print("   StepDistance only  : {}\n",summary.n_step_distance_only);
  fmt::print("   neither succeeds   : {}\n",summary.n_neither);
  fmt::print("   SizeShape no-op    : {}\n",summary.n_size_shape_noop);
  fmt::print("   StepDistance no-op : {}\n",summary.n_step_distance_noop);
  fmt::print("   StepDistance length fallback succeeds: {}\n",
             summary.n_step_distance_length_success);
  fmt::print("   CSV                 : {}\n",csv.string());
  fmt::print("   SizeShape non-success return codes:");
  for(const auto& [ierro,count] : summary.size_shape_errors){
    fmt::print(" {}:{}",ierro,count);
  }
  fmt::print("\n   StepDistance non-success return codes:");
  for(const auto& [ierro,count] : summary.step_distance_errors){
    fmt::print(" {} ({}):{}",ierro,insertion_error_name(ierro),count);
  }
  fmt::print("\n   StepDistance length non-success return codes:");
  for(const auto& [ierro,count] : summary.step_distance_length_errors){
    fmt::print(" {} ({}):{}",ierro,insertion_error_name(ierro),count);
  }
  fmt::print("\n");

  return summary;
}

template<class MFT>
ObjectiveInsertSummary dispatch_compare_objective_insertions(
    const std::string& cmd){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c,arg.v);
  MeshBase& msh = *run.msh_g;

  ObjectiveInsertSummary summary;
  bool dispatched = false;
  CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      summary = compare_objective_insertions<MFT,gdim,ideg>(cmd);
      dispatched = true;
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
  BOOST_REQUIRE(dispatched);
  return summary;
}

bool common_scan_accepts(const CommonObjectiveOptimization& result,
                         bool optimized){
  const double candidate = optimized ? result.optimized : result.initial;
  if(!std::isfinite(result.original) || !std::isfinite(candidate)
      || result.original <= 0.){
    return false;
  }
  // Match BadEntHandler(alpha = 0.00001 percent).
  return candidate <= (1. - 1.e-7)*result.original;
}

double common_scan_ratio(const CommonObjectiveOptimization& result,
                         bool optimized){
  const double candidate = optimized ? result.optimized : result.initial;
  if(!std::isfinite(result.original) || result.original <= 0.){
    return std::numeric_limits<double>::infinity();
  }
  return candidate/result.original;
}

template<class MFT, int gdim, int ideg>
CommonPScanSummary scan_common_vb_objectives(const std::string& cmd){
  static_assert(gdim == 2);
  constexpr int tdim = 2;
  constexpr int ithrd1 = 1;
  constexpr int ithrd2 = 2;

  cargHandler arg(cmd);
  MetrisRunner initial_run(arg.c,arg.v);
  Mesh<MFT>& initial_msh = static_cast<Mesh<MFT>&>(*initial_run.msh_g);
  setup_quality_mesh<MFT,gdim,ideg>(initial_msh);

  bool iinva_size = false;
  bool iinva_step = false;
  double qmin = 0., qmax = 0., qavg = 0.;
  dblAr1 lquae_size(initial_msh.nface);
  dblAr1 lquae_step(initial_msh.nface);
  getmetquamesh<MFT,QuaFun::SizeShape>(
      initial_msh,tdim,AsDeg::P1,AsDeg::P1,
      &iinva_size,&qmin,&qmax,&qavg,&lquae_size);
  getmetquamesh<MFT,QuaFun::StepDistance>(
      initial_msh,tdim,AsDeg::P1,AsDeg::P1,
      &iinva_step,&qmin,&qmax,&qavg,&lquae_step);
  BOOST_REQUIRE(!iinva_size);
  BOOST_REQUIRE(!iinva_step);

  const double length_threshold = objective_compare_length_threshold();
  const intAr2 lnoed(3,2,lnoed2[0]);
  std::set<std::pair<int,int>> seen_edges;
  std::vector<ObjectiveInsertCandidate> candidates;
  for(int iface = 0; iface < initial_msh.nface; iface++){
    if(isdeadent(iface,initial_msh.fac2poi)) continue;
    for(int iedge = 0; iedge < 3; iedge++){
      int ip1 = initial_msh.fac2poi(iface,lnoed(iedge,0));
      int ip2 = initial_msh.fac2poi(iface,lnoed(iedge,1));
      if(ip2 < ip1) std::swap(ip1,ip2);
      if(!seen_edges.emplace(ip1,ip2).second) continue;
      double sizes[2];
      const double length = getlenedg_geosz<MFT,gdim,ideg>(
          initial_msh,iface,tdim,iedge,sizes);
      if(length <= length_threshold) continue;
      candidates.push_back({iface,iedge,ip1,ip2,length,
                            lquae_size[iface],lquae_step[iface]});
    }
  }
  std::sort(candidates.begin(),candidates.end(),
            [](const ObjectiveInsertCandidate& left,
               const ObjectiveInsertCandidate& right){
              if(left.len != right.len) return left.len > right.len;
              if(left.ip1 != right.ip1) return left.ip1 < right.ip1;
              return left.ip2 < right.ip2;
            });
  const int max_candidates = objective_compare_max_candidates();
  if(max_candidates > 0 && static_cast<int>(candidates.size()) > max_candidates){
    candidates.resize(max_candidates);
  }

  CommonPScanSummary summary;
  summary.n_candidates = candidates.size();
  summary.rows.reserve(candidates.size());
  fmt::print("\n-- Common-quadrature p scan on initial a{} mesh\n",
             objective_compare_iteration());
  fmt::print("   rule                  : vertices + barycenter\n");
  fmt::print("   point optimization    : derivative-free, metric refreshed\n");
  fmt::print("   minimum child-area fraction guard: {:.3e}\n",
             objective_compare_min_child_area_fraction());
  fmt::print("   SD p=1 metric-volume barrier: rho0={:.3e}, beta={:.3e}\n",
             objective_compare_step_barrier_rho0(),
             objective_compare_step_barrier_beta());
  fmt::print("   unique long edges     : {}\n",candidates.size());

  for(std::size_t icandidate = 0; icandidate < candidates.size(); icandidate++){
    const ObjectiveInsertCandidate& candidate = candidates[icandidate];
    CommonPScanRow row;
    row.candidate = candidate;

    cargHandler candidate_arg(cmd);
    MetrisRunner candidate_run(candidate_arg.c,candidate_arg.v);
    Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*candidate_run.msh_g);
    setup_quality_mesh<MFT,gdim,ideg>(msh);

    MshCavity cav(100,100,1);
    EdgeSeed seed(msh,cav,tdim,tdim,candidate.ientt,candidate.ied);
    cav.ipins = msh.newpoint(PointType::Vertex,seed.tdimp,seed.iseed);
    int ierro = aux_bisecPointLen(msh,seed,msh.poi2bpo[cav.ipins],false,cav);
    if(ierro == 0){
      ierro = increase_cavity(msh,cav,false,ithrd1,ithrd2);
    }
    if(ierro > 0){
      row.setup_error = ierro;
      summary.n_setup_failures++;
      summary.rows.push_back(row);
      continue;
    }

    row.size_shape = optimize_common_objective_point_2d(
        msh,cav,ithrd1,CommonScanObjective::SizeShapeP1);
    row.step_p1 = optimize_common_objective_point_2d(
        msh,cav,ithrd1,CommonScanObjective::StepDistanceP1);
    row.step_p2 = optimize_common_objective_point_2d(
        msh,cav,ithrd1,CommonScanObjective::StepDistanceP2);

    if(!std::isfinite(row.size_shape.initial)
        || !std::isfinite(row.step_p1.initial)
        || !std::isfinite(row.step_p2.initial)){
      row.setup_error = -1001;
      summary.n_setup_failures++;
      summary.rows.push_back(row);
      continue;
    }

    summary.n_valid++;
    const bool size_before = common_scan_accepts(row.size_shape,false);
    const bool step_p1_before = common_scan_accepts(row.step_p1,false);
    const bool step_p2_before = common_scan_accepts(row.step_p2,false);
    const bool size_after = common_scan_accepts(row.size_shape,true);
    const bool step_p1_after = common_scan_accepts(row.step_p1,true);
    const bool step_p2_after = common_scan_accepts(row.step_p2,true);
    summary.n_size_before += size_before;
    summary.n_step_p1_before += step_p1_before;
    summary.n_step_p2_before += step_p2_before;
    summary.n_size_after += size_after;
    summary.n_step_p1_after += step_p1_after;
    summary.n_step_p2_after += step_p2_after;

    if(size_after && step_p1_after){
      summary.n_size_step_p1_both++;
    }else if(size_after){
      summary.n_size_only_vs_step_p1++;
    }else if(step_p1_after){
      summary.n_step_p1_only_vs_size++;
    }else{
      summary.n_size_step_p1_neither++;
    }
    if(step_p1_after && !step_p2_after) summary.n_step_p1_rescues_step_p2++;
    if(step_p2_after && !step_p1_after) summary.n_step_p2_only_vs_step_p1++;
    summary.rows.push_back(row);

    if((icandidate + 1) % 10 == 0 || icandidate + 1 == candidates.size()){
      fmt::print("   scanned {}/{} candidates\n",
                 icandidate + 1,candidates.size());
    }
  }

  const std::filesystem::path outdir = objective_compare_output_dir();
  std::filesystem::create_directories(outdir);
  const double min_area_guard = objective_compare_min_child_area_fraction();
  const std::string guard_suffix = min_area_guard > 0.
      ? "_guard_" + std::to_string(min_area_guard) : "";
  const double barrier_beta = objective_compare_step_barrier_beta();
  const double barrier_rho0 = objective_compare_step_barrier_rho0();
  const std::string barrier_suffix = barrier_beta > 0.
      ? fmt::format("_barrier_rho0_{:.2e}_beta_{:.2e}",
                    barrier_rho0,barrier_beta)
      : "";
  const std::filesystem::path csv = outdir /
      ("a" + std::to_string(objective_compare_iteration())
       + "_common_vb_p_scan" + guard_suffix + barrier_suffix + ".csv");
  std::ofstream fout(csv);
  BOOST_REQUIRE_MESSAGE(fout.good(),"Could not open p-scan CSV: " + csv.string());
  fout.precision(17);
  fout << "ientt,ied,ip1,ip2,length,setup_error,"
          "ss_old,ss_midpoint,ss_optimized,ss_midpoint_ratio,ss_optimized_ratio,ss_x,ss_y,ss_min_vertex_distance_ratio,ss_min_child_area_fraction,"
          "sd1_old,sd1_midpoint,sd1_optimized,sd1_midpoint_ratio,sd1_optimized_ratio,sd1_x,sd1_y,sd1_min_vertex_distance_ratio,sd1_min_child_area_fraction,sd1_min_metric_volume,sd1_old_barrier_unscaled,sd1_optimized_barrier_unscaled,"
          "sd2_old,sd2_midpoint,sd2_optimized,sd2_midpoint_ratio,sd2_optimized_ratio,sd2_x,sd2_y,sd2_min_vertex_distance_ratio,sd2_min_child_area_fraction,"
          "ss_accept_midpoint,sd1_accept_midpoint,sd2_accept_midpoint,"
          "ss_accept_optimized,sd1_accept_optimized,sd2_accept_optimized\n";
  for(const CommonPScanRow& row : summary.rows){
    const auto& c = row.candidate;
    fout << c.ientt << ',' << c.ied << ',' << c.ip1 << ',' << c.ip2 << ','
         << c.len << ',' << row.setup_error << ','
         << row.size_shape.original << ',' << row.size_shape.initial << ','
         << row.size_shape.optimized << ','
         << common_scan_ratio(row.size_shape,false) << ','
         << common_scan_ratio(row.size_shape,true) << ','
         << row.size_shape.x << ',' << row.size_shape.y << ','
         << row.size_shape.min_vertex_distance_ratio << ','
         << row.size_shape.min_child_area_fraction << ','
         << row.step_p1.original << ',' << row.step_p1.initial << ','
         << row.step_p1.optimized << ','
         << common_scan_ratio(row.step_p1,false) << ','
         << common_scan_ratio(row.step_p1,true) << ','
         << row.step_p1.x << ',' << row.step_p1.y << ','
         << row.step_p1.min_vertex_distance_ratio << ','
         << row.step_p1.min_child_area_fraction << ','
         << row.step_p1.min_metric_volume << ','
         << row.step_p1.original_barrier_unscaled << ','
         << row.step_p1.optimized_barrier_unscaled << ','
         << row.step_p2.original << ',' << row.step_p2.initial << ','
         << row.step_p2.optimized << ','
         << common_scan_ratio(row.step_p2,false) << ','
         << common_scan_ratio(row.step_p2,true) << ','
         << row.step_p2.x << ',' << row.step_p2.y << ','
         << row.step_p2.min_vertex_distance_ratio << ','
         << row.step_p2.min_child_area_fraction << ','
         << common_scan_accepts(row.size_shape,false) << ','
         << common_scan_accepts(row.step_p1,false) << ','
         << common_scan_accepts(row.step_p2,false) << ','
         << common_scan_accepts(row.size_shape,true) << ','
         << common_scan_accepts(row.step_p1,true) << ','
         << common_scan_accepts(row.step_p2,true) << '\n';
  }

  fmt::print("\n-- Common-quadrature p scan summary\n");
  fmt::print("   candidates / valid / setup failures : {} / {} / {}\n",
             summary.n_candidates,summary.n_valid,summary.n_setup_failures);
  fmt::print("   accepted at midpoint  SS / SD1 / SD2: {} / {} / {}\n",
             summary.n_size_before,summary.n_step_p1_before,
             summary.n_step_p2_before);
  fmt::print("   accepted after optimize SS / SD1 / SD2: {} / {} / {}\n",
             summary.n_size_after,summary.n_step_p1_after,
             summary.n_step_p2_after);
  fmt::print("   SizeShape & SD p=1 both              : {}\n",
             summary.n_size_step_p1_both);
  fmt::print("   SizeShape only vs SD p=1             : {}\n",
             summary.n_size_only_vs_step_p1);
  fmt::print("   SD p=1 only vs SizeShape             : {}\n",
             summary.n_step_p1_only_vs_size);
  fmt::print("   neither SizeShape nor SD p=1         : {}\n",
             summary.n_size_step_p1_neither);
  fmt::print("   SD p=1 accepts while SD p=2 rejects  : {}\n",
             summary.n_step_p1_rescues_step_p2);
  fmt::print("   SD p=2 accepts while SD p=1 rejects  : {}\n",
             summary.n_step_p2_only_vs_step_p1);
  fmt::print("   CSV                                   : {}\n",csv.string());

  auto print_mismatches = [&](bool size_only){
    std::vector<const CommonPScanRow*> mismatches;
    for(const CommonPScanRow& row : summary.rows){
      if(row.setup_error != 0) continue;
      const bool size_accept = common_scan_accepts(row.size_shape,true);
      const bool step_accept = common_scan_accepts(row.step_p1,true);
      if(size_only ? (size_accept && !step_accept)
                   : (step_accept && !size_accept)){
        mismatches.push_back(&row);
      }
    }
    std::sort(mismatches.begin(),mismatches.end(),
              [&](const CommonPScanRow* left,const CommonPScanRow* right){
                const double left_gap = std::abs(
                    common_scan_ratio(left->size_shape,true)
                    - common_scan_ratio(left->step_p1,true));
                const double right_gap = std::abs(
                    common_scan_ratio(right->size_shape,true)
                    - common_scan_ratio(right->step_p1,true));
                return left_gap > right_gap;
              });
    fmt::print("\n   {} mismatches (largest ratio gaps):\n",
               size_only ? "SizeShape-only" : "StepDistance-p1-only");
    const std::size_t nprint = std::min<std::size_t>(mismatches.size(),8);
    for(std::size_t ii = 0; ii < nprint; ii++){
      const CommonPScanRow& row = *mismatches[ii];
      fmt::print("     face {:4d} edge {} nodes {:3d}-{:3d} len {:.4f}"
                 "  SS {:+8.3f}%  SD1 {:+8.3f}%  SD2 {:+8.3f}%"
                 "  SD1 min-area {:.2e} min-rho {:.2e} min-dist {:.2e}\n",
                 row.candidate.ientt,row.candidate.ied,
                 row.candidate.ip1,row.candidate.ip2,row.candidate.len,
                 100.*(common_scan_ratio(row.size_shape,true)-1.),
                 100.*(common_scan_ratio(row.step_p1,true)-1.),
                 100.*(common_scan_ratio(row.step_p2,true)-1.),
                 row.step_p1.min_child_area_fraction,
                 row.step_p1.min_metric_volume,
                 row.step_p1.min_vertex_distance_ratio);
    }
  };
  print_mismatches(true);
  print_mismatches(false);
  return summary;
}

template<class MFT>
CommonPScanSummary dispatch_common_vb_objective_scan(const std::string& cmd){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c,arg.v);
  MeshBase& msh = *run.msh_g;
  CommonPScanSummary summary;
  bool dispatched = false;
  CT_FOR0_INC(2,2,gdim){if(gdim == msh.idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      summary = scan_common_vb_objectives<MFT,gdim,ideg>(cmd);
      dispatched = true;
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
  BOOST_REQUIRE(dispatched);
  return summary;
}

template<class MFT, int gdim, int ideg, QuaFun iquaf>
InsertAttemptResult attempt_quality_insert_with_smoothing(const std::string& cmd,
                                                          int tdim,
                                                          int ientt,
                                                          int ied,
                                                          bool smooth_initial){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
  setup_quality_mesh<MFT,gdim,ideg>(msh);

  const int ithrd1 = 1;
  const int ithrd2 = 2;

  bool iinva = false;
  double qmin = 0, qmax = 0, qavg = 0;
  dblAr1 lquae(msh.nentt(tdim));
  getmetquamesh<MFT,iquaf>(msh,tdim,AsDeg::P1,AsDeg::P1,
                           &iinva,&qmin,&qmax,&qavg,&lquae);

  BadEntHandler handler(tdim, 100., 0.00001);
  const intAr2& ent2poi = msh.ent2poi(tdim);
  handler.setCallbacks([&](int ientt_){ return lquae[ientt_]; },
                       [&](int ientt_){ return isdeadent(ientt_,ent2poi); });
  std::vector<int> sorted_ids(msh.nentt(tdim));
  std::iota(sorted_ids.begin(), sorted_ids.end(), 0);
  std::sort(sorted_ids.begin(), sorted_ids.end(),
            [&](int a, int b){ return lquae[a] > lquae[b]; });
  handler.seedFromSortedIDs(sorted_ids);

  MshCavity cav(100,100,1);
  EdgeSeed insertion_seed(msh, cav, tdim, tdim, ientt, ied);
  InsertAttemptResult result;
  if(tdim != 2){
    result.ierro = INS2D_ERR_NOQUALIMPROV;
    return result;
  }

  cav.ipins = msh.newpoint(PointType::Vertex,
                           insertion_seed.tdimp,
                           insertion_seed.iseed);
  result.ipins = cav.ipins;

  int ierro = aux_bisecPointLen(msh, insertion_seed,
                                msh.poi2bpo[cav.ipins], false, cav);
  if(ierro != 0){
    result.ierro = ierro;
    return result;
  }

  ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
  if(ierro > 0){
    result.ierro = ierro;
    return result;
  }

  if(smooth_initial){
    CavityQualityStats initial_stats =
      compute_cavity_quality_2d<MFT,iquaf>(msh, cav, ithrd1,
                                           "initial smoothing before", nullptr);

    double smooth_qsum = initial_stats.qsum1;
    double smooth_qmax = initial_stats.qmax1;
    double smooth_target_weight = 0.;

    const int nface0 = msh.nface;
    msh.set_nface(nface0 + 1);
    smoothCavity(msh, cav, handler, iquaf,
                 initial_stats.qsum1,
                 initial_stats.qmax1,
                 0.,
                 smooth_qsum,
                 smooth_qmax,
                 smooth_target_weight,
                 ithrd1, ithrd2);
    msh.fac2poi(nface0,0) = -1;
    msh.set_nface(nface0);
  }

  ierro = increase_cavity_quality(msh, cav, tdim, 5, handler, ithrd1);
  result.ncedg = cav.lcedg.get_n();
  result.ncfac = cav.lcfac.get_n();
  result.nctet = cav.lctet.get_n();

  #ifdef CAVSMOOTHING
  result.ierro = ierro;
  return result;
  #endif

  CavityQualityStats before_smoothing =
    compute_cavity_quality_2d<MFT,iquaf>(msh, cav, ithrd1,
                                         "smoothed scan before", nullptr);

  double smooth_qsum = before_smoothing.qsum1;
  double smooth_qmax = before_smoothing.qmax1;
  double smooth_target_weight = 0.;

  const int nface0 = msh.nface;
  msh.set_nface(nface0 + 1);
  smoothCavity(msh, cav, handler, iquaf,
               before_smoothing.qsum1,
               before_smoothing.qmax1,
               0.,
               smooth_qsum,
               smooth_qmax,
               smooth_target_weight,
               ithrd1, ithrd2);
  msh.fac2poi(nface0,0) = -1;
  msh.set_nface(nface0);

  result.ncedg = cav.lcedg.get_n();
  result.ncfac = cav.lcfac.get_n();
  result.nctet = cav.lctet.get_n();

  result.ierro = handler.checkSuccess(smooth_qsum, before_smoothing.qsum0)
               ? 0
               : INS2D_ERR_NOQUALIMPROV;
  return result;
}

template<class MFT, int gdim, int ideg>
DivergenceCase find_quality_insert_divergence(const std::string& cmd,
                                              bool smooth_quality,
                                              double length_worsen_pctg){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
  setup_quality_mesh<MFT,gdim,ideg>(msh);

  const int tdim = msh.get_tdim();
  BOOST_REQUIRE(tdim == 2 || tdim == 3);

  // Match MetrisRunner::adaptMesh2: the dimension-ordered variant adapts
  // CAD edges immediately before entering the first 2D quality pass.
  if constexpr(gdim == 2){
    if(tdim == 2 && msh.CAD() && msh.param->adp_line_adapt){
      adaptGeoLines<MFT>(msh);
    }
  }

  bool iinva = false;
  double qmin = 0, qmax = 0, qavg = 0;
  dblAr1 lquae(msh.nentt(tdim));

  #ifdef STEPDISTANCE
  getmetquamesh<MFT, QuaFun::StepDistance>(msh,tdim,AsDeg::P1,AsDeg::P1,
                                           &iinva,&qmin,&qmax,&qavg,&lquae);
  #else
  getmetquamesh<MFT, QuaFun::SizeShape>(msh,tdim,AsDeg::P1,AsDeg::P1,
                                        &iinva,&qmin,&qmax,&qavg,&lquae);
  #endif

  lenStat lenstat0;
  intAr2 ilned;
  dblAr1 rlned;
  getLengthEdges<MFT>(msh,tdim,-1,ilned,rlned,lenstat0);
  const double lenqua_short_max = (lenstat0.qua_short + lenstat0.qua_long) / 2;
  const bool smooth_initial = scan_smooth_initial_cavity();

  std::vector<int> sorted_ids(msh.nentt(tdim));
  std::iota(sorted_ids.begin(), sorted_ids.end(), 0);
  std::sort(sorted_ids.begin(), sorted_ids.end(),
            [&](int a, int b){ return lquae[a] > lquae[b]; });

  BadEntHandler handler( tdim, 100., 0.00001 );
  const intAr2& ent2poi = msh.ent2poi(tdim);
  handler.setCallbacks([&](int ientt){ return lquae[ientt]; },
                       [&](int ientt){ return isdeadent(ientt,ent2poi); });
  handler.seedFromSortedIDs(sorted_ids);

  const double length_threshold = std::sqrt(2.0);
  const double len_upper_bound = length_threshold;
  const int nedgl = tdim * (tdim + 1) / 2;
  const intAr2 lnoed(nedgl, 2,
                     tdim == 1 ? lnoed1[0] :
                     tdim == 2 ? lnoed2[0] :
                                 lnoed3[0]);

  for(const auto& ent_qual : handler.K){
    const int ientt = ent_qual.ientt;
    if(isdeadent(ientt, ent2poi)) continue;

    std::vector<EdgeOpCandidate> candidates;
    candidates.reserve(nedgl);

    for(int ied = 0; ied < nedgl; ied++){
      double sz[2];
      double elen = getlenedg_geosz<MFT,gdim,ideg>(msh, ientt, tdim, ied, sz);
      if(elen > len_upper_bound){
        candidates.push_back({ientt, ied, ent_qual.qentt, elen, std::log(elen)});
      }
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const EdgeOpCandidate& a, const EdgeOpCandidate& b){
                return a.dev > b.dev;
              });

    for(const EdgeOpCandidate& cand : candidates){
      InsertAttemptResult quality;
      if(smooth_quality){
        #if defined(CAVSMOOTHING)
        if(!smooth_initial){
        quality =
          attempt_insert<MFT,gdim,ideg>(cmd, tdim, cand.ientt, cand.ied,
                                        false, 0., lenqua_short_max);
        }else
        #endif
        {
          #ifdef STEPDISTANCE
          quality =
            attempt_quality_insert_with_smoothing<MFT,gdim,ideg,QuaFun::StepDistance>(
              cmd, tdim, cand.ientt, cand.ied, smooth_initial);
          #else
          quality =
            attempt_quality_insert_with_smoothing<MFT,gdim,ideg,QuaFun::SizeShape>(
              cmd, tdim, cand.ientt, cand.ied, smooth_initial);
          #endif
        }
      }else{
        quality =
          attempt_insert<MFT,gdim,ideg>(cmd, tdim, cand.ientt, cand.ied,
                                        false, 0., lenqua_short_max);
      }
      if(quality.ierro <= 0) continue;

      InsertAttemptResult length =
        attempt_insert<MFT,gdim,ideg>(cmd, tdim, cand.ientt, cand.ied,
                                      true, length_worsen_pctg,
                                      lenqua_short_max);
      if(length.ierro <= 0){
        DivergenceCase found;
        found.found = true;
        found.tdim = tdim;
        found.ientt = cand.ientt;
        found.ied = cand.ied;
        found.ip1 = ent2poi(cand.ientt, lnoed(cand.ied,0));
        found.ip2 = ent2poi(cand.ientt, lnoed(cand.ied,1));
        found.qentt = cand.qentt;
        found.len = cand.len;
        found.ierro_quality = quality.ierro;
        found.ierro_length = length.ierro;
        found.quality = quality;
        found.length = length;
        return found;
      }
    }
  }

  return DivergenceCase{};
}

template<class MFT>
DivergenceCase dispatch_find_case(const std::string& cmd,
                                  bool smooth_quality = false,
                                  double length_worsen_pctg = 0.){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  MeshBase& msh = *run.msh_g;

  DivergenceCase found;
  CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      found = find_quality_insert_divergence<MFT,gdim,ideg>(
        cmd, smooth_quality, length_worsen_pctg);
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);

  return found;
}

template<class MFT, int gdim, int ideg>
StatefulAdaptResult find_stateful_quality_adapt_divergence(const std::string& cmd,
                                                           bool save_pre_state = false){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
  setup_quality_mesh<MFT,gdim,ideg>(msh);

  const int tdim = msh.get_tdim();
  BOOST_REQUIRE(tdim == 2 || tdim == 3);

  // Match MetrisRunner::adaptMesh2: the dimension-ordered variant adapts
  // CAD edges immediately before entering the first 2D quality pass.
  if constexpr(gdim == 2){
    if(tdim == 2 && msh.CAD() && msh.param->adp_line_adapt){
      adaptGeoLines<MFT>(msh);
    }
  }

  const int ithrd1 = 1;
  const int ithrd2 = 2;
  const int ithrd3 = 3;
  const int ithrdfro = 0;
  msh.tag[ithrdfro]++;

  const double alpha = 0.00001;
  const double badX = 100;
  const double lengthThreshold = std::sqrt(2.0);

  bool iinva = false;
  double qmin = 0, qmax = 0, qavg = 0;
  dblAr1 lquae(msh.nentt(tdim));

  #ifdef STEPDISTANCE
  getmetquamesh<MFT, QuaFun::StepDistance>(msh,tdim,AsDeg::P1,AsDeg::P1,
                                           &iinva,&qmin,&qmax,&qavg,&lquae);
  #else
  getmetquamesh<MFT, QuaFun::SizeShape>(msh,tdim,AsDeg::P1,AsDeg::P1,
                                        &iinva,&qmin,&qmax,&qavg,&lquae);
  #endif

  lenStat lenstat0;
  intAr2 ilned;
  dblAr1 rlned;
  getLengthEdges<MFT>(msh,tdim,-1,ilned,rlned,lenstat0);
  const double lenqua_short_max = (lenstat0.qua_short + lenstat0.qua_long) / 2;

  std::vector<int> sorted_ids(msh.nentt(tdim));
  std::iota(sorted_ids.begin(), sorted_ids.end(), 0);
  std::sort(sorted_ids.begin(), sorted_ids.end(),
            [&](int a, int b){ return lquae[a] > lquae[b]; });

  BadEntHandler handler(tdim, badX, alpha);
  const intAr2& ent2poi = msh.ent2poi(tdim);
  const intAr2& ent2ent = msh.ent2ent(tdim);
        intAr2& ent2tag = msh.ent2tag(tdim);
  handler.setCallbacks([&](int ientt){ return lquae[ientt]; },
                       [&](int ientt){ return isdeadent(ientt,ent2poi); });
  handler.seedFromSortedIDs(sorted_ids);

  MshCavity cav(100,100,1);
  CavWrkArrs work;
  intAr1 lcaverr(CAV_ERR_NERROR);

  StatefulAdaptResult result;
  const int max_iter = stateful_adapt_max_iter();
  const int trace_interval = stateful_adapt_trace_interval();
  const int trace_start = stateful_adapt_trace_start();
  const int save_iteration = stateful_adapt_save_iteration();
  const bool require_boundary =
      stateful_adapt_require_boundary_divergence();

  while(true){

    bool didOperation = false;

    for(auto itK = handler.K.begin(); itK != handler.K.end(); itK++){

      int ientt = itK->ientt;

      const double quaent = itK->qentt;
      result.iter++;
      if(result.iter > max_iter) return result;

      if(save_iteration > 0 && result.iter == save_iteration){
        const std::string save_dir = stateful_adapt_save_dir();
        BOOST_REQUIRE(!save_dir.empty());
        std::filesystem::create_directories(save_dir);
        writeMesh(save_dir + "/mesh_MOESS_initial_a1.meshb",msh);
        msh.met.writeMetricFile(save_dir + "/met_MOESS_initial_a1.solb");
        fmt::print("Saved stateful adaptation state at iteration {} to {}\n",
                   result.iter,save_dir);
        return result;
      }

      const bool trace_iteration = trace_interval > 0
                                && result.iter >= trace_start
                                && result.iter % trace_interval == 0;
      if(trace_iteration){
        lquae.set_n(msh.nentt(tdim));
        #ifdef STEPDISTANCE
        getmetquamesh<MFT, QuaFun::StepDistance>(
            msh,tdim,AsDeg::P1,AsDeg::P1,
            &iinva,&qmin,&qmax,&qavg,&lquae);
        #else
        getmetquamesh<MFT, QuaFun::SizeShape>(
            msh,tdim,AsDeg::P1,AsDeg::P1,
            &iinva,&qmin,&qmax,&qavg,&lquae);
        #endif

        double exact_sum = 0.;
        int exact_alive = 0;
        for(int ielem = 0; ielem < msh.nentt(tdim); ielem++){
          if(isdeadent(ielem,ent2poi)) continue;
          exact_sum += lquae[ielem];
          exact_alive++;
        }

        fmt::print(
            "TRACE iter={} npoin={} slots={} alive={} K={} "
            "ins={}/{} col={}/{} smoo={}/{} "
            "handler_obj={:.16e} exact_obj={:.16e} "
            "handler_count={} handler_num={:.16e} exact_num={:.16e} "
            "qworst={:.16e} ent={}\n",
            result.iter,msh.npoin,msh.nentt(tdim),exact_alive,handler.K.size(),
            result.nSuccessInsert,result.ntryInsert,
            result.nSuccessCollapse,result.ntryCollapse,
            result.nSuccessSmoothing,result.ntrySmoothing,
            handler.getQualitySum()/handler.getQualityCount(),
            exact_sum/exact_alive,
            handler.getQualityCount(),
            handler.getQualitySum(),exact_sum,
            quaent,ientt);
      }

      const int nedgl = tdim * (tdim + 1) / 2;
      const intAr2 lnoed(nedgl, 2,
                         tdim == 1 ? lnoed1[0] :
                         tdim == 2 ? lnoed2[0] :
                                     lnoed3[0]);

      struct LocalEdgeOpCandidate{
        int ied;
        double len;
        double dev;
        bool doInsert;
      };

      std::vector<LocalEdgeOpCandidate> candidates;
      candidates.reserve(nedgl);

      const double lenLowerBound = 1.0 / lengthThreshold;
      const double lenUpperBound = lengthThreshold;

      for(int ied = 0; ied < nedgl; ied++){
        double sz[2];
        double elen = getlenedg_geosz<MFT,gdim,ideg>(msh, ientt, tdim, ied, sz);

        if(elen > lenUpperBound){
          candidates.push_back({ied, elen, std::log(elen), true});
        }else if(elen < lenLowerBound){
          candidates.push_back({ied, elen, std::log(1.0 / elen), false});
        }
      }

      std::sort(candidates.begin(), candidates.end(),
                [](const LocalEdgeOpCandidate& a, const LocalEdgeOpCandidate& b){
                  return a.dev > b.dev;
                });

      if(trace_iteration){
        fmt::print("TRACE candidates={}",candidates.size());
        for(const LocalEdgeOpCandidate& cand : candidates){
          const int ip1 = ent2poi(ientt,lnoed(cand.ied,0));
          const int ip2 = ent2poi(ientt,lnoed(cand.ied,1));
          fmt::print(" [{} edge=({}, {}) len={:.16e}]",
                     cand.doInsert ? "insert" : "collapse",
                     ip1,ip2,cand.len);
        }
        fmt::print("\n");
      }

      for(const LocalEdgeOpCandidate& cand : candidates){

        const int ied = cand.ied;

        if(cand.doInsert){

          INCVDEPTH(msh.param);

          result.ntryInsert++;

          EdgeSeed insertionSeed(msh, cav, tdim, tdim, ientt, ied);
          int ierro = insertEdge(msh, insertionSeed, lenqua_short_max, false,
                                 cav, work, lcaverr, handler, false, 0.,
                                 ithrd1, ithrd2);

          if(ierro <= 0){
            didOperation = true;
            result.nSuccessInsert++;
            msh.poicstr[cav.ipins] = false;
            break;
          }

          result.ntryInsertLength++;
          int ierro_quality = ierro;

          std::filesystem::path pre_base;
          if(save_pre_state){
            std::filesystem::create_directories(trace_dir());
            pre_base = trace_dir() / "stateful_pre_length_success";
            writeMesh(pre_base.string() + ".meshb", msh);
            msh.met.writeMetricFile(pre_base.string() + ".solb");
          }

          ierro = insertEdge(msh, insertionSeed, lenqua_short_max, false,
                             cav, work, lcaverr, handler, true, 0.,
                             ithrd1, ithrd2);

          if(ierro <= 0){
            result.nSuccessInsertLength++;
            msh.poicstr[cav.ipins] = false;
            const bool touches_boundary = cav.lcedg.get_n() > 0;
            if(touches_boundary) result.nSuccessInsertLengthBoundary++;
            else                 result.nSuccessInsertLengthInterior++;
            if(!require_boundary || touches_boundary){
              result.found = true;
              result.ientt = ientt;
              result.ied = ied;
              result.ip1 = ent2poi(ientt, lnoed(ied,0));
              result.ip2 = ent2poi(ientt, lnoed(ied,1));
              result.qentt = quaent;
              result.len = cand.len;
              result.ierro_quality = ierro_quality;
              result.ierro_length = ierro;
              result.ipins = cav.ipins;
              result.ncedg = cav.lcedg.get_n();
              result.ncfac = cav.lcfac.get_n();
              result.nctet = cav.lctet.get_n();
              if(save_pre_state){
                result.pre_mesh = pre_base.string() + ".meshb";
                result.pre_met = pre_base.string() + ".solb";
              }
              return result;
            }
            didOperation = true;
            break;
          }

        }else{

          INCVDEPTH(msh.param);

          result.ntryCollapse++;

          int ierro = collapseEdge<MFT>(msh, tdim, ientt, ied, 0,
                                        cav, work, lcaverr, handler,
                                        ithrd1, ithrd2, ithrd3);

          if(ierro == 0){
            didOperation = true;
            result.nSuccessCollapse++;
            break;
          }
        }
      }

      if(didOperation){
        handler.updateK(ientt, ent2ent, ent2tag,
                        msh.tag[ithrd1] + 1, ithrd1);
        break;
      }

      // Keep this replay synchronized with adaptMeshQuality0.
      const double quaTryThreshold = 0.01;
      if(msh.param->adp_quality_smoothing
          && quaent < quaTryThreshold){

        result.ntrySmoothing++;

        #ifdef STEPDISTANCE
        double statSmoothing = smoothElement_Ball<MFT>(msh,ientt,handler,
                                                       QuaFun::StepDistance,
                                                       ithrd1,ithrd2);
        #else
        double statSmoothing = smoothElement_Ball<MFT>(msh,ientt,handler,
                                                       QuaFun::SizeShape,
                                                       ithrd1,ithrd2);
        #endif

        if(statSmoothing > 0){
          result.nSuccessSmoothing++;
          didOperation = true;
          handler.updateK(ientt,ent2ent,ent2tag,
                          msh.tag[ithrd1]+1,ithrd1);
          break;
        }
      }
    }

    if(!didOperation) break;
  }

  return result;
}

template<class MFT>
StatefulAdaptResult dispatch_find_stateful_case(const std::string& cmd,
                                                bool save_pre_state = false){
  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  MeshBase& msh = *run.msh_g;

  StatefulAdaptResult found;
  CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      found = find_stateful_quality_adapt_divergence<MFT,gdim,ideg>(
        cmd, save_pre_state);
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);

  return found;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_exact_case_common_integration_schemes)
{
  const std::string cmd = objective_compare_input_command();
  cargHandler arg(cmd);
  MetrisRunner run(arg.c,arg.v);

  CavityIntegrationValues values;
  if(run.metricFE){
    values = dispatch_compare_exact_case_integration_schemes<MetricFieldFE>(cmd);
  }else{
    values =
        dispatch_compare_exact_case_integration_schemes<MetricFieldAnalytical>(cmd);
  }

  BOOST_TEST(values.reconnected.size_shape_one_point
             < values.original.size_shape_one_point);
  BOOST_TEST(values.reconnected.step_distance_one_point
             > values.original.step_distance_one_point);
  BOOST_TEST(values.reconnected.size_shape_vertex_barycenter
             < values.original.size_shape_vertex_barycenter);
  BOOST_TEST(values.reconnected.step_distance_vertex_barycenter
             > values.original.step_distance_vertex_barycenter);
}

BOOST_AUTO_TEST_CASE(test_sizeshape_stepdistance_initial_insert_scan)
{
#ifndef TESTQUALITYALGO
  BOOST_FAIL("test_sizeshape_stepdistance_initial_insert_scan requires TESTQUALITYALGO");
#else
  const std::string cmd = objective_compare_input_command();
  cargHandler arg(cmd);
  MetrisRunner run(arg.c,arg.v);

  ObjectiveInsertSummary summary;
  if(run.metricFE){
    summary = dispatch_compare_objective_insertions<MetricFieldFE>(cmd);
  }else{
    summary = dispatch_compare_objective_insertions<MetricFieldAnalytical>(cmd);
  }

  BOOST_TEST(summary.n_candidates ==
             summary.n_both + summary.n_size_shape_only
             + summary.n_step_distance_only + summary.n_neither);
#endif
}

BOOST_AUTO_TEST_CASE(test_stepdistance_p1_full_initial_scan)
{
#ifndef TESTQUALITYALGO
  BOOST_FAIL("test_stepdistance_p1_full_initial_scan requires TESTQUALITYALGO");
#else
  const std::string cmd = objective_compare_input_command();
  cargHandler arg(cmd);
  MetrisRunner run(arg.c,arg.v);

  CommonPScanSummary summary;
  if(run.metricFE){
    summary = dispatch_common_vb_objective_scan<MetricFieldFE>(cmd);
  }else{
    summary = dispatch_common_vb_objective_scan<MetricFieldAnalytical>(cmd);
  }

  BOOST_TEST(summary.n_candidates ==
             summary.n_valid + summary.n_setup_failures);
  BOOST_TEST(summary.n_valid ==
             summary.n_size_step_p1_both
             + summary.n_size_only_vs_step_p1
             + summary.n_step_p1_only_vs_size
             + summary.n_size_step_p1_neither);
#endif
}

BOOST_AUTO_TEST_CASE(test_quality_insert_divergence)
{
#ifndef TESTQUALITYALGO
  BOOST_FAIL("test_quality_insert_divergence requires TESTQUALITYALGO");
#else
  const std::string cmd = input_command();

  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);

  DivergenceCase found;
  if(run.metricFE){
    found = dispatch_find_case<MetricFieldFE>(cmd);
  }else{
    found = dispatch_find_case<MetricFieldAnalytical>(cmd);
  }

  if(found.found){
    fmt::print("\n-- Found quality/length insertion divergence\n");
    fmt::print("   case dir        : {}\n", case_dir());
    fmt::print("   tdim            : {}\n", found.tdim);
    fmt::print("   ientt           : {}\n", found.ientt);
    fmt::print("   ied             : {}\n", found.ied);
    fmt::print("   endpoints       : {} {}\n", found.ip1, found.ip2);
    fmt::print("   entity quality  : {:.16e}\n", found.qentt);
    fmt::print("   edge length     : {:.16e}\n", found.len);
    fmt::print("   quality ierro   : {}\n", found.ierro_quality);
    fmt::print("   length ierro    : {}\n", found.ierro_length);
    fmt::print("   quality cavity  : ipins={} nedge={} nface={} ntet={}\n",
               found.quality.ipins, found.quality.ncedg,
               found.quality.ncfac, found.quality.nctet);
    fmt::print("   length cavity   : ipins={} nedge={} nface={} ntet={}\n",
               found.length.ipins, found.length.ncedg,
               found.length.ncfac, found.length.nctet);
  }

  #ifdef CAVSMOOTHING
  if(!found.found){
    fmt::print("\n-- No case found where compiled smoothed quality insertion fails "
               "and length insertion succeeds with worsenPctg = 0\n");
    fmt::print("   case dir: {}\n", case_dir());
  }
  BOOST_TEST(true);
  #else
  BOOST_REQUIRE_MESSAGE(found.found,
                        "No long-edge case found where quality insertion fails "
                        "and length insertion succeeds with worsenPctg = 0");
  #endif
#endif
}

BOOST_AUTO_TEST_CASE(test_quality_smoothed_insert_divergence)
{
#ifndef TESTQUALITYALGO
  BOOST_FAIL("test_quality_smoothed_insert_divergence requires TESTQUALITYALGO");
#else
  const std::string cmd = input_command();

  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);

  DivergenceCase found;
  if(run.metricFE){
    found = dispatch_find_case<MetricFieldFE>(cmd, true, 0.);
  }else{
    found = dispatch_find_case<MetricFieldAnalytical>(cmd, true, 0.);
  }

  if(found.found){
    fmt::print("\n-- Found smoothed-quality/length insertion divergence\n");
    fmt::print("   case dir        : {}\n", case_dir());
    fmt::print("   tdim            : {}\n", found.tdim);
    fmt::print("   ientt           : {}\n", found.ientt);
    fmt::print("   ied             : {}\n", found.ied);
    fmt::print("   endpoints       : {} {}\n", found.ip1, found.ip2);
    fmt::print("   entity quality  : {:.16e}\n", found.qentt);
    fmt::print("   edge length     : {:.16e}\n", found.len);
    fmt::print("   smoothed quality ierro : {}\n", found.ierro_quality);
    fmt::print("   length ierro           : {}\n", found.ierro_length);
    fmt::print("   length worsenPctg      : 0\n");
    fmt::print("   smoothed quality cavity: ipins={} nedge={} nface={} ntet={}\n",
               found.quality.ipins, found.quality.ncedg,
               found.quality.ncfac, found.quality.nctet);
    fmt::print("   length cavity          : ipins={} nedge={} nface={} ntet={}\n",
               found.length.ipins, found.length.ncedg,
               found.length.ncfac, found.length.nctet);
  }else{
    fmt::print("\n-- No case found where smoothed quality insertion fails "
               "and length insertion succeeds with worsenPctg = 0\n");
    fmt::print("   case dir: {}\n", case_dir());
  }

  BOOST_TEST(true);
#endif
}

BOOST_AUTO_TEST_CASE(test_stateful_quality_adapt_divergence)
{
#ifndef TESTQUALITYALGO
  BOOST_FAIL("test_stateful_quality_adapt_divergence requires TESTQUALITYALGO");
#else
  const std::string cmd = input_command();

  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);

  StatefulAdaptResult found;
  if(run.metricFE){
    found = dispatch_find_stateful_case<MetricFieldFE>(cmd);
  }else{
    found = dispatch_find_stateful_case<MetricFieldAnalytical>(cmd);
  }

  fmt::print("\n-- Stateful quality adapt diagnostic\n");
  fmt::print("   case dir              : {}\n", case_dir());
  fmt::print("   iter                  : {}\n", found.iter);
  fmt::print("   ntrySmoothing         : {}\n", found.ntrySmoothing);
  fmt::print("   nSuccessSmoothing     : {}\n", found.nSuccessSmoothing);
  fmt::print("   ntryInsert            : {}\n", found.ntryInsert);
  fmt::print("   nSuccessInsert        : {}\n", found.nSuccessInsert);
  fmt::print("   ntryInsertLength      : {}\n", found.ntryInsertLength);
  fmt::print("   nSuccessInsertLength  : {}\n", found.nSuccessInsertLength);
  fmt::print("     interior            : {}\n",
             found.nSuccessInsertLengthInterior);
  fmt::print("     boundary            : {}\n",
             found.nSuccessInsertLengthBoundary);
  fmt::print("   ntryCollapse          : {}\n", found.ntryCollapse);
  fmt::print("   nSuccessCollapse      : {}\n", found.nSuccessCollapse);

  if(found.found){
    fmt::print("\n-- Found stateful quality/length insertion divergence\n");
    fmt::print("   ientt                 : {}\n", found.ientt);
    fmt::print("   ied                   : {}\n", found.ied);
    fmt::print("   endpoints             : {} {}\n", found.ip1, found.ip2);
    fmt::print("   entity quality        : {:.16e}\n", found.qentt);
    fmt::print("   edge length           : {:.16e}\n", found.len);
    fmt::print("   quality ierro         : {}\n", found.ierro_quality);
    fmt::print("   length ierro          : {}\n", found.ierro_length);
    fmt::print("   length cavity         : ipins={} nedge={} nface={} ntet={}\n",
               found.ipins, found.ncedg, found.ncfac, found.nctet);
  }else{
    fmt::print("\n-- No stateful case found before termination/max_iter\n");
    fmt::print("   max_iter              : {}\n", stateful_adapt_max_iter());
  }

  BOOST_TEST(true);
#endif
}

BOOST_AUTO_TEST_CASE(test_production_quality_adapt_stats)
{
#ifndef TESTQUALITYALGO
  BOOST_FAIL("test_production_quality_adapt_stats requires TESTQUALITYALGO");
#else
  const std::string cmd = input_command();

  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  run.adaptMesh2();

  BOOST_TEST(true);
#endif
}

BOOST_AUTO_TEST_CASE(test_production_runmetris_stats)
{
#ifndef TESTQUALITYALGO
  BOOST_FAIL("test_production_runmetris_stats requires TESTQUALITYALGO");
#else
  const std::string cmd = production_run_command();

  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);
  run.runMetris();

  const std::filesystem::path stats =
    runmetris_dir() / ("outputAdaptStats_MOESS_a"
                     + std::to_string(adaptive_iteration())
                     + ".txt");

  fmt::print("\n-- Production runMetris stats file: {}\n", stats.string());
  std::ifstream fin(stats);
  BOOST_REQUIRE_MESSAGE(fin.good(), "Could not open generated stats file");

  std::string line;
  int iline = 0;
  while(iline < 6 && std::getline(fin,line)){
    fmt::print("{}\n", line);
    iline++;
  }

  BOOST_TEST(true);
#endif
}

BOOST_AUTO_TEST_CASE(trace_stateful_quality_adapt_divergence_case)
{
#ifndef TESTQUALITYALGO
  BOOST_FAIL("trace_stateful_quality_adapt_divergence_case requires TESTQUALITYALGO");
#else
  const std::string cmd = input_command();

  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);

  StatefulAdaptResult found;
  if(run.metricFE){
    found = dispatch_find_stateful_case<MetricFieldFE>(cmd, true);
  }else{
    found = dispatch_find_stateful_case<MetricFieldAnalytical>(cmd, true);
  }

  if(!found.found){
    fmt::print("\n-- No stateful quality/length divergence found\n");
    BOOST_TEST(true);
    return;
  }

  const std::string cad = cad_file(case_dir());
  const std::string replay_cmd =
    input_command_from_files(found.pre_mesh, found.pre_met, cad);

  if(run.metricFE){
    dispatch_trace_exact_case<MetricFieldFE>(replay_cmd, found.ientt, found.ied);
  }else{
    dispatch_trace_exact_case<MetricFieldAnalytical>(replay_cmd,
                                                    found.ientt,
                                                    found.ied);
  }

  fmt::print("\n-- Traced stateful quality/length insertion divergence\n");
  fmt::print("   evolved mesh   : {}\n", found.pre_mesh);
  fmt::print("   evolved metric : {}\n", found.pre_met);
  fmt::print("   ientt          : {}\n", found.ientt);
  fmt::print("   ied            : {}\n", found.ied);
  fmt::print("   trace dir      : {}\n", trace_dir().string());
  fmt::print("   vizir dir      : {}\n", vizir_dir().string());
#endif
}

BOOST_AUTO_TEST_CASE(trace_quality_insert_divergence_case)
{
#ifndef TESTQUALITYALGO
  BOOST_FAIL("trace_quality_insert_divergence_case requires TESTQUALITYALGO");
#else
  const std::string cmd = input_command();

  cargHandler arg(cmd);
  MetrisRunner run(arg.c, arg.v);

  if(run.metricFE){
    dispatch_trace_exact_case<MetricFieldFE>(cmd);
  }else{
    dispatch_trace_exact_case<MetricFieldAnalytical>(cmd);
  }

  fmt::print("-- Wrote cavity trace to {}\n", trace_dir().string());
#endif
}
