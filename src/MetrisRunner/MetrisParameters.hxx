//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_METRIS_PARAMETERS__
#define __METRIS_METRIS_PARAMETERS__

#include "../Mesh/MeshFwd.hxx"
#include "../msh_anamet.hxx"
#include "../SolutionField/anasol.hxx"
#include "../metris_defaults.hxx"
#include "../metris_constants.hxx"

#include "nlohmann/json_fwd.hpp"
#include <string>
#include <cstdio>

#include "../inc_hana.hxx"

namespace Metris{

struct MetrisOptions;
enum class FEBasis;
class MetrisRunner;
class MetrisAPI;
class CADInfo;

typedef void(*anamet_proto)(const AnaMetCtx*,const double*__restrict__,double,int,double*,double*);
typedef double(*anasol_proto)(void*,const double*__restrict__,std::initializer_list<double*>);


struct MetrisParametersData;
struct MetrisParameters;

void to_json(nlohmann::json& jj, const MetrisParametersData& param);
void to_json(nlohmann::json& jj, const MetrisParameters& param);
void from_json(const nlohmann::json& jj, MetrisParametersData& param);
void from_json(const nlohmann::json& jj, MetrisParameters& param);


/*------------------------------------
2017 C++er's class reflection
Add new entries to
   METRIS_PARAMETER_OTHERFIELDS_JSON (most likely)
or METRIS_PARAMETER_OTHERFIELDS_NOJSON
depending on whether they should be serialized.
------------------------------------*/

#define METRIS_PARAMETER_FIELD1 \
    /* target degree -tardeg */ \
    FIELD(int, usrTarDeg, 1) \

// Add to this:
#define METRIS_PARAMETER_OTHERFIELDS_JSON \
    /* number of threads -nproc */ \
    FIELD(int, nproc, -1) \
    /* Scaled ccoef lower bound */ \
    FIELD(double, jtol, 1.0e-6) \
    /* Measure tolerance (similar to height) */ \
    FIELD(double, vtol, Defaults::vtol) \
    /* For adaptGeoLines: edges in target [tarlen / tolfac ; tarlen * tolfac] */ \
    /* tarlen is determined in terms of the actual curve length, to be closest */ \
    /* to 1, while having integer number of vertices on the curve. */ \
    FIELD(double, geo_lentolfac, 1.01) \
    /* Absolute tolerance for point to be on CAD edge */ \
    FIELD(double, geo_abstoledg, 1.0e-10) \
    /* ----------------- Adaptation options */ \
    FIELD(int, adp_opt_niter, 1) /* do smoothing/swapping in adaptation loop (expensive) */ \
    FIELD(int, adp_niter, 0) \
    FIELD(double, adp_unit_stop, 99.9) /* threshold to consider mesh unit and stop everything 0 - 100 */ \
    FIELD(bool, adp_line_adapt, false) \
    FIELD(double, adp_stagn_stop, 1.0e-2) /* stat threshold for stagnation (default 1e-3) */ \
    FIELD(bool, adp_smoo_len, false) /* use len-based smoothing in adaptation loop */ \
    /* Metric min/max size control */ \
    FIELD(double, hmin, 1.0e-30) \
    FIELD(double, hmax, 1.0e30) \
    FIELD(double, met_snap_tol, 1.0e-3) /* tolerance for surface snapping in intrinsic metric case */ \
    FIELD(double, anamet_dx, 0.0) \
    FIELD(double, anamet_dy, 0.0) \
    FIELD(double, anamet_dz, 0.0) \
    /* options: "curve" */ \
    /* Defaults to 0 (no curve), 1 for offsets followed by ccoef max, 2 for */ \
    /* metric-based LP, 3 for offsets followed by smoothing, 4 offsets and backtrack */ \
    /* 5 interpolation error minimization */ \
    FIELD(int, curveType, 0) \
    /* Type 0 is metric-based */ \
    /* Type 1 is interpolation error, requires -anasol be passed, possibly also */ \
    /* intp_pdeg and intp_pnorm. */ \
    FIELD(int, smoo_type, 0) \
    /* Generic integer flags */ \
    FIELD(int, iflag1, 0) \
    FIELD(int, iflag2, 0) \
    FIELD(int, iflag3, 0) \
    FIELD(double, rflag1, 0.0) \
    FIELD(double, rflag2, 0.0) \
    FIELD(double, rflag3, 0.0) \
    /* ----------------- Optimization options */ \
    /* -- Smoothing */ \
    FIELD(int, opt_niter, 5) \
    FIELD(int, opt_pnorm, 2) \
    FIELD(int, opt_power, -1) \
    FIELD(int, opt_smoo_niter, 10) \
    FIELD(double, opt_smoo_tol, 0.005) \
    /* Surface qualities weight the main term by qua_surf_w_quality */ \
    /* and the normal deviation term by qua_surf_w_normal */ \
    FIELD(double, qua_surf_wt_normal, 1) \
    FIELD(double, qua_surf_wt_quality, 1) \
    /* -- Quality (quafun_unit) */ \
    /* compute coef_det (det - 1)^powr_det + idem(tra) */ \
    FIELD(double, opt_coef_det, 1.0) \
    FIELD(double, opt_coef_tra, 1.0) \
    FIELD(int, opt_powr_det, -2) \
    FIELD(int, opt_powr_tra, 2) \
    FIELD(bool, opt_unif, false) /* experimental */ \
    /* -- Swaps */ \
    /* Maximum global swap iterations. Stops before if nothing to do */ \
    FIELD(int, opt_swap_niter, 100) \
    /* lp norm to improve: 0 is infinity (max). Governed by option -qswap-norm */ \
    /* but only in the optimization module (not adaptation) */ \
    /* Note quality norm is defined by -opt-power and such options */ \
    FIELD(int, opt_swap_pnorm, 2) \
    /* Minimum increase in quality function to go through with a swap */ \
    /* Governed by option -qswap-thres in optimization module (not adaptation) */ \
    FIELD(double, opt_swap_thres, 1.0e-16) \
    /* Whether expensive tet swaps should be done (edge -> faces) */ \
    /* These are rarely done in practice and take by far the most time, */ \
    /* seemingly 85% of total tet swapping time */ \
    FIELD(bool, opt_swap_tet_expensive, false) \
    FIELD(int, interp_err_min_algo, 1) \
    /* Metric scaling */ \
    FIELD(double, metScale, 1.0) \
    /* To use defaults, see anamet_2D and anamet_3D.cxx */ \
    FIELD(int, ianamet, -1) \
    /* See anasol.hxx. Can implement your own with same prototype */ \
    /* To use defaults, see anasol_2D and anasol_3D.cxx */ \
    FIELD(int, ianasol, -1) \
    FIELD(int, intp_pdeg, 1) /* interpolation degree */ \
    FIELD(int, intp_pnorm, 1) /* interp error norm L2 or L1 */ \
    /* ----------------- Normal deviation control */ \
    FIELD(double, nordev_tol, 0) /* must be in [0,1]: 0 means tol is that of the current cavity, and from there the tol is increased towards 1 */ \
    FIELD(double, nordev_max, 0.5) \

#define METRIS_PARAMETERS_OTHERFIELDS_NOJSON \
    /* See anamet.hxx. Can implement your own with same prototype */ \
    FIELD(bool, anaMet, false) \
    FIELD(bool, scaleMet, false) /* internal */ \
    FIELD(bool, anaSol, false) /* internal */ \
    /* Verbosity level */ \
    FIELD(int, iverb, 1) \
    /* Verbosity depth */ \
    FIELD(int, ivdepth, 1) \
    /* Full debugs (costly), wait at certain errors (debug only) */ \
    FIELD(bool, dbgfull, false) \
    FIELD(bool, interactive, false) \
    FIELD(bool, nocleanup, false) /* disable cleanup routines, easier index tracking */ \
    FIELD(FEBasis, outbasis, FEBasis::Lagrange) \
    FIELD(bool, refineConventionsInp, false) \
    FIELD(bool, refineConventionsOut, false) \


#define METRIS_PARAMETER_OTHERFIELDS \
    METRIS_PARAMETER_OTHERFIELDS_JSON \
    METRIS_PARAMETERS_OTHERFIELDS_NOJSON

// All parameter fields
#define METRIS_PARAMETER_FIELDS \
    METRIS_PARAMETER_FIELD1 \
    METRIS_PARAMETER_OTHERFIELDS \

// All parameter fields that belong in a json
#define METRIS_PARAMETER_FIELDS_JSON \
    METRIS_PARAMETER_FIELD1 \
    METRIS_PARAMETER_OTHERFIELDS_JSON \

namespace internal {
  template<typename T>
  struct FieldNameType{
    // string_view is basically a reference to a string
    // that can do non-const operations
    std::string_view name;
    using type = T;
  };

  //// Helper to get array of all field infos
  //template<typename... Args>
  //struct FieldInfoList {
  //  // sizeof...(Args) is the number of template parameters
  //  // &Args... expands to the address of each parameter
  //  // We store everything as void* since they may be FieldNameType<int>, FieldNameType<double>, etc.
  //  static constexpr std::array<const void*, sizeof...(Args)> fields = {&Args...};
  //  static constexpr size_t size = sizeof...(Args);
  //};
  // Helper to create field info
  template<typename T>
  constexpr auto make_field_info(const char* name){
    return hana::make_tuple(hana::type_c<T>, name);
  }

}

// Aggregate type, plays well with Boost:pfr
// This holds all the numeric parameters.
// Strings could also be included but:
//  - We want to compare them differently
//  - I don't expect there'll be any more in the future (or much fewer)
//  - It's good to have private setters for them
struct MetrisParametersData{
public:

  void setAnalyticalMetric(int ianamet);
  void setAnalyticalMetric(AnaMetFun anamet_ptr);
  void setAnalyticalSolution(int ianasol);
  void setAnalyticalSolution(AnaSolFun anasol_ptr);
  void setMetricScale(double sclmet);

  anamet_proto anamet_ptr = NULL;
  anasol_proto anasol_ptr = NULL;

  bool operator==(const MetrisParametersData &other) const;
  void printDifference(const MetrisParametersData &other, std::string thisName = "", FILE* logFile = stdout) const;

  #define FIELD(type, name, default_value) \
      type name{default_value};\
      static constexpr auto name##_info = internal::make_field_info<type>(#name);

  METRIS_PARAMETER_FIELDS
  #undef FIELD


  // Collection of all field infos
  static constexpr auto all_fields = hana::make_tuple(
    #define FIELD(type, name, default_value) \
      name##_info
    METRIS_PARAMETER_FIELD1
    #undef FIELD
    #define FIELD(type, name, default_value) \
      ,name##_info
    METRIS_PARAMETER_OTHERFIELDS
    #undef FIELD
  );

};

// Parameters can be set manually, or initialized by a Runner from argc/argv
// If set manually, use MetrisRunner constructor taking a MetrisParameter as input.
struct MetrisParameters : public MetrisParametersData {
public:
  friend void from_json(const nlohmann::json& jj, MetrisParameters& param);

  MetrisParameters();

  MetrisParameters(MetrisOptions &opt);

  ~MetrisParameters(){
    if(logFile_ && logFile_ != stdout && logFile_ != stderr){
      fclose(logFile_);
    }
  }

  bool operator==(const MetrisParameters &other) const;

  void checkParameters();

  void setMeshIn(std::string inp);
  void setMeshOut(std::string out);
  void setBackIn(std:: string inp);
  void setCAD(std::string cadName);
  void setMetricFile(std::string metName);
  void setLogFile(std::string fname);


  void printDifference(const MetrisParameters &other, std::string thisName = "") const;

public:
  const std::string& outmFileName;
  std::string outmPrefix;
  const std::string& cadFileName;
  const std::string& backFileName;
  const std::string& metFileName;
  const std::string& logFileName;
  const std::string& meshFileName;
  FILE* const & logFile;

private:
  std::string outmFileName_;
public:
  bool main_in_prefix{false}; // whether main output goes in prefix (default no)
private:

  bool wrtMesh{false};
  std::string meshFileName_;

  std::string cadFileName_;
  std::string backFileName_;
  std::string metFileName_;

  FILE* logFile_{stdout};
  std::string logFileName_;

  bool inpMet{false};
  bool inpBack{false};
  bool inpCAD{false};
  bool inpMesh{false};

  friend class MetrisRunner;
  friend class MetrisAPI;
  friend class MeshBase;
  friend class MeshBack;
  template<class MFT>
  friend class Mesh;
  friend class CADInfo;

};

} // End namespace
#endif
