//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../metris_options.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../metris_defaults.hxx"
#include "../metris_constants.hxx"
#include "../io_libmeshb.hxx"
#include "../utils/mprintf.hxx"

#include "nlohmann/json.hpp"
#include <string>

namespace Metris{

bool MetrisParametersData::operator==(const MetrisParametersData &other) const{
  #define FIELD(type, name, default_value) \
  if(name != other.name) return false;
  METRIS_PARAMETER_FIELDS_JSON
  #undef FIELD
  return true;
}

bool MetrisParameters::operator==(const MetrisParameters &other) const{
  if(!MetrisParametersData::operator==(other)) return false;

  // Check filename component of file names match, ignore path.
  std::string inpCAD = this->cadFileName.substr(this->cadFileName.find_last_of("/\\") + 1);
  std::string othCAD = other.cadFileName.substr(other.cadFileName.find_last_of("/\\") + 1);
  if(inpCAD != othCAD) return false;

  std::string inpBack = this->backFileName.substr(this->backFileName.find_last_of("/\\") + 1);
  std::string othBack = other.backFileName.substr(other.backFileName.find_last_of("/\\") + 1);
  if(inpBack != othBack) return false;

  std::string inpMesh = this->meshFileName.substr(this->meshFileName.find_last_of("/\\") + 1);
  std::string othMesh = other.meshFileName.substr(other.meshFileName.find_last_of("/\\") + 1);
  if(inpMesh != othMesh) return false;

  std::string inpMet = this->metFileName.substr(this->metFileName.find_last_of("/\\") + 1);
  std::string othMet = other.metFileName.substr(other.metFileName.find_last_of("/\\") + 1);
  if(inpMet != othMet) return false;

  return true;
}

void MetrisParametersData::printDifference(const MetrisParametersData &other, std::string thisName, FILE* logFile) const{

  #define FIELD(type, name, default_value) \
  if(name != other.name){\
    if(!thisName.empty()) fmt::print(logFile,"-- {} differs from {}: {} -> {}\n", #name, thisName, this->name, other.name);\
    else fmt::print(logFile,"-- {} differs: {} -> {}\n",  #name, this->name, other.name);\
  }
  METRIS_PARAMETER_FIELDS_JSON
  #undef FIELD
}

void MetrisParameters::printDifference(const MetrisParameters &other, std::string thisName) const {

  MetrisParametersData::printDifference(other, thisName, logFile_);
  // Check filename component of file names match, ignore path.
  std::string inpCAD = this->cadFileName.substr(this->cadFileName.find_last_of("/\\") + 1);
  std::string othCAD = other.cadFileName.substr(other.cadFileName.find_last_of("/\\") + 1);
  if(inpCAD != othCAD) fmt::print("-- CAD differs: {} != {}\n", inpCAD, othCAD);

  std::string inpBack = this->backFileName.substr(this->backFileName.find_last_of("/\\") + 1);
  std::string othBack = other.backFileName.substr(other.backFileName.find_last_of("/\\") + 1);
  if(inpBack != othBack) fmt::print("-- Back Mesh differs: {} != {}\n", inpBack, othBack);

  std::string inpMesh = this->meshFileName.substr(this->meshFileName.find_last_of("/\\") + 1);
  std::string othMesh = other.meshFileName.substr(other.meshFileName.find_last_of("/\\") + 1);
  if(inpMesh != othMesh) fmt::print("-- Input Mesh differs: {} != {}\n", inpMesh, othMesh);

  std::string inpMet = this->metFileName.substr(this->metFileName.find_last_of("/\\") + 1);
  std::string othMet = other.metFileName.substr(other.metFileName.find_last_of("/\\") + 1);
  if(inpMet != othMet) fmt::print("-- Input Metric differs: {} != {}\n", inpMet, othMet);

}

void to_json(nlohmann::json& jj, const MetrisParametersData& param) {

  jj = nlohmann::json{
    #define FIELD(type, name, default_value) {#name, param.name}
    METRIS_PARAMETER_FIELD1
    #undef FIELD
    #define FIELD(type, name, default_value) ,{#name, param.name}
    METRIS_PARAMETER_OTHERFIELDS_JSON
    #undef FIELD
  };

}

void from_json(const nlohmann::json& jj, MetrisParametersData& param) {

  #define FIELD(type, name, default_value) \
  if(jj.contains(#name)) jj.at(#name).get_to(param.name);
  METRIS_PARAMETER_FIELDS_JSON
  #undef FIELD

}

void to_json(nlohmann::json& jj, const MetrisParameters& param) {
  jj = nlohmann::json{{"cadFileName", param.cadFileName}
                     ,{"backFileName", param.backFileName}
                     ,{"metFileName", param.metFileName}
                     ,{"meshFileName", param.meshFileName}
                     };
  nlohmann::json jj2;
  to_json(jj2, static_cast<const MetrisParametersData&>(param));
  jj.update(jj2); // merge
}

void from_json(const nlohmann::json& jj, MetrisParameters& param) {
  from_json(jj, static_cast<MetrisParametersData&>(param));
  jj.at("cadFileName").get_to(param.cadFileName_);
  jj.at("backFileName").get_to(param.backFileName_);
  jj.at("metFileName").get_to(param.metFileName_);
  jj.at("meshFileName").get_to(param.meshFileName_);
}



MetrisParameters::MetrisParameters() :
  outmFileName(outmFileName_),
  cadFileName(cadFileName_), backFileName(backFileName_),
  metFileName(metFileName_), logFileName(logFileName_),
  meshFileName(meshFileName_), logFile(logFile_)
{
}

MetrisParameters::MetrisParameters(const MetrisParameters &other) :
  MetrisParametersData(other),
  // Public references → bind to THIS object's private strings
  outmFileName(outmFileName_),
  outmPrefix(other.outmPrefix),
  cadFileName(cadFileName_),
  backFileName(backFileName_),
  metFileName(metFileName_),
  logFileName(logFileName_),
  meshFileName(meshFileName_),
  logFile(logFile_),
  // Private data — declaration order
  outmFileName_(other.outmFileName_),
  main_in_prefix(other.main_in_prefix),
  wrtMesh(other.wrtMesh),
  meshFileName_(other.meshFileName_),
  cadFileName_(other.cadFileName_),
  backFileName_(other.backFileName_),
  metFileName_(other.metFileName_),
  logFile_(other.logFile_),
  logFileOwner_(false),       // copy does NOT own the FILE*
  logFileName_(other.logFileName_),
  inpMet(other.inpMet),
  inpBack(other.inpBack),
  inpCAD(other.inpCAD),
  inpMesh(other.inpMesh)
{
}

MetrisParameters::MetrisParameters(MetrisOptions &opt) : MetrisParameters(){

  if(opt.count("help")) {
    fmt::print("Flag --help:\n");
    std::cout << opt.s << std::endl;
    exit(0);
  }

  // These need to be first, or subsequent prints will be ill-defined.
  if(opt.count("verb")){
    iverb = opt.m["verb"].template as<int>();
  }
  if(opt.count("vdepth")){
    ivdepth = opt.m["vdepth"].template as<int>();
  }
  if(opt.count("log")){
    setLogFile(opt.m["log"].template as<std::string>());
  }

  GETVDEPTH(this);

  if(opt.count("refine-conventions-inp")){
    refineConventionsInp = true;
  }
  if(opt.count("refine-conventions-out")){
    refineConventionsOut = true;
  }


  if(opt.count("opt-unif")){
    opt_unif = true;
    CPRINTF1("-- Set opt-unif\n");
  }

  if(opt.count("in")){
    setMeshIn(opt.m["in"].template as<std::string>());
    CPRINTF1("-- Read input mesh name {}\n", meshFileName.c_str());
  }

  if(opt.count("prefix")){
    outmPrefix = opt.m["prefix"].template as<std::string>();
    CPRINTF1("-- File prefix: {}\n", outmPrefix.c_str());
  }
  if(opt.count("main-in-prefix")){
    main_in_prefix = true;
  }

  if(opt.count("bez")){
    outbasis = FEBasis::Bezier;
    CPRINTF1("-- Bézier output basis\n");
  }


  if(opt.count("out")) {
    setMeshOut(opt.m["out"].template as<std::string>());
    CPRINTF1("-- Read output file name {}.\n", outmFileName.c_str());
  }else{
    CPRINTF1("# Output mesh file name not set. Use --out or -o <filename>.\n");
    CPRINTF1("# Running but skipping mesh output.\n");
  }

  // usrMaxDeg is the very maximum the user is allowing for storage. It is hard bounded by the constant METRIS_MAX_DEG
  // usrTarDeg is the minimum degree the user wants.
  if(opt.count("tardeg")){
    usrTarDeg = opt.m["tardeg"].template as<int>();
  }

  if(opt.count("nproc")){
    nproc = opt.m["nproc"].template as<int>();
    CPRINTF1("-- Running with nproc = {} \n",nproc);
  }

  if(opt.count("cad")){
    setCAD(opt.m["cad"].template as<std::string>());
  }

  if(opt.count("dbgfull")){
    CPRINTF1("-- Full debugs activated\n");
    dbgfull = true;
  }
  if(opt.count("interactive")){
    CPRINTF1("-- Wait calls activated\n");
    interactive = true;
  }
  if(opt.count("nocleanup")){
    CPRINTF1("-- Cleanup calls deactivated\n");
    nocleanup = true;
  }

  if(opt.count("back")){
    setBackIn(opt.m["back"].template as<std::string>());
    CPRINTF1(" - Read back mesh name {}\n", backFileName.c_str());
  }

  if(opt.count("met")){
    setMetricFile(opt.m["met"].template as<std::string>());
  }

  if(opt.count("anamet")){
    setAnalyticalMetric(opt.m["anamet"].template as<int>());
    CPRINTF1("Using analytical metric {} \n", ianamet);
  }

  if(opt.count("anasol")){
    setAnalyticalSolution(opt.m["anasol"].template as<int>());
    CPRINTF1("Using analytical metric {} \n", ianamet);
  }

  if(opt.count("intp-pdeg")){
    intp_pdeg = opt.m["intp-pdeg"].template as<int>();
  }
  if(opt.count("intp-pnorm")){
    intp_pnorm = opt.m["intp-pnorm"].template as<int>();
  }

  if(opt.count("nordev-tol")){
    nordev_tol = opt.m["nordev-tol"].template as<double>();
  }
  if(opt.count("nordev-max")){
    nordev_max = opt.m["nordev-max"].template as<double>();
  }

  if(opt.count("progressiveAdapt")){
    progressiveAdapt = opt.m["progressiveAdapt"].template as<bool>();
  }

  if(opt.count("MOESS_adapt_it")){
    MOESS_adapt_it = opt.m["MOESS_adapt_it"].template as<int>();
  }

  if(opt.count("sclmet")){
    setMetricScale(opt.m["sclmet"].template as<double>());
  }

  if(opt.count("adapt")){
    adp_niter = opt.m["adapt"].template as<int>();
  }
  if(opt.count("adp-opt-niter")){
    adp_opt_niter = opt.m["adp-opt-niter"].as<int>();
  }
  if(opt.count("adp-unit-stop")){
    adp_unit_stop = opt.m["adp-unit-stop"].as<double>();
  }
  if(opt.count("adp-stagn-stop")){
    adp_stagn_stop = opt.m["adp-stagn-stop"].as<double>();
  }
  if(opt.count("adp-smoo-len")){
    adp_smoo_len = true;
  }

  if(opt.count("do-line-adp")){
    adp_line_adapt = true;
  }




  if(opt.count("curve")){
    curveType = opt.m["curve"].template as<int>();
    //METRIS_ENFORCE(curveType == 1 || curveType == 2 || curveType == 3);
  }

  if(opt.count("smoo-type")){
    smoo_type = opt.m["smoo-type"].template as<int>();
    //METRIS_ENFORCE(curveType == 1 || curveType == 2 || curveType == 3);
  }

  if(opt.count("jtol")){
    jtol = opt.m["jtol"].as<double>();
  }

  if(opt.count("vtol")){
    vtol = opt.m["vtol"].as<double>();
  }


  if(opt.count("geo-lentolfac")){
    geo_lentolfac = opt.m["geo-lentolfac"].as<double>();
  }

  if(opt.count("geo-abstoledg")){
    geo_abstoledg = opt.m["geo-abstoledg"].as<double>();
  }



  if(opt.count("hmin")){
    hmin = opt.m["hmin"].as<double>();
  }
  if(opt.count("hmax")){
    hmax = opt.m["hmax"].as<double>();
  }

  if(opt.count("mdx")){
    anamet_dx = opt.m["mdx"].as<double>();
  }
  if(opt.count("mdy")){
    anamet_dy = opt.m["mdy"].as<double>();
  }
  if(opt.count("mdz")){
    anamet_dz = opt.m["mdz"].as<double>();
  }


  if(opt.count("met-snap-tol")){
    met_snap_tol = opt.m["met-snap-tol"].as<double>();
  }

  if(opt.count("adp-quality-smoothing")){
    adp_quality_smoothing = true;
  }

  if(opt.count("opt-niter")){
    opt_niter = opt.m["opt-niter"].as<int>();
  }
  if(opt.count("opt-pnorm")){
    opt_pnorm = opt.m["opt-pnorm"].as<int>();
  }
  if(opt.count("opt-power")){
    opt_power = opt.m["opt-power"].as<int>();
  }
  if(opt.count("objective-p")){
    objective_p = opt.m["objective-p"].as<double>();
  }
  if(opt.count("objective-quadrature-order")){
    objective_quadrature_order
        = opt.m["objective-quadrature-order"].as<int>();
  }
  if(opt.count("step-distance-p")){
    const double compatibility_value
        = opt.m["step-distance-p"].as<double>();
    METRIS_ENFORCE_MSG(
        !opt.count("objective-p") || compatibility_value == objective_p,
        "Conflicting --objective-p ({}) and deprecated --step-distance-p "
        "({}) values",
        objective_p,compatibility_value);
    objective_p = compatibility_value;
  }
  if(opt.count("step-distance-regularization")){
    step_distance_regularization =
        opt.m["step-distance-regularization"].as<double>();
  }
  if(opt.count("step-distance-shape-volume")){
    step_distance_shape_volume = true;
  }
  if(opt.count("step-distance-cavity-target-average")){
    step_distance_cavity_target_average = true;
  }
  if(opt.count("step-distance-barrier-rho0")){
    step_distance_barrier_rho0 =
        opt.m["step-distance-barrier-rho0"].as<double>();
  }
  if(opt.count("step-distance-barrier-beta")){
    step_distance_barrier_beta =
        opt.m["step-distance-barrier-beta"].as<double>();
  }
  if(opt.count("opt-smoo-niter")){
    opt_smoo_niter = opt.m["opt-smoo-niter"].as<int>();
  }
  if(opt.count("opt-smoo-tol")){
    opt_smoo_tol = opt.m["opt-smoo-tol"].as<double>();
  }

  if(opt.count("opt-swap-pnorm")){
    opt_swap_pnorm = opt.m["opt-swap-pnorm"].as<int>();
  }
  if(opt.count("opt-swap-niter")){
    opt_swap_niter = opt.m["opt-swap-niter"].as<int>();
  }
  if(opt.count("opt-swap-tet-expensive")){
    opt_swap_tet_expensive = true;
  }


  if(opt.count("qua-surf-wt-quality")){
    qua_surf_wt_quality = opt.m["qua-surf-wt-quality"].as<double>();
  }

  if(opt.count("qua-surf-wt-normal")){
    qua_surf_wt_normal = opt.m["qua-surf-wt-normal"].as<double>();
  }

  if(opt.count("iflag1")){
    iflag1 = opt.m["iflag1"].as<int>();
  }
  if(opt.count("iflag2")){
    iflag2 = opt.m["iflag2"].as<int>();
  }
  if(opt.count("iflag3")){
    iflag3 = opt.m["iflag3"].as<int>();
  }

  if(opt.count("rflag1")){
    rflag1 = opt.m["rflag1"].as<double>();
  }
  if(opt.count("rflag2")){
    rflag2 = opt.m["rflag2"].as<double>();
  }
  if(opt.count("rflag3")){
    rflag3 = opt.m["rflag3"].as<double>();
  }

  if(opt.count("interp-err-min-algo")){
    interp_err_min_algo = opt.m["interp-err-min-algo"].as<int>();
    METRIS_ENFORCE_MSG(interp_err_min_algo == 0 || interp_err_min_algo == 1,
      "interp_err_min_algo has to be 0 (Newton) or 1 (DIRECT)");
  }
}

void MetrisParameters::checkParameters(){
  METRIS_ENFORCE_MSG(abs(opt_power) == 1, "opt_power can be set to -1 or 1");
  METRIS_ENFORCE_MSG(objective_p >= 1.0,
                     "objective_p must be greater than or equal to 1");
  METRIS_ENFORCE_MSG(
      objective_quadrature_order == -1
          || objective_quadrature_order == 0
          || (objective_quadrature_order >= 2
              && objective_quadrature_order <= 5),
      "objective-quadrature-order must be -1 (automatic), 0, 2, 3, 4, or 5");
  METRIS_ENFORCE_MSG(step_distance_regularization > 0.0,
                     "step_distance_regularization must be positive");
  METRIS_ENFORCE_MSG(
      !(step_distance_shape_volume && step_distance_cavity_target_average),
      "step-distance-shape-volume and "
      "step-distance-cavity-target-average are distinct variants");
  METRIS_ENFORCE_MSG(step_distance_barrier_rho0 >= 0.0,
                     "step_distance_barrier_rho0 must be nonnegative");
  METRIS_ENFORCE_MSG(step_distance_barrier_beta >= 0.0,
                     "step_distance_barrier_beta must be nonnegative");
  METRIS_ENFORCE_MSG(opt_smoo_tol >= 0.0 && opt_smoo_tol < 1.0,
                     "opt-smoo-tol must be in [0,1)");
  METRIS_ENFORCE(geo_abstoledg >= 0.0);
  METRIS_ENFORCE(geo_lentolfac >= 1.0);
  METRIS_ENFORCE_MSG(usrTarDeg >= 1, "Degree < 1 provided through tardeg.");
  METRIS_ENFORCE_MSG(usrTarDeg <= METRIS_MAX_DEG, "Opt -tardeg > METRIS_MAX_DEG = {}",METRIS_MAX_DEG);
}

void MetrisParameters::setMeshIn(std::string meshName){
  inpMesh = true;
  meshFileName_ = correctExtension_meshb(meshName);
}

void MetrisParameters::setBackIn(std::string backName){
  inpBack = true;
  backFileName_ = correctExtension_meshb(backName);
}

void MetrisParameters::setCAD(std::string CADname){
  inpCAD = true;
  cadFileName_ = correctExtension_egads(CADname);
}

void MetrisParameters::setMetricFile(std::string metName){
  inpMet = true;
  metFileName_ = correctExtension_solb(metName);
}

void MetrisParameters::setMeshOut(std::string out){
  wrtMesh = true;
  outmFileName_ = correctExtension_meshb(out);
}


void MetrisParametersData::setAnalyticalMetric(int ianamet){
  anaMet = true;
  this->ianamet = ianamet;
}

void MetrisParametersData::setAnalyticalMetric(AnaMetFun anamet_ptr){
  anaMet = true;
  this->anamet_ptr = (anamet_proto) anamet_ptr;
}

void MetrisParametersData::setAnalyticalSolution(int ianasol){
  anaSol = true;
  this->ianasol = ianasol;
}

void MetrisParametersData::setAnalyticalSolution(AnaSolFun anasol_ptr){
  anaSol = true;
  this->anasol_ptr = (anasol_proto) anasol_ptr;
}


void MetrisParametersData::setMetricScale(double sclmet){
  scaleMet = true;
  this->metScale = sclmet;
}


void MetrisParameters::setLogFile(std::string fname){
  if(logFileOwner_ && logFile_ && logFile_ != stdout && logFile_ != stderr){
    fclose(logFile_);
  }
  logFileName_ = fname;
  if(fname == "stdout"){
    logFile_ = stdout;
    logFileOwner_ = false;
  }else if(fname == "stderr"){
    logFile_ = stderr;
    logFileOwner_ = false;
  }else{
    logFile_ = fopen(fname.c_str(), "w");
    METRIS_ENFORCE_MSG(logFile_ != NULL,"Error opening log file {}", fname);
    logFileOwner_ = true;
  }
}

} // End namespace
