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


void to_json(nlohmann::json& jj, const MetrisParameters& param) {
  jj = nlohmann::json{{"cadFileName", param.cadFileName}
                     ,{"backFileName", param.backFileName}
                     ,{"metFileName", param.metFileName}
                     ,{"meshFileName", param.meshFileName}
                     ,{"usrTarDeg", param.usrTarDeg}
                     ,{"nproc", param.nproc}
                     ,{"jtol", param.jtol}
                     ,{"vtol", param.vtol}
                     ,{"geo_lentolfac", param.geo_lentolfac}
                     ,{"geo_abstoledg", param.geo_abstoledg}
                     ,{"adp_opt_niter", param.adp_opt_niter}
                     ,{"adp_niter", param.adp_niter}
                     ,{"adp_unit_stop", param.adp_unit_stop}
                     ,{"adp_line_adapt", param.adp_line_adapt}
                     ,{"adp_stagn_stop", param.adp_stagn_stop}
                     ,{"adp_smoo_len", param.adp_smoo_len}
                     ,{"hmin", param.hmin}
                     ,{"hmax", param.hmax}
                     ,{"met_snap_tol", param.met_snap_tol}
                     ,{"anamet_dx", param.anamet_dx}
                     ,{"anamet_dy", param.anamet_dy}
                     ,{"anamet_dz", param.anamet_dz}
                     ,{"curveType", param.curveType}
                     ,{"smoo_type", param.smoo_type}
                     ,{"nocleanup", param.nocleanup}
                     ,{"iflag1", param.iflag1}
                     ,{"iflag2", param.iflag2}
                     ,{"iflag3", param.iflag3}
                     ,{"rflag1", param.rflag1}
                     ,{"rflag2", param.rflag2}
                     ,{"rflag3", param.rflag3}
                     ,{"opt_niter", param.opt_niter}
                     ,{"opt_pnorm", param.opt_pnorm}
                     ,{"opt_power", param.opt_power}
                     ,{"opt_smoo_niter", param.opt_smoo_niter}
                     ,{"opt_smoo_tol", param.opt_smoo_tol}
                     ,{"qua_surf_wt_normal", param.qua_surf_wt_normal}
                     ,{"qua_surf_wt_quality", param.qua_surf_wt_quality}
                     ,{"opt_coef_det", param.opt_coef_det}
                     ,{"opt_coef_tra", param.opt_coef_tra}
                     ,{"opt_powr_det", param.opt_powr_det}
                     ,{"opt_powr_tra", param.opt_powr_tra}
                     ,{"opt_unif", param.opt_unif}
                     ,{"opt_swap_niter", param.opt_swap_niter}
                     ,{"opt_swap_pnorm", param.opt_swap_pnorm}
                     ,{"opt_swap_thres", param.opt_swap_thres}
                     ,{"opt_swap_tet_expensive", param.opt_swap_tet_expensive}
                     ,{"interp_err_min_algo", param.interp_err_min_algo}
                     ,{"metScale", param.metScale}
                     ,{"ianamet", param.ianamet}
                     ,{"intp_pdeg", param.intp_pdeg}
                     ,{"intp_pnorm", param.intp_pnorm}
                     };
}

void from_json(const nlohmann::json& jj, MetrisParameters& param) {
  jj.at("cadFileName").get_to(param.cadFileName);
  jj.at("backFileName").get_to(param.backFileName);
  jj.at("metFileName").get_to(param.metFileName);
  jj.at("meshFileName").get_to(param.meshFileName);
  jj.at("usrTarDeg").get_to(param.usrTarDeg);
  jj.at("nproc").get_to(param.nproc);
  jj.at("jtol").get_to(param.jtol);
  jj.at("vtol").get_to(param.vtol);
  jj.at("geo_lentolfac").get_to(param.geo_lentolfac);
  jj.at("geo_abstoledg").get_to(param.geo_abstoledg);
  jj.at("adp_opt_niter").get_to(param.adp_opt_niter);
  jj.at("adp_niter").get_to(param.adp_niter);
  jj.at("adp_unit_stop").get_to(param.adp_unit_stop);
  jj.at("adp_line_adapt").get_to(param.adp_line_adapt);
  jj.at("adp_stagn_stop").get_to(param.adp_stagn_stop);
  jj.at("adp_smoo_len").get_to(param.adp_smoo_len);
  jj.at("hmin").get_to(param.hmin);
  jj.at("hmax").get_to(param.hmax);
  jj.at("met_snap_tol").get_to(param.met_snap_tol);
  jj.at("anamet_dx").get_to(param.anamet_dx);
  jj.at("anamet_dy").get_to(param.anamet_dy);
  jj.at("anamet_dz").get_to(param.anamet_dz);
  jj.at("curveType").get_to(param.curveType);
  jj.at("smoo_type").get_to(param.smoo_type);
  jj.at("nocleanup").get_to(param.nocleanup);
  jj.at("iflag1").get_to(param.iflag1);
  jj.at("iflag2").get_to(param.iflag2);
  jj.at("iflag3").get_to(param.iflag3);
  jj.at("rflag1").get_to(param.rflag1);
  jj.at("rflag2").get_to(param.rflag2);
  jj.at("rflag3").get_to(param.rflag3);
  jj.at("opt_niter").get_to(param.opt_niter);
  jj.at("opt_pnorm").get_to(param.opt_pnorm);
  jj.at("opt_power").get_to(param.opt_power);
  jj.at("opt_smoo_niter").get_to(param.opt_smoo_niter);
  jj.at("opt_smoo_tol").get_to(param.opt_smoo_tol);
  jj.at("qua_surf_wt_normal").get_to(param.qua_surf_wt_normal);
  jj.at("qua_surf_wt_quality").get_to(param.qua_surf_wt_quality);
  jj.at("opt_coef_det").get_to(param.opt_coef_det);
  jj.at("opt_coef_tra").get_to(param.opt_coef_tra);
  jj.at("opt_powr_det").get_to(param.opt_powr_det);
  jj.at("opt_powr_tra").get_to(param.opt_powr_tra);
  jj.at("opt_unif").get_to(param.opt_unif);
  jj.at("opt_swap_niter").get_to(param.opt_swap_niter);
  jj.at("opt_swap_pnorm").get_to(param.opt_swap_pnorm);
  jj.at("opt_swap_thres").get_to(param.opt_swap_thres);
  jj.at("opt_swap_tet_expensive").get_to(param.opt_swap_tet_expensive);
  jj.at("interp_err_min_algo").get_to(param.interp_err_min_algo);
  jj.at("metScale").get_to(param.metScale);
  jj.at("ianamet").get_to(param.ianamet);
  jj.at("intp_pdeg").get_to(param.intp_pdeg);
  jj.at("intp_pnorm").get_to(param.intp_pnorm);
}

MetrisParameters::MetrisParameters(){
  usrTarDeg = 1;
  nproc     = -1;

  jtol = Defaults::jtol;
  vtol = Defaults::vtol;

  geo_lentolfac = Defaults::geo_lentolfac;
  geo_abstoledg = Defaults::geo_abstoledg;

  anamet_ptr= NULL;
  ianamet   = -1;
  metScale  = 1;
  hmin = 1.0e-30;
  hmax = 1.0e30;
  anamet_dx = anamet_dy = anamet_dz = 0.0;
  met_snap_tol = Defaults::met_snap_tol;

  anasol_ptr= NULL;
  ianasol   = -1;
  anaSol    = false;

  main_in_prefix = false;

  adp_niter     = 0;
  adp_unit_stop = 99.9;
  adp_stagn_stop = Defaults::adp_stagn_stop;
  adp_opt_niter = 1;
  adp_line_adapt = false;
  adp_smoo_len   = false;

  // 0 is none, default
  // 3 is offsets followed by smoothing 
  // 4 is offsets then backtrack and stop there 
  curveType = 0; 
  smoo_type = 0;

  opt_unif = false;

  iverb     = 0;
  ivdepth   = 0;

  interactive = dbgfull = nocleanup = false;

  refineConventionsInp = refineConventionsOut = false;

  opt_pnorm = Defaults::opt_pnorm;
  opt_power = Defaults::opt_power;
  opt_niter = Defaults::opt_niter;
  opt_smoo_niter = Defaults::opt_smoo_niter;
  opt_smoo_tol   = Defaults::opt_smoo_tol;

  opt_swap_pnorm = Defaults::opt_swap_pnorm;
  opt_swap_niter = Defaults::opt_swap_niter;
  opt_swap_thres = Defaults::opt_swap_thres;
  opt_swap_tet_expensive = false;


  qua_surf_wt_normal  = Defaults::qua_surf_wt_normal;
  qua_surf_wt_quality = Defaults::qua_surf_wt_quality;


  opt_coef_det = 1.0;
  opt_powr_det = -2;

  opt_coef_tra = 1.0;
  opt_powr_tra =  2;

  intp_pdeg  = 1;
  intp_pnorm = 1;

  // Private members (internal use)
  wrtMesh   = false;
  inpMesh   = false;
  inpBack   = false;
  inpCAD    = false;
  inpMet    = false;
  anaMet    = false;
  scaleMet  = false;
  outbasis  = FEBasis::Lagrange;

  logFile = stdout;

  iflag1 = iflag2 = iflag3 = 0;
  rflag1 = rflag2 = rflag3 = 0.0;
  interp_err_min_algo = 1; // 1 for Newton, 0 for DIRECT
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
    inpBack = true;
    backFileName = correctExtension_meshb(opt.m["back"].template as<std::string>());
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
  if(opt.count("adp-stat-stop")){
    adp_stagn_stop = opt.m["adp-stat-stop"].as<double>();
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

  if(opt.count("opt-niter")){
    opt_niter = opt.m["opt-niter"].as<int>();
  }
  if(opt.count("opt-pnorm")){
    opt_pnorm = opt.m["opt-pnorm"].as<int>();
  }
  if(opt.count("opt-power")){
    opt_power = opt.m["opt-power"].as<int>();
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
  METRIS_ENFORCE(geo_abstoledg >= 0.0);
  METRIS_ENFORCE(geo_lentolfac >= 1.0);
  METRIS_ENFORCE_MSG(usrTarDeg >= 1, "Degree < 1 provided through tardeg.");
  METRIS_ENFORCE_MSG(usrTarDeg <= METRIS_MAX_DEG, "Opt -tardeg > METRIS_MAX_DEG = {}",METRIS_MAX_DEG);
}

void MetrisParameters::setMeshIn(std::string meshName){
  inpMesh = true;
  meshFileName = correctExtension_meshb(meshName);
}


void MetrisParameters::setCAD(std::string CADname){
  inpCAD = true;
  cadFileName = correctExtension_egads(CADname);
}

void MetrisParameters::setMetricFile(std::string metName){
  inpMet = true;
  metFileName = correctExtension_solb(metName);
}

void MetrisParameters::setMeshOut(std::string out){
  wrtMesh = true;
  outmFileName = correctExtension_meshb(out);
}


void MetrisParameters::setAnalyticalMetric(int ianamet){
  anaMet = true;
  this->ianamet = ianamet;
}

void MetrisParameters::setAnalyticalMetric(AnaMetFun anamet_ptr){
  anaMet = true;
  this->anamet_ptr = (anamet_proto) anamet_ptr;
}

void MetrisParameters::setAnalyticalSolution(int ianasol){
  anaSol = true;
  this->ianasol = ianasol;
}

void MetrisParameters::setAnalyticalSolution(AnaSolFun anasol_ptr){
  anaSol = true;
  this->anasol_ptr = (anasol_proto) anasol_ptr;
}


void MetrisParameters::setMetricScale(double sclmet){
  scaleMet = true;
  this->metScale = sclmet;
}


void MetrisParameters::setLogFile(std::string fname){
  logFileName = fname;
  if(fname == "stdout"){
    logFile = stdout;
  }else if(fname == "stderr"){
    logFile = stderr;
  }else{
    logFile = fopen(fname.c_str(), "w");
    METRIS_ENFORCE_MSG(logFile != NULL,"Error opening log file {}", fname);
  }
}

} // End namespace
