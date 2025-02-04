//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../metris_options.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../metris_defaults.hxx"
#include "../metris_constants.hxx"
#include "../io_libmeshb.hxx"
#include <string>

namespace Metris{


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


  adp_niter     = 0;
  adp_opt_niter = 1;
  adp_line_adapt = true;

  // 3 is offsets followed by smoothing 
  // 4 is offsets then backtrack and stop there 
  curveType = 3; 

  opt_unif = false;

  iverb     = 2;
  ivdepth   = 0;
  dbgfull   = false;
  refineConventionsInp = false;
  refineConventionsOut = false;

  opt_pnorm = Defaults::opt_pnorm;
  opt_power = Defaults::opt_power;
  opt_niter = Defaults::opt_niter;
  opt_smoo_niter = Defaults::opt_smoo_niter;
  opt_swap_pnorm = Defaults::opt_swap_niter;


  opt_coef_det = 1.0;
  opt_powr_det = -2;

  opt_coef_tra = 1.0;
  opt_powr_tra =  2;

  // Private members (internal use)
  wrtMesh   = false;
  inpMesh   = false;
  inpBack   = false;
  inpCAD    = false;
  inpMet    = false;
  anaMet    = false;
  scaleMet  = false;
  outbasis  = FEBasis::Lagrange;

  iflag1 = iflag2 = iflag3 = 0;
}

MetrisParameters::MetrisParameters(MetrisOptions &opt) : MetrisParameters(){

  if(opt.count("refine-conventions-inp")){
    refineConventionsInp = true;
  }
  if(opt.count("refine-conventions-out")){
    refineConventionsOut = true;
  }
  
  if(opt.count("help")) { 
    std::cout << "Cf MetrisOptions class" <<"\n";
    exit(1);
  }
  if(opt.count("verb")){
    iverb = opt.m["verb"].template as<int>();
  }
  if(opt.count("vdepth")){
    ivdepth = opt.m["vdepth"].template as<int>();
  }

  if(opt.count("opt-unif")){
    opt_unif = true;
    if(iverb >= 1) std::cout << "-- Set opt-unif \n";
  }

  if(opt.count("in")){
    inpMesh = true;
    meshFileName = correctExtension_meshb(opt.m["in"].template as<std::string>());
    if(iverb >= 1) std::cout << "-- Read input mesh name " << meshFileName << "\n";
  }

  if(opt.count("prefix")){
    outmPrefix = opt.m["prefix"].template as<std::string>();
    if(iverb >= 1) std::cout << "-- File prefix: " << outmPrefix << "\n";
  }

  if(opt.count("bez")){
    outbasis = FEBasis::Bezier;
    if(iverb >= 1) std::cout << "-- Bézier output basis\n";
  }


  if(opt.count("out")) { 
    setMeshOut(opt.m["out"].template as<std::string>());
    if(iverb >= 1) std::cout << "-- Read output file name " << outmFileName << ".\n";
  }else{
    if(iverb >= 1) std::cout << "# Output mesh file name not set. Use --out or -o <filename>.\n";
    if(iverb >= 1) std::cout << "# Running but skipping mesh output."<<"\n";
  }

  // usrMaxDeg is the very maximum the user is allowing for storage. It is hard bounded by the constant METRIS_MAX_DEG
  // usrTarDeg is the minimum degree the user wants. 
  if(opt.count("tardeg")){  
    usrTarDeg = opt.m["tardeg"].template as<int>();
    METRIS_ENFORCE_MSG(usrTarDeg >= 1, "Degree < 1 provided through tardeg.");
    METRIS_ENFORCE_MSG(usrTarDeg <= METRIS_MAX_DEG, "Opt -tardeg > METRIS_MAX_DEG = "<<METRIS_MAX_DEG);
  }

  if(opt.count("nproc")){
    nproc = opt.m["nproc"].template as<int>();
    if(iverb >= 1) printf("-- Running with nproc = %d \n",nproc);
  }

  if(opt.count("cad")){
    inpCAD = true;
    cadFileName = correctExtension_egads(opt.m["cad"].template as<std::string>());
  }

  if(opt.count("dbgfull")){
    if(iverb >= 1) printf("-- Full debugs activated\n");
    dbgfull = true;
  }

  if(opt.count("back")){
    inpBack = true;
    backFileName = correctExtension_meshb(opt.m["back"].template as<std::string>());
    if(iverb >= 1) std::cout<<" - Read back mesh name "<<backFileName<<"\n";
  }

  if(opt.count("met")){
    inpMet = true;
    metFileName = correctExtension_solb(opt.m["met"].template as<std::string>());
  }

  if(opt.count("anamet")){
    anaMet  = true;
    ianamet = opt.m["anamet"].template as<int>();
    if(iverb >= 1) printf("Using analytical metric %d \n", ianamet);
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
  if(opt.count("no-line-adp")){
    adp_line_adapt = false; 
  }

  


  if(opt.count("curve")){
    curveType = opt.m["curve"].template as<int>();
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
    METRIS_ENFORCE(geo_lentolfac >= 1.0);
  }

  if(opt.count("geo-abstoledg")){
    geo_abstoledg = opt.m["geo-abstoledg"].as<double>();
    METRIS_ENFORCE(geo_abstoledg >= 0.0);
  }
  
  

  if(opt.count("hmin")){
    hmin = opt.m["hmin"].as<double>();
  }
  if(opt.count("hmax")){
    hmax = opt.m["hmax"].as<double>();
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

  if(opt.count("opt-swap-pnorm")){
    opt_swap_pnorm = opt.m["opt-swap-pnorm"].as<int>();
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


void MetrisParameters::setMetricScale(double sclmet){
  scaleMet = true;
  this->metScale = sclmet;
}


} // End namespace
