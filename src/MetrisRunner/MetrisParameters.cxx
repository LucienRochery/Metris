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

  adaptIter = 0;

  curveType = 0;

  opt_unif = false;

  iverb     = 1;
  dbgfull   = false;
  refineConventions = false;

  opt_pnorm = Defaults::qopt_pnorm;
  opt_power = Defaults::qopt_power;
  opt_niter = Defaults::qopt_niter;
  opt_smoo_niter = Defaults::qopt_smoo_niter;
  opt_swap_pnorm = Defaults::qopt_swap_niter;


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
}

MetrisParameters::MetrisParameters(MetrisOptions &opt) : MetrisParameters(){

  if(opt.count("refine-conventions")){
    refineConventions = true;
  }
  
  if(opt.count("help")) { 
    std::cout << "Cf MetrisOptions class" <<"\n";
    exit(1);
  }
  if(opt.count("verb")){
    iverb = opt.m["verb"].template as<int>();
  }else{
    iverb = 1;
  }

  if(opt.count("qopt-unif")){
    opt_unif = true;
    if(iverb >= 1) std::cout << "-- Set qopt-unif \n";
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
    if(opt.count("sclmet")){
      scaleMet = true;
      double scl = opt.m["sclmet"].template as<double>();
      if(iverb >= 1) printf("Analytical metric scaling by %f \n",scl);
      metScale = scl;
    }
  }

  if(opt.count("adapt")){
    adaptIter = opt.m["adapt"].template as<int>();
  }


  if(opt.count("curve")){
    curveType = opt.m["curve"].template as<int>();
    METRIS_ENFORCE(curveType == 1 || curveType == 2);
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

  if(opt.count("qopt-niter")){
    opt_niter = opt.m["qopt-niter"].as<int>();
  }
  if(opt.count("qopt-pnorm")){
    opt_pnorm = opt.m["qopt-pnorm"].as<int>();
  }
  if(opt.count("qopt-power")){
    opt_power = opt.m["qopt-power"].as<int>();
  }
  if(opt.count("qopt-smoo-niter")){
    opt_smoo_niter = opt.m["qopt-smoo-niter"].as<int>();
  }

  if(opt.count("qopt-swap-pnorm")){
    opt_swap_pnorm = opt.m["qopt-swap-pnorm"].as<int>();
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
