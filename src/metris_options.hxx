//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_OPTIONS__
#define __METRIS_OPTIONS__


#include <boost/program_options.hpp>
#include "metris_defaults.hxx"
#include "aux_exceptions.hxx"


namespace Metris{

// Note: implicit_value: value if specified but no argument 
//       default_value: value if not specified at all
  
struct MetrisOptions{

  MetrisOptions(): s("Metris"){ 
    namespace po = boost::program_options;
    s.add_options()
      ("help"   , "Print help message") // a bool parameter
      ("verb"   , po::value<int>(), 
        "Verbosity level 0-3. 0: minimal prints. 1: large steps. 2: 1+. 3: full debug.")
      ("in"     , po::value<std::string>(), "Input mesh file "     )
      ("cad"    , po::value<std::string>(), "Input CAD file "      )
      ("met"    , po::value<std::string>(), "Input metric file "   )
      ("back"   , po::value<std::string>(), "Input back mesh file ")
      ("out"    , po::value<std::string>(), "Output mesh file "    )
      ("prefix", po::value<std::string>(), "Output prefix, default ./"    )
      ("bez"    , "Output format Bézier (default Lagrange)")
      ("adapt"  , po::value<int>()
                  ->default_value(0)->implicit_value(-1),"Adapt mesh")
      ("jtol", po::value<double>(), "Scaled Jacobian control coefficient minimum")
      ("vtol", po::value<double>(), "Flatness tolerance")
      ("curve" , po::value<int>(), "Apply metric-based smoothing."
          " Type 1: Offsets followed by ccoef maximization."
          " Type 2: metric-based LP.")
      ("tardeg" , po::value<int>(), "Target mesh degree"   )
      ("nosort" , "Disable Hilbert reordering"   )
      ("dbgfull" , "Enable expensive debugs"   )
      ("nproc"  , po::value<int>(), "Maximum number of CPU cores for multi-threading"   )
      ("anamet" , po::value<int>(), "Analytical metric index, see src/anamet.hxx for options" )
      ("sclmet" , po::value<double>(), "Analytical metric scaling"   ) 
      ("maxmem" , po::value<double>(), "Maximum memory to use in Mib"   )
      ("ratmem" , po::value<double>(), "Maximum memory to use as proportion (0 to 1) of total memory" )
      ("hmin" , po::value<double>(), "Minimum metric size"   )
      ("hmax" , po::value<double>(), "Maximum metric size"   )
      ("adp-opt-niter", po::value<int>(), "Flag for smoothing to unstuck adaptation (expensive)")
      ("opt-unif" , "Shape preserving uniformization");


    s.add_options()
      ("refine-conventions", "Adopt Refine conventions for VerticesOnGeometricX");

    // ----------------- Adaptation options  
    s.add_options()
      ("geo-lentolfac", po::value<double>()->default_value(Defaults::geo_lentolfac),
        "Tolerance factor for geometric edge length in adaptGeoLines")
      ("geo-abstoledg", po::value<double>(), "Absolute distance tolerance "
        "such that point is considered on CAD edge");

    // ----------------- Optimization options  
    s.add_options()  
      ("opt-niter" , po::value<int>(),
                    "Apply <x> itertions quality-based optim after of adaptation")
      ("opt-pnorm", po::value<int>(),
                    "Optimization pnorm parameter. Compute Q^(opt-power) in norm pnorm")
      ("opt-power", po::value<int>(),
                    "Optimization power parameter. If power == -1, Q ~ det / tra. "
                    "Otherwise ~ tra / det.")
      ("opt-smoo-niter", po::value<int>(),
                    "Inner optimization loop global smoothing iterations")
      ("opt-swap-pnorm", po::value<int>(),
                    "Optimization pnorm parameter (default same as smoothing). "
                    "Compute Q^(opt-power) in norm pnorm")
      ("opt-swap-thres", po::value<double>(),
                    "Quality p-norm (over shell) increase threshold for swaps.")
      ("opt-swap-niter", po::value<int>(),
                    "Inner optimization loop global swapping iterations");


  }

  MetrisOptions(int argc, char **argv): MetrisOptions() {
      parse(argc, argv);
  }

  boost::program_options::options_description s;
  boost::program_options::variables_map m;

  void parse(int argc, char **argv){
    namespace po = boost::program_options;
    po::store(
      po::command_line_parser(argc,argv)
      .options(s)
      .style(po::command_line_style::unix_style |
             po::command_line_style::allow_long_disguise)
      .allow_unregistered()
      .run(),
      m);
    po::notify(m);
  }

  int count(std::string str){
    return m.count(str);
  }
  
};


struct cargHandler{
  cargHandler() : margv(256){
    v = (char**) malloc(margv*sizeof(char*));
    c = 0;
  }
  cargHandler(std::string cmd) : cargHandler() {
    setArgs(cmd);
  }
  ~cargHandler(){
    while(c --> 1){ // Goes to zero operator :)
      free(v[c]);
    }
    free(v);
  }

  void setArgs(std::string cmd){
    while(c --> 1){ // Goes to zero operator :)
      free(v[c]);
    }

    // First argument ignored
    c = 1;
    std::string str;
    std::stringstream cmd_(cmd);
    while(std::getline(cmd_, str, ' ')){
      int n = str.length();
      if(c >= margv) METRIS_THROW_MSG(SMemExcept(), "256 options? if legitimate, increase margv");
      v[c] = (char *) malloc((n+1)*sizeof(char));
      strncpy(v[c],str.c_str(),n+1);
      c++;
    }
  }

  int c;
  char **v;
  const int margv;
};


} // End namespace

#endif