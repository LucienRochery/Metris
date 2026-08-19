//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_OPTIONS__
#define __METRIS_OPTIONS__


#include <boost/program_options.hpp>
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
      ("vdepth"   , po::value<int>(),
        "Verbosity depth 0+. Call stack / loop depth prints")
      ("in"     , po::value<std::string>(), "Input mesh file "     )
      ("cad"    , po::value<std::string>(), "Input CAD file "      )
      ("met"    , po::value<std::string>(), "Input metric file "   )
      ("back"   , po::value<std::string>(), "Input back mesh file ")
      ("out"    , po::value<std::string>(), "Output mesh file "    )
      ("prefix" , po::value<std::string>(), "Output prefix, default ./"    )
      ("main-in-prefix", "Prefix applies to main output (default no)")
      ("bez"    , "Output format Bézier (default Lagrange)")
      ("jtol", po::value<double>(), "Scaled Jacobian control coefficient minimum")
      ("vtol", po::value<double>(), "Flatness tolerance")
      ("met-snap-tol", po::value<double>(), "Surface metric snapping tolerance")
      ("curve" , po::value<int>(), "Metric-based smoothing type."
          " Type 1: Offsets followed by ccoef maximization."
          " Type 2: metric-based LP.")
      ("smoo-type", po::value<int>(), "Smoothing type: 0 (default) = metric-based"
                                      " 1 = interpolation error based (use -anasol)")
      ("tardeg" , po::value<int>(), "Target mesh degree"   )
      ("nosort" , "Disable Hilbert reordering"   )
      ("dbgfull"    , "Enable expensive debugs"   )
      ("interactive", "Enable wait() calls (debug)"   )
      ("nocleanup", "Disable cleanup (debug)"   )
      ("nproc"  , po::value<int>(), "Maximum number of CPU cores for multi-threading"   );

    s.add_options()
      ("log", po::value<std::string>(),
       "Log file name, default stdout");

    s.add_options()
      ("refine-conventions-inp", "Adopt Refine conventions for VerticesOnGeometricX");
    s.add_options()
      ("refine-conventions-out", "Adopt Refine conventions for VerticesOnGeometricX");

    // ----------------- Metric and solution options
    s.add_options()
      ("anamet" , po::value<int>(),
        "Analytical metric index, see src/anamet.hxx for options")
      ("sclmet" , po::value<double>(),
        "Analytical metric scaling")
      ("hmin" , po::value<double>(), "Minimum metric size"   )
      ("hmax" , po::value<double>(), "Maximum metric size"   )
      ("mdx" , po::value<double>(), "Analytical metric x offset"   )
      ("mdy" , po::value<double>(), "Analytical metric y offset"   )
      ("mdz" , po::value<double>(), "Analytical metric z offset"   );

    s.add_options()
      ("anasol" , po::value<int>(),
        "Analytical solution index, see src/anasol.hxx for options")
      ("intp-pdeg" , po::value<int>(),
        "Solution interpolation degree <= mesh degree.")
      ("intp-pnorm" , po::value<int>(),
        "Interpolation error norm 1 or 2.");

    // ----------------- Adaptation options
    s.add_options()
      ("adapt"  , po::value<int>(),
        "Adaptation iterations")
      ("adp-unit-stop", po::value<double>(),
        "Percent unit edges to stop adaptation, default 99.9%")
      ("adp-stagn-stop", po::value<double>(),
        "Stat value (work / entities) threshold to stop adaptation. Default 1e-3")
      ("adp-opt-niter", po::value<int>(),
        "Smoothing in adaptation: -1 unlimited, N > 0 number of iter")
      ("do-line-adp",
        "Use adaptGeoLines (not very robust if boundary very coarse)")
      ("adp-smoo-len",
        "Use length-based smoothing in adaptation loop")
      ("opt-unif" ,
        "Shape preserving uniformization")
      ("geo-lentolfac", po::value<double>(),
        "Tolerance factor for geometric edge length in adaptGeoLines")
      ("geo-abstoledg", po::value<double>(),
        "Absolute distance tolerance such that point is considered on CAD edge");

    // ----------------- Optimization options
    s.add_options()
      ("opt-niter" , po::value<int>(),
                    "Apply <x> itertions quality-based optim after of adaptation")
      ("opt-pnorm", po::value<int>(),
                    "Optimization pnorm parameter. Compute Q^(opt-power) in norm pnorm")
      ("opt-power", po::value<int>(),
                    "Optimization power parameter. If power == -1, Q ~ det / tra. "
                    "Otherwise ~ tra / det.")
      ("objective-p", po::value<double>(),
                    "Pointwise objective exponent p (>= 1, default 1).")
      ("objective-quadrature-order", po::value<int>(),
                    "Objective quadrature order. By default, order 4 is used "
                    "for triangles and order 3 for tetrahedra. Order -1 "
                    "selects this automatic behavior explicitly; order 0 "
                    "selects the historical vertex-barycenter rule; orders "
                    "2 through 5 select positive simplex rules of those "
                    "degrees.")
      ("step-distance-p", po::value<double>(),
                    "Deprecated alias for --objective-p.")
      ("step-distance-regularization", po::value<double>(),
                    "Smooth StepDistance norm regularization epsilon (> 0).")
      ("step-distance-shape-volume",
                    "Use the volume-stiffened shape/volume SPD distance. "
                    "Requires objective-p >= 1. "
                    "This path uses the frozen geometric volume factor and "
                    "does not use the collapse barrier.")
      ("step-distance-cavity-target-average",
                    "Use the arithmetic mean of unweighted reference-space "
                    "elemental StepDistance integrals. Cavity replacements "
                    "must strictly improve the mesh-wide mean. This path "
                    "does not use the collapse barrier.")
      ("step-distance-cavity-global-tolerance", po::value<double>(),
                    "Deprecated compatibility option; ignored by the strict "
                    "global-improvement CavityTargetAverage criterion.")
      ("step-distance-cavity-global-gain-fraction", po::value<double>(),
                    "Deprecated compatibility option; ignored by the strict "
                    "global-improvement CavityTargetAverage criterion.")
      ("step-distance-barrier-rho0", po::value<double>(),
                    "Metric-volume threshold rho0 for the StepDistance collapse barrier.")
      ("step-distance-barrier-beta", po::value<double>(),
                    "Nonnegative StepDistance collapse-barrier coefficient beta.")
      ("opt-smoo-niter", po::value<int>(),
                    "Inner optimization loop global smoothing iterations")
      ("opt-smoo-tol", po::value<double>(),
                    "Quality improvement tolerance (absolute for qual in range 0-1)"
                    " to freeze vertex smoothing. Default 0.005")
      ("opt-swap-pnorm", po::value<int>(),
                    "Optimization pnorm parameter (default same as smoothing). "
                    "Compute Q^(opt-power) in norm pnorm")
      ("opt-swap-thres", po::value<double>(),
                    "Quality p-norm (over shell) increase threshold for swaps.")
      ("opt-swap-niter", po::value<int>(),
                    "Inner optimization loop global swapping iterations")
      ("opt-swap-tet-expensive",
                    "Do expensive tet swaps (~ 6x slowdown on 3D swaps)")
      ("qua-surf-wt-quality", po::value<double>(),
                    "Weight of raw quality in surface quality.")
      ("qua-surf-wt-normal", po::value<double>(),
                    "Weight of normal deviation in surface quality.");


    // ----------------- Normal deviation options
    s.add_options()
      ("nordev-tol" , po::value<double>(),
                    "range in [0,1]. Loosen nordev tol obtained from the current cav. 0 keeps one from cav")
      ("nordev-max", po::value<double>(),
                    "range in [0,1]. Maximum nordev allowed");

    // ----------------- Progressive adaptation
    s.add_options()
      ("progressiveAdapt" , po::value<bool>(),
                    "For qual-based algo, must be used if target metric too different from intial mesh");

    // ----------------- Progressive adaptation
    s.add_options()
      ("MOESS_adapt_it" , po::value<int>(),
                    "To easily append this number to any outputs and do not overwrite them");

                    // ----------------- Generic flags. Used for quick debugging
    s.add_options()
      ("iflag1", po::value<int>(), "Generic integer flag")
      ("iflag2", po::value<int>(), "Generic integer flag")
      ("iflag3", po::value<int>(), "Generic integer flag")
      ("rflag1", po::value<double>(), "Generic real flag")
      ("rflag2", po::value<double>(), "Generic real flag")
      ("rflag3", po::value<double>(), "Generic real flag")
      ("interp-err-min-algo", po::value<int>(),
        "Interpolation error minimization algo: 0 for Newton, 1 for DIRECT");
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
      METRIS_ENFORCE_MSG(c < margv,">= {}  options? if legitimate, increase margv",margv)
      METRIS_ENFORCE(n >= 0);
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
