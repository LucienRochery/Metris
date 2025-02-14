//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_METRIS_PARAMETERS__
#define __METRIS_METRIS_PARAMETERS__

#include "../Mesh/MeshFwd.hxx"
#include "../msh_anamet.hxx"
#include <string>


namespace Metris{

struct MetrisOptions;
enum class FEBasis;
class MetrisRunner;
class MetrisAPI;
class CADInfo;

// Parameters can be set manually, or initialized by a Runner from argc/argv 
// If set manually, use MetrisRunner constructor taking a MetrisParameter as input. 
struct MetrisParameters{

  MetrisParameters();

  MetrisParameters(MetrisOptions &opt);

  void setMeshIn(std::string inp);
  void setMeshOut(std::string out);
  void setAnalyticalMetric(int ianamet);
  void setAnalyticalMetric(AnaMetFun anamet_ptr); 
  void setMetricScale(double sclmet);

  std::string outmFileName;
  std::string outmPrefix;
  bool main_in_prefix; // whether main output goes in prefix (default no)

  std::string cadFileName;

  std::string backFileName;

  std::string metFileName;

  // Target degree as specified by -tardeg
  int usrTarDeg;

  // Number of cores for multi-threading
  int nproc;

  // Minimum admissible scaled control coefficient 
  double jtol;

  // Measure tolerance (similar to height)
  double vtol;

  // For adaptGeoLines: edges in target [tarlen / tolfac ; tarlen * tolfac]
  // tarlen is determined in terms of the actual curve length, to be closest 
  // to 1, while having integer number of vertices on the curve. 
  double geo_lentolfac;
  // Absolute tolerance for point to be on CAD edge
  double geo_abstoledg;


  // ----------------- Adaptation options  
  int adp_opt_niter; // do smoothing/swapping in adaptation loop (expensive)
  int adp_niter;
  bool adp_line_adapt;
  // ----------------- END Adaptation options  

  // Metric min/max size control
  double hmin, hmax;

  // options: "curve"
  // Defaults to 0 (no curve), 1 for offsets followed by ccoef max, 2 for 
  // metric-based LP, 3 for offsets followed by smoothing, 4 offsets and backtrack
  int curveType;



  
  // Verbosity level
  int iverb;
  // Verboity depth
  int ivdepth;

  // Full debugs (costly), wait at certain errors (debug only)
  bool dbgfull, interactive;

  // Generic integer flags
  int iflag1, iflag2, iflag3;

  bool refineConventionsInp;
  bool refineConventionsOut;
  
  // ----------------- Optimization options  
  // -- Smoothing 
  int opt_niter;
  int opt_pnorm;
  int opt_power;
  int opt_smoo_niter;

  // -- Quality (quafun_unit)
  // compute coef_det (det - 1)^powr_det + idem(tra)
  double opt_coef_det, opt_coef_tra;
  int    opt_powr_det, opt_powr_tra;

  // experimental 
  bool opt_unif;

  // -- Swaps  
  // Maximum global swap iterations. Stops before if nothing to do. 
  int opt_swap_niter; 
  // lp norm to improve: 0 is infinity (max). Governed by option -qswap-norm
  // but only in the optimization module (not adaptation). 
  // Note quality norm is defined by -opt-power and such options.
  int opt_swap_pnorm;
  // Minimum increase in quality function to go through with a swap. Governed by
  // option -qswap-thres in optimization module (not adaptation). 
  double opt_swap_thres;
  // ----------------- END Optimization options  


  FEBasis outbasis;

  friend class MetrisRunner;
  friend class MetrisAPI;
  friend class MeshBase;
  friend class MeshBack;
  template<class MFT>
  friend class Mesh;
  friend class CADInfo;
public:  
  bool anaMet;
  bool inpMet;
  bool inpBack;
  bool inpCAD;
  bool inpMesh;

  bool wrtMesh;
  std::string meshFileName;

  // Metric scaling 
  bool scaleMet;
  double metScale; 

  // See anamet.hxx. Can implement your own with same prototype
  anamet_proto anamet_ptr;
  // To use defaults, see anamet_2D and anamet_3D.cxx 
  int ianamet; 

};

} // End namespace
#endif
