//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_METRIS_RUNNER__
#define __METRIS_METRIS_RUNNER__

#include "MetrisParameters.hxx"

#include "../Mesh/MeshBack.hxx"
#include "../metris_options.hxx"
#include "../utils/aux_MinMaxAvg.hxx"

#include "fmt/format.h"
#include "nlohmann/json_fwd.hpp"

#include <fstream>
// #include <string>

namespace Metris{

/*
MetrisRunner coordinates front and back, provides high-level functionality
Basically the API as a whole.
See src/main_adap.cxx or bunit/common_setup.hxx for useage
The MetrisAPI class only handles interfacing (data) i.e. gets, sets and file IO
*/

class MetrisAPI;
class MeshStat;

class MetrisRunner{
public:
  MetrisRunner();

  // Initialize from argc, argv identical to CLI call; see metris_options.hxx
  MetrisRunner(int argc, char** argv);

  // Initialize from API objects. These are destroyed by constructor.
  // They can be NULL if the parameters specify input mesh
  MetrisRunner(MetrisAPI *data_front, MetrisAPI *data_back, MetrisParameters &p);
  MetrisRunner(MetrisAPI *data_front, MetrisParameters &p) : MetrisRunner(data_front,NULL,p) { }

  void constructorCommon(MetrisAPI *data_front, MetrisAPI *data_back);

  ~MetrisRunner();


  // -- Primary

  // Calls the other functions in a certain order, default call.
  void runMetris();
  /* Similar to the above, but approaches slowly to the target metric,
     taking into account the scale ratio between current metric and target
  */
  template <class MetricFieldType>
  void runMetrisProgressive();

  // -tardeg <d> Mesh goes to degree d while conserving geometry
  int degElevate();

  // -adapt Adapt to the metric field:
  //   -anamet <int> (cf src/anamet.hxx and anamet2D.cxx anamet3D.cxx)
  //   -met <fname.sol(b)>
  void adaptMesh();
  void adaptMesh2();

  double optimMesh();

  // -smooth <int> High-order smoothing, only -s 4 is vaguely functional and no Jacobian correction yet
  void curveMesh();

  // Print mesh statistics
  void statMesh(int tdim = 0, MeshStat *stat = NULL);

  // -out <fname(.mesh(b))>
  void writeOutputs();

public:
  bool metricFE;

  // The order of these two MUST NOT be changed ! (construction order)
  MetrisOptions opt;
  MetrisParameters param_;
  MetrisParameters * const param; // for INCVDEPTH convenience

  MeshBase *msh_g;
  MeshBack bak;


private:
  friend class MetrisAPI;

  template<class MetricFieldType>
  void iniMetris(MetrisAPI *data_front,MetrisAPI *data_back);

  template<class MetricFieldType>
  void degElevate0();

  template<class MetricFieldType, int gdim, int ideg>
  void adaptMesh0(int tdim);

  template<class MetricFieldType, int gdim, int ideg>
  void adaptMeshQuality0(int tdim);

  template<class MetricFieldType, int gdim, int ideg>
  double optimMesh0();

  template<class MetricFieldType>
  void curveMesh0();

  template<class MetricFieldType>
  void statMesh0(int tdim, MeshStat *stat);

  template<class MetricFieldType>
  void writeOutputs0();

  template<class MetricFieldType>
  void set_array_debugids();

  int nbpo0;
  bool objectiveLineAdapted = false;

  std::fstream foutputAdaptStats;

#ifdef OUTPUTTIMEANDUNITINFO
  bool printUnit;
  std::fstream foutputTimeUnit;
#endif
  //// If an API has been initialized with this runner, we need to let it hard-copy
  //// when we free. Only one is allowed, we can make this an array in the future
  //// if necessary.
  //MetrisAPI *hookedAPI;
  //void moveAPI();
  //// These are used to check the runner didn't change state.
  //// Otherwise the user could have some data change under their feet
  //int onhook_npoin, onhook_nedge, onhook_nface, onhook_nelem;
  //// These are used to convert back if needed
  //FEBasis onhook_mshbasis, onhook_metbasis;
  //MetSpace onhook_metspace;
};




} // End namespace



#endif
