//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_METRIS_RUNNER__
#define __METRIS_METRIS_RUNNER__

#include "../Mesh/MeshBack.hxx"
#include "../metris_options.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"


namespace Metris{

/*
MetrisRunner coordinates front and back, provides high-level functionality 
Basically the API as a whole. 
See src/main_adap.cxx or bunit/common_setup.hxx for useage
The MetrisAPI class only handles interfacing (data) i.e. gets, sets and file IO
*/

class MetrisAPI;

class MetrisRunner{
public:
  MetrisRunner();

  // Initialize from argc, argv identical to CLI call; see metris_options.hxx
  MetrisRunner(int argc, char** argv, bool isilent = false); 

  // Initialize from API objects. These are destroyed by constructor. 
  MetrisRunner(MetrisAPI *data_front, MetrisAPI *data_back, MetrisParameters &p);
  MetrisRunner(MetrisAPI *data_front, MetrisParameters &p) : MetrisRunner(data_front,NULL,p) { }

  void constructorCommon(MetrisAPI *data_front, MetrisAPI *data_back);

  ~MetrisRunner();

  // -- Primary 
  // -tardeg <d> Mesh goes to degree d while conserving geometry
  int degElevate();

  // -adapt Adapt to the metric field:
  //   -anamet <int> (cf src/anamet.hxx and anamet2D.cxx anamet3D.cxx)
  //   -met <fname.sol(b)> 
  void adaptMesh();

  double optimMesh();

  // -smooth <int> High-order smoothing, only -s 4 is vaguely functional and no Jacobian correction yet 
  void curveMesh();

  // Print mesh statistics 
  void statMesh();

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
  void adaptMesh0();
  
  template<class MetricFieldType, int gdim, int ideg>
  double optimMesh0();
  
  template<class MetricFieldType>
  void curveMesh0();

  template<class MetricFieldType>
  void statMesh0();

  template<class MetricFieldType>
  void writeOutputs0();

  template<class MetricFieldType>
  void set_array_debugids();
  
  int nbpo0;

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