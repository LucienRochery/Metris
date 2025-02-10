//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MESH_BACK__
#define __METRIS_MESH_BACK__


#include "MeshMetric.hxx"


namespace Metris{

class MetrisAPI;
class MetricFieldFE;

class MeshBack : public MeshMetric<MetricFieldFE>{

public:


  // Store geometric deviation for edges
  // geodev = 1 - abs(dtprd) at nodes 
  // dtprd between CAD and discrete tangent 
  dblAr1 edg2dev;
  // Geometric deviation for faces 
  dblAr1 fac2dev;

        dblAr1&  ent2dev(int tdimn);
  const dblAr1&  ent2dev(int tdimn) const;


  //// Store outgoing normals cone based on CAD for each edge.
  //// This is an upgrade to the geodev because geodev is unsigned, this is 
  //// essentially signed geodev. 
  //// We store in 1:idim the cone center axis and at idim + 1 the dtprd tolerance
  //dblAr2 edg2con;


  MeshClass meshClass() const override { return MeshClass::MeshBack; }

	MeshBack() : MeshMetric<MetricFieldFE>(10,10,10,10,10) {};
	~MeshBack(){}

  //template<class MetricFieldType>
  //void initialize(MetrisAPI *data, Mesh<MetricFieldType> *msh, )

  void initialize(MetrisAPI *data, 
    MetrisParameters &param);

	void copyConstants(const MeshBase &msh);

  void set_nedge(int nedge, bool skipallocf = false) override; 
  void set_nface(int nedge, bool skipallocf = false) override; 

public:
  int newpoitopo(int tdimn, int ientt = -1);
  
	void readConstants(int64_t libIdx, int usrMinDeg);
  void readConstants(const MetrisAPI &data, int usrMinDeg);

public:
	void getEnttMemCosts(int *memCostPpoi, int *memCostPbpo, int *memCostPedg, 
		                 int *memCostPfac, int *memCostPelt);

	//Estimation of the number of elements in the mesh, as prescribed by the metric
	double getMetComplexity();


};


} // End namespace

#endif
