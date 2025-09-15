//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_MESH__
#define __METRIS_MESH__

#include "../types.hxx"
#include "../MetricField/MetricField.hxx"
#include "MeshMetric.hxx"
#include "MeshBack.hxx"


//template<class MetricFieldType> class MetrisRunner<MetricFieldType>;
namespace Metris{


template<class MetricFieldType>
class Mesh : public MeshMetric<MetricFieldType>{
public:

  MeshClass meshClass() const override { return MeshClass::Mesh; }

	Mesh() : MeshMetric<MetricFieldType>() {
    bak = NULL;
  }
  // poi2bak gives MeshBack germs (elements) of dimension same as the point. 
  // if the point is corner, the "element" is a corner in the back mesh.
	intAr1 poi2bak;
  MeshBack *bak;

  /* ----- Back mesh interpolation ----- 
  Compute metric at mesh vertex ipoin and update met, poi2bak
  Common args:
    - iseed/tdim: FRONT element close to coop of dimension tdim
    - iref: reference of elements to localize in. -1 unconstrained
    - algnd: for boundary localization. In case tdim == 1, give tangent. 
    This works both in 2D and 3D, unlike normal. 
    In case tdim == 2, give normal. */

  // If ipoin is properly initialized (has poi2ent and t/(u,v)), use this:
  int interpMetBack(int ipoin);
  // As previously but provide the point seed explicitly. 
  int interpMetBack(int ipoin, int ipseed);
  // Otherwise supply element seed and compute algnd outside:
  int interpMetBack(int ipoin, int tdim, int iseed, int iref, 
                    const double* algnd);
  // When the seed is a point, not an element. The seed point should be dim
  // <= dim of ipoi0. Algnd should be provided if needed (tdim < msh.idim)
  int interpMetBack00(int ipoi0, int tdim, int iref, int ipseed,
                      const double*__restrict__ algnd);

private:
  int interpMetBack0(int ipoi0, 
                     int tdim, int iseed, 
                     int iref, 
                     const double*__restrict__ algnd);



public:
  // Delete dead elements and minimize array sizes 
  void cleanup();


// -- Internal
public:
  int newpoitopo(PointType ptype, int tdimn, int ientt = -1);

  // Called from MetrisRunner
  void initialize(MetrisAPI *data, MeshBack &bak, MetrisParameters &param);

  void set_npoin(int npoin, bool skipallocf = false) override;
  void set_nentt(int tdimn, int nentt, bool skipallocf = false) override;

protected:
  void initializeCommon(MetrisAPI *data, MeshBack &bak, MetrisParameters &param);
};

} // End namespace

#endif