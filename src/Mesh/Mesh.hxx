//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_MESH__
#define __METRIS_MESH__

#include "../types.hxx"
#include "../Mesh/MeshBack.hxx"


//template<class MetricFieldType> class MetrisRunner<MetricFieldType>;
namespace Metris{

class MeshBack;

template<class T>
class Mesh;


template<class MetricFieldType>
class Mesh : public MeshMetric<MetricFieldType>{
public:

  MeshClass meshClass() const override { return MeshClass::Mesh; }

	Mesh() : MeshMetric<MetricFieldType>(1,1,1,1,1) {
    bak = NULL;
  }
  // poi2bak gives MeshBack germs (elements) of dimension given by second index
  // +1. So poi2bak(ipoin,tdim-1) gives a tdim germ for ipoin. 
	intAr2 poi2bak;
  MeshBack *bak;

  /* ----- Back mesh interpolation ----- 
  Compute metric at mesh vertex ipoin and update met, poi2bak
  Common args:
    - iseed/tdim: FRONT element close to coop of dimension tdim
    - iref: reference of elements to localize in. -1 unconstrained
    - algnd: for boundary localization. In case tdim == 1, give tangent. 
    This works both in 2D and 3D, unlike normal. 
    In case tdim == 2, give normal. */

  int interpMetBack(int ipoin, int tdim, int iseed, int iref, 
                    const double* algnd);

private:
  int interpMetBack0(int ipoi0, 
                     int tdim, int iseed, 
                     int iref, const double*__restrict__ algnd,
                     int*__restrict__ ieleb,
                     double*__restrict__ barb);

public:
  // Delete dead elements and minimize array sizes 
  void cleanup();


  // -- Internal or deprecated: 
public:

  int newpoitopo(int tdimn, int ientt = -1);

  // Called from MetrisRunner
  void initialize(MetrisAPI *data, MeshBack &bak, 
    MetrisParameters &param);

  
  void set_npoin(int npoin, bool skipallocf = false) override;
  void set_nentt(int tdimn, int nentt, bool skipallocf = false) override;


protected:
  void initializeCommon(MetrisAPI *data, MeshBack &bak, 
    MetrisParameters &param);
};

} // End namespace

#endif