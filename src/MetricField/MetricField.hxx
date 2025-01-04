//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_METRIC_FIELD__
#define __METRIS_METRIC_FIELD__

#include "../metris_constants.hxx"
#include "../aux_exceptions.hxx"
#include "../linalg/explogmet.hxx"
#include "../low_eval.hxx"
#include "../msh_anamet.hxx"

#include <boost/pool/pool_alloc.hpp>



// curiously reoccuring template pattern
// -> template arg is the base class
// inherit from classes which implement the 0() 
// "API class" -> 

namespace Metris{


class MeshBase;
class MetricFieldAnalytical;

template<class T> class MeshMetric;

class MetricFieldFE{

public:
	MetricFieldFE() = delete;
	MetricFieldFE(MeshBase &msh_);


	//void setdeg(int ideg_){ideg = ideg_;}
	//void setnentt(int nentt_){nentt = nentt_;}

	int getnnmet() const;

	void allocate();

	inline const double *operator[](int i) const {return rfld[i];}
	inline double *operator[](int i){return rfld[i];}


  inline const double& operator()(int i, int j) const {return rfld(i,j);}
  inline double& operator()(int i, int j){return rfld(i,j);}


	MetricFieldFE &operator=(const MetricFieldFE& inp);

	MetSpace getSpace()const{return ispace;}
	void setSpace(MetSpace ispacn){
		METRIS_ASSERT(ispacn != MetSpace::Undefined);
		if(ispacn == MetSpace::Log) setLog();
		else												setExp();
	}
	
  void forceSpaceFlag(MetSpace newsp){ispace = newsp;}
	
  FEBasis getBasis()const{return ibasis;}

	void killBasis(){
		std::cout<<"## Warning killBasis() called\n";
		ibasis = FEBasis::Undefined;
	}

	void setBasis(FEBasis ibasn){
    if(ibasn == getBasis()) return;
    METRIS_ASSERT(ibasn == FEBasis::Lagrange || ibasn == FEBasis::Bezier);

    MetSpace ispac0 = getSpace();
    setSpace(MetSpace::Log);
    
		if(ibasn == FEBasis::Lagrange) setLagrange();
		else                           setBezier();

    setSpace(ispac0);
	}

  // To be used in initialization
  void forceBasisFlag(FEBasis ibasn){ibasis = ibasn;}

	void readMetricFile(int64_t libIdx);
	void writeMetricFile(std::string outname, MetSpace outspac = MetSpace::Exp);

	template<int gdim, int tdim, int mshdeg, int tardeg>
	void getMetNodes(int ientt, double *metnod) const;

  void correctMetric();

	void normalize(double coeff);

	/*
	Evaluation functions. If you know: 
	- barycentrics: call getMetBary
	- physical coordinates: call getMetPhys
	- both: call getMetFullinfo, will dispatch to appropriate depending on Analytical or FE. 
	Inputs:
	- msh: the mesh attached to this metric field
	- ientt, tdimn: tetra, triangle, edge...
	- ieleg: highest topo dim of msh, guess for localization
					 typically host element or, if back, check poi2bak
	*/

	// Derivatives barycentric of physical as specified in idiff, see metris_constants enums
	// DifVar::Phys or DifVar::Bary (or DifVar::None = 0). 
	void getMetBary(AsDeg asdmet,
                  DifVar idiff, MetSpace tarspac, 
                  const int*__restrict__ ent2pol, 
                  int tdimn,  const double* __restrict__ bary, 
                  double*__restrict__ metl, double*__restrict__ dmet);

	void getMetPhys(AsDeg metdeg, DifVar idiff, MetSpace tarspac, 
                  int *ieleg, 
		              const double* __restrict__ coop, 
		              double*__restrict__ metl, 
                  double*__restrict__ dmet, int ithread = 0) ;

	void getMetFullinfo(AsDeg asdmet,
                      DifVar idiff, MetSpace tarspac, 
                      const int*__restrict__ ent2pol, 
                      int tdimn, 
											const double* __restrict__ bary, 
                      const double* __restrict__ coop, 
											double*__restrict__ metl, 
                      double*__restrict__ dmet) ;

	//// Differentiated, see low_eval_d, for instance with respect to control point position. 
	//void getMetBary_d(const MeshBase &msh, int ientt, int tdimn, DifVar idiff, 
	//	                double* __restrict__ bary, double*__restrict__ metl, double*__restrict__ dmet);

public: 
	dblAr2 rfld;

protected:
	MetSpace ispace;
	FEBasis ibasis;
	//int ideg;
	
	//int nnmet,nentt;
	//const int* npoin; // Attention: the pointer itself is not const, it is the pointeE which is const (what we want)
	//const intAr2 &ent2dof;
	MeshBase &msh; // Not const because of localization and work array usage

	void setLagrange();//(int nwork, double *rwork);
	void setBezier();//(int nwork, double *rwork);

	void setLog();
	void setExp();

  boost::fast_pool_allocator<double> rwrkmem;

  int nspace_miss;

private:
	//template<int gdim, int metdeg>
	//void getMetFullinfo0(const MeshBase &msh, int ientt, int tdimn, DifVar idiff, 
	//										double* __restrict__ bary, double* __restrict__ coop, double*__restrict__ metl, double*__restrict__ dmet);
	template<int gdim, int ideg>
	void getMetBary0(DifVar idiff, MetSpace tarspac, const int*__restrict__ ent2pol, 
                   int tdimn, const double*__restrict__ bary, 
                   double*__restrict__ metl, double*__restrict__ dmet) ;
	template<int gdim, int ideg>
	void getMetPhys0(DifVar idiff, MetSpace tarspac,  int *ieleg, 
		               const double*__restrict__ coop,
		                     double*__restrict__ metl, double*__restrict__ dmet, int ithread = 0) ;
};


/*
operator[] deals with the FE field in all cases
*/

class MetricFieldAnalytical : public MetricFieldFE{

public:
	MetricFieldAnalytical() = delete;
	MetricFieldAnalytical(MeshMetric<MetricFieldAnalytical> &msh_);

  void setAnalyticalMetric(int ianamet_);
	void setAnalyticalMetric(anamet_proto fptr);


	MetricFieldAnalytical &operator=(const MetricFieldAnalytical& inp);

	void normalize(double coeff);

	template<int gdim, int tdim, int mshdeg, int tardeg>
	void getMetNodes(int ientt, double *metnod) const;


  // AsDeg just to ensure compatible prototypes, unused option
  void getMetBary(AsDeg asdmet,
                  DifVar idiff, MetSpace tarspac, 
                  const int*__restrict__ ent2pol, 
                  int tdimn,  const double* __restrict__ bary, 
                  double*__restrict__ metl, double*__restrict__ dmet);

  void getMetPhys(AsDeg metdeg, DifVar idiff, MetSpace tarspac, int *ieleg, 
                  const double* __restrict__ coop, 
                  double*__restrict__ metl, 
                  double*__restrict__ dmet, int ithread = 0) ;

  void getMetFullinfo(AsDeg asdmet,
                      DifVar idiff, MetSpace tarspac, 
                      const int*__restrict__ ent2pol, 
                      int tdimn, 
                      const double* __restrict__ bary, 
                      const double* __restrict__ coop, 
                      double*__restrict__ metl, 
                      double*__restrict__ dmet) ;

protected:
	int ianamet;
  void (*anamet)(void* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet);
  double scale;
  //const MeshMetric<MetricFieldAnalytical> &msh;
  
	// Physical derivatives
	template<int gdim, int ideg>
	void getMetBary0(DifVar idiff, MetSpace tarspac, const int*__restrict__ ent2poi, 
		               int tdimn,  const double*__restrict__ bary,  
                   double*__restrict__ metl, double*__restrict__ dmet) ;
	template<int gdim>
	void getMetPhys0(DifVar idiff, MetSpace tarspac, 
		               const double* __restrict__ coop, 
		                     double*__restrict__ metl, double*__restrict__ dmet) ;
};
















} // End namespace





































#endif