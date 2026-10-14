//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
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
#include <vector>


namespace Metris{

class MeshBase;
struct MetrisParameters;

enum class MetricClass{None, MetricFieldFE, MetricFieldAnalytical};
enum class MetricFieldLayout{Continuous, Discontinuous};


class MetricFieldFE{


public:
  virtual MetricClass metricClass() const { return MetricClass::MetricFieldFE; }
	MetricFieldFE() = delete;
	MetricFieldFE(MeshBase &msh_);
  // Copies retain the supporting mesh and own their coefficient storage.
  MetricFieldFE(const MetricFieldFE &inp);
  virtual ~MetricFieldFE() = default;


	//void setdeg(int ideg_){ideg = ideg_;}
	//void setnentt(int nentt_){nentt = nentt_;}

	int getnnmet() const;

	void allocate();

  MetricFieldLayout getLayout() const { return layout; }
  int getDegree() const;
  int getCoefficientCount() const;

  // Validate the supported DG representation against its supporting mesh.
  void validateDiscontinuous() const;
  void requireContinuous(const char *operation) const {
    METRIS_ENFORCE_MSG(layout == MetricFieldLayout::Continuous,
                      "{} requires a continuous metric layout", operation);
  }

  // Initialize independent P1 Lagrange log-metric coefficients to zero (M = I).
  // This replaces the metric representation without changing mesh connectivity.
  void initializeDiscontinuous(int degree = 1);
  int getElementCoefficientIndex(int ientt, int localNode) const;
  const double *elementCoefficient(int ientt, int localNode) const;
  double *elementCoefficient(int ientt, int localNode);

  // Evaluate L = log(M) on a full-dimensional owning entity. Bary derivatives
  // use xi_i = bary[i+1]; Phys derivatives use physical Cartesian coordinates.
  // dlogm[direction*getnnmet()+component] stores the requested first derivatives.
  void getLogMetOnEntity(int ientt, DifVar idiff, const double *bary,
                         double *logm, double *dlogm) const;

  // Evaluate M and derivatives in the same coordinates as getLogMetOnEntity.
  void getMetOnEntity(int ientt, DifVar idiff, const double *bary,
                     double *met, double *dmet) const;

	inline const double *operator[](int i) const {
    requireContinuous("Point-indexed metric access"); return rfld[i];
  }
	inline double *operator[](int i){
    requireContinuous("Point-indexed metric access"); return rfld[i];
  }


  inline const double& operator()(int i, int j) const {
    requireContinuous("Point-indexed metric access"); return rfld(i,j);
  }
  inline double& operator()(int i, int j){
    requireContinuous("Point-indexed metric access"); return rfld(i,j);
  }


	MetricFieldFE &operator=(const MetricFieldFE& inp);

	MetSpace getSpace()const{return ispace;}
  // Never set iforce to true
	void setSpace(MetSpace ispacn, bool iforce = false);
	
  void forceSpaceFlag(MetSpace newsp){
    if(newsp != ispace) requireContinuous("Changing the metric space flag");
    ispace = newsp;
  }
	
  FEBasis getBasis()const{return ibasis;}

	void killBasis(){
    requireContinuous("Clearing the metric basis");
		std::cout<<"## Warning killBasis() called\n";
		ibasis = FEBasis::Undefined;
	}

	void setBasis(FEBasis ibasn){
    validateDiscontinuous();
    if(ibasn == getBasis()) return;
    if(layout == MetricFieldLayout::Discontinuous){
      METRIS_ENFORCE_MSG(ibasn == FEBasis::Lagrange || ibasn == FEBasis::Bezier,
                        "DG metrics support only P1 Lagrange or Bernstein basis");
      // P1 Lagrange and Bernstein coefficients coincide.
      ibasis = ibasn;
      return;
    }
    METRIS_ASSERT(ibasn == FEBasis::Lagrange || ibasn == FEBasis::Bezier);

    MetSpace ispac0 = getSpace();
    setSpace(MetSpace::Log);
    
		if(ibasn == FEBasis::Lagrange) setLagrange();
		else                           setBezier();

    setSpace(ispac0);
	}

  // To be used in initialization
  void forceBasisFlag(FEBasis ibasn){
    if(ibasn != ibasis) requireContinuous("Changing the metric basis flag");
    ibasis = ibasn;
  }

	void readMetricFile(std::string outname);
  // Optional map: source point ID minus one -> current zero-based point ID.
  void readDiscontinuousMetricFile(const std::string &filename,
                                   const std::vector<int> &sourceToPoint = {});
  void writeDiscontinuousMetricFiles(const std::string &meshFilename,
                                     const std::string &metFilename) const;
	void writeMetricFile(std::string outname, bool iprefix = true);

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
                  int tdimn,  const double*__restrict__ bary, 
                  double*__restrict__ metl, double*__restrict__ dmet);

  #if 0
	void getMetPhys(AsDeg metdeg, DifVar idiff, MetSpace tarspac, 
                  int *ieleg, 
		              const double*__restrict__ coop, 
		              double*__restrict__ metl, 
                  double*__restrict__ dmet, int ithread) ;
  #endif

	void getMetFullinfo(AsDeg asdmet,
                      DifVar idiff, MetSpace tarspac, 
                      const int*__restrict__ ent2pol, 
                      int tdimn, 
											const double*__restrict__ bary, 
                      const double*__restrict__ coop, 
											double*__restrict__ metl, 
                      double*__restrict__ dmet) ;

	//// Differentiated, see low_eval_d, for instance with respect to control point position. 
	//void getMetBary_d(const MeshBase &msh, int ientt, int tdimn, DifVar idiff, 
	//	                double*__restrict__ bary, double*__restrict__ metl, double*__restrict__ dmet);

public: 
	// Point-indexed storage for the continuous layout; inactive in the DG layout.
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
  MetricFieldLayout layout = MetricFieldLayout::Continuous;
  int metricDegree = 0;
  intAr2 elementToCoefficient;
  dblAr2 localCoefficients;

  template<int gdim, int metdeg, int geomdeg>
  void getLogMetOnEntity0(int ientt, DifVar idiff, const double *bary,
                          double *logm, double *dlogm) const;

	//template<int gdim, int metdeg>
	//void getMetFullinfo0(const MeshBase &msh, int ientt, int tdimn, DifVar idiff, 
	//										double*__restrict__ bary, double*__restrict__ coop, double*__restrict__ metl, double*__restrict__ dmet);
	template<int gdim, int ideg>
	void getMetBary0(DifVar idiff, MetSpace tarspac, const int*__restrict__ ent2pol, 
                   int tdimn, const double*__restrict__ bary, 
                   double*__restrict__ metl, double*__restrict__ dmet) ;
	template<int gdim, int ideg>
	void getMetPhys0(DifVar idiff, MetSpace tarspac,  int *ieleg, 
		               const double*__restrict__ coop,
		                     double*__restrict__ metl, double*__restrict__ dmet) ;
};


/*
operator[] deals with the FE field in all cases
*/

class MetricFieldAnalytical : public MetricFieldFE{

public:
  friend class MeshBack;
  
  MetricClass metricClass() const override { return MetricClass::MetricFieldAnalytical; }

	MetricFieldAnalytical() = delete;
	MetricFieldAnalytical(MeshBase &msh_);
  virtual ~MetricFieldAnalytical() = default;

  void setAnalyticalMetric(const MetrisParameters& param);
protected:
  void setAnalyticalMetric(int ianamet_);
	void setAnalyticalMetric(anamet_proto fptr);

public:
	MetricFieldAnalytical &operator=(const MetricFieldAnalytical& inp);

	void normalize(double coeff);

	template<int gdim, int tdim, int mshdeg, int tardeg>
	void getMetNodes(int ientt, double *metnod) const;


  // AsDeg just to ensure compatible prototypes, unused option
  void getMetBary(AsDeg asdmet,
                  DifVar idiff, MetSpace tarspac, 
                  const int*__restrict__ ent2pol, 
                  int tdimn,  const double*__restrict__ bary, 
                  double*__restrict__ metl, double*__restrict__ dmet);

  void getMetPhys(DifVar idiff, MetSpace tarspac,
                  const double*__restrict__ coop, 
                  double*__restrict__ metl, 
                  double*__restrict__ dmet) const ;

  void getMetFullinfo(AsDeg asdmet,
                      DifVar idiff, MetSpace tarspac, 
                      const int*__restrict__ ent2pol, 
                      int tdimn, 
                      const double*__restrict__ bary, 
                      const double*__restrict__ coop, 
                      double*__restrict__ metl, 
                      double*__restrict__ dmet) ;

protected:
	int ianamet;
  void (*anamet)(const AnaMetCtx*, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet);
  double scale;
  AnaMetCtx ctx;
  //const MeshMetric<MetricFieldAnalytical> &msh;
  
	// Physical derivatives
	template<int gdim, int ideg>
	void getMetBary0(DifVar idiff, MetSpace tarspac, const int*__restrict__ ent2poi, 
		               int tdimn,  const double*__restrict__ bary,  
                   double*__restrict__ metl, double*__restrict__ dmet) ;
	template<int gdim>
	void getMetPhys0(DifVar idiff, MetSpace tarspac, 
		               const double*__restrict__ coop, 
		                     double*__restrict__ metl, 
                         double*__restrict__ dmet) const;
};
















} // End namespace





































#endif
