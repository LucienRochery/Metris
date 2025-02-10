//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __MSH_CAVITY__
#define __MSH_CAVITY__

#include "../Mesh/MeshFwd.hxx"
#include "../types.hxx"


namespace Metris{


enum Cavity_Errors {CAV_NOERR = 0, 
                    CAV_ERR_NOBPO = 1,
                    CAV_ERR_TDIMN = 2,
                    CAV_ERR_DUPEDG = 3,
                    CAV_ERR_INTEDG = 4,
                    CAV_ERR_FLATFAC = 5,
                    CAV_ERR_NEGFAC = 6,
                    CAV_ERR_DUPEDG2 = 7,
                    CAV_ERR_QMAXNEC = 8,
                    CAV_ERR_QMAXIFF = 9,
                    CAV_ERR_QFACNEG = 10,
                    CAV_ERR_CADFAR = 11,
                    CAV_ERR_INCORRECTIBLE = 12,
                    CAV_ERR_DRYFAIL1 = 13,
                    CAV_ERR_DRYFAIL2 = 14,
                    CAV_ERR_LINETOPO = 15,
                    CAV_ERR_GEODEVLIN = 16,
                    CAV_ERR_FLATEDG = 17,
                    CAV_ERR_NERROR = 18
                    };


// Minimum print level for cavity
#define METRIS_CAV_PRTLEV 4


class MshCavity{
public:
	MshCavity(){
		//mctet = mcfac = mcedg = 0;
		//nctet = ncfac = ncedg = 0;
		//lctet = lcfac = lcedg = NULL;
	}
  MshCavity(int mctet_, int mcfac_, int mcedg_) 
    : lctet(mctet_), lcfac(mcfac_), lcedg(mcedg_) {}
  
	MshCavity(int mctet_, int mcfac_, int mcedg_, 
		        int nbuff, int *lbuff){

		int tmem_ = mctet_ + mcfac_ + mcedg_;
		if(tmem_ > nbuff)
			METRIS_THROW_MSG(DMemExcept(),"PROVIDE LARGER CAVITY BUFFER")

		int ovft = (nbuff - tmem_) / 2;
		int ovff = (nbuff - tmem_ - ovft) / 2;
		int ovfe = (nbuff - tmem_ - ovff - ovft);


		int mctet = mctet_ + ovft;
		int mcfac = mcfac_ + ovff;
		int mcedg = mcedg_ + ovfe;

		//nctet = ncfac = ncedg = 0;
		lctet.set_buffer(mctet,&(lbuff[0]));
		lcfac.set_buffer(mcfac,&(lbuff[mctet]));
		lcedg.set_buffer(mcedg,&(lbuff[mctet+mcfac]));

    //lctet = &(lbuff[0]);
    //lcfac = &(lbuff[mctet]);
    //lcedg = &(lbuff[mctet+mcfac]);
		nrempts = 0;
		iremcor = -1;
		ipins   = -1;

    nrmal  = NULL;
	}
	void reset(){
		//nctet = 0;
		//ncedg = 0;
		//ncfac = 0;
		ipins =-1;
		nrempts = 0;
		iremcor = -1;
    lcedg.set_n(0);
    lcfac.set_n(0);
    lctet.set_n(0);
    nrmal = NULL;
	}
	~MshCavity(){
		//lctet = NULL;
		//lcfac = NULL;
		//lcedg = NULL;
	}

  intAr1& lcent(int tdimn){
    METRIS_ASSERT(tdimn >= 1 && tdimn <= 3);
    if(tdimn == 1){
      return lcedg;
    }else if (tdimn == 2){
      return lcfac;
    }else{
      return lctet;
    }
  }

  const intAr1& lcent(int tdimn) const {
    METRIS_ASSERT(tdimn >= 1 && tdimn <= 3);
    if(tdimn == 1){
      return lcedg;
    }else if (tdimn == 2){
      return lcfac;
    }else{
      return lctet;
    }
  }

  template<int tdimn>
  constexpr intAr1& lcent(){
    static_assert(tdimn >= 1 && tdimn <= 3);
    if constexpr(tdimn == 1){
      return lcedg;
    }else if (tdimn == 2){
      return lcfac;
    }else{
      return lctet;
    }
  }

  template<int tdimn>
  constexpr const intAr1& lcent() const {
    static_assert(tdimn >= 1 && tdimn <= 3);
    if constexpr(tdimn == 1){
      return lcedg;
    }else if (tdimn == 2){
      return lcfac;
    }else{
      return lctet;
    }
  }
  /* User set data */
  // Usage examples bunit/face_cavityX.cxx 
  
  // m = storage size
  // n = element count 
  // l = element list
	//int  mctet, mcfac, mcedg;
	//int  nctet, ncfac, ncedg;
	//int *lctet,*lcfac,*lcedg;
  intAr1 lctet,lcfac,lcedg;

  // Point to be inserted. ibpoi must be set. If needed, use newedgvirtual/newfacvirtual to give a surface ref 
  // of the lowest topo dim, for a new point. See e.g. bunit/face_cavityX.cxx 
  // Metric at ipins must be computed and naturally stored in msh.met[ipins] if quality is requested, i.e. 
  // if options qmax_nec, qmax_suf or qmax_iff are set > 0. 
	int ipins;

  // Surface normal at ipins. Can be NULL if not surface. 2D doesnt need either. 
  double *nrmal;

  /* End user set data */


  // Internal use
	// Store removed points, whether a corner is removed and if so which one (one at the most)
	int nrempts, iremcor;
};

struct CavOprOpt{
	// If a partial initial cavity is supplied, this should be set to 1.
	// Otherwise, if one is confident in the initial cavity, this should
	// not be set if behaviour as close to initially intended is desired. 
	bool allow_topological_correction;


	// If confident in one's cavity, this setting should be set to true. 
	// It will no doubt provide a large speedup to the first stage. 
	bool skip_topo_checks;

	// Only concerns volume or manifold surface points. 
	bool allow_remove_points;

	// A corner (triple+ pt) is only removed if one is reinserting the same corner. 
	// This defaults to false to avoid accidents. 
	bool allow_remove_corners;

	// How many iterations of cavity extension for geometric reasons are allowed
	int max_increase_cav_geo;

  // Geometric deviation tolerance for edges expressed in normalized dotprod
  double geodev1;

  // In this mode, created elements are immediately tested for P1 validity and 
  // the whole operation is rejected at the first failure. Thought for collapses 
  // where several reinsertions are tested. 
  // It may be better (not tested) to disable this in order to know about all
  // bad entities when using cavity correction after failure. 
  // If cavity correction is disabled (max_increase_cav_geo == 0), this should always be true 
  bool fast_reject;

  // In dryrun mode, modifications are rejected, but qualities are computed
  // This is for collapses or swaps. 
  bool dryrun;


  // Quality = conformity error (higher is worse)
  // qmax is max quality error across tdim == msh.get_tdim()
  // elements and options only in effect if options > 0. 
  // Strongest prevails: iff > (nec | suf)
  //                     nec + suf == iff.
  // If qmax is above value, run rejected. 
  double qmax_nec;
  // If qmax is below value, run validated regardless of dryrun
  double qmax_suf;
  // Combines the two: run validated iff qmax < qmax_iff
  double qmax_iff;
  // Other remarks on quality: 
  // To delay metric interpolation at new vertices and save time, we use 
  // the metric as a P1 field to compute qualities.
  // However this still means the metric at ipins must be computed !!

  // Power: > 0 = (tra / det)^|power| 
  // Negative is more stable but less discriminant 
  // This can be solved by using higher pnorm (or power). 
  int qpower;
  int qpnorm;


	CavOprOpt():allow_topological_correction(true)
             ,skip_topo_checks(true)
	           ,allow_remove_points(true)
	           ,allow_remove_corners(false)
	           ,max_increase_cav_geo(0)
             ,geodev1(0.5)
             ,fast_reject(false)
             ,dryrun(false)
             ,qmax_nec(-1.0)
             ,qmax_suf(-1.0)
             ,qmax_iff(-1.0)
             ,qpower(-1)
             ,qpnorm(2)
             {}
};

// Cavity operator returns
struct CavOprInfo{
  // Anisotropic qualities start, end (for swaps, not implemented yet)
  double qmax_ini,qavg_ini;
  double qmax_end,qavg_end;
  bool done; // flags whether change was done (dryrun) ; different from an error
};

// Don't worry about this, simply declare one and reuse for all cavity calls. 
// Automatic reallocation is done, this is 10x faster than using boost pool allocators, if less elegant.
struct CavWrkArrs{
  intAr2 lbad;

  dblAr2 lmeas;

  intAr1 edtyp;
  intAr1 edent;

  intAr1 lfcco;

  // Store normals for each connex component of the cavity.
  dblAr2 lnorf;

  int tagf0;

  CavWrkArrs(){
    lbad.allocate(100,2);

    edtyp.allocate(100);
    edent.allocate(100);

    lfcco.allocate(100);

    lnorf.allocate(10,3);

    tagf0 = -1;
  }

  //void set_mmeas(int mmeas_){
  //  if(mmeas_ <= mmeas) return;
  //  mmeas = mmeas_;
  //  lmeas.allocate(2,mmeas); 
  //}

  //void set_mfcco(int mfcco_){
  //  if(mfcco_ <= mfcco) return;
  //  mfcco = mfcco_;
  //  lfcco.allocate(mfcco);
  //}

  //void set_nedex(int nedex){
  //  METRIS_ASSERT(nedex >= 0);
  //  if(nedex >= medex){
  //    medex *= 1.5;
  //    edtyp.allocate(medex);
  //    edent.allocate(medex);
  //  }
  //  edtyp.set_n(nedex);
  //  edent.set_n(nedex);
  //}
  //void check_medex(int nedex){
  //  if(nedex >= medex){
  //    medex *= 1.5;
  //    edtyp.allocate(medex);
  //    edent.allocate(medex);
  //  }
  //}

};



template<class MetricFieldType, int ideg>
int cavity_operator(Mesh<MetricFieldType> &msh ,
                   MshCavity &cav,
                   CavOprOpt &opts,
                    CavWrkArrs &work,
                    CavOprInfo &info,
                   int ithread = 0);

// Check at most one corner
// Check edges attached to faces attached to tets
// Check no int points if norempts option set
template<class MetricFieldType>
int check_cavity_topo(Mesh<MetricFieldType> &msh, MshCavity &cav, 
                      CavOprOpt &opts,//RoutineWorkMemory<int> &iwrk,
                      int ithread = 0);

template<class MetricFieldType, int ideg>
int update_cavity(Mesh<MetricFieldType> &msh, const MshCavity &cav, const CavWrkArrs &work, 
                 int npoi0, int nedg0, int nfac0, int nele0, int ithread = 0);


// The boundary of the line cavity is simply the set of end points.
// There can be arbitrarily many in a non-manifold mesh if a corner 
// point is included in the line cavity. 
template<class MetricFieldType, int ideg>
int reconnect_lincav(Mesh<MetricFieldType> &msh, const MshCavity& cav, CavOprOpt &opts, 
                     int ithread = 0);//,int mnwedg, int *nnwedg,int lnwedg[]);


// The boundary here is a set of edges. They will be reconnected to 
// ipins. 
// Note: if the initial cavity includes edges, these have been reconnected
// and the result is in lnewed
template<class MetricFieldType, int ideg>
int reconnect_faccav(Mesh<MetricFieldType> &msh, const MshCavity& cav, CavOprOpt &opts, CavWrkArrs &work,
                     int nedg0, double *qumin, int ithread = 0);
                     //,int nnwedg,const int* lnwedg,int mnwfac,int *nnwfac,int lnwfac[]);

template<class MetricFieldType, int ideg>
int crenewfa(Mesh<MetricFieldType> &msh, const MshCavity& cav, CavOprOpt &opts, CavWrkArrs &work,
              int ient1, int ied, int iref ,
              bool check_val, bool check_qua,
              int nedg0, int nfac0, double* qmax, int ithread);

// Triangles, likewise. 
template<class MetricFieldType, int ideg>
int reconnect_tetcav(Mesh<MetricFieldType> &msh, const MshCavity& cav, CavOprOpt &opts,  
                     double *qumin, int ithread = 0);//,int nnwfac,const int* lnwfac,
//	                                                              int mnwtet,int *nnwtet,int lnwtet[]);

// After first reconnection round, check if cavity is valid or can be made valid with quick corrections
// If not, tag facets supporting invalid elements as "bad", proposing them for extension.
// Note: if ideg = 1, this leads to the usual cavity correction by extension through faces yielding
// invalid elements.
template<class MetricFieldType,int ideg>
int correct_cavity_fast(Mesh<MetricFieldType> &msh, 
                        MshCavity &cav, CavOprOpt &opts, 
                        int npoi0, int nedg0, int nfac0, int nele0,
                        intAr2& lbad, CavWrkArrs &work, int ithread);
template<class MetricFieldType, int gdim, int ideg>
int correct_cavity_fast0(Mesh<MetricFieldType> &msh, 
                         MshCavity &cav, CavOprOpt &opts, 
                         int npoi0, int nedg0, int nfac0, int nele0,
                         intAr2& lbad,CavWrkArrs &work,int ithread);





}// End namespace
#endif