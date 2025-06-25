//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MESH_BASE__
#define __METRIS_MESH_BASE__


#include "../types.hxx"
#include "../metris_constants.hxx"
#include "../Mesh/CADInfo.hxx"
#include "../aux_topo.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include <boost/pool/poolfwd.hpp>
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>


namespace Metris{

class MetricFieldFE;
class MetricFieldAnalytical;
class MetrisAPI; 
class MeshBase;

/* N Integer Boundary Information: stride of bpo2ibi */                   
const int nibi = 4; 

/* N Real Boundary Information: stride of bpo2rbi */
// 0-2 U-V coordinates
// Normals for edge / corner? 
const int nrbi = 2; 

enum class MeshClass{MeshBase, MeshMetric, MeshBack, Mesh};
enum class MetricClass;



class MeshBase{
public: 
  friend class MetricFieldFE;
  friend class MetricFieldAnalytical;

	/* Entity types and definitions:
	- Vertices (poi) and boundary points (bpo): 
		- Points (generic ipoin) index into:
				coord,poi2ref,poi2tag[0],poi2bpo
			and are indexed from 
			  all ent2poi arrays
			Points are uniquely defined by their coordinates. 
		- Boundary points (generic ibpoi) index into:
		    bpo2ibi,bpo2rbi
		  and are indexed from 
		    poi2bpo,fac2bpo,edg2bpo
			Boundary points are uniquely defined by their (u,v) found in bpo2rbi. 
			They point to ipoin, iref through bpo2ibi. This is not unique. 
			Distinguishing between two ibpois is required by edges and triangles to compute 
			new (u,v) coordinates on e.g. an edge. 
			This is to handle the case of periodic mappings (cylinder) or even simply the (u,v) from two
			CAD faces meeting at an edge (and we also need the edge (u,1-u)). 
			To obtain the appropriate ibpoi from a triangle, use the hash table fac2bpo. 
	- Navigating point <-> boundary point links:
   - poi2bpo[ipoin] = ibpoi attached to lowest-dimensional CAD entity
	  - bpo2ibi(ibpoi,0) = ipoin 
	  - bpo2ibi(ibpoi,1) = tdimn: 2 face int    -| this may be replaced by an ego** array    
	                               1 edge int     | though perhaps the type should be kept 
	                               0 node/corner  | without indirection
	  - bpo2ibi(ibpoi,2) = ientt (mesh ent idx) -|  -> CAD ref can be fetched from (edg|fac)2ref
	  - bpo2ibi(ibpoi,3) = next ibpoi sharing same ipoin
	- Tetrahedra:
	  - tet2poi(ielem,i) = ipoin
	- Triangles: surface only. 
		- tri2poi(iface,i) = ipoin
		- When we need the (u,v), get (ibpo1, ibpo2, ibpo3) in hash table fac2bpo
		- This is better than storing directly the ibpoi because most uses of triangles need to avoid the indirection
			and are better off with tri2poi(iface,i) than bpo2ibi[tri2bpo[iface][i]][0]
		- Only triangle edges that neighbour different ref neighbours or CAD edges need this information
	- Edges: idem

	- Everything is 0-based including refs. 

	-poi2ent:  lowest-dim entity attached. Check poi2bpo: 
	  - if ibpoi < 0, volume point, poi2ent is a tetra
	  - if ibpoi == 2, surface interior point, poi2ent is a triangle. Get tetras using fac2tet.
		- if ibpoi == 1, edge interior point, poi2ent is an edge. Get triangle using edg2fac. 
		- if ibpoi == 0, corner point. Should practically never happen but poi2ent can be -1 in this case. 
			Otherwise, it points to an edge. 
    -> UPDATE: Now intAr2, poi2ent(ipoin,1) gives topo dimn
	- tet2ftg (tet to face tag): 0 if no face attached, 1 otherwise. To speed up ball checks in cavity.
	*/
	intAr2  tet2poi; 
	intAr1  tet2ref;
	intAr2r tet2tag;
	bolAr1  tet2ftg; 
  intAr2  dom2tag; // tetra domain tags
  // Return global face at tet face if exists, -1 otherwise
  int tet2fac(int ielem, int ifal);
  const int &nelem = nelem_;
  int ndomn; // number of domain refs (some may be empty)

	intAr2  fac2poi; 
	intAr1  fac2ref;
	intAr2r fac2tag;
  // Return global edge at face edge if exists, -1 otherwise
  int fac2edg(int iface, int iedl);
  const int &nface = nface_;

	intAr2  edg2poi; 
	intAr1  edg2ref;
	intAr2r edg2tag; 
  const int &nedge = nedge_;

	intAr2  poi2ent;
	intAr2r poi2tag; 
	dblAr2  coord;
  const int &npoin = npoin_;

  // Boundary link 
  intAr1  poi2bpo; 
  intAr2  bpo2ibi; 
  intAr2r bpo2tag;
  dblAr2  bpo2rbi;
  const int &nbpoi = nbpoi_;

  // point to element boundary point
  // Seek a bpo for ipoin of dim tdim that matches either ientt or, if not, 
  // then iref either can be -1 (or both, then return first dim matching ent)
  int poi2ebp(int ipoin, int tdim, int ientt, int iref) const;

  // point to dimension element
  // Give an entity of dimension tdim that contains ipoin and has ref iref 
  // if iref >= 0.
  int poi2del(int ipoin, int tdim, int iref) const;



	int tag[METRIS_MAXTAGS];

  // tet -> tet neighbours in tet2tet
  // edg -> edg neighbours in edg2edg
  // tet -> edg and fac -> edg neighbours in edgHshTab:
  //    simply search for ip1, ip2 in hash tab, get iedge
  // tet -> fac neighbours in facHshTab
  //    idem, search ip1,ip2,ip3 in tab, get iface
  // fac -> tet in fac2tet (holds both)
  // edg -> fac in edg2fac (holds one). There can be arbitrarily many faces around a given edge. 
  // edg -> tet:
  //   if ever needed: edg -> fac through edg2fac then fac2tet then shell. Or if HO, use HO node
  // then poi2ent. 
	intAr2 edg2edg, // 2xnedge ; >= 0 edge neighbour ; -1 no nei; < -1 nm edge neighbour
	       fac2fac, // 3xnface ; idem
	       fac2tet, // 2xnface ; [][0] first adjacent tet (e.g. bdry fac) [][1] second (internal face)
	       tet2tet; // 4xnelem ; >= 0 neighbour ; -1 no neighbour ; <-1 invalid
	intAr1 edg2fac;
	// Geometric edges (init as read from file)
	// Store ip1,ip2,iedge. This is used to point to an edge from a triangle or tetrahedron
  HshTab_I2I edgHshTab; 
  // Similarly, used to point to a triangle from a tetrahedron. 
  HshTab_I3I facHshTab;

	int curdeg,strdeg; 
  int idim;


	// This is intended to be kept throughout execution and 
	// passed to e.g. the cavity operator as needed. 
	// Nothing inside of this allocator is kept valid from
	// one call to another. 
	boost::fast_pool_allocator<int>    iwrkmem;
	boost::fast_pool_allocator<double> rwrkmem;

	// [x,y,z][min,max]
	double bb[3][2];

	// Some flags to set and communicate with lower level routines
	// without having to add arguments. 
	int idbg[10];


  CADInfo CAD;
  intAr2 cfa2tag, ced2tag, cno2tag;

  virtual MeshClass meshClass() const { return MeshClass::MeshBase; }
  virtual MetricClass metricClass() const;

	MeshBase();

  virtual ~MeshBase(){}


  const int &mbpoi = mbpoi_;
  const int &mpoin = mpoin_;
  const int &medge = medge_;
  const int &mface = mface_;
  const int &melem = melem_;

	/* START INIT */

public:
//protected:
  void readConstants(int64_t libIdx, int usrMinDeg);
	void readConstants(const MetrisAPI &data, int usrMinDeg);
public:
  void copyConstants(const MeshBase &msh);

//protected:
	unsigned long long int  getMemCost(); // in bytes
	void zeroArrays();

	void readMeshFile(int64_t libmeshbIdx, int ithread);
  // This destroys the data 
  void readMeshData(MetrisAPI &data);

	void iniNeighbours();
  // Returns nbpo0 s.t. ibpoi <= nbpo0 needs no projection. 
	int iniBdryPoints(int ithread);
  void iniCADLink(int nbpo0);

  /* END INIT */


public:

	MeshBase &operator=(const MeshBase &msh);

	void setBasis(FEBasis ibasis_);
	FEBasis getBasis() const {return ibasis;}
  // To be used in initialization
  void forceBasisFlag(FEBasis ibasn){ibasis = ibasn;}

  // Return highest topo dim in mesh
  int get_tdim() const;

  // Dimension generic helpers
              int  nentt(int tdimn) const;
              //int &nentt(int tdimn);
              int  mentt(int tdimn) const;
              //int &mentt(int tdimn);
              int  nnode(int tdimn) const;
        intAr2&  ent2poi(int tdimn);
  const intAr2&  ent2poi(int tdimn) const;
  template<int tdimn>       intAr2& ent2poi();
  template<int tdimn> const intAr2& ent2poi() const;
        intAr1&  ent2ref(int tdimn);
  const intAr1&  ent2ref(int tdimn) const;
        intAr2r& ent2tag(int tdimn);
  const intAr2r& ent2tag(int tdimn) const;
  template<int tdimn>       intAr2r& ent2tag();
  template<int tdimn> const intAr2r& ent2tag() const;
        intAr2&  ent2ent(int tdimn);
  const intAr2&  ent2ent(int tdimn) const;
  template<int tdimn>       intAr2& ent2ent();
  template<int tdimn> const intAr2& ent2ent() const;

  // points to ced2tag, cfa2tag or dom2tag
        intAr2r& ref2tag(int tdimn);
  const intAr2r& ref2tag(int tdimn) const;

  // return edgHshTab (tdimn = 1) or facHshTab reference 
  template<int tdimn> typename 
              std::conditional<tdimn==1,HshTab_I2I,HshTab_I3I>::type & hshTab();


	// Flag whether being on an edge or triangle makes us an ibpoi
	// Instead of checking everywhere if msh.nface > 0 || msh.nelem > 0
	// And in case we want to change this, keep it in one place. 
	bool isboundary_edges()const{return idim >= 2;}
	bool isboundary_faces()const{return idim >= 3;}
  bool isboundary_tdim(int tdimn)const{return tdimn == 0 ? true : 
                                              tdimn == 1 ? isboundary_edges() :
                                              tdimn == 2 ? isboundary_faces() : false;}

  // Lowest point topological dimension
  int getpoitdim(int ipoin) const;

  // Get ipoin local index in entity of dim tdimn
  template<int ideg> 
  int getverent(int ientt, int tdimn, int ipoin) const;
  int getverent(int ientt, int tdimn, int ipoin) const;

  // Get ipoin local index in edge, face, tetra
  template<int ideg> 
  int getveredg(int ientt, int ipoin) const;
  template<int ideg> 
  int getverfac(int ientt, int ipoin) const;
  template<int ideg> 
  int getvertet(int ientt, int ipoin) const;

  // Mostly for internal use, compute the localization alignment direction
  // (edge tangent, face normal) of a point with CAD link and poi2ent initialized.
  void get_algnd(int ipoin, double *algnd);

protected:
  MeshArray1D<intAr1> iwork;
  MeshArray1D<dblAr1> rwork;
  intAr1 iwork_lock, rwork_lock;

public:
  // Request rwork/iwork sized nn
  dblWrkAr1 get_rwork(int nn = 0);
  intWrkAr1 get_iwork(int nn = 0);
  // No need to call these.
  template<typename T>
  void free_work(int ii){
    static_assert(std::is_same<T,double>::value || std::is_same<T,int>::value);
    if(std::is_same<T,double>::value){
      rwork_lock[ii] = 0;
    }else if(std::is_same<T,int>::value){
      iwork_lock[ii] = 0;
    }
  }


  // Flag skipallocf determines whether main data (found in file) should be 
  // allocated or not. This is for initialization from API, to use std::move. 
  virtual void set_nbpoi(int nbpoi);
  virtual void set_npoin(int npoin, bool skipallocf = false);
  virtual void set_nedge(int nedge, bool skipallocf = false);
  virtual void set_nface(int nface, bool skipallocf = false);
  virtual void set_nelem(int nelem, bool skipallocf = false);
  virtual void set_nentt(int tdimn, int nentt, bool skipallocf = false);


//protected:
public:
  // This should only be called from the top levels (Mesh and MeshBack)
  // Otherwise some auxiliary data structs may not be properly set. 
	int newpoitopo(int tdimn, int ientt = -1);
  friend void debugInveval(std::string meshName_, MeshBase &msh, int tdim,  int* ent2pol, double *coop);

public:
	// Create new face by copying from tetrahedron
	template <int ideg>
	void newfactopo(int ielem, int ifael, int iref = -1, int iele2 = -1);
  // Create virtual face or edge
  // Used in the cavity operator to give boundary information to a new point.
  int newfacvirtual(int iref);
  int newedgvirtual(int iref);


	template <int ideg>
	void newedgtopo(int iface, int iedfa, int iref = -1);

	int newbpotopo(int ipoin, int tdim, int ientt = -1);
  void killpoint(int ipoin);

  // Remove all tagged entities from ipoin
  void rembpotag(int ipoin, int ithread);

	// Check if local edge iedl of iface is a global edge, and return index (or < 0)
	int facedg2glo(int iface, int iedl) const;
	// Check if local face ifal of ielem is a global face, and return index (or < 0)
	int tetfac2glo(int ielem, int ifal) const;
	// Check if local edge iedl of ielem is a global edge, and return index (or < 0)
	int tetedg2glo(int ielem, int iedl) const;

  // Return computed geodev (see below)
  double get_geodev(int tdim) const {return geodev[tdim-1];}


  MetrisParameters* param;


protected:
  FEBasis ibasis;

	void setLagrange();
	void setBezier();

  // Main init routine: pass a NULL data to read from param instead. 
  // Note this shouldn't be called manually, use a MetrisRunner. 
  void initialize(MetrisAPI *data, MetrisParameters &param);
  void iniFromFile(std::string fname, int usrTarDeg);
  void iniFromData(MetrisAPI &data, int usrTarDeg);

  int npoin_,nbpoi_,nedge_,nface_,nelem_;
  int mpoin_,mbpoi_,medge_,mface_,melem_;

  // Store maximum deviation between surface directions and element directions:
  // 1: edges, tangent 
  // 2: faces, normal
  // Deviation is 1 - abs(dtprd) 
  double geodev[2];  // Also in back mesh... 
};



template<typename T>
class WorkArray1D{
public:
  WorkArray1D(MeshBase& msh_, int ilock_, MeshArray1D<T>& array_) 
    : msh(msh_), ilock(ilock_), array(array_) {}

  ~WorkArray1D(){
    msh.free_work<T>(ilock);
  }

  ALWAYS_INLINE T &operator[](const int &ii){
    return array[ii];
  }
  ALWAYS_INLINE const T &operator[](const int &ii) const {
    return array[ii];
  }

  void set_n(int nn){array.set_n(nn);}
  int get_n() const {return array.get_n();}
  void stack(T val){array.stack(val);}
  T pop(){return array.pop();}
  bool allocate(int m){return array.allocate(m);}

  MeshArray1D<T>& get_array(){return array;}
  const MeshArray1D<T>& get_array() const {return array;}

private:
  int ilock, itype;
  MeshArray1D<T>& array;
  MeshBase& msh;
};


} // End namespace

#endif
