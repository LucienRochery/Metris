//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MESH_BASE__
#define __METRIS_MESH_BASE__


#include "../types.hxx"
#include "../metris_constants.hxx"
#include "../Mesh/CADInfo.hxx"
#include "../aux_topo.hxx"

#include <boost/pool/poolfwd.hpp>
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>


namespace Metris{

class MetricFieldFE;
class MetricFieldAnalytical;
class MetrisAPI; 
struct MetrisParameters; 



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
  int tet2fac(int ielem, int ifal);

	intAr2  fac2poi; 
	intAr1  fac2ref;
	intAr2r fac2tag;
  // Return global edge at face edge if exists, -1 otherwise
  int fac2edg(int iface, int iedl);

	intAr2  edg2poi; 
	intAr1  edg2ref;
	intAr2r edg2tag; 

	intAr2  poi2ent;
	intAr2r poi2tag; 
	dblAr2  coord;
  // Boundary link 
  intAr1  poi2bpo; 
  intAr2  bpo2ibi; 
  intAr2r bpo2tag;
  dblAr2  bpo2rbi;

  // Seek a bpo for ipoin of dim tdim that matches either ientt or, if not, then iref
  // either can be -1 (or both, then return first dim matching ent)
  int poi2ebp(int ipoin, int tdim, int ientt, int iref) const;

	// Work arrays: to be freely used by any routine
  intAr1 iwork;
  dblAr1 rwork;



	int tag[METRIS_MAXTAGS];

	// Walk linked list instead. Can use getbpos in low_topo. 
	//// Indexed by unique pairs (iface,ivert) or (iedge,ivert)
	//// Yielding 
	//absl::flat_hash_map<std::tuple<int,int>, int> fac2bpo;
	//absl::flat_hash_map<std::tuple<int,int>, int> edg2bpo;

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
  HshTabInt2 edgHshTab; 
  // Similarly, used to point to a triangle from a tetrahedron. 
  HshTabInt3 facHshTab;

	int curdeg,strdeg; 
  int idim;


  // Different subclasses will need different values. 
  // For now it's only MeshBack that needs 10.
  const int nipwk, niewk, nifwk, nitwk;
  const int nrpwk;
	

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

	MeshBase() = delete;
	MeshBase(int nipwk_=1, int niewk_=1, int nifwk_=1, int nitwk_=1, int nrpwk_=1);
	virtual ~MeshBase(){}

  const int &nbpoi = nbpoi_;
  const int &npoin = npoin_;
  const int &nedge = nedge_;
  const int &nface = nface_;
  const int &nelem = nelem_;

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
  void copyConstants(const MeshBase &msh, int MAX_DEG = METRIS_MAX_DEG);

//protected:
	unsigned long long int  getMemCost(); // in bytes
	void zeroArrays();

	void readMeshFile(int64_t libmeshbIdx, int ithread = 0);
  // This destroys the data 
  void readMeshData(MetrisAPI &data);

	void iniNeighbours();
	void iniBdryPoints(int ithread = 0);
  void iniCADLink(int nbpo0);

  /* END INIT */


public:

	MeshBase &operator=(const MeshBase &msh);

	void setBasis(FEBasis ibasis_);
	FEBasis getBasis() const {return ibasis;}
  // To be used in initialization
  void forceBasisFlag(FEBasis ibasn){ibasis = ibasn;}

  int get_tdim() const{
    if(nelem > 0) return 3;
    else if (nface > 0) return 2;
    return 1;
  }

  int get_gdim() const{return idim;}

  void set_gdim(int idim){
    METRIS_ASSERT(idim == 2 || idim == 3);
    this->idim = idim;
  }


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
  template<int tdimn>       intAr2& ent2tag();
  template<int tdimn> const intAr2& ent2tag() const;
        intAr2&  ent2ent(int tdimn);
  const intAr2&  ent2ent(int tdimn) const;
  template<int tdimn>       intAr2& ent2ent();
  template<int tdimn> const intAr2& ent2ent() const;

  // return edgHshTab (tdimn = 1) or facHshTab reference 
  template<int tdimn> typename 
              std::conditional<tdimn==1,HshTabInt2,HshTabInt3>::type & hshTab();


	// Flag whether being on an edge or triangle makes us an ibpoi
	// Instead of checking everywhere if msh.nface > 0 || msh.nelem > 0
	// And in case we want to change this, keep it in one place. 
	bool isboundary_edges()const{return nface > 0 || nelem > 0 || idim > 1;}
	bool isboundary_faces()const{return nelem > 0 || idim > 2;}
  bool isboundary_tdimn(int tdimn)const{return tdimn == 0 ? true : 
                                               tdimn == 1 ? isboundary_edges() :
                                               tdimn == 2 ? isboundary_faces() : false;}

  // Lowest point topological dimension
  int getpoitdim(int ipoin) const;

  template<int ideg> 
  int getverent(int ientt, int tdimn, int ipoin);
  int getverent(int ientt, int tdimn, int ipoin);


  // Whenever some nbpoi, npoin, nedge, nface, nelem is incremented, call these guys beforehand
  // They will reallocate if needed, and throw an exception if impossible. 
  // Flag skipallocf determines whether main data (found in file) should be 
  // allocated or not. This is for initialization from API, to use std::move. 
  virtual void set_nbpoi(int nbpoi);
  virtual void set_npoin(int npoin, bool skipallocf = false);
  virtual void set_nedge(int nedge, bool skipallocf = false);
  virtual void set_nface(int nface, bool skipallocf = false);
  virtual void set_nelem(int nelem, bool skipallocf = false);

  virtual void set_nentt(int tdimn, int nentt, bool skipallocf = false);


protected:
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
	template<int tdim>
	int newbpotopo(int ipoin, int ientt = -1);

  void killpoint(int ipoin);

  // Remove all tagged entities from ipoin
  void rembpotag(int ipoin, int ithread = 0);

	// Check if local edge iedl of iface is a global edge, and return index (or < 0)
	int facedg2glo(int iface, int iedl) const;
	// Check if local face ifal of ielem is a global face, and return index (or < 0)
	int tetfac2glo(int ielem, int ifal) const;
	// Check if local edge iedl of ielem is a global edge, and return index (or < 0)
	int tetedg2glo(int ielem, int iedl) const;

  // Return computed geodev (see below)
  double get_geodev(int tdim) const {return geodev[tdim-1];}
  //// return surface normal for surface iref 
  //// averaged
  //// bad for periodic but good enough for first approximation
  //void getpoinormal(int ipoin, int iref, double *nrmal);


  MetrisParameters* param;


protected:
  FEBasis ibasis;


	void getEnttMemCosts(int *memCostPpoi, int *memCostPbpo, int *memCostPedg, int *memCostPfac, int *memCostPelt) const;
	// From maximum point count, deduce maximum other entity count
	void setMpoiToMent();

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


} // End namespace

#endif
