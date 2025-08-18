//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MESH_BASE__
#define __METRIS_MESH_BASE__


#include "CADInfo.hxx"
#include "WorkArray.hxx"

#include "../types.hxx"
#include "../metris_constants.hxx"
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
class TagArray;

/* N Integer Boundary Information: stride of bpo2ibi */                   
const int nibi = 4; 

/* N Real Boundary Information: stride of bpo2rbi */
// 0-2 U-V coordinates
// Normals for edge / corner? 
const int nrbi = 2; 

enum class MeshClass{MeshBase, MeshMetric, MeshBack, Mesh};
enum class MetricClass;

// For tag and work arrays. Only positive values can be used, 0 is "not tagged".
enum class MeshSize{
   Untracked = 1 // internal use
  ,Point = 2
  ,Edge = 3
  ,Face = 4
  ,Tetra = 5
  ,BPoint = 6
  ,Domain = 7 // tetra domain
  ,CADNode = 8
  ,CADEdge = 9
  ,CADFace = 10
  ,nTrackedType = 11
};

class MeshBase{
public: 
  friend class MetricFieldFE;
  friend class MetricFieldAnalytical;
  template<typename T>
  friend class WorkArray1D;
  friend class TagArray;

  // ---- Connectivity arrays:
	intAr2  edg2poi, fac2poi, tet2poi; // (ientt,inode) -> ipoin
	int curdeg,strdeg; // current mesh degree / storage degree (internal)
  int nnode(int tdim) const; // can also use getnnod1/2/3 or getnnode which are constexpr, see ho_constants.hxx
  // - Entity counts:
  const int &nedge = nedge_;
  const int &nface = nface_;
  const int &nelem = nelem_;
  int nentt(int tdim) const;
  // - Dimension-generic helpers:
        intAr2&  ent2poi(int tdim);
  const intAr2&  ent2poi(int tdim) const;
  template<int tdim>       intAr2& ent2poi();
  template<int tdim> const intAr2& ent2poi() const;


  // ---- Vertices
  const int &npoin = npoin_;
	dblAr2 coord;
  int idim; // geometric dimension
	intAr2  poi2ent; // 2xnpoin: 0 = lowest-dim element attached, 1 = dimension of element
  int getpoitdim(int ipoin) const;  // Point topological dimension (lowest dim attached element)
  bolAr1  poicstr; // is point constrained (internal)
	double bb[3][2]; // bounding box: [x,y,z][min,max]
  // - ipoin -> local
  template<int ideg>  int getveredg(int ientt, int ipoin) const;
  template<int ideg>  int getverfac(int ientt, int ipoin) const;
  template<int ideg>  int getvertet(int ientt, int ipoin) const;
  template<int ideg>  int getverent(int ientt, int tdim, int ipoin) const;
                      int getverent(int ientt, int tdim, int ipoin) const;


  // ---- Neighbour arrays. i-th entry is element opposite i-th vertex (edges included)
	intAr2 edg2edg, // 2xnedge: >= 0 edge neighbour ; -1 no nei; < -1 nm edge neighbour
	       fac2fac, // 3xnface: idem
	       fac2tet, // 2xnface: [*][0] first adjacent tet (e.g. bdry fac) [*][1] second (internal face)
	       tet2tet; // 4xnelem: >= 0 neighbour ; -1 no neighbour ; <-1 invalid
	intAr1 edg2fac; //   nedge: any adjacent face
  // - Dimension-generic helpers:
        intAr2&  ent2ent(int tdim);
  const intAr2&  ent2ent(int tdim) const;
  template<int tdim>       intAr2& ent2ent();
  template<int tdim> const intAr2& ent2ent() const;
  // - Lower-dimensional neighbours:
  int fac2edg(int iface, int iedl); // Face local edge to global edge (-1 if none)
  int tet2fac(int ielem, int ifal); // Tet local face to global face (-1 if none)


  // ---- Boundary information
  CADInfo CAD; // CAD object
	intAr1 edg2ref, fac2ref, tet2ref; // Surface/domain refs, map to CAD.cad2fac/edg. 
  int ndomn; // number of tetra refs. Others are found in CAD.ncadfa, CAD.ncaded
  bool is_nonmanifold() const {return !is_manifold;} // is the domain manifold?
  // - Dimension-generic helpers:
        intAr1&  ent2ref(int tdim);
  const intAr1&  ent2ref(int tdim) const;
  // - Internal use (interpMetBack)
  // Compute localization alignment direction (edge tangent, face normal)
  void get_algnd(int ipoin, double *algnd);


  // ---- Vertex boundary links: ipoin -> (u,v) or t CAD coordinates.
  // Each boundary vertex has a linked list of t / (u,v) coordinates and seed information. 
  // The entries in a given linked list are in ascending order of tdim. 
  // To handle periodic patches/curves:
  // - a dimension 0 point has an ibpoi per edge
  // - a dimension <= 1 point has an ibpoi per face
  // This holds regardless of whether the surface is in fact periodic. In these cases, a vertex can have 
  // two distinct t or (u,v) coordinates for a single CAD edge/face ref. 
  intAr1  poi2bpo; // ipoin -> ibpoi (-1 if interior)
  intAr2  bpo2ibi; // ibpoi -> (ipoin, tdim, ientt, next ibpoi)
  dblAr2  bpo2rbi; // ibpoi -> t / (u,v) coordinates
  const int &nbpoi = nbpoi_;
  // - Helpers
  // Seek a ibpoi for ipoin of dim tdim that matches either ientt or, if not, then iref.
  // Either (but not both) can be -1. 
  int poi2ebp(int ipoin, int tdim, int ientt, int iref) const; // POInt to Element Boundary Point


  // ---- Utilities
  // - Work arrays can be manipulated like regular MeshArray1D. This tries to avoid unnecessary allocations. 
  // The memory is automatically "freed" (not deallocated but allowed to be reused) when the objects are destroyed.
  // Shortcut names are intWrkAr1 and dblWrkAr1
  template<typename T>
  WorkArray1D<T> get_work(int nn = 0);
  template<typename T>
  WorkArray1D<T> get_work(MeshSize nentt);

  // Proxies for calling get_work with <int> or <double>
  intWrkAr1 get_iwork(int nn);
  dblWrkAr1 get_rwork(int nn);
  intWrkAr1 get_iwork(MeshSize size);
  dblWrkAr1 get_rwork(MeshSize size);


  TagArray get_tagarray(MeshSize itype);

  int tag[METRIS_MAXTAGS];
	intAr2r tet2tag;
  intAr2  dom2tag; // tetra domain tags
	intAr2r fac2tag;
	intAr2r edg2tag; 
	intAr2r poi2tag; 
  intAr2r bpo2tag;
  intAr2 cfa2tag, ced2tag, cno2tag;

        intAr2r& ent2tag(int tdimn);
  const intAr2r& ent2tag(int tdimn) const;
  template<int tdimn>       intAr2r& ent2tag();
  template<int tdimn> const intAr2r& ent2tag() const;
  // points to ced2tag, cfa2tag or dom2tag
        intAr2r& ref2tag(int tdimn);
  const intAr2r& ref2tag(int tdimn) const;


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
	// Geometric edges (init as read from file)
	// Store ip1,ip2,iedge. This is used to point to an edge from a triangle or tetrahedron
  HshTab_I2I edgHshTab; 
  // Similarly, used to point to a triangle from a tetrahedron. 
  HshTab_I3I facHshTab;




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


  // return edgHshTab (tdim = 1) or facHshTab reference 
  template<int tdim> typename 
              std::conditional<tdim==1,HshTab_I2I,HshTab_I3I>::type & hshTab();


	// Flag whether being on an edge or triangle makes us an ibpoi
	// Instead of checking everywhere if msh.nface > 0 || msh.nelem > 0
	// And in case we want to change this, keep it in one place. 
	bool isboundary_edges()const{return idim >= 2;}
	bool isboundary_faces()const{return idim >= 3;}
  bool isboundary_tdim(int tdim)const{return tdim == 0 ? true : 
                                             tdim == 1 ? isboundary_edges() :
                                             tdim == 2 ? isboundary_faces() : false;}


protected:
  MeshArray1D<intAr1> iwork;
  MeshArray1D<dblAr1> rwork;
  intAr1 iwork_lock, rwork_lock;

  // Track reference to underlying MeshArray1D
  // of WorkArrays that track mesh sizes. 
  // This is so we can reallocate them. 
  MeshArray1D<intWrkAr1*> iwork_Point, 
                          iwork_Edge,
                          iwork_Face,
                          iwork_Tetra,
                          iwork_BPoint;
  MeshArray1D<dblWrkAr1*> rwork_Point, 
                          rwork_Edge,
                          rwork_Face,
                          rwork_Tetra,
                          rwork_BPoint;  
  // Returns arrays iwork_Point, etc. 
  MeshArray1D<intWrkAr1*> &get_iwork_tracked(MeshSize itype);
  MeshArray1D<dblWrkAr1*> &get_rwork_tracked(MeshSize itype);
  template<typename T> 
  auto& get_work_tracked(MeshSize itype){
    static_assert(std::is_same_v<T,int> || std::is_same_v<T,double>);
    if constexpr(std::is_same_v<T,int>) return get_iwork_tracked(itype);
    else                                return get_rwork_tracked(itype);
  }


  // Count how many currently tracked, so we can reset the
  // arrays (i|r)work_Point, etc. 
  int n_iwork_tracked[(int) MeshSize::nTrackedType];
  int n_rwork_tracked[(int) MeshSize::nTrackedType];

  template<typename T> 
  constexpr auto n_work_tracked(){
    if constexpr (std::is_same_v<T,int>) return n_iwork_tracked;
    else                                 return n_rwork_tracked;
  }
public: // for debug purposes only
  // We are made to give it a unique name, otherwise the 
  // different scopes make this const version invisible.
  template<typename T> 
  constexpr const int* debug_n_work_tracked() const {
    if constexpr (std::is_same_v<T,int>) return n_iwork_tracked;
    else                                 return n_rwork_tracked;
  }
  template<typename T> 
  auto& debug_get_work_tracked(MeshSize itype){
    static_assert(std::is_same_v<T,int> || std::is_same_v<T,double>);
    if constexpr(std::is_same_v<T,int>) return get_iwork_tracked(itype);
    else                                return get_rwork_tracked(itype);
  }

protected:
  template<typename T> 
  MeshArray1D<MeshArray1D<T>>& lwork_map(){
    static_assert(std::is_same<T,double>::value || std::is_same<T,int>::value);
    if constexpr(std::is_same<T, double>::value){
      return rwork;
    }else{
      return iwork;
    }
  }

  template<typename T> 
  intAr1& lwork_lock_map(){
    static_assert(std::is_same<T,double>::value || std::is_same<T,int>::value);
    if constexpr(std::is_same<T, double>::value){
      return rwork_lock;
    }else{
      return iwork_lock;
    }
  }

  // For untracked work
  template<typename T>
  void free_work(int ii);

  // For work arrays tracking mesh sizes
  template<typename T>
  void free_work(int ii, WorkArray1D<T>* obj);

  void free_tag(int ii){
    tagarr_locks[ii] = 0;
  }

  void update_tracked_work_arrays(MeshSize itype, int mentt, int nentt);

protected:
  MeshArray1D<intAr1> tagarrs;
  intAr1 tagarr_locks;
  intAr1 itags;

public:

  // Flag skipallocf determines whether main data (found in file) should be 
  // allocated or not. This is for initialization from API, to use std::move. 
  virtual void set_nbpoi(int nbpoi);
  virtual void set_npoin(int npoin, bool skipallocf = false);
  virtual void set_nedge(int nedge, bool skipallocf = false);
  virtual void set_nface(int nface, bool skipallocf = false);
  virtual void set_nelem(int nelem, bool skipallocf = false);
  virtual void set_nentt(int tdim, int nentt, bool skipallocf = false);


//protected:
public:
  // This should only be called from the top levels (Mesh and MeshBack)
  // Otherwise some auxiliary data structs may not be properly set. 
	int newpoitopo(int tdim, int ientt = -1);
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

  // Store whether the mesh is non-manifold or not 
  bool is_manifold;
public:
  MeshArray1D<bool> isperiodic_face;
  int nperiodic_face; // Number of periodic faces

// --- DEBUG functions
public:
  // These are only used in a unit test:
  const intAr1 &debug_get_iwork_lock() const {return iwork_lock;}
  const intAr1 &debug_get_rwork_lock() const {return rwork_lock;}
  const MeshArray1D<intAr1> &debug_get_iwork() const {return iwork;}
  const MeshArray1D<dblAr1> &debug_get_rwork() const {return rwork;}
  template<typename T>
  const auto &debug_get_work() const {
    if constexpr (std::is_same_v<T,int>) return iwork;
    else                                 return rwork;
  }
  template<typename T>
  const auto &debug_get_work_lock() const {
    if constexpr (std::is_same_v<T,int>) return iwork_lock;
    else                                 return rwork_lock;
  }
  intAr1 &debug_get_iwork_lock() {return iwork_lock;}
  intAr1 &debug_get_rwork_lock() {return rwork_lock;}
  MeshArray1D<intAr1> &debug_get_iwork() {return iwork;}
  MeshArray1D<dblAr1> &debug_get_rwork() {return rwork;}
  template<typename T>
  auto &debug_get_work(){
    if constexpr (std::is_same_v<T,int>) return iwork;
    else                                 return rwork;
  }
  template<typename T>
  auto &debug_get_work_lock(){
    if constexpr (std::is_same_v<T,int>) return iwork_lock;
    else                                 return rwork_lock;
  }
private:
  void get_nMeshSize(MeshSize itype, int* nn, int* mm);
};



class TagArray{
public:
  friend class MeshBase;

  TagArray(MeshBase& msh_, int ilock_, intAr1& array_, int& itag_) 
    : ilock(ilock_), array(array_), msh(msh_), itag(itag_), itag1(itag_){
      #ifndef NDEBUG
      itag0 = itag;
      #endif
      if(itag > INT_MAX - 1000){
        itag = 1;
        array.fill(0);
      }
    }

  ~TagArray(){
    itag = MAX(itag1, itag);
    msh.free_tag(ilock);
  }

  // increment tag
  void operator++(){
    itag++;
    METRIS_ASSERT_MSG(itag < itag0 + 1000,"Single tag limits need to be revised");
  }
  // Check element is tagged
  bool is_tagged(int ientt, int ivalue) const { return array[ientt] >= itag + ivalue; }
  bool is_tagged(int ientt) const { return is_tagged(ientt,0); }

  int get_tag(int ientt) const { return array[ientt] - itag; }

  // Set tag value of element
  void tag(int ientt, int value){
    array[ientt] = itag + value; 
    itag1 = MAX(itag1, itag + value);
  }

  void tag(int ientt){ tag(ientt, 0); }

  // Get entity tag. Should not be used often. 
  ALWAYS_INLINE int &operator[](const int &ii){ return array[ii]; }
  ALWAYS_INLINE const int &operator[](const int &ii) const { return array[ii]; }

private:
  int ilock;
  intAr1& array;
  MeshBase& msh;
  int &itag, itag1;
  #ifndef NDEBUG
  int itag0;
  #endif

  int get_tag() const{ return itag; }
};


} // End namespace

#endif
