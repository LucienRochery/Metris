//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_INITIALIZER__
#define __METRIS_INITIALIZER__


#include "../types.hxx"
#include "../ho_constants.hxx"
#include "../metris_constants.hxx"
#include "../Mesh/CADInfo.hxx"
#include "../Mesh/MeshFwd.hxx"

namespace Metris{

class MetrisRunner;

/*
The MetrisAPI class implements a "file in memory" to be passed to a MetrisRunner
together with some MetrisParameters to initialize Metris. 
Back/Front mesh distinction is done simply by creating two MetrisAPI objects. 
When read by a MetrisRunner, a MetrisAPI object is invalidated. 
However, it can be refilled using MetrisAPI::initialize(MetrisRunner &run). 
*/

class MetrisAPI{
public:
  friend class MeshBase;
  friend class MeshBack;
  friend class MetrisRunner;

  MetrisAPI();
  // -- Setter mode
  // Later call setX() functions. 
  // Call setElementsOrdering() to set ordering from defaults or custom. 
  // gpoe are VertexOnGeometricX entries 
  MetrisAPI(int idim, int ideg, 
            int ncorn, int ngpoe, int ngpof, 
            int npoin, bool imet, 
            int nedge, int nface, int nelem, 
            FEBasis mshbasis, FEBasis metbasis, MetSpace metspace);
  // Pass everything in. Default ordering. 
  // Disabled because of orderings. 
  //MetrisAPI(int idim, int ideg, 
  //          int ncorn, int ngpoe, int ngpof, 
  //          int npoin, int npmet, 
  //          int nedge, int nface, int nelem,
  //          FEBasis mshbasis, FEBasis metbasis, MetSpace metspace,
  //          const int *lcorn, 
  //          const int *lgpoe, const double *rgpoe, 
  //          const int *lgpof, const double *rgpof, 
  //          const double *coord, const double *metfld, 
  //          int iordering, 
  //          const int *edg2poi, const int *edg2ref, 
  //          const int *fac2poi, const int *fac2ref, 
  //          const int *tet2poi, const int *tet2ref);

  // -- Getter mode
  // Constructor destroys the input MetrisRunner. 
  MetrisAPI(MetrisRunner &run); 


  // Equivalently call this at any point (overwrite) -> can be used to 
  // "resurrect" a MetrisAPI
  void initialize(MetrisRunner &run); 
  // Note the following can be called twice. The arrays can be reallocated over
  // and the flags reset. 
  // This can be helpful if information "trickles in" 
  // EGADS_model can be NULL
  void initialize(int idim, int ideg, 
                  int ncorn, int ngpoe, int ngpof, 
                  int npoin, bool imet, 
                  int nedge, int nface, int nelem, 
                  FEBasis mshbasis, FEBasis metbasis, MetSpace metspace);
  // Can set only the flags, and later setNCorner etc. 
  // EGADS_model can be NULL
  void initialize(int idim, int ideg, bool imet,
                  FEBasis mshbasis, FEBasis metbasis, MetSpace metspace);



  void copyFlags(MetrisAPI *into) const;

  // Change degree, allocate to desired. 
  // IF LESS THAN CURRENT, WILL TRUNCATE!
  void setDegree(int tardeg);


  //~MetrisAPI(); 

  void getConstants(int* idim, int* ideg,
                    int* ncorn, int* npoin, bool* imet, 
                    int* nedge, int* nface, int* nelem,
                    FEBasis* mshbasis, FEBasis* metbasis) const;

  /* Vertices & metrics */
  void setNPoints(int npoin_);
  int  getNPoints() const {return npoin;}

  void setCoord(int ipoin, const double *coord);
  void setCoord(int ipoi1, int ipoi2, const double *coord);
  void setCoord(dblAr2 &&coord);

  void getCoord(int ipoin, double *coord) const;
  void getCoord(int ipoi1, int ipoi2, double *coord) const;

  void setMetric(int ipoin, const double *met);
  void setMetric(int ipoi1, int ipoi2, const double *met);
  void setMetric(dblAr2 &&metfld);

  void getMetric(int ipoin, double *met) const;
  void getMetric(int ipoi1, int ipoi2, double *met) const;


  /* Elements */
  // Call prior to any setElement() or program will not work as intended. 
  // ordering: array dim nnode x (tdimn + 1). nnode distinct entries 
  //           for all inode, sum ordering(inode,:) = ideg 
  //        -> if NULL then Metris internal ordering 
  // Note, only for ideg > 1
  void setElementsOrdering(int tdimn, const int *ordering);
  void setElementsOrdering(int iordering);

  void setNEdges(int nedge_);
  //int  getNEdges() const;
  void setNFaces(int nface_);
  //int  getNFaces() const;
  void setNTetrahedra(int nelem_);
  //int  getNTetrahedra() const;

  void setElement(int tdimn, int ielem, const int* lnode, int iref);
  void setElement(int tdimn, int iele1, int iele2, const int* lnode, const int* lref);
  void setElement(int tdimn, intAr2 &&ent2poi, intAr1 &&ent2ref);

  void getElement(int tdimn, int ielem, int* lnode, int* iref) const;
  void getElement(int tdimn, int iele1, int iele2, int* lnode, int* lref) const;
  void getElementRef(int tdimn, int ielem, int *iref) const;


  /* Corners (geometric nodes) */
  void setNCorners(int ncorn_);
  //int  getNCorners() const;

  void setCorner(int icorn, int ipoin);
  void setCorner(int icor1, int icor2, const int* lpoin);

  void getCorner(int icorn, int *ipoin) const;
  void getCorner(int icor1, int icor2, int* lpoin) const;
  void copyCorners(MetrisAPI *into) const;


  /* Geometric edges */
  void setNVerticesOnGeometricEdges(int ngpoe_);
  //int  getNVerticesOnGeometricEdges() const;

  void setVerticesOnGeometricEdges(int ientry, const int *lgpoe, const double *rgpoe);
  void setVerticesOnGeometricEdges(int ientr1, int ientr2, const int *lgpoe, const double *rgpoe);

  void getVerticesOnGeometricEdges(int ientry, int *lgpoe, double *rgpoe) const;
  void getVerticesOnGeometricEdges(int ientr1, int ientr2, int *lgpoe, double *rgpoe) const;

  void copyVerticesOnGeometricEdges(MetrisAPI *into) const;

  /* Geometric faces */
  void setNVerticesOnGeometricTriangles(int ngpof_);
  //int  getNVerticesOnGeometricTriangles() const;

  void setVerticesOnGeometricTriangles(int ientry, const int *lgpof, const double *rgpof);
  void setVerticesOnGeometricTriangles(int ientr1, int ientr2, const int *lgpof, const double *rgpof);

  void getVerticesOnGeometricTriangles(int ientry, int *lgpof, double *rgpof) const;
  void getVerticesOnGeometricTriangles(int ientr1, int ientr2, int *lgpof, double *rgpof) const;

  void copyVerticesOnGeometricTriangles(MetrisAPI *into) const;

  // Set the EGADS model
  // If context is NULL, then a new context is created, and the model is hard
  // copied.
  // Otherwise, it is assumed the context owning the model will live at least 
  // as long as this API object; the original model is referenced. 
  void setCADModel(ego EGADS_context, ego EGADS_model);
  void copyCAD(MetrisAPI *into) const;


private:

  void free();

  FEBasis mshbasis, metbasis; 
  MetSpace metspace; 

  int idim, ideg;
  bool imet; 
  int npoin, nedge, nface, nelem;

  dblAr2 coord;
  dblAr2 metfld;

  intAr2 edg2poi, fac2poi, tet2poi;
  intAr1 edg2ref, fac2ref, tet2ref; 

public: // Unfortunate 
  int ncorn, ngpoe, ngpof; 
  intAr1 lcorn;
  intAr2 lgpoe, lgpof;
  dblAr2 rgpoe, rgpof; 
private:
  bool flagsInit;
  CADInfo CAD_; 
public:
  const CADInfo& CAD;
private:
  int usrord[3][tetnpps[METRIS_MAX_DEG]];

  //MetrisRunner *run;
};




} // end namespace
#endif