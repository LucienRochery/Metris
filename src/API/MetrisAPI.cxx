//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "MetrisAPI.hxx"
#include "../MetrisRunner/MetrisRunner.hxx"
#include "../Mesh/Mesh.hxx"
#include "../Boundary/msh_inisurf.hxx"


namespace Metris{

//MetrisAPI::~MetrisAPI(){
////  if(run != NULL) run->moveAPI(); 
////  run = NULL; 
//}


MetrisAPI::MetrisAPI() : flagsInit(false), CAD(CAD_) {}

MetrisAPI::MetrisAPI(int idim, int ideg, int ncorn, int ngpoe, int ngpof, 
                     int npoin, bool imet, int nedge, int nface, int nelem,
                     FEBasis mshbasis, FEBasis metbasis, MetSpace metspace) : 
  flagsInit(false), CAD(CAD_)
{

  initialize(idim,ideg,ncorn,ngpoe,ngpof,
             npoin,imet,nedge,nface,nelem,mshbasis,metbasis,metspace);

  //run = NULL;

  this->idim = idim;
  this->ideg = ideg;
  this->npoin = npoin;
  this->nedge = nedge;
  this->nface = nface;
  this->nelem = nelem;
  this->imet  = imet ;
  this->ncorn = ncorn;
  this->ngpoe = ngpoe; 
  this->ngpof = ngpof; 

  this->mshbasis = mshbasis;
  this->metbasis = metbasis;
  this->metspace = metspace;

  coord.allocate(npoin,idim);
  coord.set_n(npoin); 

  if(imet) metfld.allocate(npoin,(idim*(idim+1))/2);
  if(imet) metfld.set_n(npoin);


  edg2poi.allocate(nedge,edgnpps[ideg]);
  edg2poi.set_n(nedge);
  edg2ref.allocate(nedge);
  edg2ref.set_n(nedge);

  fac2poi.allocate(nface,facnpps[ideg]);
  fac2poi.set_n(nface);
  fac2ref.allocate(nface);
  fac2ref.set_n(nface);

  tet2poi.allocate(nelem,tetnpps[ideg]);
  tet2poi.set_n(nelem);
  tet2ref.allocate(nelem);
  tet2ref.set_n(nelem);

  lcorn.allocate(ncorn);
  lcorn.set_n(ncorn);  
  lgpoe.allocate(ngpoe,2);
  lgpoe.set_n(ngpoe);  
  rgpoe.allocate(ngpoe,2);
  rgpoe.set_n(ngpoe);  
  lgpof.allocate(ngpof,2);
  lgpof.set_n(ngpof);  
  rgpof.allocate(ngpof,3);
  rgpof.set_n(ngpof);  


  for(int tdimn = 1; tdimn <= 3; tdimn++) usrord[tdimn-1][0] = -1;
}

MetrisAPI::MetrisAPI(MetrisRunner &run) :  flagsInit(true), CAD_(), CAD(CAD_){
  initialize(run); 
}

void MetrisAPI::free(){
  coord.free();
  metfld.free();

  edg2poi.free();
  edg2ref.free();

  fac2poi.free();
  fac2ref.free();

  tet2poi.free();
  tet2ref.free();

  lcorn.free();
  lgpoe.free();
  rgpoe.free();
  lgpof.free();
  rgpof.free();

  flagsInit = false;
}

void MetrisAPI::setCADModel(ego EGADS_context, ego EGADS_model){
  METRIS_ENFORCE(EGADS_model != NULL);
  CAD_.setModel(EGADS_context, EGADS_model);
}


void MetrisAPI::initialize(MetrisRunner &run){

  flagsInit = true;

  for(int tdimn = 1; tdimn <= 3; tdimn++) usrord[tdimn-1][0] = -1;

  idim = run.msh_g->idim;
  ideg = run.msh_g->curdeg;
  if(ideg != run.msh_g->strdeg) METRIS_THROW_MSG(TODOExcept(), 
    "Implement resizing arrays in API curdeg = "
    <<run.msh_g->curdeg<<" strdeg = "<<run.msh_g->strdeg);

  //if(run.hookedAPI != NULL)METRIS_THROW_MSG(TODOExcept(), 
  //    "Either avoid initializing several MetrisAPI using single Runner,\
  //     or implement hookedAPI as array")

  //this->run = &run;

  //run.hookedAPI = this;
  

  ncorn = 0;
  //run.onhook_npoin = npoin = run.msh_g->npoin;
  //run.onhook_nedge = nedge = run.msh_g->nedge;
  //run.onhook_nface = nface = run.msh_g->nface;
  //run.onhook_nelem = nelem = run.msh_g->nelem;
  //                   npmet = run.msh_g->npoin;
  npoin = run.msh_g->npoin;
  nedge = run.msh_g->nedge;
  nface = run.msh_g->nface;
  nelem = run.msh_g->nelem;
  imet  = true; 

  int mcorn = run.msh_g->nbpoi;
  int mgpoe = run.msh_g->nbpoi;
  int mgpof = run.msh_g->nbpoi;
  lcorn.allocate(mcorn);
  lgpoe.allocate(mgpoe,2);
  rgpoe.allocate(mgpoe,2);
  lgpof.allocate(mgpof,2);
  rgpof.allocate(mgpof,3);
  intAr1 dum;
  genOnGeometricEntLists(*(run.msh_g), lcorn, dum ,
                                       lgpoe, rgpoe,
                                       lgpof, rgpof);
  ncorn = lcorn.get_n();
  ngpoe = lgpoe.get_n();
  ngpof = lgpof.get_n();
  
  coord.free();
  coord   = std::move(run.msh_g->coord);
  edg2poi.free();
  edg2poi = std::move(run.msh_g->edg2poi);
  fac2poi.free();
  fac2poi = std::move(run.msh_g->fac2poi);
  tet2poi.free();
  tet2poi = std::move(run.msh_g->tet2poi);

  edg2ref.free();
  edg2ref = std::move(run.msh_g->edg2ref);
  fac2ref.free();
  fac2ref = std::move(run.msh_g->fac2ref);
  tet2ref.free();
  tet2ref = std::move(run.msh_g->tet2ref);

  mshbasis = run.msh_g->getBasis();
  //run.onhook_mshbasis = mshbasis = run.msh_g->getBasis();

  if(run.metricFE){
    Mesh<MetricFieldFE> *msh = (Mesh<MetricFieldFE>*) run.msh_g;
    metfld.free();
    metfld = std::move(msh->met.rfld);
    //run.onhook_metbasis = metbasis = msh->met.getBasis();
    //run.onhook_metspace = metspace = msh->met.getSpace();
    metbasis = msh->met.getBasis();
    metspace = msh->met.getSpace();
  }else{
    Mesh<MetricFieldAnalytical> *msh = (Mesh<MetricFieldAnalytical>*) run.msh_g;
    metfld.free();
    metfld = std::move(msh->met.rfld);
    //run.onhook_metbasis = metbasis = msh->met.getBasis();
    //run.onhook_metspace = metspace = msh->met.getSpace();
    metbasis = msh->met.getBasis();
    metspace = msh->met.getSpace();
  }

  CAD_ = run.msh_g->CAD; 
  //run.~MetrisRunner();
}


void MetrisAPI::initialize(int idim, int ideg, int ncorn, int ngpoe, int ngpof, 
                           int npoin, bool imet, int nedge, int nface, int nelem,
                           FEBasis mshbasis, FEBasis metbasis, MetSpace metspace){

  initialize(idim,ideg,imet,mshbasis,metbasis,metspace);

  setNPoints(npoin);

  setNEdges(nedge);
  setNFaces(nface);
  setNTetrahedra(nelem);

  setNCorners(ncorn);
  setNVerticesOnGeometricEdges(ngpoe);
  setNVerticesOnGeometricTriangles(ngpof);

  for(int tdimn = 1; tdimn <= 3; tdimn++) usrord[tdimn-1][0] = -1;
}

void MetrisAPI::initialize(int idim, int ideg, bool imet,
                  FEBasis mshbasis, FEBasis metbasis, MetSpace metspace){
  flagsInit = true;
  this->idim = idim;
  this->ideg = ideg;
  this->imet  = imet ;
  this->mshbasis = mshbasis;
  this->metbasis = metbasis;
  this->metspace = metspace;
}



void MetrisAPI::getConstants(int* idim, int* ideg,
                             int* ncorn, int* npoin, bool* imet, 
                             int* nedge, int* nface, int* nelem,
                             FEBasis* mshbasis, FEBasis* metbasis) const{
  *idim = this->idim;
  *ideg = this->ideg;
  *ncorn = this->ncorn;
  *npoin = this->npoin;
  *imet = this->imet ;
  *nedge = this->nedge;
  *nface = this->nface;
  *nelem = this->nelem;
  *mshbasis = this->mshbasis;
  *metbasis = this->metbasis;
}


void MetrisAPI::setDegree(int tardeg){
  METRIS_ENFORCE(tardeg >= 1 && tardeg < METRIS_MAX_DEG);
  if(tardeg < ideg) METRIS_ENFORCE(tardeg == 1)

  intAr1 lpoin;
  bool ialloc = false;

  int ndead[3] = {0};
  for(int tdim = 1; tdim <= 3; tdim++){
    int nentt = tdim == 1 ? nedge :
                tdim == 2 ? nface : 
                            nelem;
    if(nentt <= 0) continue;

    intAr2 &ent2poi = tdim == 1 ? edg2poi :
                      tdim == 2 ? fac2poi : 
                                  tet2poi;
    int nnode = tdim == 1 ? edgnpps[tardeg] :
                tdim == 2 ? facnpps[tardeg] : 
                            tetnpps[tardeg];


    if(tardeg < ideg){
      METRIS_ASSERT(nnode == tdim + 1);
      if(!ialloc){
        // lpoin will serve as the mapping but also tagging array 
        // We need to do this in two phases to avoid a sort -> naturally 
        // guarantee iponw <= ipoin by setting iponw in loop over ipoin.
        ialloc = true;
        lpoin.allocate(npoin);
        lpoin.set_n(npoin);
        lpoin.fill(-2);
      }
      // Points into same memory region but lets us access with diff stride
      intAr2 ent2pon(nentt, tdim + 1, ent2poi[0]);
      ent2pon.set_n(nentt);
      // We can skip first because they won't move 
      for(int ientt = 1; ientt < nentt; ientt++){
        if(isdeadent(ientt, ent2poi)){
          ndead[tdim-1]++;
        }
        for(int ii = 0; ii < tdim + 1; ii++){
          int ipoin = ent2poi(ientt,ii);
          ent2pon(ientt,ii) = ipoin;
          lpoin[ipoin] = -1;
        }
      }// for ientt
      ent2poi.set_stride(tdim+1);
    }else{//if ideg 
      ent2poi.allocate(nentt,nnode);
    }//elseif ideg 

    if(ndead[tdim-1] > 0) 
      printf("# Warning: %d dead tdim %d elements in API\n",
             ndead[tdim-1], tdim);
  }//for tdim
  ideg = tardeg;


  // Reorder vertices 
  if(ialloc){

    int nponw = 0;
    for(int ipoin = 0; ipoin < npoin; ipoin++){
      if(lpoin[ipoin] == -2) continue; // never seen
      if(lpoin[ipoin] >=  0) continue; // already set 
      lpoin[ipoin] = nponw;
      nponw++;
    }

    const int nnmet = (idim * (idim + 1))/2;

    //// Inverse mapping     
    //intAr1 lponw(nponw);
    //lponw.set_n(nponw);
    //for(int ipoin = 0; ipoin < npoin; ipoin++){
    //  int iponw = lpoin[ipoin];
    //  if(iponw < 0) continue;
    //  lponw[iponw] = ipoin;
    //}

    for(int tdim = 1; tdim <= 3; tdim++){
      int nentt = tdim == 1 ? nedge :
                  tdim == 2 ? nface : 
                              nelem;
      intAr2 &ent2poi = tdim == 1 ? edg2poi :
                        tdim == 2 ? fac2poi : 
                                    tet2poi;
      //int nnode = tdim == 1 ? edgnpps[ideg] :
      //            tdim == 2 ? facnpps[ideg] : 
      //                        tetnpps[ideg];

      for(int ientt = 0; ientt < nentt; ientt++){
        for(int ii = 0; ii < tdim + 1; ii++){
          int ipoin = ent2poi(ientt,ii);
          int iponw = lpoin[ipoin];
          ent2poi(ientt,ii) = iponw;
        }
      }// for ientt
    }

    for(int ipoin = 0; ipoin < npoin; ipoin++){
      int iponw = lpoin[ipoin];
      METRIS_ASSERT(iponw != -1);
      if(iponw == -2) continue;
      // iponw is always <= ipoin
      // Hence there is no risk of overwriting.

      for(int ii = 0; ii < idim; ii++) coord(iponw,ii) = coord(ipoin,ii);

      if(!imet) continue;
      for(int ii = 0; ii < nnmet; ii++) metfld(iponw,ii) = metfld(ipoin,ii);
    }

    for(int icorn = 0; icorn < ncorn; icorn++){
      int ipoin = lcorn[icorn];
      int iponw = lpoin[ipoin];
      lcorn[icorn] = iponw;
    }

    for(int igpoe = 0; igpoe < ngpoe; igpoe++){
      int ipoin = lgpoe(igpoe,0);
      int iponw = lpoin[ipoin];
      lgpoe(igpoe,0) = iponw;
    }

    for(int igpof = 0; igpof < ngpof; igpof++){
      int ipoin = lgpof(igpof,0);
      int iponw = lpoin[ipoin];
      lgpof(igpof,0) = iponw;
    }

    printf("## DEBUG: MetrisAPI cleanup() npoin %d -> %d\n",npoin,nponw);
    npoin = nponw;

  }

  return;
}



// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Setters

void MetrisAPI::setNPoints(int npoin_){
  // Only wrong program flow would lead to this, so assert 
  METRIS_ASSERT(flagsInit);
  // Input dependent 
  METRIS_ENFORCE(npoin_ >= 0);

  npoin = npoin_;
  coord.allocate(npoin,idim);
  coord.set_n(npoin); 

  if(imet) metfld.allocate(npoin,(idim*(idim+1))/2);
  if(imet) metfld.set_n(npoin);
}

void MetrisAPI::setNEdges(int nedge_){
  // Only wrong program flow would lead to this, so assert 
  METRIS_ASSERT(flagsInit);
  // Input dependent 
  METRIS_ENFORCE(nedge_ >= 0);

  nedge = nedge_;
  edg2poi.allocate(nedge,edgnpps[ideg]);
  edg2poi.set_n(nedge);
  edg2ref.allocate(nedge);
  edg2ref.set_n(nedge);
}
void MetrisAPI::setNFaces(int nface_){
  // Only wrong program flow would lead to this, so assert 
  METRIS_ASSERT(flagsInit);
  // Input dependent 
  METRIS_ENFORCE(nface_ >= 0);

  nface = nface_;
  fac2poi.allocate(nface,facnpps[ideg]);
  fac2poi.set_n(nface);
  fac2ref.allocate(nface);
  fac2ref.set_n(nface);
}
void MetrisAPI::setNTetrahedra(int nelem_){
  // Only wrong program flow would lead to this, so assert 
  METRIS_ASSERT(flagsInit);
  // Input dependent 
  METRIS_ENFORCE(nelem_ >= 0);

  nelem = nelem_;
  tet2poi.allocate(nelem,tetnpps[ideg]);
  tet2poi.set_n(nelem);
  tet2ref.allocate(nelem);
  tet2ref.set_n(nelem);
}
void MetrisAPI::setNCorners(int ncorn_){
  // Only wrong program flow would lead to this, so assert 
  METRIS_ASSERT(flagsInit);
  // Input dependent 
  METRIS_ENFORCE(ncorn_ >= 0);

  ncorn = ncorn_;
  lcorn.allocate(ncorn_);
  lcorn.set_n(ncorn_);  
}
void MetrisAPI::setNVerticesOnGeometricEdges(int ngpoe_){
  // Only wrong program flow would lead to this, so assert 
  METRIS_ASSERT(flagsInit);
  // Input dependent 
  METRIS_ENFORCE(ngpoe_ >= 0);

  ngpoe = ngpoe_;
  lgpoe.allocate(ngpoe_,2);
  lgpoe.set_n(ngpoe_);  
  rgpoe.allocate(ngpoe_,2);
  rgpoe.set_n(ngpoe_);  
}
void MetrisAPI::setNVerticesOnGeometricTriangles(int ngpof_){
  // Only wrong program flow would lead to this, so assert 
  METRIS_ASSERT(flagsInit);
  // Input dependent 
  METRIS_ENFORCE(ngpof_ >= 0);

  ngpof = ngpof_;
  lgpof.allocate(ngpof_,2);
  lgpof.set_n(ngpof_);  
  rgpof.allocate(ngpof_,3);
  rgpof.set_n(ngpof_);  
}


//// -----------------------------------------------------------------------------
//// -----------------------------------------------------------------------------
//// Getters
//
//int MetrisAPI::getNPoints() const {
//  return npoin;
//}
//int MetrisAPI::getNEdges() const {
//  return nedge;
//}
//int MetrisAPI::getNFaces() const {
//  return nface;
//}
//int MetrisAPI::getNTetrahedra() const {
//  return nelem;
//}
//int MetrisAPI::getNCorners() const {
//  return ncorn; 
//}
//int MetrisAPI::getNVerticesOnGeometricEdges() const {
//  return ngpoe;
//}
//int MetrisAPI::getNVerticesOnGeometricTriangles() const {
//  return ngpof;
//}

// -----------------------------------------------------------------------------
// Copy into other API 

void MetrisAPI::copyFlags(MetrisAPI *into) const{
  METRIS_ENFORCE(flagsInit);
  into->initialize(idim,ideg,imet,mshbasis,metbasis,metspace);
}

void MetrisAPI::copyCorners(MetrisAPI *into) const{
  into->setNCorners(ncorn);
  if(ncorn == 0) return;
  into->setCorner(0,ncorn,&lcorn[0]);
}

void MetrisAPI::copyVerticesOnGeometricEdges(MetrisAPI *into) const{
  into->setNVerticesOnGeometricEdges(ngpoe);
  if(ngpoe == 0) return;
  into->setVerticesOnGeometricEdges(0,ngpoe,lgpoe[0],rgpoe[0]);
}

void MetrisAPI::copyVerticesOnGeometricTriangles(MetrisAPI *into) const{
  into->setNVerticesOnGeometricTriangles(ngpof);
  if(ngpof == 0) return;
  into->setVerticesOnGeometricTriangles(0,ngpof,lgpof[0],rgpof[0]);
}

void MetrisAPI::copyCAD(MetrisAPI *into) const{
  into->CAD_ = CAD;
}





// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Coordinates

void MetrisAPI::setCoord(int ipoin, const double *coord){
  METRIS_ASSERT(ipoin >= 0 && ipoin < npoin);
  for(int ii = 0; ii < idim; ii++) this->coord(ipoin,ii) = coord[ii];
}

void MetrisAPI::setCoord(int ipoi1, int ipoi2, const double *coord){
  METRIS_ASSERT(ipoi1 >= 0 && ipoi1 < npoin);
  METRIS_ASSERT(ipoi2 >= 0 && ipoi2 <= npoin);
  METRIS_ASSERT(ipoi2 > ipoi1);
  for(int ipoin = ipoi1; ipoin < ipoi2; ipoin++){
    for(int ii = 0; ii < idim; ii++) 
      this->coord(ipoin,ii) = coord[idim*(ipoin - ipoi1) + ii];
  }
}

void MetrisAPI::setCoord(dblAr2 &&coord){
  this->coord = std::move(coord); 
}


void MetrisAPI::getCoord(int ipoin, double *coord) const {
  METRIS_ASSERT(ipoin >= 0 && ipoin < npoin);
  for(int ii = 0; ii < idim; ii++) coord[ii] = this->coord(ipoin,ii);
}

void MetrisAPI::getCoord(int ipoi1, int ipoi2, double *coord) const {
  METRIS_ASSERT(ipoi1 >= 0 && ipoi1 < npoin);
  METRIS_ASSERT(ipoi2 >= 0 && ipoi2 <= npoin);
  METRIS_ASSERT(ipoi2 > ipoi1);
  for(int ipoin = ipoi1; ipoin < ipoi2; ipoin++){
    for(int ii = 0; ii < idim; ii++) 
      coord[idim*(ipoin - ipoi1) + ii] = this->coord(ipoin,ii);
  }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Metric field

void MetrisAPI::setMetric(int ipoin, const double *metfld){
  if(!imet) return;
  METRIS_ASSERT_MSG(ipoin >= 0 && ipoin < npoin,
    "ipoin = "<<ipoin<<" npoin = "<<npoin);
  int nnmet = (idim*(idim+1))/2;
  for(int ii = 0; ii < nnmet; ii++) this->metfld(ipoin,ii) = metfld[ii];
}

void MetrisAPI::setMetric(int ipoi1, int ipoi2, const double *metfld){
  if(!imet) return;
  METRIS_ASSERT(ipoi1 >= 0 && ipoi1 < npoin);
  METRIS_ASSERT(ipoi2 >= 0 && ipoi2 <= npoin);
  METRIS_ASSERT(ipoi2 > ipoi1);
  int nnmet = (idim*(idim+1))/2;
  for(int ipoin = ipoi1; ipoin < ipoi2; ipoin++){
    for(int ii = 0; ii < nnmet; ii++) 
      this->metfld(ipoin,ii) = metfld[nnmet*(ipoin - ipoi1) + ii];
  }
}

void MetrisAPI::setMetric(dblAr2 &&metfld){
  this->metfld = std::move(metfld); 
}




void MetrisAPI::getMetric(int ipoin, double *metfld) const {
  METRIS_ASSERT(ipoin >= 0 && ipoin < npoin);
  int nnmet = (idim*(idim+1))/2;
  for(int ii = 0; ii < nnmet; ii++) metfld[ii] = this->metfld(ipoin,ii);
}

void MetrisAPI::getMetric(int ipoi1, int ipoi2, double *metfld) const {
  if(!imet) return;
  METRIS_ASSERT(ipoi1 >= 0 && ipoi1 < npoin);
  METRIS_ASSERT(ipoi2 >= 0 && ipoi2 < npoin);
  METRIS_ASSERT(ipoi2 > ipoi1);
  int nnmet = (idim*(idim+1))/2;
  for(int ipoin = ipoi1; ipoin < ipoi2; ipoin++){
    for(int ii = 0; ii < nnmet; ii++) 
      metfld[nnmet*(ipoin - ipoi1) + ii] = this->metfld(ipoin,ii);
  }
}




// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Elements

void MetrisAPI::setElementsOrdering(int tdimn, const int *ordering){
  METRIS_ASSERT(ordering != NULL);

  int nnode = tdimn == 1 ? edgnpps[ideg] : 
              tdimn == 2 ? facnpps[ideg] : tetnpps[ideg];

  for(int ii = 0; ii < nnode; ii++){
    int inode = mul2nod(tdimn,&ordering[(tdimn+1) * ii]);
    METRIS_ENFORCE_MSG(usrord[tdimn-1][inode] == -1, "Ordering repeats itself inode = "<<ii<<
      " previous = "<<usrord[tdimn-1][inode]);
    usrord[tdimn-1][inode] = ii;
  }
}

void MetrisAPI::setElementsOrdering(int iordering){
  METRIS_THROW_MSG(TODOExcept(), 
    "Default orderings not implemented. "
    "Either don't call setElementsOrdering, or provide explicitely. Given"<<iordering);
}

void MetrisAPI::setElement(int tdimn, int ielem, const int *ent2pol, int iref){
  METRIS_ASSERT(tdimn >= 1 && tdimn <= 3);
  intAr2 &ent2poi = tdimn == 1 ? edg2poi : 
                    tdimn == 2 ? fac2poi : tet2poi;
  intAr1 &ent2ref = tdimn == 1 ? edg2ref : 
                    tdimn == 2 ? fac2ref : tet2ref;
  int nnode = tdimn == 1 ? edgnpps[ideg] : 
              tdimn == 2 ? facnpps[ideg] : tetnpps[ideg];

  if(usrord[tdimn-1][0] == -1){
    for(int ii = 0; ii < nnode; ii++) ent2poi(ielem,ii) = ent2pol[ii];
  }else{
    for(int ii = 0; ii < nnode; ii++) ent2poi[ielem][usrord[tdimn-1][ii]] = ent2pol[ii];
  }
  ent2ref[ielem] = iref;
}


void MetrisAPI::setElement(int tdimn, int iele1, int iele2, const int *ent2pol, const int *lref){
  METRIS_ASSERT(tdimn >= 1 && tdimn <= 3);
  intAr2 &ent2poi = tdimn == 1 ? edg2poi : 
                    tdimn == 2 ? fac2poi : tet2poi;
  intAr1 &ent2ref = tdimn == 1 ? edg2ref : 
                    tdimn == 2 ? fac2ref : tet2ref;
  int nnode = tdimn == 1 ? edgnpps[ideg] : 
              tdimn == 2 ? facnpps[ideg] : tetnpps[ideg];

  if(usrord[tdimn-1][0] == -1){
    for(int ielem = iele1; ielem < iele2; ielem++){
      for(int ii = 0; ii < nnode; ii++) ent2poi(ielem,ii) = ent2pol[nnode*(ielem-iele1) + ii];
      ent2ref[ielem] = lref[ielem-iele1];
    }
  }else{
    for(int ielem = iele1; ielem < iele2; ielem++){
      for(int ii = 0; ii < nnode; ii++) ent2poi[ielem][usrord[tdimn-1][ii]] = ent2pol[nnode*(ielem-iele1) + ii];
      ent2ref[ielem] = lref[ielem-iele1];
    }
  }
}

void MetrisAPI::setElement(int tdimn, intAr2 &&ent2poi, intAr1 &&ent2ref){
  intAr2 &ent2poi_ = tdimn == 1 ? edg2poi : 
                     tdimn == 2 ? fac2poi : tet2poi;
  intAr1 &ent2ref_ = tdimn == 1 ? edg2ref : 
                     tdimn == 2 ? fac2ref : tet2ref;
  ent2poi_ = std::move(ent2poi); 
  ent2ref_ = std::move(ent2ref); 
}


void MetrisAPI::getElement(int tdimn, int ielem, int *ent2pol, int *iref) const {
  METRIS_ASSERT(tdimn >= 1 && tdimn <= 3);
  const intAr2 &ent2poi = tdimn == 1 ? edg2poi : 
                          tdimn == 2 ? fac2poi : tet2poi;
  const intAr1 &ent2ref = tdimn == 1 ? edg2ref : 
                          tdimn == 2 ? fac2ref : tet2ref;
  int nnode = tdimn == 1 ? edgnpps[ideg] : 
              tdimn == 2 ? facnpps[ideg] : tetnpps[ideg];

  *iref = ent2ref[ielem];

  if(ent2pol == NULL) return;

  if(usrord[tdimn-1][0] == -1){
    for(int ii = 0; ii < nnode; ii++) ent2pol[ii] = ent2poi(ielem,ii);
  }else{
    for(int ii = 0; ii < nnode; ii++) ent2pol[ii] = ent2poi[ielem][usrord[tdimn-1][ii]];
  }
}

void MetrisAPI::getElement(int tdimn, int iele1, int iele2, int *ent2pol, int *lref) const {
  METRIS_ASSERT(tdimn >= 1 && tdimn <= 3);
  const intAr2 &ent2poi = tdimn == 1 ? edg2poi : 
                          tdimn == 2 ? fac2poi : tet2poi;
  const intAr1 &ent2ref = tdimn == 1 ? edg2ref : 
                          tdimn == 2 ? fac2ref : tet2ref;
  int nnode = tdimn == 1 ? edgnpps[ideg] : 
              tdimn == 2 ? facnpps[ideg] : tetnpps[ideg];

  for(int ielem = iele1; ielem < iele2; ielem++){
    lref[ielem-iele1] = ent2ref[ielem];
  }

  if(ent2pol == NULL) return;

  if(usrord[tdimn-1][0] == -1){
    for(int ielem = iele1; ielem < iele2; ielem++){
      for(int ii = 0; ii < nnode; ii++){
        ent2pol[nnode*(ielem-iele1) + ii] = ent2poi(ielem,ii);
      }
    }
  }else{
    for(int ielem = iele1; ielem < iele2; ielem++){
      for(int ii = 0; ii < nnode; ii++){
        ent2pol[nnode*(ielem-iele1) + ii] = ent2poi[ielem][usrord[tdimn-1][ii]];
      }
    }
  }
}

void MetrisAPI::getElementRef(int tdimn, int ielem, int *iref) const{
  METRIS_ASSERT(tdimn >= 1 && tdimn <= 3);
  if(tdimn == 1){
    *iref = edg2ref[ielem];
  }else if(tdimn == 2){
    *iref = fac2ref[ielem];
  }else{
    *iref = tet2ref[ielem];
  }
  return;
}



// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Boundary info 

void MetrisAPI::setCorner(int icorn, int ipcor){
  if(ncorn == 0) return;
  METRIS_ASSERT(ipcor >= 0 && ipcor < ncorn);
  lcorn[icorn] = ipcor; 
}

void MetrisAPI::setCorner(int icor1, int icor2, const int *lpcor){
  if(ncorn == 0) return;
  METRIS_ASSERT(icor1 >= 0 && icor1 < ncorn);
  METRIS_ASSERT(icor2 >= 0 && icor2 <= ncorn);
  METRIS_ASSERT(icor2 > icor1);
  for(int icorn = icor1; icorn < icor2; icorn++){
    lcorn[icorn] = lpcor[icorn - icor1];
  }
}

void MetrisAPI::setVerticesOnGeometricEdges(int ientry, const int *lgpoe, const double *rgpoe){
  METRIS_ASSERT(ientry < ngpoe && ientry >= 0);
  for(int ii = 0; ii < 2; ii++)
    this->lgpoe[ientry][ii] = lgpoe[ii];

  for(int ii = 0; ii < 2; ii++)
    this->rgpoe[ientry][ii] = rgpoe[ii];
}

void MetrisAPI::setVerticesOnGeometricEdges(int ientr1, int ientr2, const int *lgpoe, const double *rgpoe){
  if(ngpoe == 0) return; 
  METRIS_ASSERT(ientr1 >= 0 && ientr1 < ngpoe && ientr1 < ientr2);
  METRIS_ASSERT(ientr2 >= 0 && ientr2 <= ngpoe);
  for(int ientry = ientr1; ientry < ientr2; ientry++){
    for(int ii = 0; ii < 2; ii++)
      this->lgpoe[ientry][ii] = lgpoe[2*(ientry - ientr1) + ii];

    for(int ii = 0; ii < 2; ii++)
      this->rgpoe[ientry][ii] = rgpoe[2*(ientry - ientr1) + ii];
  }
}

void MetrisAPI::setVerticesOnGeometricTriangles(int ientry, const int *lgpof, const double *rgpof){
  METRIS_ASSERT(ientry < ngpof && ientry >= 0);
  for(int ii = 0; ii < 2; ii++)
    this->lgpof[ientry][ii] = lgpof[ii];

  for(int ii = 0; ii < 3; ii++)
    this->rgpof[ientry][ii] = rgpof[ii];
}

void MetrisAPI::setVerticesOnGeometricTriangles(int ientr1, int ientr2, const int *lgpof, const double *rgpof){
  if(ngpof == 0) return; 
  METRIS_ASSERT(ientr1 >= 0 && ientr1 < ngpof && ientr1 < ientr2);
  METRIS_ASSERT(ientr2 >= 0 && ientr2 <= ngpof);
  for(int ientry = ientr1; ientry < ientr2; ientry++){
    for(int ii = 0; ii < 2; ii++)
      this->lgpof[ientry][ii] = lgpof[2*(ientry - ientr1) + ii];

    for(int ii = 0; ii < 3; ii++)
      this->rgpof[ientry][ii] = rgpof[3*(ientry - ientr1) + ii];
  }
}

// -----------------------------------------------------------------------------
// Getters

void MetrisAPI::getCorner(int icorn, int *ipcor) const {
  if(ncorn == 0) return;
  METRIS_ASSERT(icorn >= 0 && icorn < ncorn);
  *ipcor = lcorn[icorn]; 
}

void MetrisAPI::getCorner(int icor1, int icor2, int *lpcor) const {
  if(ncorn == 0) return;
  METRIS_ASSERT(icor1 >= 0 && icor1 < ncorn);
  METRIS_ASSERT(icor2 >= 0 && icor2 <= ncorn);
  METRIS_ASSERT(icor2 > icor1);
  for(int icorn = icor1; icorn < icor2; icorn++){
    lpcor[icorn - icor1] = lcorn[icorn];
  }
}

void MetrisAPI::getVerticesOnGeometricEdges(int ientry, int *lgpoe, double *rgpoe) const {
  METRIS_ASSERT(ientry < ngpoe && ientry >= 0);
  for(int ii = 0; ii < 2; ii++)
    lgpoe[ii] = this->lgpoe[ientry][ii];

  for(int ii = 0; ii < 2; ii++)
    rgpoe[ii] = this->rgpoe[ientry][ii];
}

void MetrisAPI::getVerticesOnGeometricEdges(int ientr1, int ientr2, int *lgpoe, double *rgpoe) const {
  if(ngpoe == 0) return; 
  METRIS_ASSERT(ientr1 >= 0 && ientr1 < ngpoe && ientr1 < ientr2);
  METRIS_ASSERT(ientr2 >= 0 && ientr2 <= ngpoe);
  for(int ientry = ientr1; ientry < ientr2; ientry++){
    for(int ii = 0; ii < 2; ii++)
      lgpoe[2*(ientry - ientr1) + ii] = this->lgpoe[ientry][ii];

    for(int ii = 0; ii < 2; ii++)
      rgpoe[2*(ientry - ientr1) + ii] = this->rgpoe[ientry][ii];
  }
}

void MetrisAPI::getVerticesOnGeometricTriangles(int ientry, int *lgpof, double *rgpof) const {
  METRIS_ASSERT(ientry < ngpof && ientry >= 0);
  for(int ii = 0; ii < 2; ii++)
    lgpof[ii] = this->lgpof[ientry][ii];

  for(int ii = 0; ii < 3; ii++)
    rgpof[ii] = this->rgpof[ientry][ii];
}

void MetrisAPI::getVerticesOnGeometricTriangles(int ientr1, int ientr2, int *lgpof, double *rgpof) const {
  if(ngpof == 0) return; 
  METRIS_ASSERT(ientr1 >= 0 && ientr1 < ngpof && ientr1 < ientr2);
  METRIS_ASSERT(ientr2 >= 0 && ientr2 <= ngpof);
  for(int ientry = ientr1; ientry < ientr2; ientry++){
    for(int ii = 0; ii < 2; ii++)
      lgpof[2*(ientry - ientr1) + ii] = this->lgpof[ientry][ii];

    for(int ii = 0; ii < 3; ii++)
      rgpof[3*(ientry - ientr1) + ii] = this->rgpof[ientry][ii];
  }
}


} // end namespace
