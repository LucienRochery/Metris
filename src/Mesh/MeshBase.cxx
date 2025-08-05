//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../Mesh/MeshBase.hxx"

#include "../MetricField/MetricField.hxx"


#include "../MetrisRunner/MetrisParameters.hxx"
#include "../aux_exceptions.hxx"
#include "../metris_constants.hxx"
#include "../ho_constants.hxx"
#include "../utils/CT_loop.hxx"
#include "../msh_lag2bez.hxx"
#include "../low_geo/normal.hxx"
#include "../linalg/det.hxx"
#include "../low_geo/misc.hxx"
#include "../utils/mprintf.hxx"

namespace Metris{

MetricClass MeshBase::metricClass() const { return MetricClass::None; }
  
MeshBase::MeshBase(){
  param  = NULL;
  nbpoi_ = (mbpoi_ = 0);
  npoin_ = (mpoin_ = 0);
  nedge_ = (medge_ = 0);
  nface_ = (mface_ = 0);
  nelem_ = (melem_ = 0);
  ibasis  = FEBasis::Undefined; 
  curdeg  = 0;
  strdeg  = 0; 
  //hasbak  = false;
  //ianamet = -1;
  idim    = 0;
}


// Mostly for internal use, compute the localization alignment direction
// (edge tangent, face normal) of a point with CAD link and poi2ent initialized.
void MeshBase::get_algnd(int ipoin, double* algnd){
  GETVDEPTH(param);
  int tdim = getpoitdim(ipoin);
  if(tdim == idim) return;

  int ientt = poi2ent(ipoin,0);
  int iref  = ent2ref(tdim)[ientt];
  bool done = false;
  if(CAD()){
    done = true; 
    int ibpoi = poi2bpo[ipoin];
    METRIS_ASSERT(bpo2ibi(ibpoi,1) == tdim);

    double result[18];
    ego obj = tdim == 1 ? CAD.cad2edg[iref] : CAD.cad2fac[iref];
    int ierro = EG_evaluate(obj, bpo2rbi[ibpoi], result);
    METRIS_ASSERT(ierro == 0); // EGADS errors are not expected 

    if(tdim == 1){
      for(int ii = 0; ii < idim; ii++) algnd[ii] = result[3+ii];
    }else{
      vecprod(&result[3],&result[6],algnd);
      if(normalize_vec<3>(algnd)){
        CPRINTF1("# ZERO NORMAL IN get_algnd")
        done = false;
      }
    }
  }

  if(!done){
    if(tdim == 1){
      int ipoi1 = edg2poi(ientt,0);
      int ipoi2 = edg2poi(ientt,1);
      for(int ii = 0; ii < idim; ii++)
        algnd[ii] = coord(ipoi1,ii) - coord(ipoi2,ii);
    }else{
      getnorfacP1(fac2poi[ientt],coord,algnd);
    }
  }


}


int MeshBase::get_tdim() const{
  if(nelem > 0) return 3;
  else if (nface > 0) return 2;
  return 1;
}

int MeshBase::nentt(int tdimn) const {
  switch(tdimn){
  case(1):
    return nedge;
  case(2):
    return nface;
  case(3):
    return nelem;
  default:
    METRIS_THROW_MSG(WArgExcept(),"nentt tdimn not in range = "<<tdimn);
  }
}


int MeshBase::nnode(int tdimn) const {
  switch(tdimn){
  case(1):
    return getnnod1(curdeg);
  case(2):
    return getnnod2(curdeg);
  case(3):
    return getnnod3(curdeg);
  default:
    METRIS_THROW_MSG(WArgExcept(),"nnode tdimn not in range = "<<tdimn);
  }
}


intAr2& MeshBase::ent2poi(int tdimn){
  switch(tdimn){
  case(1):
    return edg2poi;
  case(2):
    return fac2poi;
  case(3):
    return tet2poi;
  default:
    METRIS_THROW_MSG(WArgExcept(),"ent2poi (1) tdimn not in range = "<<tdimn);
  }
}

const intAr2& MeshBase::ent2poi(int tdimn) const{
  switch(tdimn){
  case(1):
    return edg2poi;
  case(2):
    return fac2poi;
  case(3):
    return tet2poi;
  default:
    METRIS_THROW_MSG(WArgExcept(),"ent2poi (2) tdimn not in range = "<<tdimn);
  }
}


template<int tdimn>
      intAr2& MeshBase::ent2poi(){
  static_assert(tdimn >= 1 && tdimn <= 3);
  if constexpr(tdimn == 1){
    return edg2poi;
  }else if(tdimn == 2){
    return fac2poi;
  }else{
    return tet2poi;
  }
}
template<int tdimn>
const intAr2& MeshBase::ent2poi() const{
  static_assert(tdimn >= 1 && tdimn <= 3);
  if constexpr(tdimn == 1){
    return edg2poi;
  }else if(tdimn == 2){
    return fac2poi;
  }else{
    return tet2poi;
  }
}

template intAr2& MeshBase::ent2poi<1>();
template intAr2& MeshBase::ent2poi<2>();
template intAr2& MeshBase::ent2poi<3>();
template const intAr2& MeshBase::ent2poi<1>()const;
template const intAr2& MeshBase::ent2poi<2>()const;
template const intAr2& MeshBase::ent2poi<3>()const;


intAr1& MeshBase::ent2ref(int tdimn){
  switch(tdimn){
  case(1):
    return edg2ref;
  case(2):
    return fac2ref;
  case(3):
    return tet2ref;
  default:
    METRIS_THROW_MSG(WArgExcept(),"ent2ref (1) tdimn not in range = "<<tdimn);
  }
}

const intAr1& MeshBase::ent2ref(int tdimn) const{
  switch(tdimn){
  case(1):
    return edg2ref;
  case(2):
    return fac2ref;
  case(3):
    return tet2ref;
  default:
    METRIS_THROW_MSG(WArgExcept(),"ent2ref (2) tdimn not in range = "<<tdimn);
  }
}


intAr2& MeshBase::ent2ent(int tdimn){
  switch(tdimn){
  case(1):
    return edg2edg;
  case(2):
    return fac2fac;
  case(3):
    return tet2tet;
  default:
    METRIS_THROW_MSG(WArgExcept(),"ent2ent (1) tdimn not in range = "<<tdimn);
  }
}

const intAr2& MeshBase::ent2ent(int tdimn) const{
  switch(tdimn){
  case(1):
    return edg2edg;
  case(2):
    return fac2fac;
  case(3):
    return tet2tet;
  default:
    METRIS_THROW_MSG(WArgExcept(),"ent2ent (2) tdimn not in range = "<<tdimn);
  }
}


template<int tdimn>
      intAr2& MeshBase::ent2ent(){
  static_assert(tdimn >= 1 && tdimn <= 3);
  if constexpr(tdimn == 1){
    return edg2edg;
  }else if(tdimn == 2){
    return fac2fac;
  }else{
    return tet2tet;
  }
}
template<int tdimn>
const intAr2& MeshBase::ent2ent() const{
  static_assert(tdimn >= 1 && tdimn <= 3);
  if constexpr(tdimn == 1){
    return edg2edg;
  }else if(tdimn == 2){
    return fac2fac;
  }else{
    return tet2tet;
  }
}

template intAr2& MeshBase::ent2ent<1>();
template intAr2& MeshBase::ent2ent<2>();
template intAr2& MeshBase::ent2ent<3>();
template const intAr2& MeshBase::ent2ent<1>()const;
template const intAr2& MeshBase::ent2ent<2>()const;
template const intAr2& MeshBase::ent2ent<3>()const;




template<int tdimn> typename std::conditional<tdimn==1,HshTab_I2I,HshTab_I3I>::type &
MeshBase::hshTab(){
  if constexpr(tdimn == 1) return edgHshTab;
  else                     return facHshTab;
}
template HshTab_I2I& MeshBase::hshTab<1>();
template HshTab_I3I& MeshBase::hshTab<2>();



MeshBase& MeshBase::operator=(const MeshBase &inp){

  strdeg = MAX(inp.strdeg,strdeg);
  curdeg = inp.curdeg;
  ibasis = inp.ibasis;
  idim   = inp.idim;

  is_manifold = inp.is_manifold;

  set_npoin(inp.npoin);
  set_nedge(inp.nedge);
  set_nface(inp.nface);
  set_nelem(inp.nelem);
  set_nbpoi(inp.nbpoi);

  if(inp.npoin > 0) inp.coord.copyTo(coord,inp.npoin);
  if(inp.nbpoi > 0) inp.bpo2ibi.copyTo(bpo2ibi,inp.nbpoi);
  if(inp.nbpoi > 0) inp.bpo2rbi.copyTo(bpo2rbi,inp.nbpoi);
  if(inp.npoin > 0) inp.poi2bpo.copyTo(poi2bpo,inp.npoin);
  if(inp.npoin > 0) inp.poi2ent.copyTo(poi2ent,inp.npoin);
  if(inp.npoin > 0) inp.poicstr.copyTo(poicstr,inp.npoin);
                                                           
  if(inp.nedge > 0) inp.edg2edg.copyTo(edg2edg,inp.nedge);
  if(inp.nedge > 0) inp.edg2fac.copyTo(edg2fac,inp.nedge);
  if(inp.nface > 0) inp.fac2fac.copyTo(fac2fac,inp.nface);
  if(inp.nface > 0 
  && inp.nelem > 0) inp.fac2tet.copyTo(fac2tet,inp.nface);
  if(inp.nelem > 0) inp.tet2tet.copyTo(tet2tet,inp.nelem);
                                                           
  if(inp.nedge > 0) inp.edg2ref.copyTo(edg2ref,inp.nedge);
  if(inp.nface > 0) inp.fac2ref.copyTo(fac2ref,inp.nface);
  if(inp.nelem > 0) inp.tet2ref.copyTo(tet2ref,inp.nelem);
                                                           
  //if(inp.nelem > 0) inp.tet2ftg.copyTo(tet2ftg,inp.nelem);

  for(int iedge = 0; iedge < inp.nedge; iedge++){
    for(int i = 0; i < getnnod1(curdeg); i++)
      edg2poi(iedge,i) = inp.edg2poi(iedge,i);
  }
  for(int iface = 0; iface < inp.nface; iface++){
    for(int i = 0; i < getnnod2(curdeg); i++)
      fac2poi(iface,i) = inp.fac2poi(iface,i);
  }
  for(int ielem = 0; ielem < inp.nelem; ielem++){
    for(int i = 0; i < getnnod3(curdeg); i++)
      tet2poi(ielem,i) = inp.tet2poi(ielem,i);
  }

  edgHshTab = inp.edgHshTab;
  facHshTab = inp.facHshTab;

  for(int i = 0; i < idim; i++){
    bb[i][0] = inp.bb[i][0];
    bb[i][1] = inp.bb[i][1];
  }

  this->CAD = inp.CAD; 
  
  cfa2tag.set_n(0);
  bool ine1 = cfa2tag.allocate(METRIS_MAXTAGS, CAD.ncadfa);
  ced2tag.set_n(0);
  bool ine2 = ced2tag.allocate(METRIS_MAXTAGS, CAD.ncaded);
  cno2tag.set_n(0);
  bool ine3 = cno2tag.allocate(METRIS_MAXTAGS, CAD.ncadno);

  cfa2tag.set_n(METRIS_MAXTAGS);
  ced2tag.set_n(METRIS_MAXTAGS);
  cno2tag.set_n(METRIS_MAXTAGS);

  if(ine1) cfa2tag.fill(METRIS_MAXTAGS, CAD.ncadfa,0);
  if(ine2) ced2tag.fill(METRIS_MAXTAGS, CAD.ncaded,0);
  if(ine3) cno2tag.fill(METRIS_MAXTAGS, CAD.ncadno,0);

  ndomn = inp.ndomn;
  dom2tag.allocate(METRIS_MAXTAGS, ndomn);
  dom2tag.set_n(METRIS_MAXTAGS);
  dom2tag.fill(0);

  bpo2tag.fill(METRIS_MAXTAGS,nbpoi,0);
  poi2tag.fill(METRIS_MAXTAGS,npoin,0);
  edg2tag.fill(METRIS_MAXTAGS,nedge,0);
  fac2tag.fill(METRIS_MAXTAGS,nface,0);
  if(idim >= 3) tet2tag.fill(METRIS_MAXTAGS,nelem,0);

  for(int itag = 0; itag < METRIS_MAXTAGS; itag++) tag[itag] = 0;



  return *this;
}




// Update work arrays currently locked to entity type point
#define UPDATE_WORK_TYPE(DATATYPE, MESHTYPE, NAME)\
{\
  intAr1 &lwork_lock = lwork_lock_map<DATATYPE>();\
  MeshArray1D<MeshArray1D<DATATYPE>> &lwork = lwork_map<DATATYPE>();\
  int nwork = lwork.get_n();\
  for(int iarray = 0; iarray < nwork; iarray++){\
    int itype = lwork_lock[iarray];\
    if(itype != (int)MeshSize::MESHTYPE) continue;\
    lwork[iarray].allocate(m##NAME);\
    lwork[iarray].set_n(n##NAME);\
  }\
}

void MeshBase::set_nbpoi(int nbpoi){
  METRIS_ASSERT(Defaults::mem_growfac > 1); 

  nbpoi_ = nbpoi;
  if(nbpoi > mbpoi_) mbpoi_ = MAX(nbpoi, mbpoi_*Defaults::mem_growfac);

  bpo2ibi.allocate(mbpoi, nibi);
  bpo2ibi.set_n(nbpoi);

  bpo2rbi.allocate(mbpoi, nrbi);
  bpo2rbi.set_n(nbpoi);

  UPDATE_WORK_TYPE(int, BPoint, bpoi);
  UPDATE_WORK_TYPE(double, BPoint, bpoi);

  // Update tag arrays currently locked to entity type BPoint
  int ntag = tagarrs.get_n();
  for(int ii = 0; ii < ntag; ii++){
    int itype = tagarr_locks[ii];
    if(itype != (int)MeshSize::BPoint) continue;
    tagarrs[ii].allocate(mbpoi);
    int nn0 = tagarrs[ii].get_n();
    tagarrs[ii].set_n(nbpoi);
    int itag = itags[ii];
    for(int jj = nn0; jj < nbpoi; jj++) tagarrs[ii][jj] = itag - 1;
  }

  bpo2tag.allocate(METRIS_MAXTAGS, mbpoi);
  bpo2tag.set_n(METRIS_MAXTAGS);
}


void MeshBase::set_npoin(int npoin, bool skipallocf){
  METRIS_ASSERT(Defaults::mem_growfac > 1); 

  npoin_ = npoin;
  if(npoin > mpoin_) mpoin_ = MAX(npoin, mpoin_*Defaults::mem_growfac);


  poi2ent.allocate(mpoin, 2);
  poi2ent.set_n(npoin);

  poi2bpo.allocate(mpoin);
  poi2bpo.set_n(npoin);

  poi2tag.allocate(METRIS_MAXTAGS, mpoin);
  poi2tag.set_n(METRIS_MAXTAGS);

  UPDATE_WORK_TYPE(int, Point, poin);
  UPDATE_WORK_TYPE(double, Point, poin);


  // Update tag arrays currently locked to entity type point
  int ntag = tagarrs.get_n();
  for(int ii = 0; ii < ntag; ii++){
    int itype = tagarr_locks[ii];
    if(itype != (int)MeshSize::Point) continue;
    tagarrs[ii].allocate(mpoin);
    int nn0 = tagarrs[ii].get_n();
    tagarrs[ii].set_n(npoin);
    int itag = itags[ii];
    for(int jj = nn0; jj < npoin; jj++) tagarrs[ii][jj] = itag - 1;
  }

  

  poicstr.allocate(mpoin);
  poicstr.set_n(npoin);

  if(skipallocf) return;

  coord.allocate(mpoin, idim);
  coord.set_n(npoin);

}

void MeshBase::set_nedge(int nedge, bool skipallocf){
  METRIS_ASSERT(Defaults::mem_growfac > 1); 

  nedge_ = nedge;
  if(nedge > medge_) medge_ = MAX(nedge, medge_*Defaults::mem_growfac);

  edg2tag.allocate(METRIS_MAXTAGS, medge);
  edg2tag.set_n(METRIS_MAXTAGS);


  UPDATE_WORK_TYPE(int, Edge, edge);
  UPDATE_WORK_TYPE(double, Edge, edge);
  
  // Update tag arrays currently locked to entity type Edge
  int ntag = tagarrs.get_n();
  for(int ii = 0; ii < ntag; ii++){
    int itype = tagarr_locks[ii];
    if(itype != (int)MeshSize::Edge) continue;
    tagarrs[ii].allocate(medge);
    int nn0 = tagarrs[ii].get_n();
    tagarrs[ii].set_n(nedge);
    int itag = itags[ii];
    for(int jj = nn0; jj < nedge; jj++) tagarrs[ii][jj] = itag - 1;
  }

  edg2edg.allocate(medge, 2);
  edg2edg.set_n(nedge);

  edg2fac.allocate(medge);
  edg2fac.set_n(nedge);
  
  edgHshTab.reserve(nedge);

  if(skipallocf) return;

  edg2poi.allocate(medge, getnnod1(strdeg)); 
  edg2poi.set_n(nedge); 

  edg2ref.allocate(medge);
  edg2ref.set_n(nedge);

}

void MeshBase::set_nface(int nface, bool skipallocf){
  METRIS_ASSERT(Defaults::mem_growfac > 1); 

  nface_ = nface;
  if(nface > mface_) mface_ = MAX(nface, mface_*Defaults::mem_growfac);

  fac2tag.allocate(METRIS_MAXTAGS, mface);
  fac2tag.set_n(METRIS_MAXTAGS); 


  UPDATE_WORK_TYPE(int, Face, face);
  UPDATE_WORK_TYPE(double, Face, face);

  // Update tag arrays currently locked to entity type Face
  int ntag = tagarrs.get_n();
  for(int ii = 0; ii < ntag; ii++){
    int itype = tagarr_locks[ii];
    if(itype != (int)MeshSize::Face) continue;
    tagarrs[ii].allocate(mface);
    int nn0 = tagarrs[ii].get_n();
    tagarrs[ii].set_n(nface);
    int itag = itags[ii];
    for(int jj = nn0; jj < nface; jj++) tagarrs[ii][jj] = itag - 1;
  }


  fac2fac.allocate(mface, 3);
  fac2fac.set_n(nface);

  if(idim >= 3){
    fac2tet.allocate(mface, 2);
    fac2tet.set_n(nface);
    facHshTab.reserve(mface);
  }

  if(skipallocf) return; 

  fac2poi.allocate(mface, getnnod2(strdeg));
  fac2poi.set_n(nface);

  fac2ref.allocate(mface);
  fac2ref.set_n(nface);

}


void MeshBase::set_nelem(int nelem, bool skipallocf){
  METRIS_ASSERT(Defaults::mem_growfac > 1); 

  if(nelem <= 0) return;

  nelem_ = nelem;
  if(nelem > melem_) melem_ = MAX(nelem, melem_*Defaults::mem_growfac);


  UPDATE_WORK_TYPE(int, Tetra, elem);
  UPDATE_WORK_TYPE(double, Tetra, elem);

  //tet2ftg.allocate(melem);
  //tet2ftg.set_n(nelem);
  tet2tag.allocate(METRIS_MAXTAGS, melem);
  tet2tag.set_n(METRIS_MAXTAGS);

  tet2tet.allocate(melem, 4);
  tet2tet.set_n(nelem);

  // Update tag arrays currently locked to entity type Tetra
  int ntag = tagarrs.get_n();
  for(int ii = 0; ii < ntag; ii++){
    int itype = tagarr_locks[ii];
    if(itype != (int)MeshSize::Tetra) continue;
    tagarrs[ii].allocate(melem);
    int nn0 = tagarrs[ii].get_n();
    tagarrs[ii].set_n(nelem);
    int itag = itags[ii];
    for(int jj = nn0; jj < nelem; jj++) tagarrs[ii][jj] = itag - 1;
  }


  if(skipallocf) return;

  tet2poi.allocate(melem, getnnod3(strdeg));
  tet2poi.set_n(nelem);

  tet2ref.allocate(melem);
  tet2ref.set_n(nelem);


}

void MeshBase::set_nentt(int tdimn, int nentt, bool skipallocf){
  switch(tdimn){
  case(-1):
    set_nbpoi(nentt);
    break;
  case(0):
    set_npoin(nentt,skipallocf);
    break;
  case(1):
    set_nedge(nentt,skipallocf);
    break;
  case(2):
    set_nface(nentt,skipallocf);
    break;
  case(3):
    set_nelem(nentt,skipallocf);
    break;
  }
}

#undef UPDATE_WORK_TYPE



//// This'll average between periodic sides in that case
//void MeshBase::getpoinormal(int ipoin, int iref, double *norpoi){
//  METRIS_ASSERT(idim == 3);
//  METRIS_ASSERT(isboundary_faces());
//
//  int ibpoi = msh.poi2bpo[ipoin];
//  METRIS_ASSERT(ibpoi >= 0);
//
//  if(EGADS_context != NULL) METRIS_THROW_MSG(TODOExcept(),
//    "getpoinormal with egads context")
//
//  for(int ii = 0; ii < 3; ii++) norpoi[ii] = 0;
//  double norfa[3];
//  double nrmtot = 0;
//
//  int ibpo2 = ibpoi;
//  do{
//    int ityp = msh.bpo2ibi(ibpo2,1);
//    if(ityp == 2){
//      int iface = msh.bpo2ibi(ibpo2,2);
//      int iref2 = msh.fac2ref[iface];
//      if(iref2 == iref){  
//        getnorfacP1(fac2poi[iface],coord,norfa);
//        for(int ii = 0; ii < 3; ii++) norpoi[ii] += norfa[ii];
//        nrmtot += sqrt(getnrml2<3>(norfa));
//      }
//    }
//    ibpo2 = msh.bpo2ibi(ibpo2,3);
//  }while(ibpo2 >= 0 && ibpo2 != ibpoi);
//
//  for(int ii = 0; ii < 3; ii++) norpoi[ii] /= nrmtot;
//}


void MeshBase::get_nMeshSize(MeshSize itype, int* nn, int* mm){
  int nentt = -1, mentt = -1;
  switch(itype){
    case MeshSize::Point:
      nentt = npoin;
      mentt = mpoin;
      break;
    case MeshSize::Edge:
      nentt = nedge;
      mentt = medge;
      break;
    case MeshSize::Face:
      nentt = nface;
      mentt = mface;
      break;
    case MeshSize::Tetra:
      nentt = nelem;
      mentt = melem;
      break;
    case MeshSize::BPoint:
      nentt = nbpoi;
      mentt = mbpoi;
      break;
    case MeshSize::Domain:
      nentt = ndomn;
      mentt = ndomn;
      break;
    case MeshSize::CADNode:
      nentt = CAD.ncadno;
      mentt = CAD.ncadno;
      break;
    case MeshSize::CADEdge:
      nentt = CAD.ncaded;
      mentt = CAD.ncaded;
      break;
    case MeshSize::CADFace:
      nentt = CAD.ncadfa;
      mentt = CAD.ncadfa;
      break;
    default:
      METRIS_THROW_MSG(WArgExcept(), "get_nMeshSize: invalid entity type = "<< (int) itype);
  }
  *nn = nentt;
  *mm = mentt;
}

// To be used both by get_iwork and get_tagarray
// Returns the index of the locked array
template<typename T>
inline int get_locked_array(int nn, MeshArray1D<T> &array_pool, intAr1 &array_locks){
  int npool = array_pool.get_n();
  int ibest = -1, nbest = -1;
  for(int iarray = 0; iarray < npool; iarray++){
    if(array_locks[iarray] > 0) continue;
    int array_size = array_pool[iarray].size();
    if(array_size < nn) continue;

    if(array_size < nbest || nbest < 0){
      ibest = iarray;
      nbest = array_size;
    }
  }

  return ibest;
}


template<typename T>
WorkArray1D<T> MeshBase::get_work(int nn){
  int iarray = get_locked_array<MeshArray1D<T>>(nn, lwork_map<T>(), lwork_lock_map<T>());
  if(iarray < 0){
    iarray = lwork_map<T>().get_n();
    lwork_map<T>().inc_n();
    lwork_lock_map<T>().inc_n();
    lwork_map<T>()[iarray].allocate(nn);
  }
  lwork_map<T>()[iarray].set_n(nn);
  lwork_lock_map<T>()[iarray] = (int) MeshSize::Untracked;
  return WorkArray1D<T>(*this, iarray, lwork_map<T>()[iarray]);
}

template intWrkAr1 MeshBase::get_work<int   >(int nn);
template dblWrkAr1 MeshBase::get_work<double>(int nn);


template<typename T>
WorkArray1D<T> MeshBase::get_work(MeshSize itype){
  int nentt = -1, mentt = -1;
  get_nMeshSize(itype, &nentt, &mentt);

  int iarray = get_locked_array<MeshArray1D<T>>(mentt, lwork_map<T>(), lwork_lock_map<T>());
  if(iarray < 0){
    iarray = lwork_map<T>().get_n();
    lwork_map<T>().inc_n();
    lwork_lock_map<T>().inc_n();
    lwork_map<T>()[iarray].allocate(mentt);
  }
  lwork_map<T>()[iarray].set_n(nentt);
  lwork_lock_map<T>()[iarray] = (int) itype;
  return WorkArray1D<T>(*this, iarray, lwork_map<T>()[iarray]);
}

template intWrkAr1 MeshBase::get_work<int   >(MeshSize itype);
template dblWrkAr1 MeshBase::get_work<double>(MeshSize itype);


intWrkAr1 MeshBase::get_iwork(int nn){
  return get_work<int>(nn);
}
dblWrkAr1 MeshBase::get_rwork(int nn){
  return get_work<double>(nn);
}
intWrkAr1 MeshBase::get_iwork(MeshSize size){
  return get_work<int>(size);
}
dblWrkAr1 MeshBase::get_rwork(MeshSize size){
  return get_work<double>(size);
}


TagArray MeshBase::get_tagarray(MeshSize itype){
  METRIS_ASSERT((int) itype > 0);

  int nentt = -1, mentt = -1;
  get_nMeshSize(itype, &nentt, &mentt);
  
  int iarray = get_locked_array<intAr1>(mentt, tagarrs, tagarr_locks);
  
  printf("Got iarray = %d \n",iarray);
  if(iarray < 0){
    iarray = tagarrs.get_n();

    tagarrs.inc_n();
    tagarr_locks.inc_n();
    itags.stack(0);

    tagarrs[iarray].allocate(mentt);
    tagarrs[iarray].set_n(nentt);
    tagarrs[iarray].fill(0);

    printf("Allocated iarray %d nmem %d size %d check : %d %d \n",iarray,mentt,nentt,
           tagarrs[iarray].size(),tagarrs[iarray].get_n());
  }

  itags[iarray]++;
  tagarr_locks[iarray] = (int) itype;
  tagarrs[iarray].set_n(nentt); // could be originally larger than needed

  return TagArray(*this, iarray, tagarrs[iarray], itags[iarray]);
}







void MeshBase::setBasis(FEBasis ibasis_){
  if(ibasis_ == FEBasis::Lagrange){
    setLagrange();
  }else if(ibasis_ == FEBasis::Bezier){
    setBezier();
  }else{
    METRIS_THROW_MSG(WArgExcept(), "Invalid basis for coordinates (MeshBase)")
  }
}

void MeshBase::setLagrange(){
  if(this->ibasis == FEBasis::Lagrange) return; 

  CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){
    if(ideg == this->curdeg){
      if(this->idim == 2){
        setFieldLagrange<ideg,2>(*this,this->coord);
      }else{
        setFieldLagrange<ideg,3>(*this,this->coord);
      }
    }
  }CT_FOR1(ideg);

  this->ibasis = FEBasis::Lagrange;
  return;
}

void MeshBase::setBezier(){
  if(this->ibasis == FEBasis::Bezier) return;

  CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){
    if(ideg == this->curdeg){
      if(this->idim == 2){
        setFieldBezier<ideg,2>(*this,this->coord);
      }else{
        setFieldBezier<ideg,3>(*this,this->coord);
      }
    }
  }CT_FOR1(ideg);

  this->ibasis = FEBasis::Bezier;
  return;
}




intAr2r& MeshBase::ent2tag(int tdimn){
  switch(tdimn){
  case(1):
    return edg2tag;
  case(2):
    return fac2tag;
  case(3):
    return tet2tag;
  default:
    METRIS_THROW_MSG(WArgExcept(),"ent2tag (1) tdimn not in range = "<<tdimn);
  }
}

const intAr2r& MeshBase::ent2tag(int tdimn) const{
  switch(tdimn){
  case(1):
    return edg2tag;
  case(2):
    return fac2tag;
  case(3):
    return tet2tag;
  default:
    METRIS_THROW_MSG(WArgExcept(),"ent2tag (1) tdimn not in range = "<<tdimn);
  }
}


template<int tdimn>
      intAr2r& MeshBase::ent2tag(){
  static_assert(tdimn >= 1 && tdimn <= 3);
  if constexpr(tdimn == 1){
    return edg2tag;
  }else if(tdimn == 2){
    return fac2tag;
  }else{
    return tet2tag;
  }
}
template<int tdimn>
const intAr2r& MeshBase::ent2tag() const{
  static_assert(tdimn >= 1 && tdimn <= 3);
  if constexpr(tdimn == 1){
    return edg2tag;
  }else if(tdimn == 2){
    return fac2tag;
  }else{
    return tet2tag;
  }
}

template intAr2r& MeshBase::ent2tag<1>();
template intAr2r& MeshBase::ent2tag<2>();
template intAr2r& MeshBase::ent2tag<3>();
template const intAr2r& MeshBase::ent2tag<1>()const;
template const intAr2r& MeshBase::ent2tag<2>()const;
template const intAr2r& MeshBase::ent2tag<3>()const;

  // points to ced2tag, cfa2tag or dom2tag
intAr2r& MeshBase::ref2tag(int tdimn){
  METRIS_ASSERT(tdimn >= 1 && tdimn <= 3);
  if(tdimn == 1){
    return ced2tag;
  }else if(tdimn == 2){
    return cfa2tag;
  }else{
    return dom2tag;
  }
}

const intAr2r& MeshBase::ref2tag(int tdimn) const{
  METRIS_ASSERT(tdimn >= 1 && tdimn <= 3);
  if(tdimn == 1){
    return ced2tag;
  }else if(tdimn == 2){
    return cfa2tag;
  }else{
    return dom2tag;
  }
}








} // End namespace




