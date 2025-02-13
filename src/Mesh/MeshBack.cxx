//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Mesh/MeshBack.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../ho_constants.hxx"
#include "../aux_utils.hxx"
#include "../aux_timer.hxx"
#include "../low_eval.hxx"
#include "../low_geo.hxx"
#include "../mprintf.hxx"
#include "../msh_intrinsicmet.hxx"
#include "../io_libmeshb.hxx"
#include "../MetricField/msh_checkmet.hxx"
#include "../API/MetrisAPI.hxx"
#include "../linalg/det.hxx"


namespace Metris{


int MeshBack::newpoitopo(int tdimn, int ientt){
  return MeshBase::newpoitopo(tdimn, ientt);
}

void MeshBack::copyConstants(const MeshBase &msh){
  MeshMetric<MetricFieldFE>::copyConstants(msh);
  mpoin_ = npoin;
  mbpoi_ = nbpoi;
  medge_ = nedge;
  mface_ = nface;
  melem_ = nelem;
  int mbpo_guess = edgnpps[strdeg]*nedge + facnpps[strdeg]*nface;
  if(nbpoi == 0) mbpoi_ = mbpo_guess;
}
void MeshBack::readConstants(int64_t libIdx, int usrMinDeg){
  MeshMetric<MetricFieldFE>::readConstants(libIdx,usrMinDeg);
  mpoin_ = npoin;
  mbpoi_ = nbpoi;
  medge_ = nedge;
  mface_ = nface;
  melem_ = nelem;
  int mbpo_guess = edgnpps[strdeg]*nedge + facnpps[strdeg]*nface;
  if(nbpoi == 0) mbpoi_ = mbpo_guess;
}
void MeshBack::readConstants(const MetrisAPI &data, int usrMinDeg){
  MeshMetric<MetricFieldFE>::readConstants(data,usrMinDeg);
  mpoin_ = npoin;
  mbpoi_ = nbpoi;
  medge_ = nedge;
  mface_ = nface;
  melem_ = nelem;
  int mbpo_guess = edgnpps[strdeg]*nedge + facnpps[strdeg]*nface;
  if(nbpoi == 0) mbpoi_ = mbpo_guess;
}

dblAr1& MeshBack::ent2dev(int tdimn){
  switch(tdimn){
  case(1):
    return edg2dev;
  case(2):
    return fac2dev;
  default:
    METRIS_THROW_MSG(WArgExcept(),"tdimn not in range = "<<tdimn);
  }
}

const dblAr1& MeshBack::ent2dev(int tdimn) const{
  switch(tdimn){
  case(1):
    return edg2dev;
  case(2):
    return fac2dev;
  default:
    METRIS_THROW_MSG(WArgExcept(),"tdimn not in range = "<<tdimn);
  }
}


void MeshBack::set_nedge(int nedge, bool skip_allocf){
  MeshMetric<MetricFieldFE>::set_nedge(nedge, skip_allocf);

  // Allocate edg2dev
  edg2dev.allocate(medge);
  edg2dev.set_n(nedge);

  //edg2con.allocate(medge,idim + 1);
  //edg2con.set_n(nedge);
}

void MeshBack::set_nface(int nface, bool skip_allocf){
  MeshMetric<MetricFieldFE>::set_nface(nface, skip_allocf);
  
  // Allocate fac2dev
  fac2dev.allocate(mface);
  fac2dev.set_n(nface);
}

void MeshBack::initialize(MetrisAPI *data, 
  MetrisParameters &param){
  MeshBase::initialize(data,param);
  setBasis(FEBasis::Bezier);

  GETVDEPTH((*this));

  if(DOPRINTS2()) writeMesh("mshDATAback", *this);

  met.forceBasisFlag(FEBasis::Undefined);

  // If no input file: 
  if(!param.inpMet && data == NULL || (data != NULL && !data->imet)){
    // If analytical:
    if(param.anamet_ptr != NULL || param.ianamet >= 1){
      auto anamet = param.anamet_ptr;
      if(param.ianamet >= 0) anamet = (idim == 2 ? __ANAMET2D[param.ianamet-1] : __ANAMET3D[param.ianamet-1]);

      met.forceBasisFlag(FEBasis::Lagrange);
      met.forceSpaceFlag(MetSpace::Exp);
      int nnmet = (idim * (idim+1))/2;
      met.rfld.allocate(npoin,nnmet);
      met.rfld.set_n(npoin);
      FEBasis ibas0 = ibasis;
      setBasis(FEBasis::Lagrange);
      for(int ipoin = 0; ipoin < npoin; ipoin++){
        anamet(NULL, coord[ipoin], param.metScale, 0, met[ipoin], NULL);
      } 
      setBasis(ibas0);

    // Else intrinsic:
    }else{
      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == curdeg){
        if(DOPRINTS1()) std::cout << "(back)  - Compute intrinsic metric field \n";
        double t0, t1;

        if(DOPRINTS1()) t0 = get_wall_time();
        getMetMesh<MetricFieldFE,ideg>(param,*this);
        if(DOPRINTS1()) t1 = get_wall_time();
        if(DOPRINTS1()) std::cout << "(back)  - Done time = "<<t1-t0<<std::endl;
      }}CT_FOR1(ideg);
    }

    if(DOPRINTS2()) met.writeMetricFile("backmet.solb");
    #ifndef NDEBUG
      checkMet(*this);
    #endif

  }else if(param.inpMet){
    
    if(data != NULL)METRIS_ENFORCE_MSG(!data->imet, "Metric specified both in data and file");

    met.readMetricFile(param.metFileName);

  }else if(data != NULL && data->imet){

    met.rfld = std::move(data->metfld); 
    //int nnmet = (idim * (idim + 1)) / 2;
    //for(int ipoin = 0; ipoin < npoin; ipoin++){
    //  for(int ii = 0; ii < nnmet; ii++){
    //    met(ipoin,ii) = data->metfld[ipoin][ii];
    //  }
    //}
    met.forceBasisFlag(data->metbasis);
    met.forceSpaceFlag(data->metspace);

    if(DOPRINTS2()) met.writeMetricFile("metDATA");

  }else{
    METRIS_THROW_MSG(WArgExcept(), "No metric info for back ?? ");
  }


  // Correct metric if needed 
  //met.correctMetric();

  met.setSpace(MetSpace::Log);
  met.setBasis(FEBasis::Bezier);

  if(param.scaleMet){
    CPRINTF1("-- Back scaling metric by %15.7e\n", param.metScale);
    met.normalize(param.metScale);
  }


  #if 0
  // Compute edg2con 
  if(this->CAD()){

    for(int iedge = 0; iedge < nedge; iedge++){
      if(isdeadent(iedge,edg2poi)) continue;
      int iref = edg2ref[iedge];
      ego obj = CAD.cad2edg[iref];

      for(int ii = 0; ii < 2; ii++){
        int ipoin = edg2poi(iedge,ii);

        // Get CAD tangent 
        int ibpoi = poi2bpo[ipoin];
        METRIS_ASSERT(ibpoi >= 0);
        ibpoi = getref2bpo(*this,ibpoi,iref,1);
        double result[18];
        int ierro = EG_evaluate(obj, bpo2rbi[ibpoi], result);
        METRIS_ENFORCE(ierro == 0);
        double *du = &result[3];

        double norpoi[3];
        ierro = getnorpoiCAD(msh, ipoin, edgorient, norpoi);



      }// for int ii 

    }// for int iedge 

  }else{
    printf("## WARNING: No CAD: no normal control in loc\n");
  }
  #endif

  // Compute geodev 
  
  // For edges, the tangent deviation is computed at the vertices. 
  // check edges only belong to one loop
  geodev[0] = -1;
  geodev[1] = -1;
  int imax[2] = {-1};
  if(this->CAD()){

    for(int tdim = 1; tdim <= 2; tdim++){

      if(tdim == 1 && !isboundary_edges()){
        CPRINTF1("-- Skipping edges: not linked to CAD\n");
        continue;
      }
      if(tdim == 2 && !isboundary_faces()){
        CPRINTF1("-- Skipping face: not linked to CAD\n");
        continue;
      }

      int nentt = this->nentt(tdim);
      const intAr2& ent2poi = this->ent2poi(tdim);
      const intAr1& ent2ref = this->ent2ref(tdim);
      dblAr1& ent2dev = this->ent2dev(tdim);

      int nCADsing = 0;

      tag[0]++;
      for(int ientt = 0; ientt < nentt; ientt++){
        INCVDEPTH((*this));
        if(isdeadent(ientt,ent2poi)) continue;
        int iref = ent2ref[ientt];
        ego obj = tdim == 1 ? CAD.cad2edg[iref] : CAD.cad2fac[iref];
        int mtype = obj->mtype;


        ent2dev[ientt] = -1;

        for(int ii = 0; ii < tdim + 1; ii++){
          int ipoin = ent2poi(ientt,ii);

          // Get CAD tangent 
          int ibpoi = poi2ebp(ipoin,tdim,-1,iref);
          METRIS_ASSERT(ibpoi >= 0);

          double result[18];
          int ierro = EG_evaluate(obj, bpo2rbi[ibpoi], result);
          METRIS_ENFORCE(ierro == 0);
          double *dirCAD;

          double buf[3];
          if(tdim == 1){ 
            // If edge, simply take the tangent. 
            dirCAD = &result[3];
          }else{
            // Otherwise, compute normal
            vecprod(&result[3], &result[6], buf);
            dirCAD = buf;
          }

          CPRINTF1("ientt %d iref %d ipoin = %d ibpoi = %d : ",
                               ientt,iref,ipoin,ibpoi);
          if(DOPRINTS1()) intAr1(nibi,bpo2ibi[ibpoi]).print();

          // Get elt direction (tangent if dim = 1, normal otherwise)
          double dirent[3], dum[3];
          double bary[3] = {0}; // not typo
          bary[ii] = 1;

          if(tdim == 1){
            CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == curdeg){
              if(idim == 2){
                eval1<2,ideg>(coord, ent2poi[ientt],
                              getBasis(), DifVar::Bary, DifVar::None,
                              bary, dum, dirent, NULL);
              }else{
                eval1<3,ideg>(coord, ent2poi[ientt],
                              getBasis(), DifVar::Bary, DifVar::None,
                              bary, dum, dirent, NULL);
              }
            }}CT_FOR1(ideg); 
          }else{
            double jmat[2][3];
            CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == curdeg){
              if(idim == 3){
                eval2<3,ideg>(coord, ent2poi[ientt],
                              getBasis(), DifVar::Bary, DifVar::None,
                              bary, dum, jmat[0], NULL);
              }else{
                METRIS_THROW_MSG(TODOExcept(),"How are we in idim = "<<idim<<" in face bdry case?")
              }
            }}CT_FOR1(ideg); 
            vecprod(jmat[0], jmat[1], dirent);
          }

          // Normalize both 
          int iCADsing = 0;
          if(idim == 2){
            iCADsing = normalize_vec<2>(dirCAD);
            METRIS_ENFORCE(normalize_vec<2>(dirent) == 0); 
          }else{
            iCADsing = normalize_vec<3>(dirCAD);
            METRIS_ENFORCE(normalize_vec<3>(dirent) == 0); 
          }

          if(iCADsing){
            if(poi2tag(0,ipoin) >= tag[0]) continue;
            poi2tag(0,ipoin) = tag[0];
            CPRINTF1("## CAD normal singular at point %d \n",ipoin);
            nCADsing++;
            continue;
          }


          double dtprd = idim == 2 ? getprdl2<2>(dirCAD,dirent)
                                   : getprdl2<3>(dirCAD,dirent);

          // Force dtprd to be positive. We can do this because entities have been 
          // oriented, so there's nothing to check here. 
          dtprd = abs(dtprd);

          if(DOPRINTS1()){
            CPRINTF1(" - ientt %d ipoin %d mtype %d dtprd %f dirCAD = ",ientt,ipoin,mtype,dtprd);
            dblAr1(idim,dirCAD).print();
            CPRINTF1(" - dirent = ");
            dblAr1(idim,dirent).print();
            CPRINTF1(" - (u,v) = ");
            dblAr1(nrbi,bpo2rbi[ibpoi]).print();
            if(1 - abs(dtprd) >= 0.9){
              CPRINTF1("## LARGE GEODEV \n");
              for(int ibpo0 = poi2bpo[ipoin]; ibpo0 >= 0; ibpo0 = bpo2ibi(ibpo0,3)){
                CPRINTF1("   - ibpoi %d : ",ibpo0);
                intAr1(nibi,bpo2ibi[ibpo0]).print();
                MPRINTF(" (u,v) = %f %f \n",bpo2rbi(ibpo0,0),bpo2rbi(ibpo0,1));
              }
            }
          }

          METRIS_ASSERT_MSG(dtprd >= 1.0e-16,"zero dtprd = "<<dtprd);


          double dev = 1 - abs(dtprd);

          ent2dev[ientt] = MAX(ent2dev[ientt], dev);

          if(dev >= geodev[tdim-1]){
            imax[tdim-1] = ientt;
            geodev[tdim-1] = dev;
          }

        }// for int ii 
        CPRINTF1(" - ientt %d final dev = %15.7e \n",ientt,ent2dev[ientt]);

      }
      CPRINTF1("-- %d singular CAD normals / tangents at topo dim %d entities\n",nCADsing,tdim);
    }

  }// if CAD()
  else{
    geodev[0] = 1;
    geodev[1] = 1;
  }
  
  CPRINTF1("-- Computed max geodev, edges: %15.7e at %d faces %15.7e at %d\n",
                               geodev[0], imax[0], geodev[1], imax[1]);



}


double MeshBack::getMetComplexity(){

  double volM = getDomainVolume();

  //std::cout<<"Domain volume = " << vol0 << "\n";// << " aniso = "<<volM<<"\n";
  //// 1) If an element is unit, its volume under the metric is vK0
  //// Thus the volume of a unit mesh is Nelem * vK0 
  //// 2) The volume under cM, v_{cM}, of the mesh is c^{1/2} v_M
  //// Conclusion) c such that the mesh is unit under cM verifies 
  //// c^{1/2} v_M = Nelem * vK_0
  //double coeff = (nelem * Constants::vK0[gdim] / volM)
  //             * (nelem * Constants::vK0[gdim] / volM);

  //std::cout<<"Metric normalization factor = "<<coeff<<" \n ";

  //met.normalize(coeff);

  // If an element is unit, its volume under the metric is vK0
  // Thus we require volM / vK0 unit elements to fill the domain.
  double comp = volM / Constants::vK0[idim];

  std::cout<<"Estimated tar nelem = "<<comp<<" current = "<<nelem<<"\n";

  return comp;

  //double ret = -1;
  //CT_FOR0_INC(2,3,idim_){if(idim_ == this->idim){
  //  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){
  //    if(ideg == this->curdeg) ret = getMetComplexity0<idim_,ideg>();
  //  }CT_FOR1(ideg);
  //}}CT_FOR1(idim_);
  //return ret;
}

//void MeshBack::allocate(){
//  setMpoiToMent();
//  MeshBase::allocate();
//
//  this->met.allocate();
//}

} // End namespace
