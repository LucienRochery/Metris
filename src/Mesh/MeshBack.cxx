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
#include "../msh_intrinsicmet.hxx"
#include "../io_libmeshb.hxx"
#include "../MetricField/msh_checkmet.hxx"
#include "../API/MetrisAPI.hxx"


namespace Metris{


int MeshBack::newpoitopo(int tdimn, int ientt){
  return MeshBase::newpoitopo(tdimn, ientt);
}

void MeshBack::copyConstants(const MeshBase &msh, int MAX_DEG){
  MeshMetric<MetricFieldFE>::copyConstants(msh,MAX_DEG);
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

void MeshBack::set_nedge(int nedge, bool skip_allocf){
  MeshMetric<MetricFieldFE>::set_nedge(nedge, skip_allocf);

  // Allocate edg2dev
  edg2dev.allocate(medge);
  edg2dev.set_n(nedge);

  //edg2con.allocate(medge,idim + 1);
  //edg2con.set_n(nedge);
}


void MeshBack::initialize(MetrisAPI *data, 
  #ifdef NDEBUG 
  const 
  #endif
  MetrisParameters &param){
  MeshBase::initialize(data,param);
  setBasis(FEBasis::Bezier);

  if(param.iverb >= 1) writeMesh("mshDATAback", *this);

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
        if(param.iverb >= 1) std::cout << "(back)  - Compute intrinsic metric field \n";
        double t0, t1;

        if(param.iverb >= 1) t0 = get_wall_time();
        getMetMesh<MetricFieldFE,ideg>(param,*this);
        if(param.iverb >= 1) t1 = get_wall_time();
        if(param.iverb >= 1) std::cout << "(back)  - Done time = "<<t1-t0<<std::endl;
      }}CT_FOR1(ideg);
    }

    if(param.iverb >= 3) met.writeMetricFile("backmet.solb");
    #ifndef NDEBUG
      checkMet(*this);
    #endif

  }else if(param.inpMet){
    
    if(data != NULL)METRIS_ENFORCE_MSG(!data->imet, "Metric specified both in data and file");

    int metdim;
    int64_t libIdxMet = MetrisOpenMeshFile<GmfRead>(param.metFileName.c_str(), &metdim);

    METRIS_ENFORCE_MSG(idim == metdim, "Back and metric dimensions must agree");

    met.readMetricFile(libIdxMet);
    GmfCloseMesh(libIdxMet);

  }else if(data != NULL && data->imet){

    met.rfld = std::move(data->metfld); 
    //int nnmet = (idim * (idim + 1)) / 2;
    //for(int ipoin = 0; ipoin < npoin; ipoin++){
    //  for(int ii = 0; ii < nnmet; ii++){
    //    met[ipoin][ii] = data->metfld[ipoin][ii];
    //  }
    //}
    met.forceBasisFlag(data->metbasis);
    met.forceSpaceFlag(data->metspace);

    if(param.iverb >= 1) met.writeMetricFile("metDATA", MetSpace::Exp);

  }else{
    METRIS_THROW_MSG(WArgExcept(), "No metric info for back ?? ");
  }

  // Correct metric if needed 
  //met.correctMetric();

  met.setSpace(MetSpace::Log);
  met.setBasis(FEBasis::Bezier);

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
  int imax = -1;
  if(this->CAD()){

    for(int iedge = 0; iedge < nedge; iedge++){
      if(isdeadent(iedge,edg2poi)) continue;
      int iref = edg2ref[iedge];
      ego obj = CAD.cad2edg[iref];

      edg2dev[iedge] = -1;

      for(int ii = 0; ii < 2; ii++){
        int ipoin = edg2poi(iedge,ii);

        // Get CAD tangent 
        int ibpoi = poi2bpo[ipoin];
        METRIS_ASSERT(ibpoi >= 0);

        ibpoi = getref2bpo(*this,ibpoi,iref,1);
        METRIS_ASSERT(bpo2ibi(ibpoi,0) == ipoin);
        double result[18];
        int ierro = EG_evaluate(obj, bpo2rbi[ibpoi], result);
        METRIS_ENFORCE(ierro == 0);
        double *tanCAD = &result[3];

        for(int ii = 0; ii < 18; ii++) printf(" %d : %15.7e \n",ii,result[ii]);

        if(param.iverb >= 3) printf("iedge %d iref %d ipoin = %d ibpoi = %d : ",
                             iedge,iref,ipoin,ibpoi);
        if(param.iverb >= 3) intAr1(nibi,bpo2ibi[ibpoi]).print();

        // Get elt tangent
        double tanedg[3], dum[3];
        double bary[2] = {1.0 - ii, (double)ii};
        CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == curdeg){
          if(idim == 2){
            eval1<2,ideg>(coord, edg2poi[iedge],
                          getBasis(), DifVar::Bary, DifVar::None,
                          bary, dum, tanedg, NULL);
          }else{
            eval1<3,ideg>(coord, edg2poi[iedge],
                          getBasis(), DifVar::Bary, DifVar::None,
                          bary, dum, tanedg, NULL);
          }
        }}CT_FOR1(ideg); 

        // Normalize both 
        double normCAD, normedg; 
        if(idim == 2){
          normCAD = getnrml2<2>(tanCAD);
          normedg = getnrml2<2>(tanedg);
        }else{
          normCAD = getnrml2<3>(tanCAD);
          normedg = getnrml2<3>(tanedg);
        }
        METRIS_ENFORCE(normCAD >= 1.0e-32);
        METRIS_ENFORCE(normedg >= 1.0e-32);

        for(int ii = 0; ii < idim; ii++) tanCAD[ii] /= sqrt(normCAD);
        for(int ii = 0; ii < idim; ii++) tanedg[ii] /= sqrt(normedg);

        double dtprd = idim == 2 ? getprdl2<2>(tanCAD,tanedg)
                                 : getprdl2<3>(tanCAD,tanedg);

        // Force dtprd to be positive. We can do this because edges have been 
        // oriented, so there's nothing to check here.
        dtprd = abs(dtprd);

        if(param.iverb >= 3){
          printf("  - iedge %d ipoin %d dtprd %f tanCAD = ",ipoin,iedge,dtprd);
          dblAr1(idim,tanCAD).print();
          printf("  - tanedg = ");
          dblAr1(idim,tanedg).print();
        }

        METRIS_ASSERT_MSG(dtprd >= 1.0e-16,"zero dtprd = "<<dtprd);

        double dev = 1 - abs(dtprd);

        edg2dev[iedge] = MAX(edg2dev[iedge], dev);


        if(dev >= geodev[0]){
          imax = iedge;
          geodev[0] = dev;
        }

      }// for int ii 
      if(param.iverb >= 3) printf(" - iedge %d final dev = %15.7e \n",iedge,edg2dev[iedge]);
    }

  }// if msh.CAD()
  else{
    geodev[0] = 1;
  }
  
  if(param.iverb >= 2) printf("-- Computed max tandev edges, got %15.7e at %d\n",
                               geodev[0], imax);


  // Compute maximum normal deviation for localization tolerance 
  if(isboundary_faces()) METRIS_THROW_MSG(TODOExcept(), 
                     "Implement geodev for triangles in 3D: \n - Add fac2dev \n" 
                     " - Don't forget override set_nface()")
  geodev[1] = -1;

}



void MeshBack::getEnttMemCosts(int *memCostPpoi, int *memCostPbpo, int *memCostPedg, int *memCostPfac, int *memCostPelt){
  
	MeshBase::getEnttMemCosts(memCostPpoi,memCostPbpo,memCostPedg,memCostPfac,memCostPelt);

  int memCostPdbl = sizeof(double);

  int nnmet = (this->idim*(this->idim + 1))/2;

  *memCostPpoi += nnmet*memCostPdbl ;/* met     */
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
