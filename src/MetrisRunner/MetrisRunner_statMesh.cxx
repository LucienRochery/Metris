//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "MetrisRunner.hxx"
#include "MeshStat.hxx"

#include "../Mesh/Mesh.hxx"

#include "../msh_lenedg.hxx"
#include "../aux_histogram.hxx"
#include "../quality/msh_metqua.hxx"
#include "../low_geo/ccoef.hxx"
#include "../utils/mprintf.hxx"

#include "../SolutionField/SolutionField.hxx"
#include "../SolutionField/minInterpError.hxx"
#include "../SolutionField/interpError.hxx"

namespace Metris{

void MetrisRunner::statMesh(int tdim,MeshStat* stat){
  if(this->metricFE){
    statMesh0<MetricFieldFE>(tdim, stat);
  }else{
    statMesh0<MetricFieldAnalytical>(tdim, stat);
  }
}

template<class MFT>
void MetrisRunner::statMesh0(int tdim, MeshStat* stat_){


  Mesh<MFT> &msh = *( (Mesh<MFT>*) msh_g );
  GETVDEPTH(msh.param);


  if(tdim <= 0) tdim = msh.get_tdim();

  if(!DOPRINTS1() && stat_ == NULL) return;

  MeshStat stat;

  msh.cleanup();
  msh.met.setSpace(MetSpace::Exp);

  static const dblAr1 qbnds = {1.0e-8, 0.1};
  static const dblAr1 jbnds = {msh.param->jtol, 1};
  static const dblAr1 lenbds = {1.0/sqrt(2.0), sqrt(2.0)};

  // Length
  intAr2 ilned;
  dblAr1 rlned;
  for(int tdim_ = 1; tdim_ <= tdim; tdim_++){
    lenStat dummy;
    getLengthEdges<MFT>(msh,tdim_ ,-1,ilned,rlned,dummy,LenTyp::GeoSiz);
    stat.setLength(tdim_,rlned);
    if(DOPRINTS2()) print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length ("+std::to_string(tdim_)+"D)");
  }

  // Quality P1
  dblAr1 rquel;
  double qua_surf_wt_normal = msh.param->qua_surf_wt_normal;
  msh.param->qua_surf_wt_normal = 0;
  for(AsDeg asdeg : {AsDeg::P1, AsDeg::Pk}){
    if(asdeg == AsDeg::Pk && msh.curdeg == 1) continue;
    for(int tdim_ = 2; tdim_ <= tdim; tdim_++){
      bool iinva;
      double qmin, qmax, qavg;
      getmetquamesh<MFT>(msh,tdim_,asdeg,asdeg,
                        &iinva,&qmin,&qmax,&qavg,&rquel);
      stat.setQuality(tdim_,rquel,asdeg);
      std::string degstr = asdeg == AsDeg::P1 ? "P1" : "Pk";
      if(DOPRINTS2()) print_histogram(msh,rquel,IntrpTyp::Geometric,qbnds,"q","Quality "+degstr+" ("+std::to_string(tdim_)+"D)");
    }
  }
  msh.param->qua_surf_wt_normal = qua_surf_wt_normal;



  #if METRIS_MAX_DEG >= 2
  if(param->anaSol && msh.idim == 2 && DOPRINTS1()){
    SolutionFieldAnalytical sol(msh);
    sol.setAnalyticalSolution(param->ianasol);
    double errGlo = -1;
    //CT_FOR0_INC(2,3,idim){if(idim == msh.idim){
    CT_FOR0_INC(1,2,ideg){if(ideg == msh.curdeg){
    CT_FOR0_INC(1,2,pnorm){if(pnorm == param->intp_pnorm){
    CT_FOR0_INC(1,2,pdeg){if(pdeg == param->intp_pdeg){
      errGlo = interpErrGlo<2,ideg,pdeg,pnorm,true>(sol);
    }}CT_FOR1(pdeg);
    }}CT_FOR1(pnorm);
    }}CT_FOR1(ideg);
    //}}CT_FOR1(idim);
    CPRINTF1("-- Analytical solution interpolation error pnorm {} pdeg {}: {}\n",
             param->intp_pdeg, param->intp_pnorm, errGlo);
  }
  #endif

  if(DOPRINTS1()) stat.print("",param->logFile);
  if(stat_) *stat_ = stat;
}

//#define BOOST_PP_LOCAL_MACRO(n)
template void MetrisRunner::statMesh0<MetricFieldAnalytical>(int tdim, MeshStat *stat);
template void MetrisRunner::statMesh0<MetricFieldFE        >(int tdim, MeshStat *stat);
//#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
//#include BOOST_PP_LOCAL_ITERATE()


}