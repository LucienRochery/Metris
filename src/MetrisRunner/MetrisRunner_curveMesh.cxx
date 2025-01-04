//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "MetrisRunner.hxx"

#include "../LPopt/msh_maxccoef.hxx"
#include "../Mesh/Mesh.hxx"

#include "../aux_histogram.hxx"
#include "../aux_timer.hxx"
#include "../io_libmeshb.hxx"
#include "../quality/msh_metqua.hxx"
#include "../BezierOffsets/msh_curve_offsets.hxx"

namespace Metris{

void MetrisRunner::curveMesh(){
  if(this->metricFE){
    curveMesh0<MetricFieldFE>();
  }else{
    curveMesh0<MetricFieldAnalytical>();
  }
}
template<class MFT>
void MetrisRunner::curveMesh0(){

  Mesh<MFT> &msh = static_cast<Mesh<MFT>&>(*msh_g);

  //int itype = opt.m["smooth"].as<int>();
  int itype = param.curveType;
  int iverb = param.iverb;

  if(iverb >= 1) printf("-- Metric based smoothing requested type = %d \n",itype);


  msh.met.setSpace(MetSpace::Log);


  double qmin, qmax, qavg;
  bool iinva;
  dblAr1 lquae, dum = {0.1, 0.9};
  if(iverb >= 1){
    if(msh.idim == 2){
      getmetquamesh<MFT,2,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    }else{
      getmetquamesh<MFT,3,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    }
    print_histogram(lquae,IntrpTyp::Geometric,dum,"q","Element quality");
  }


  double t1 = get_wall_time();
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){
    if(ideg == msh.curdeg){
      if(itype == 1){
        if(iverb >= 1) printf(" - using curveMeshOffsets\n");
        if constexpr(ideg > 1){
          if(msh.idim == 2){
            curveMeshOffsets<MFT,2,ideg>(msh);
          }else{
            curveMeshOffsets<MFT,3,ideg>(msh);
          }
        }
      }else if(itype == 2){
        if(iverb >= 1) printf(" - using maximizeMetCcoef\n");
        if constexpr(ideg > 1){
          if(msh.idim == 2){
            maximizeMetCcoef<MFT,2,2,ideg>(msh, OptDoF::HO, LPMethod::IPM, LPLib::alglib);
          }else{
            maximizeMetCcoef<MFT,3,3,ideg>(msh, OptDoF::HO, LPMethod::IPM, LPLib::alglib);
          }
        }
      }
    }
  }CT_FOR1(ideg);
  double t2 = get_wall_time();
  if(iverb >= 1) printf("-- Curving end time = %f\n",t2-t1);


  if(iverb >= 2){
    writeMesh("smooth_end.meshb",msh);
    int ideg = msh.curdeg;
    msh.curdeg = 1;
    writeMesh("smooth_end.qua.meshb",msh);
    msh.met.writeMetricFile("smooth_end.solb");
    msh.curdeg = ideg;

    writeField("smooth_end.qua.solb",msh,SolTyp::P0Elt,lquae);
  }

  if(iverb >= 1){
    if(msh.idim == 2){
      getmetquamesh<MFT,2,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    }else{
      getmetquamesh<MFT,3,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    }
    print_histogram(lquae,IntrpTyp::Geometric,dum,"q","Element quality");
  }


}
template void MetrisRunner::curveMesh0<MetricFieldFE>();
template void MetrisRunner::curveMesh0<MetricFieldAnalytical>();


} // end namespace 
