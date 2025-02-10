//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../MetrisRunner/MetrisRunner.hxx"
#include "../Mesh/Mesh.hxx"
#include "../adapt/msh_swap2D.hxx"
#include "../adapt/msh_reinsert_flat.hxx"
#include "../aux_utils.hxx"
#include "../aux_timer.hxx"
#include "../quality/msh_metqua.hxx"
#include "../io_libmeshb.hxx"
#include "../smoothing/msh_smooball.hxx"
#include "../msh_checktopo.hxx"
#include "../aux_histogram.hxx"
#include "../msh_lenedg.hxx"
#include "../mprintf.hxx"

namespace Metris{


double MetrisRunner::optimMesh(){
  double stat = 0;
  CT_FOR0_INC(2,3,gdim){if(gdim == msh_g->idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh_g->curdeg == ideg){
      if(this->metricFE){
        stat = optimMesh0<MetricFieldFE        ,gdim,ideg>();
      }else{
        stat = optimMesh0<MetricFieldAnalytical,gdim,ideg>();
      }
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
  return stat;
}



template<class MFT, int gdim, int ideg>
double MetrisRunner::optimMesh0(){
  GETVDEPTH((*this));
  Mesh<MFT> &msh = static_cast<Mesh<MFT>&>(*msh_g);
  if(param_.opt_niter == 0) return 0;

  double t01 = get_wall_time();

  if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(), 
                                   "Implement tetras in cavity and optimMesh0");


  CPRINTF2("-- START optimMesh\n");

  // the higher this is, the more permissive 
  //const int qpower = -1;
  //const int qpnorm = 2;

  msh.met.setSpace(MetSpace::Log);
  msh.setBasis(FEBasis::Lagrange);

  double qmin, qmax, qavg;
  //double qmax_suf;
  bool iinva;

  //#ifndef NDEBUG
  //check_topo(msh);
  //#endif

  const int miter = param_.opt_niter;
  swapOptions swapOpt(*(msh.param));

  if(msh.param->dbgfull) check_topo(msh);

  intAr2 ilned;
  dblAr1 rlned;
  dblAr1 lenbds = {1.0/sqrt(2), sqrt(2)};
  dblAr1 lquae,dum = {0.1, 0.9};
  double t0,t1;

  if(DOPRINTS1()){
    getLengthEdges(msh,ilned,rlned);
    print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
  
    getmetquamesh<MFT,gdim,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    print_histogram(msh,lquae,IntrpTyp::Geometric,dum,"q","Element quality (As Pk)");

    getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,NULL);
    CPRINTF1(" - Quality min = %15.7e \n",qmin);
    CPRINTF1("           max = %15.7e \n",qmax);
    CPRINTF1("           avg = %15.7e \n",qavg);
  }

  double stat0 = 1;
  for(int niter = 1; niter <= miter && stat0 >= 0.01; niter++){
    INCVDEPTH(msh);
    stat0 = 0;

    //qmax_suf = qavg * MAX(10 / (niter * 1.0), 1.0);
    ////qmax_suf = 1.0 - (niter - 1) / (double) miter; 

    double stat;


    //if(msh.param->opt_unif){
    //  t0 = get_wall_time();
    //  double stat = rebalanceMesh<MFT,gdim>(msh);
    //  stat0 = MAX(stat, stat0); 
    //  t1 = get_wall_time();
    //  if(DOPRINTS1()){
    //    if(DOPRINTS2()) writeMesh("debug_rebalance_glo"+ std::to_string(niter)+".meshb",msh);
    //    if(DOPRINTS2()) msh.met.writeMetricFile("debug_rebalance_glo"+ std::to_string(niter)+".solb");    
    //    CPRINTF1("------------------------------------------------------------\n");
    //    CPRINTF1("- iteration %d rebalance stat = %f time = %f \n",niter,stat,t1-t0);
    //    CPRINTF1("------------------------------------------------------------\n");
    //  }
    //}




    //// 4. Smoothing (heuristic) -> fast but bad; improve
    //t0 = get_wall_time();
    //stat = smoothInterior_Ball<MFT>(msh,QuaFun::Unit);
    //stat0 = MAX(stat, stat0); 
    //t1 = get_wall_time();
    //if(DOPRINTS1()){
    //  if(DOPRINTS2()) writeMesh("debug_unit_glo"+ std::to_string(niter)+".meshb",msh);
    //  if(DOPRINTS2()) msh.met.writeMetricFile("debug_unit_glo"+ std::to_string(niter)+".solb");    
    //  CPRINTF1("------------------------------------------------------------\n");
    //  CPRINTF1("- iteration %d smooth unit stat = %f time = %f \n",niter,stat,t1-t0);
    //  CPRINTF1("------------------------------------------------------------\n");
    //  getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,NULL);
    //  CPRINTF1(" - Quality min = %15.7e \n",qmin);
    //  CPRINTF1("           max = %15.7e \n",qmax);
    //  CPRINTF1("           avg = %15.7e \n",qavg);
    //}

    //if(msh.param->dbgfull) check_topo(msh);


    // 4. Smoothing (iterative)

    if(msh.idim == msh.get_tdim()){

      t0 = get_wall_time();
      stat = smoothInterior_Ball<MFT>(msh,QuaFun::Distortion);
      stat0 = MAX(stat, stat0); 
      t1 = get_wall_time();
      if(DOPRINTS1()){
        if(DOPRINTS2()) writeMesh("v2_smooth_opt"+ std::to_string(niter)+".meshb",msh);
        if(DOPRINTS2()) msh.met.writeMetricFile("v2_smooth_opt"+ std::to_string(niter)+".solb");    
        CPRINTF1("------------------------------------------------------------\n");
        CPRINTF1("- iteration %d smooth ball stat = %f time = %f \n",niter,stat,t1-t0);
        CPRINTF1("------------------------------------------------------------\n");
        getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,NULL);
        CPRINTF1(" - Quality min = %15.7e \n",qmin);
        CPRINTF1("           max = %15.7e \n",qmax);
        CPRINTF1("           avg = %15.7e \n",qavg);
      }
      if(msh.param->dbgfull) check_topo(msh);

    }else{
      CPRINTF1("## Smoothing disabled in case gdim = %d tdim = %d \n", msh.idim, msh.get_tdim());
    }


    #if 0
    //getmetquamesh<MFT,gdim,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    //print_histogram(msh,lquae,IntrpTyp::Geometric,dum,"q","Element quality (As Pk)");
    // 5. Smoothing  2
    t0 = get_wall_time();
    //smoothInterior_Full<MFT>(msh);
    //smoothInterior_Full_TAO<MFT>(msh); 
    smoothInterior_Full_custom<MFT>(msh,&stat);
    stat0 = MAX(stat0,stat);
    t1 = get_wall_time();
    if(DOPRINTS1()){
      if(DOPRINTS2()) writeMesh("debug_smooth_glo"+ std::to_string(niter)+".meshb",msh);
      if(DOPRINTS2()) msh.met.writeMetricFile("debug_smooth_glo"+ std::to_string(niter)+".solb");    
      CPRINTF1("------------------------------------------------------------\n");
      CPRINTF1("- iteration %d smooth full stat = %15.7e time = %f \n",niter,stat,t1-t0);
      CPRINTF1("------------------------------------------------------------\n");
      getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,NULL);
      CPRINTF1(" - Quality min = %15.7e \n",qmin);
      CPRINTF1("           max = %15.7e \n",qmax);
      CPRINTF1("           avg = %15.7e \n",qavg);
    }
    #endif

    if(msh.idim == msh.get_tdim()){
      t0 = get_wall_time();
      int noper = reinsertFlat<MFT,gdim,ideg>(msh);
      t1 = get_wall_time();
      msh.cleanup();
      stat  = noper / (double) msh.nface; 
      stat0 = MAX(stat0, stat);
      if(DOPRINTS1()){
        if(DOPRINTS2() && noper >= 0) writeMesh("v2_flat_opt"+ std::to_string(niter)+".meshb",msh);
        if(DOPRINTS2() && noper >= 0) msh.met.writeMetricFile("v2_flat_opt"+ std::to_string(niter)+".solb");
        CPRINTF1("------------------------------------------------------------\n");
        CPRINTF1("- iteration %d flat collapse noper = %d stat = %f time = %f \n",niter,noper,stat,t1-t0);
        CPRINTF1("------------------------------------------------------------\n");
      }

      if(msh.param->dbgfull) check_topo(msh);

    }else{
      CPRINTF1("## reinsertFlat disabled in case gdim = %d tdim = %d \n", msh.idim, msh.get_tdim());
    }
    //getmetquamesh<MFT,gdim,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    //print_histogram(msh,lquae,IntrpTyp::Geometric,dum,"q","Element quality (As Pk)");
    
    // ALWAYS SWAP LAST ! 
    // 2. Swaps
    t0 = get_wall_time();
    int nswap;
    stat  = swap2D<MFT,gdim,ideg>(msh, swapOpt, &nswap);
    stat0 = MAX(stat0,stat);
    t1 = get_wall_time();
    if(msh.param->dbgfull) check_topo(msh);
    msh.cleanup();
    if(msh.param->dbgfull) check_topo(msh);
    if(DOPRINTS1()){
      if(DOPRINTS2()) writeMesh("v2_swap_opt"+ std::to_string(niter)+".meshb",msh);
      if(DOPRINTS2()) msh.met.writeMetricFile("v2_swap_opt"+ std::to_string(niter)+".solb");
      CPRINTF1("------------------------------------------------------------\n");
      CPRINTF1("- iteration %d swaps stat = %f time = %f \n",niter,stat,t1-t0);
      CPRINTF1("------------------------------------------------------------\n");
      getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,NULL);
      CPRINTF1(" - Quality min = %15.7e \n",qmin);
      CPRINTF1("           max = %15.7e \n",qmax);
      CPRINTF1("           avg = %15.7e \n",qavg);
    }

    CPRINTF1("-- END OPTIM LOOP %d / %d, op stat = %15.7e \n",niter,miter,stat0);
    if(stat0 < 1.0e-8){
      CPRINTF1("------------------------------------------------------------\n");
      CPRINTF1(" - low stat = %e break \n",stat0);
      break;
    }

  }

  //if(param_.iverb >= 1){
  //  getmetquamesh<MFT,gdim,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
  //  print_histogram(msh,lquae,IntrpTyp::Geometric,dum,"q","Element quality (As Pk)");

  //  getLengthEdges(msh,ilned,rlned);
  //  print_histogram(rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
  //}

  double t11 = get_wall_time();

  CPRINTF1("-- OptimMesh end runtime = %fs \n",t11-t01);
  return stat0;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template double MetrisRunner::optimMesh0<MetricFieldFE        ,2,n>();\
template double MetrisRunner::optimMesh0<MetricFieldFE        ,3,n>();\
template double MetrisRunner::optimMesh0<MetricFieldAnalytical,2,n>();\
template double MetrisRunner::optimMesh0<MetricFieldAnalytical,3,n>();
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

}//end namespace
