//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "MetrisRunner.hxx"

#include "../adapt/msh_collapse.hxx"
#include "../adapt/msh_swap2D.hxx"
#include "../adapt/msh_insert2D.hxx"

#include "../Mesh/Mesh.hxx"
#include "../aux_utils.hxx"
#include "../aux_timer.hxx"
#include "../quality/msh_metqua.hxx"
#include "../io_libmeshb.hxx"
#include "../adapt/msh_lineadapt.hxx"
#include "../smoothing/msh_smooball.hxx"
#include "../msh_checktopo.hxx"
#include "../aux_histogram.hxx"
#include "../CT_loop.hxx"
#include "../msh_lenedg.hxx"
#include "../linalg/det.hxx"

namespace Metris{


void MetrisRunner::adaptMesh(){
  CT_FOR0_INC(2,3,gdim){if(gdim == msh_g->idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh_g->curdeg == ideg){
      if(this->metricFE){
        adaptMesh0<MetricFieldFE        ,gdim,ideg>();
      }else{
        adaptMesh0<MetricFieldAnalytical,gdim,ideg>();
      }
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
}


// Profiling is attrocious if the template parameters are unrolled within the function
template<class MFT,int gdim,int ideg>
void MetrisRunner::adaptMesh0(){
  Mesh<MFT> &msh = static_cast<Mesh<MFT>&>(*msh_g);

  //printf("## DEBUG calling only reinsertFlat ! \n");
  //writeMesh("reinsertFlat0",msh);
  //reinsertFlat<MFT,gdim,ideg>(msh, true);
  //writeMesh("reinsertFlat1",msh);
  //return;

  // Make it an option 
  const double minstat = 1.0e-12;
  const int miter = param.adp_niter;
  const int iverb = param.iverb;

  if(msh.param->adp_niter == 0)  return;  
  if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(), "Implement tetras in cavity and adaptMesh0");

  msh.cleanup();

  if(param.iverb >= 1) printf("-- START adaptMesh with miter = %d \n",miter);
  if(param.iverb >= 1){
    writeMesh("debug_adapt_inp", msh);
    msh.met.writeMetricFile("debug_adapt_inp");
    if(iverb >= 2) writeBackLinks("debug_adapt_inp_poi2bak", msh);
    double rat = msh.npoin / (double) msh.nface;
    dblAr1 pdens(msh.npoin);
    //dblAr1 fdens(msh.nface);
    //fdens.set_n(msh.nface);
    pdens.fill(0);
    //for(int iface = 0; iface < msh.nface; iface++){
    //  if(isdeadent(iface,msh.fac2poi)) continue;
    //  double meas = getmeasentP1<2,2>(msh.face2poi[iface],
    //                                  msh.coord);
    //  for(int ii = 0; ii < 3; ii++){
    //    int ipoin = msh.fac2poi(iface,ii);
    //    double det = detsym<2>(msh.met[ipoin]);
    //    pdens[ipoin] += sqrt(det) / meas;
    //  }

    //}
    pdens.set_n(msh.npoin);
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      double det = detsym<2>(msh.met[ipoin]);
      pdens[ipoin] = sqrt(det) / sqrt(3) * 2 * rat;
    }
    writeField("debug_adapt_inp.dens.solb",
                              msh,SolTyp::CG,pdens);
  }

  msh.met.setSpace(MetSpace::Log);
  msh.setBasis(FEBasis::Lagrange);

  double qmin, qmax, qavg;
  double qmax_suf;
  bool iinva;

  if(msh.param->dbgfull) check_topo(msh);

  double pct_unit;

  intAr2 ilned;
  ilned.set_n(0);
  dblAr1 rlned;
  dblAr1 lenbds = {1.0/sqrt(2), sqrt(2)};
  if(msh.param->iverb >= 1){
    pct_unit = getLengthEdges<MFT>(msh,ilned,rlned);
    print_histogram(rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
  }
  double t0,t1;

  if(msh.CAD()){

    t0 = get_wall_time();
    adaptGeoLines<MFT>(msh,0);
    t1 = get_wall_time();
    if(iverb >= 1) printf(" - adaptGeoLines time = %fs \n",t1-t0);
    if(iverb >= 2) writeMesh("v2_geolines_adp",msh);
    if(iverb >= 2) msh.met.writeMetricFile("v2_geolines_adp");
    //wait();

    pct_unit = getLengthEdges_Bdry<MFT>(msh,ilned,rlned);
    print_histogram(rlned,IntrpTyp::Linear,lenbds,"l","Edge length (bdry)");

    //wait();
    //pct_unit = getLengthEdges<MFT>(msh,ilned,rlned,LenTyp::Quad);
    //print_histogram(rlned,IntrpTyp::Linear,lenbds,"l","Edge length (bdry, quad)");


    swap2D<MFT,gdim,ideg>(msh, Defaults::swapOptAdapt);

    //if(iverb >= 2) writeMesh("debug_swap_pre.meshb",msh);
    //if(iverb >= 2) msh.met.writeMetricFile("debug_swap_pre.solb");
  }


  if(msh.param->dbgfull) check_topo(msh);
  
  // Will never exceed this 
  const int miter_max = 100;

  dblAr1 lquae(msh.nface);

  double tinsert = 0, tcollapse = 0, tswap = 0, tsmooth = 0;
  double ttotal = get_wall_time();

  int iopt_niter = 0;
  double stat0 = 1;
  //for(int niter = 1; niter <= miter || ( miter < 0 && niter <= miter_max
  //                                    && stat0 > 0.1); niter++){
  for(int niter = 1; niter <= miter || (miter < 0 && niter < miter_max); niter++){
    stat0 = 0;
    //if(niter == miter - 1){
    //  iverb++;
    //  printf("last iter iverb++\n");
    //  wait();
    //}


    qmax_suf = qavg * MAX(10 / (niter * 1.0), 1.0);
    //qmax_suf = 1.0 - (niter - 1) / (double) miter; 
    //qmax_suf = 0;
    double stat;

    // 1. Collapse short edges
    t0 = get_wall_time();
    stat  = collapseShortEdges<MFT,gdim,ideg>(msh, qmax_suf);
    stat0 = MAX(stat0,stat);
    t1 = get_wall_time();
    tcollapse += t1-t0;

    if(iverb >= 1){
      if(iverb >= 2) writeMesh("v2_collapse_adp"+ std::to_string(niter)+".meshb",msh);
      if(iverb >= 2) msh.met.writeMetricFile("v2_collapse_adp"+ std::to_string(niter)+".solb");
      if(iverb >= 2) writeBackLinks("v2_collapse_adp_poi2bak" + std::to_string(niter), msh);
      getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
      printf(" - Quality min = %15.7e \n",qmin);
      printf("           max = %15.7e \n",qmax);
      printf("           avg = %15.7e \n",qavg);
      if(iverb >= 2) writeField("v2_collapsequa_adp"+ std::to_string(niter)+".solb",
                                msh,SolTyp::P0Elt,lquae);
      printf("------------------------------------------------------------\n");
      printf("\n - iteration %d collapse stat = %f time = %f \n\n",niter,stat,t1-t0);
      printf("------------------------------------------------------------\n");
    }
    if(msh.param->dbgfull) check_topo(msh);
    // 2. Swaps

    if(niter%2 == 0 || qmax_suf < 0.5){
      t0 = get_wall_time();
      stat  = swap2D<MFT,gdim,ideg>(msh, Defaults::swapOptAdapt);
      stat0 = MAX(stat0,stat);
      t1 = get_wall_time();
      tswap += t1-t0;

      if(iverb >= 1){
        if(iverb >= 2) writeMesh("v2_swap_adp"+ std::to_string(niter)+".meshb",msh);
        if(iverb >= 2) msh.met.writeMetricFile("v2_swap_adp"+ std::to_string(niter)+".solb");
        if(iverb >= 2) writeBackLinks("v2_swap_adp_poi2bak" + std::to_string(niter), msh);
        getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
        printf(" - Quality min = %15.7e \n",qmin);
        printf("           max = %15.7e \n",qmax);
        printf("           avg = %15.7e \n",qavg);
        if(iverb >= 2) writeField("v2_swapqua_adp"+ std::to_string(niter)+".solb",
                                 msh,SolTyp::P0Elt,lquae);
        printf("------------------------------------------------------------\n");
        printf("\n - iteration %d swaps stat = %f time = %f \n\n",niter,stat,t1-t0);
        printf("------------------------------------------------------------\n");
      }
    }
    if(msh.param->dbgfull) check_topo(msh);


    // 3. Insert on long edges 

    t0 = get_wall_time();
    stat  = insertLongEdges<MFT,gdim,ideg>(msh);
    stat0 = MAX(stat0,stat);
    t1 = get_wall_time();
    tinsert += t1-t0;

    if(msh.param->dbgfull) check_topo(msh);
    msh.cleanup();
    if(msh.param->dbgfull) check_topo(msh);

    if(iverb >= 1){
      if(iverb >= 2) writeMesh("v2_insert_adp"+ std::to_string(niter)+".meshb",msh);
      if(iverb >= 2) msh.met.writeMetricFile("v2_insert_adp"+std::to_string(niter)+".solb");
      if(iverb >= 2) writeBackLinks("v2_insert_adp_poi2bak" + std::to_string(niter), msh);
      //getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
      //printf(" - Quality min = %15.7e \n",qmin);
      //printf("           max = %15.7e \n",qmax);
      //printf("           avg = %15.7e \n",qavg);
      //if(iverb >= 2) writeField("debug_insert_glo_qua"+ std::to_string(niter)+".solb",
      //                          msh.get_tdim(),SolTyp::P0Elt,lquae);
      printf("------------------------------------------------------------\n");
      printf("\n - iteration %d insertions stat = %f time = %f \n\n",niter,stat,t1-t0);
      printf("------------------------------------------------------------\n");
    }


    //if(msh.param->opt_unif){
    //  if(iverb >= 1) t0 = get_wall_time();
    //  stat = rebalanceMesh<MFT,gdim>(msh);
    //  stat0 = MAX(stat0,stat);
    //  if(iverb >= 1) t1 = get_wall_time();
    //  if(iverb >= 1){
    //    if(iverb >= 2) writeMesh("debug_rebalance_glo"+ std::to_string(niter)+".meshb",msh);
    //    if(iverb >= 2) msh.met.writeMetricFile("debug_rebalance_glo"+ std::to_string(niter)+".solb");    
    //    printf("------------------------------------------------------------\n");
    //    printf("\n - iteration %d rebalance stat = %f time = %f \n\n",niter,stat,t1-t0);
    //    printf("------------------------------------------------------------\n");
    //  }
    //}


    if(msh.param->opt_unif){
      // 4. Smoothing (heuristic) -> fast but bad; improve
      if(iverb >= 1) t0 = get_wall_time();
      double stat = smoothInterior_Ball<MFT>(msh,QuaFun::Unit);
      stat0 = MAX(stat, stat0); 
      if(iverb >= 1) t1 = get_wall_time();
      if(iverb >= 1){
        if(iverb >= 2) writeMesh("v2_unif_adp"+ std::to_string(niter)+".meshb",msh);
        if(iverb >= 2) msh.met.writeMetricFile("v2_unif_adp"+ std::to_string(niter)+".solb");    
        if(iverb >= 2) writeBackLinks("v2_unif_adp_poi2bak" + std::to_string(niter), msh);
        printf("------------------------------------------------------------\n");
        printf("\n - iteration %d unif ball stat = %f time = %f \n\n",niter,stat,t1-t0);
        printf("------------------------------------------------------------\n");
        getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,NULL);
        printf(" - Quality min = %15.7e \n",qmin);
        printf("           max = %15.7e \n",qmax);
        printf("           avg = %15.7e \n",qavg);
      }
    }



    
    //#ifdef USE_PETSC
    //// 4. Smoothing 
    //t0 = get_wall_time();
    ////smoothInterior_Ball<MFT>(msh,iverb-1,0);
    //PetscErrorCode ierro = smoothInterior_Full_custom<MFT>(msh,0);
    //if(ierro != PETSC_SUCCESS){
    //  printf("PetscErrorCode %d \n",ierro);
    //}
    //#endif

    //t1 = get_wall_time();
    //tsmooth += t1-t0;
    //if(iverb >= 2) writeMesh("debug_smooth_glo"+ std::to_string(niter)+".meshb",msh);
    //if(iverb >= 2) msh.met.writeMetricFile("debug_smooth_glo"+ std::to_string(niter)+".solb");
    //printf("------------------------------------------------------------\n");
    //printf("\n - iteration %d smooth stat = %f time = %f \n\n",niter,stat,t1-t0);
    //printf("------------------------------------------------------------\n");
    //#ifndef NDEBUG
    //check_topo(msh);
    //#endif


    pct_unit = getLengthEdges(msh,ilned,rlned);
    if(iverb >= 1) 
      printf("-- END ADAPT LOOP %d / %d  %f edges unit, op stat = %f \n",niter,miter,
            pct_unit,stat0);
    if(pct_unit >= 0.999){
      if(iverb >= 1) printf("------------------------------------------------------------\n");
      if(iverb >= 1) printf(" - 99.9%% edges unit exit\n");
      break;
    }
    if(stat0 < 1.0e-3 && miter < 0 || stat0 < minstat){
      if(iverb >= 1) printf("------------------------------------------------------------\n");
      if(iverb >= 1) printf(" - low stat = %e break or optimize\n",stat0);
      if(niter >= miter -1) break;
      if(msh.param->opt_niter > 0 && 
        (msh.param->adp_opt_niter < iopt_niter || msh.param->adp_opt_niter < 0)
         && !msh.param->opt_unif){
        iopt_niter++;
        stat = optimMesh();
        if(iverb >= 2){
          writeMesh("v2_optim_adp" + std::to_string(iopt_niter), msh);
          msh.met.writeMetricFile("v2_optim_adp" + std::to_string(iopt_niter));
          if(iverb >= 2) writeBackLinks("v2_optim_adp_poi2bak" + std::to_string(niter), msh);
        }
        if(stat < minstat){
          if(iverb >= 1) printf(" - low optim stat %e break\n",stat);
          break;
        }
      }else{
        break;
      }
    }
  }


  ttotal = get_wall_time() - ttotal;

  msh.cleanup();

  if(iverb >= 1){
    printf("-- Adaptation end total time = %f \n",ttotal);
    printf(" - insertion time = %f \n",tinsert);
    printf(" -  collapse time = %f \n",tcollapse);
    printf(" -      swap time = %f \n",tswap);
    printf(" - smoothing time = %f \n",tsmooth);
  }

  if(iverb >= 1 && msh.param->opt_niter > 0 || iverb >= 2){
    writeMesh("adapt_end.meshb",msh);
    msh.met.writeMetricFile("adapt_end.solb");
    if(iverb >= 2){
      getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
      writeField("adapt_end.qua.solb",msh,SolTyp::P0Elt,lquae);
    }
  }


  if(msh.param->iverb >= 1) statMesh();

  //pct_unit = getLengthEdges(msh,ilned,rlned);
  //print_histogram(rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
  //writeEdgesLengths(msh,"edgelen", ilned, rlned);
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void MetrisRunner::adaptMesh0<MetricFieldFE        ,2,n>();\
template void MetrisRunner::adaptMesh0<MetricFieldFE        ,3,n>();\
template void MetrisRunner::adaptMesh0<MetricFieldAnalytical,2,n>();\
template void MetrisRunner::adaptMesh0<MetricFieldAnalytical,3,n>();
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

}//end namespace
