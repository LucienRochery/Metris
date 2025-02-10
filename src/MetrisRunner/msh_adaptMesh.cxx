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
#include "../mprintf.hxx"


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
  GETVDEPTH((*this));
  Mesh<MFT> &msh = static_cast<Mesh<MFT>&>(*msh_g);


  // Make it an option 
  const double minstat = 1.0e-12;
  const int miter = param_.adp_niter;

  if(msh.param->adp_niter == 0)  return;  
  if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(), "Implement tetras in cavity and adaptMesh0");

  msh.cleanup();

  CPRINTF1("-- START adaptMesh with miter = %d \n",miter);
  if(DOPRINTS1()){
    writeMesh("debug_adapt_inp", msh);
    msh.met.writeMetricFile("debug_adapt_inp");
    if(DOPRINTS2()) writeBackLinks("debug_adapt_inp_poi2bak", msh);
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

  if(msh.param->dbgfull) check_topo(msh,1);

  double pct_unit;
  int nswap, ninser, ncoll;

  intAr2 ilned;
  ilned.set_n(0);
  dblAr1 rlned;
  dblAr1 lenbds = {1.0/sqrt(2), sqrt(2)};
  if(DOPRINTS1()){
    pct_unit = getLengthEdges<MFT>(msh,ilned,rlned);
    print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
  }
  double t0,t1;

  // This is the common thread for all routines. Tagged elements are ignored
  const int ithrdfro = 0;
  msh.tag[ithrdfro]++;

  if(msh.CAD() && msh.param->adp_line_adapt){

    t0 = get_wall_time();
    adaptGeoLines<MFT>(msh,0);
    t1 = get_wall_time();
    CPRINTF1(" - adaptGeoLines time = %fs \n",t1-t0);
    if(DOPRINTS2()) writeMesh("v2_geolines_adp",msh);
    if(DOPRINTS2()) msh.met.writeMetricFile("v2_geolines_adp");

    if(DOPRINTS1()){
      pct_unit = getLengthEdges_Bdry<MFT>(msh,ilned,rlned);
      print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length (bdry)");
    }

    //#ifndef NDEBUG
    //  printf("wait after geolines\n");
    //  wait();
    //#endif

    //wait();
    //pct_unit = getLengthEdges<MFT>(msh,ilned,rlned,LenTyp::Quad);
    //print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length (bdry, quad)");


    swap2D<MFT,gdim,ideg>(msh, Defaults::swapOptAdapt, &nswap, ithrdfro, 1);

    //if(DOPRINTS2()) writeMesh("debug_swap_pre.meshb",msh);
    //if(DOPRINTS2()) msh.met.writeMetricFile("debug_swap_pre.solb");
  }


  if(msh.param->dbgfull) check_topo(msh,1);
  

  // Will never exceed this 
  const int miter_max = 100;

  dblAr1 lquae(msh.nface);

  double tinsert = 0, tcollapse = 0, tswap = 0, tsmooth = 0;
  double ttotal = get_wall_time();


  int iopt_niter = 0;
  double stat0 = 1;
  //for(int niter = 1; niter <= miter || ( miter < 0 && niter <= miter_max
  //                                    && stat0 > 0.1); niter++){
  //double stat_prev = stat0;
  msh.tag[ithrdfro]++;
  for(int niter = 1; niter <= miter || (miter < 0 && niter < miter_max); niter++){
    stat0 = 0;
    //if(niter == miter - 1){
    //  iverb++;
    //  printf("last iter iverb++\n");
    //  wait();
    //}

    double tloop0 = get_wall_time();


    qmax_suf = qavg * MAX(10 / (niter * 1.0), 1.0);
    //qmax_suf = 1.0 - (niter - 1) / (double) miter; 
    //qmax_suf = 0;
    double stat;

    // 1. Collapse short edges
    t0 = get_wall_time();
    stat  = collapseShortEdges<MFT,gdim,ideg>(msh, qmax_suf, &ncoll, ithrdfro, 1, 2);
    stat0 = MAX(stat0,stat);
    t1 = get_wall_time();
    tcollapse += t1-t0;

    if(DOPRINTS2()){
      writeMesh("v2_collapse_adp"+ std::to_string(niter)+".meshb",msh);
      msh.met.writeMetricFile("v2_collapse_adp"+ std::to_string(niter)+".solb");
      writeBackLinks("v2_collapse_adp_poi2bak" + std::to_string(niter), msh);
    }

    if(DOPRINTS2()){
      getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
      CPRINTF1(" - Quality min = %15.7e \n",qmin);
      CPRINTF1(" max = %15.7e \n",qmax);
      CPRINTF1(" avg = %15.7e \n",qavg);
      if(DOPRINTS2()) writeField("v2_collapsequa_adp"+ std::to_string(niter)+".solb",
                                msh,SolTyp::P0Elt,lquae);
      CPRINTF2("------------------------------------------------------------\n");
      CPRINTF2("- iteration %d collapse stat = %f time = %f \n",niter,stat,t1-t0);
      CPRINTF2("------------------------------------------------------------\n");
    }
    if(msh.param->dbgfull) check_topo(msh,1);
    // 2. Swaps

    if(niter%2 == 0 || qmax_suf < 0.5){
      t0 = get_wall_time();
      stat  = swap2D<MFT,gdim,ideg>(msh, Defaults::swapOptAdapt, &nswap, ithrdfro, 1);
      stat0 = MAX(stat0,stat);
      t1 = get_wall_time();
      tswap += t1-t0;

      if(DOPRINTS2()){
        writeMesh("v2_swap_adp"+ std::to_string(niter)+".meshb",msh);
        msh.met.writeMetricFile("v2_swap_adp"+ std::to_string(niter)+".solb");
        writeBackLinks("v2_swap_adp_poi2bak" + std::to_string(niter), msh);
      }
      if(DOPRINTS2()){
        getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
        CPRINTF2(" - Quality min = %15.7e \n",qmin);
        CPRINTF2(" max = %15.7e \n",qmax);
        CPRINTF2(" avg = %15.7e \n",qavg);
        if(DOPRINTS2()) writeField("v2_swapqua_adp"+ std::to_string(niter)+".solb",
                                 msh,SolTyp::P0Elt,lquae);
        CPRINTF2("------------------------------------------------------------\n");
        CPRINTF2("- iteration %d swaps stat = %f time = %f \n",niter,stat,t1-t0);
        CPRINTF2("------------------------------------------------------------\n");
      }
    }
    if(msh.param->dbgfull) check_topo(msh,1);

    // 3. Insert on long edges 

    t0 = get_wall_time();
    stat  = insertLongEdges<MFT,gdim,ideg>(msh, &ninser, ithrdfro, 1, 2);
    stat0 = MAX(stat0,stat);
    t1 = get_wall_time();
    tinsert += t1-t0;

    if(msh.param->dbgfull) check_topo(msh,1);
    msh.cleanup();
    if(msh.param->dbgfull) check_topo(msh,1);

    if(DOPRINTS2()){
      writeMesh("v2_insert_adp"+ std::to_string(niter)+".meshb",msh);
      msh.met.writeMetricFile("v2_insert_adp"+std::to_string(niter)+".solb");
      writeBackLinks("v2_insert_adp_poi2bak" + std::to_string(niter), msh);
    }
    if(DOPRINTS2()){
      //getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
      //printf(" - Quality min = %15.7e \n",qmin);
      //printf("           max = %15.7e \n",qmax);
      //printf("           avg = %15.7e \n",qavg);
      //if(DOPRINTS2()) writeField("debug_insert_glo_qua"+ std::to_string(niter)+".solb",
      //                          msh.get_tdim(),SolTyp::P0Elt,lquae);
      CPRINTF2("------------------------------------------------------------\n");
      CPRINTF2("- iteration %d insertions stat = %f time = %f \n",niter,stat,t1-t0);
      CPRINTF2("------------------------------------------------------------\n");
    }


    //if(msh.param->opt_unif){
    //  if(iverb >= 1) t0 = get_wall_time();
    //  stat = rebalanceMesh<MFT,gdim>(msh);
    //  stat0 = MAX(stat0,stat);
    //  if(iverb >= 1) t1 = get_wall_time();
    //  if(iverb >= 1){
    //    if(DOPRINTS2()) writeMesh("debug_rebalance_glo"+ std::to_string(niter)+".meshb",msh);
    //    if(DOPRINTS2()) msh.met.writeMetricFile("debug_rebalance_glo"+ std::to_string(niter)+".solb");    
    //    printf("------------------------------------------------------------\n");
    //    printf("\n - iteration %d rebalance stat = %f time = %f \n",niter,stat,t1-t0);
    //    printf("------------------------------------------------------------\n");
    //  }
    //}


    if(msh.param->opt_unif){
      // 4. Smoothing (heuristic) -> fast but bad; improve
      t0 = get_wall_time();
      double stat = smoothInterior_Ball<MFT>(msh,QuaFun::Unit);
      stat0 = MAX(stat, stat0); 
      t1 = get_wall_time();
      if(DOPRINTS2()) writeMesh("v2_unif_adp"+ std::to_string(niter)+".meshb",msh);
      if(DOPRINTS2()) msh.met.writeMetricFile("v2_unif_adp"+ std::to_string(niter)+".solb");    
      if(DOPRINTS2()) writeBackLinks("v2_unif_adp_poi2bak" + std::to_string(niter), msh);
      if(DOPRINTS2()){
        CPRINTF2("------------------------------------------------------------\n");
        CPRINTF2("- iteration %d unif ball stat = %f time = %f \n",niter,stat,t1-t0);
        CPRINTF2("------------------------------------------------------------\n");
        getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,NULL);
        CPRINTF2(" - Quality min = %15.7e \n",qmin);
        CPRINTF2(" max = %15.7e \n",qmax);
        CPRINTF2(" avg = %15.7e \n",qavg);
      }
    }




    pct_unit = 100*getLengthEdges(msh,ilned,rlned);
    int ndigit = ceil(log10((double)msh.npoin
                        + (double)msh.nelem
                        + (double)msh.nface
                        + (double)msh.nedge) ) + 1;

    double tloop1 = get_wall_time();

    //int ncheck = 0;
    //for(int iface = 0; iface < msh.nface; iface++){
    //  if(isdeadent(iface,msh.fac2poi)) continue;
    //  if(msh.fac2tag(ithrdfro,iface) >= msh.tag[ithrdfro]) continue;
    //  ncheck++;
    //}

    //double tloop2 = get_wall_time();
    //printf("## DEBUG tloop empty %f ncehck %d / %d \n",tloop2-tloop1,ncheck,msh.nface);
    //wait();


    std::string fmt = 
    "%s-- Adp loop %3d / %3d time %fs " 
    "%" + std::to_string(ndigit) + "d inser "
    "%" + std::to_string(ndigit) + "d coll "
    "%" + std::to_string(ndigit) + "d swap, "
    + "%7.2fpct unit, op stat = %f \n";
    // "-- Adp loop %3d / %3d inser %d coll %d swap %d, %fpct unit, op stat = %f \n"
    if(DOPRINTS1())printf(fmt.c_str(), spaces_string__, 
             niter,miter, tloop1 - tloop0, ninser,ncoll,nswap, pct_unit,stat0);
    if(pct_unit >= 99.9){
      CPRINTF1("------------------------------------------------------------\n");
      CPRINTF1("- 99.9%% edges unit exit\n");
      break;
    }

    bool stagn = stat0 < 1.0e-3
              || stat0 < minstat;
             // || abs(stat0 - stat_prev) < 1.0e-6 ;// This last criterion just catches cycles

    //stat_prev = stat0;

    if(stagn){
      CPRINTF1(" - low stat = %e break or optimize\n",stat0);
      if(niter >= miter -1) break;
      if(msh.param->opt_niter > 0 && 
        (iopt_niter < msh.param->adp_opt_niter|| msh.param->adp_opt_niter < 0)
         && !msh.param->opt_unif){ 
        iopt_niter++;
        double tsmo0 = get_wall_time();
        stat = optimMesh();
        double tsmo1 = get_wall_time();
        tsmooth += tsmo1 - tsmo0;
        msh.tag[ithrdfro]++;
        if(DOPRINTS2()){
          writeMesh("v2_optim_adp" + std::to_string(iopt_niter), msh);
          msh.met.writeMetricFile("v2_optim_adp" + std::to_string(iopt_niter));
          writeBackLinks("v2_optim_adp_poi2bak" + std::to_string(niter), msh);
        }
        if(stat < minstat){
          CPRINTF1(" - low optim stat %e break\n",stat);
          break;
        }
      }else{
        break;
      }
    }
  }


  ttotal = get_wall_time() - ttotal;

  msh.cleanup();

  CPRINTF1("-- Adaptation end total time = %f \n",ttotal);
  CPRINTF1(" - insertion time = %f \n",tinsert);
  CPRINTF1(" -  collapse time = %f \n",tcollapse);
  CPRINTF1(" -      swap time = %f \n",tswap);
  CPRINTF1(" - smoothing time = %f \n",tsmooth);

  if(DOPRINTS1() && DOPRINTS2()){
    writeMesh("adapt_end.meshb",msh);
    msh.met.writeMetricFile("adapt_end.solb");
    if(DOPRINTS2()){
      getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
      writeField("adapt_end.qua.solb",msh,SolTyp::P0Elt,lquae);
    }
  }


  if(DOPRINTS1()) statMesh();

  //pct_unit = getLengthEdges(msh,ilned,rlned);
  //print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
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
