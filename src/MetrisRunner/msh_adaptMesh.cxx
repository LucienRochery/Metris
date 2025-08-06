//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "MetrisRunner.hxx"

#include "../adapt/msh_collapse.hxx"
#include "../adapt/msh_swap.hxx"
#include "../adapt/msh_insert.hxx"
#include "../adapt/msh_reinsert_flat.hxx"

#include "../Mesh/Mesh.hxx"
#include "../utils/aux_misc.hxx"
#include "../utils/aux_timer.hxx"
#include "../quality/msh_metqua.hxx"
#include "../io_libmeshb.hxx"
#include "../adapt/msh_lineadapt.hxx"
#include "../smoothing/msh_smooball.hxx"
#include "../msh_checktopo.hxx"
#include "../aux_histogram.hxx"
#include "../utils/CT_loop.hxx"
#include "../msh_lenedg.hxx"
#include "../linalg/det.hxx"
#include "../utils/mprintf.hxx"

#include "../low_geo/misc.hxx"


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
  if(param->adp_niter == 0) return;

  GETVDEPTH(this->param);
  Mesh<MFT> &msh = static_cast<Mesh<MFT>&>(*msh_g);


  //METRIS_THROW(TODOExcept());

  // Make it an option 
  const double minstat = 1.0e-12;
  const int miter = param_.adp_niter;

  msh.cleanup();

  CPRINTF1("-- START adaptMesh with miter = {} \n",miter);
  if(DOPRINTS1()){
    writeMesh("debug_adapt_inp", msh);
    msh.met.writeMetricFile("debug_adapt_inp");
  }

  msh.met.setSpace(MetSpace::Exp);
  msh.setBasis(FEBasis::Lagrange);

  double qmin, qmax, qavg;
  double qmax_suf;
  bool iinva;

  #ifndef NDEBUG
  check_topo(msh,1);
  #endif

  int nswap, ninser, ncoll;

  intAr2 ilned;
  ilned.set_n(0);
  dblAr1 rlned;
  dblAr1 lenbds = {1.0/sqrt(2), sqrt(2)};
  lenStat lenstat;
  if(DOPRINTS1()){
    getLengthEdges<MFT>(msh,msh.get_tdim(),-1,ilned,rlned,lenstat);
    print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
  }
  double t0,t1;

  // This is the common thread for all routines. Tagged elements are ignored
  const int ithrdfro = 0;
  const int ithrd1 = 1;
  const int ithrd2 = 2;
  const int ithrd3 = 3;
  msh.tag[ithrdfro]++;

  if(msh.CAD() && msh.param->adp_line_adapt){

    t0 = get_wall_time();
    adaptGeoLines<MFT>(msh);
    t1 = get_wall_time();
    CPRINTF1(" - adaptGeoLines time = {:.2e}ss \n",t1-t0);
    if(DOPRINTS2()) writeMesh("v2_geolines_adp",msh);
    if(DOPRINTS2()) msh.met.writeMetricFile("v2_geolines_adp");

    if(DOPRINTS1()){
      getLengthEdges<MFT>(msh,1,-1,ilned,rlned,lenstat);
      print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length (lines)");
    }
    
    swapMesh<MFT,gdim,ideg>(msh, Defaults::swapOptAdapt, &nswap, ithrdfro, ithrd1, ithrd2);

    #ifndef NDEBUG
    check_topo(msh,1);
    #endif
  }



  

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
    double tloop0 = get_wall_time();


    qmax_suf = qavg * MAX(10 / (niter * 1.0), 1.0);
    //qmax_suf = 1.0 - (niter - 1) / (double) miter; 
    //qmax_suf = 0;
    double stat;
    // 1. Collapse short edges
    t0 = get_wall_time();
    stat  = collapseShortEdges<MFT,gdim,ideg>(msh, qmax_suf, &ncoll, ithrdfro, ithrd1, ithrd2, ithrd3);
    stat0 = MAX(stat0,stat);
    t1 = get_wall_time();
    tcollapse += t1-t0;

    if(DOPRINTS2()){
      writeMesh("v2_collapse_adp"+ std::to_string(niter)+".meshb",msh);
      msh.met.writeMetricFile("v2_collapse_adp"+ std::to_string(niter)+".solb");
      writeBackLinks("v2_collapse_adp_poi2bak" + std::to_string(niter), msh);
    }

    if(DOPRINTS2()){
      getLengthEdges(msh,msh.get_tdim(),-1,ilned,rlned,lenstat);
      CPRINTF1(" - Length qua short = {}\n",lenstat.qua_short);
      CPRINTF1(" -            long  = {}\n",lenstat.qua_long);
      if(DOPRINTS3()) print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");

      getmetquamesh<MFT>(msh,msh.get_tdim(),AsDeg::P1,AsDeg::P1,
                         &iinva,&qmin,&qmax,&qavg,&lquae);
      CPRINTF1(" - Quality min = {:15.7e} \n",qmin);
      CPRINTF1("           max = {:15.7e} \n",qmax);
      CPRINTF1("           avg = {:15.7e} \n",qavg);
      writeField("v2_collapsequa_adp"+ std::to_string(niter)+".solb",
                                msh,SolTyp::P0Elt,lquae);
      CPRINTF2("------------------------------------------------------------\n");
      CPRINTF2("- iteration {} collapse stat = {:.2e} time = {:.2e}s \n",niter,stat,t1-t0);
      CPRINTF2("------------------------------------------------------------\n");
    }
    #ifndef NDEBUG
    check_topo(msh,1);
    #endif

    // 2. Swaps

    if(niter%2 == 0 || qmax_suf < 0.5){
      t0 = get_wall_time();
      stat  = swapMesh<MFT,gdim,ideg>(msh, Defaults::swapOptAdapt, &nswap, ithrdfro, ithrd1, ithrd2);
      stat0 = MAX(stat0,stat);
      t1 = get_wall_time();
      tswap += t1-t0;

      if(DOPRINTS2()){
        writeMesh("v2_swap_adp"+ std::to_string(niter)+".meshb",msh);
        msh.met.writeMetricFile("v2_swap_adp"+ std::to_string(niter)+".solb");
        writeBackLinks("v2_swap_adp_poi2bak" + std::to_string(niter), msh);
      }
      if(DOPRINTS2()){
        getLengthEdges(msh,msh.get_tdim(),-1,ilned,rlned,lenstat);
        CPRINTF1(" - Length qua short = {}\n",lenstat.qua_short);
        CPRINTF1(" -            long  = {}\n",lenstat.qua_long);
        if(DOPRINTS3()) print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
        getmetquamesh<MFT>(msh,msh.get_tdim(),AsDeg::P1,AsDeg::P1,
                           &iinva,&qmin,&qmax,&qavg,&lquae);
        CPRINTF2(" - Quality min = {:15.7e} \n",qmin);
        CPRINTF2("           max = {:15.7e} \n",qmax);
        CPRINTF2("           avg = {:15.7e} \n",qavg);
        if(DOPRINTS2()) writeField("v2_swapqua_adp"+ std::to_string(niter)+".solb",
                                 msh,SolTyp::P0Elt,lquae);
        CPRINTF2("------------------------------------------------------------\n");
        CPRINTF2("- iteration {} swaps stat = {:.2e} time = {:.2e}s \n",niter,stat,t1-t0);
        CPRINTF2("------------------------------------------------------------\n");
      }
    }
    #ifndef NDEBUG
    check_topo(msh,1);
    #endif

    // 3. Insert on long edges 


    t0 = get_wall_time();
    stat  = insertLongEdges<MFT,gdim,ideg>(msh, &ninser,ithrd1, ithrd2);
    stat0 = MAX(stat0,stat);
    t1 = get_wall_time();
    tinsert += t1-t0;

    msh.cleanup();
    #ifndef NDEBUG
    check_topo(msh,1);
    #endif

    if(DOPRINTS2()){
      writeMesh("v2_insert_adp"+ std::to_string(niter)+".meshb",msh);
      msh.met.writeMetricFile("v2_insert_adp"+std::to_string(niter)+".solb");
      writeBackLinks("v2_insert_adp_poi2bak" + std::to_string(niter), msh);
    }
    if(DOPRINTS2()){
      getLengthEdges(msh,msh.get_tdim(),-1,ilned,rlned,lenstat);
      CPRINTF1(" - Length qua short = {}\n",lenstat.qua_short);
      CPRINTF1(" -            long  = {}\n",lenstat.qua_long);
      if(DOPRINTS3()) print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
      CPRINTF2("------------------------------------------------------------\n");
      CPRINTF2("- iteration {} insertions stat = {:.2e} time = {:.2e}s \n",niter,stat,t1-t0);
      CPRINTF2("------------------------------------------------------------\n");
    }

    if(msh.param->opt_unif){
      // 4. Smoothing (heuristic) -> fast but bad; improve
      t0 = get_wall_time();
      double stat = smoothInterior_Ball<MFT>(msh,QuaFun::Unit,ithrd1,ithrd2);
      stat0 = MAX(stat, stat0); 
      t1 = get_wall_time();
      if(DOPRINTS2()) writeMesh("v2_unif_adp"+ std::to_string(niter)+".meshb",msh);
      if(DOPRINTS2()) msh.met.writeMetricFile("v2_unif_adp"+ std::to_string(niter)+".solb");    
      if(DOPRINTS2()) writeBackLinks("v2_unif_adp_poi2bak" + std::to_string(niter), msh);
      if(DOPRINTS2()){
        CPRINTF2("------------------------------------------------------------\n");
        CPRINTF2("- iteration {} unif ball stat = {:.2e} time = {:.2e}s \n",niter,stat,t1-t0);
        CPRINTF2("------------------------------------------------------------\n");
        getmetquamesh<MFT>(msh,msh.get_tdim(),AsDeg::P1,AsDeg::P1,
                           &iinva,&qmin,&qmax,&qavg,&lquae);
        CPRINTF2(" - Quality min = {:15.7e} \n",qmin);
        CPRINTF2("           max = {:15.7e} \n",qmax);
        CPRINTF2("           avg = {:15.7e} \n",qavg);
      }
    }



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
        CPRINTF1("- iteration {} flat collapse noper = {} stat = {:.2e} time = {:.2e}s \n",niter,noper,stat,t1-t0);
        CPRINTF1("------------------------------------------------------------\n");
      }

      #ifndef NDEBUG
      check_topo(msh,0);
      #endif

    }else{
      CPRINTF1("## reinsertFlat disabled in case gdim = {} tdim = {} \n", msh.idim, msh.get_tdim());
    }
    

    getLengthEdges(msh,msh.get_tdim(),-1,ilned,rlned,lenstat);
    int ndigit = ceil(log10((double)msh.npoin
                        + (double)msh.nelem
                        + (double)msh.nface
                        + (double)msh.nedge) ) + 1;

    double tloop1 = get_wall_time();

    std::string fmt = 
    "{}-- Adp loop {:3} / {:3} time {:.2e}ss " 
    "{:" + std::to_string(ndigit) + "} inser "
    "{:" + std::to_string(ndigit) + "} coll "
    "{:" + std::to_string(ndigit) + "} swap, "
    + "{:7.2f}% unit, op stat = {:.2e} \n";
    // "-- Adp loop {:3} / {:3} inser {} coll {} swap {}, {}% unit, op stat = {:.2e} \n"
    //CPRINTF1(fmt.c_str(), niter, miter, tloop1 - tloop0,
    //         msh.npoin, ninser, ncoll, nswap,
    //         100*lenstat.prop_unit, stat0);
    if(DOPRINTS1()) fmt::print(LOGFILE__, fmt.c_str(), spaces_string__, 
             niter,miter, tloop1 - tloop0, ninser,ncoll,nswap, 100*lenstat.prop_unit,stat0);

    //if(niter == 1){
    //  printf("## DEBUG SET MAX PRINTS\n");
    //  wait();
    //  msh.param->iverb = 5;
    //  msh.param->ivdepth = 5;
    //}
    if(lenstat.prop_unit*100 >= msh.param->adp_unit_stop){
      CPRINTF1("------------------------------------------------------------\n");
      CPRINTF1("- {:7.2f}% edges unit exit threshold = {:7.2f}\n",100*lenstat.prop_unit,
                msh.param->adp_unit_stop);
      break;
    }

    //bool dosmoo_adp = (iopt_niter < msh.param->adp_opt_niter|| msh.param->adp_opt_niter < 0)
    //                    && !msh.param->opt_unif;
    bool stagn = stat0 < msh.param->adp_stagn_stop
              || stat0 < minstat;
              //|| stat0 < 5.0e-2 && dosmoo_adp;
             // || abs(stat0 - stat_prev) < 1.0e-6 ;// This last criterion just catches cycles

    //stat_prev = stat0;

    if(stagn){
      CPRINTF1(" - low stat = {:.2e} break or optimize\n",stat0);
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
          CPRINTF1(" - low optim stat {:.2e} break\n",stat);
          break;
        }

        // If we continue, unconstrain the points now
        // We notice that the boundary in 3D, and the whole mesh in 2D 
        // improves by unconstraining the points, but not the interior in 3D.
        // So, regardless of dimension, we unconstrain only dim <= 2 points.
        for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
          if(msh.poi2ent(ipoin,0) < 0) continue;
          if(msh.getpoitdim(ipoin) > 2) continue;
          msh.poicstr[ipoin] = false;
        }
        //msh.poicstr.fill(false);
      }else{
        break;
      }
    }

  }// for niter


  ttotal = get_wall_time() - ttotal;

  msh.cleanup();

  CPRINTF1("-- Adaptation end total time = {:.2e}s \n",ttotal);
  CPRINTF1(" - insertion time = {:.2e}s \n",tinsert);
  CPRINTF1(" -  collapse time = {:.2e}s \n",tcollapse);
  CPRINTF1(" -      swap time = {:.2e}s \n",tswap);
  CPRINTF1(" - smoothing time = {:.2e}s \n",tsmooth);

  if(DOPRINTS1() || DOPRINTS3()){
    writeMesh("adapt_end.meshb",msh);
    msh.met.writeMetricFile("adapt_end.solb");
    if(DOPRINTS3()){
      getmetquamesh<MFT>(msh,msh.get_tdim(),AsDeg::P1,AsDeg::P1,
                         &iinva,&qmin,&qmax,&qavg,&lquae);
      writeField("adapt_end.qua.solb",msh,SolTyp::P0Elt,lquae);
    }
  }

  if(DOPRINTS1()) statMesh();

  #if 0
  printf("## DEBUG REMOVE THIS\n");
  double voltot = msh.getDomainVolume();
  int npvol = 0;
  for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
    if(msh.poi2ent(ipoin,0) < 0) continue;
    npvol ++;
  }
  printf("-- Point density = {} ; vol = {} np = {}\n",voltot/npvol, voltot, npvol);
  const double pi = 3.141592653589793238462643383279502884;
  //double dens0 = msh.get_tdim() == 2 ? pi / 4 : 2*pi/3;
  double dens0 = msh.get_tdim() == 2 ? pi / 4 : 1.0/sqrt(2);
  printf(" - expected = {} err = {}\n",dens0,abs(dens0 - voltot/npvol));
  wait();
  #endif
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void MetrisRunner::adaptMesh0<MetricFieldFE        ,2,n>();\
template void MetrisRunner::adaptMesh0<MetricFieldFE        ,3,n>();\
template void MetrisRunner::adaptMesh0<MetricFieldAnalytical,2,n>();\
template void MetrisRunner::adaptMesh0<MetricFieldAnalytical,3,n>();
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

}//end namespace
