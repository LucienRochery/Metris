//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "MetrisRunner.hxx"

#include "../Adaptation/msh_collapse.hxx"
#include "../Adaptation/msh_swap.hxx"
#include "../Adaptation/Insertion/msh_insert.hxx"
#include "../Adaptation/msh_reinsert_flat.hxx"
#include "../Adaptation/msh_lineadapt.hxx"

#include "../Mesh/Mesh.hxx"
#include "../quality/msh_metqua.hxx"
#include "../io_libmeshb.hxx"
#include "../smoothing/msh_smooball.hxx"
#include "../smoothing/msh_smoolen.hxx"
#include "../msh_checktopo.hxx"
#include "../aux_histogram.hxx"
#include "../msh_lenedg.hxx"
#include "../linalg/det.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/CT_loop.hxx"
#include "../utils/aux_misc.hxx"
#include "../utils/aux_timer.hxx"

#include "../low_geo/misc.hxx"

#include "../cavity/msh_cavity.hxx"
#include "../Adaptation/Insertion/low_insert.hxx"
#include "../Adaptation/Insertion/EdgeSeed.hxx"
#include "../low_geo/lenedg.hxx"
#include "../Adaptation/low_collapse.hxx"
#include "../aux_badEntHandler.hxx"

#ifdef DIAGNOSIS_QUALALGO
#include <iostream>
#include <fstream>
#endif

namespace Metris{


void MetrisRunner::adaptMesh(){
  CT_FOR0_INC(2,3,gdim){if(gdim == msh_g->idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh_g->curdeg == ideg){
      if(this->metricFE){
        adaptMesh0<MetricFieldFE        ,gdim,ideg>(msh_g->get_tdim());
      }else{
        adaptMesh0<MetricFieldAnalytical,gdim,ideg>(msh_g->get_tdim());
      }
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
}


void MetrisRunner::adaptMesh2(){
  CT_FOR0_INC(2,3,gdim){if(gdim == msh_g->idim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh_g->curdeg == ideg){
      #ifdef TESTQUALITYALGO
      if(this->metricFE){
        adaptMeshQuality0<MetricFieldFE        ,gdim,ideg>(msh_g->get_tdim());
      }else{
        adaptMeshQuality0<MetricFieldAnalytical,gdim,ideg>(msh_g->get_tdim());
      }
      #else
      for(int tdim = 1; tdim <= msh_g->get_tdim(); tdim++){
        GETVDEPTH(this->param);
        if(this->metricFE){
          adaptMesh0<MetricFieldFE        ,gdim,ideg>(tdim);
        }else{
          adaptMesh0<MetricFieldAnalytical,gdim,ideg>(tdim);
        }
        //CPRINTF1("\n\n");
        //if(tdim == 2){
        //  printf("\n\n## DEBUG BREAK AT TDIM = 2\n\n\n");
        //  wait();
        //  break;
        //}
      }
      #endif
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
}


// Profiling is attrocious if the template parameters are unrolled within the function
template<class MFT,int gdim,int ideg>
void MetrisRunner::adaptMesh0(int tdim){

#ifdef TESTQUALITYALGO
  METRIS_ASSERT_MSG(0, "TESTQUALITYALGO must be undefined to use adaptMesh0");
#else
  if(tdim > gdim) return;

  if(param->adp_niter == 0) return;

  GETVDEPTH(this->param);

  Mesh<MFT> &msh = static_cast<Mesh<MFT>&>(*msh_g);

  //METRIS_THROW(TODOExcept());

  // Make it an option
  const double minstat = 1.0e-12;
  const int miter = param_.adp_niter;
  int irestart = 0;

  msh.cleanup();

  CPRINTF1("-- START adaptMesh tdim = {} with miter = {} \n",tdim,miter);
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
    getLengthEdges<MFT>(msh,tdim,-1,ilned,rlned,lenstat);
    print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
  }
  double t0,t1;

  // This is the common thread for all routines. Tagged elements are ignored
  const int ithrdfro = 0;
  const int ithrd1 = 1;
  const int ithrd2 = 2;
  const int ithrd3 = 3;
  msh.tag[ithrdfro]++;

  if(msh.CAD() && msh.param->adp_line_adapt && tdim == 1){

    t0 = get_cpu_time();
    adaptGeoLines<MFT>(msh);
    t1 = get_cpu_time();
    CPRINTF1(" - adaptGeoLines time = {:.2e}s \n",t1-t0);
    if(DOPRINTS2()) writeMesh("v2_geolines_adp",msh);
    if(DOPRINTS2()) msh.met.writeMetricFile("v2_geolines_adp");


    //fmt::print("## WAIT HERE \n");
    //wait();

    if(DOPRINTS1()){
      getLengthEdges<MFT>(msh,1,-1,ilned,rlned,lenstat);
      print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length (lines)");
    }

    swapMesh<MFT,gdim,ideg>(msh, Defaults::swapOptAdapt, &nswap, ithrdfro, ithrd1, ithrd2);

    #ifndef NDEBUG
    check_topo(msh,1);
    #endif


    if(tdim == 1){
      if(DOPRINTS1()) statMesh(tdim);
      return;
    }

  }

  // Will never exceed this
  const int miter_max = 100;

  dblAr1 lquae(msh.nface);

  double tinsert = 0, tcollapse = 0, tswap = 0, tsmooth = 0;
  double ttotal = get_cpu_time();

  int iopt_niter = 0;
  double stat0 = 1;
  //for(int niter = 1; niter <= miter || ( miter < 0 && niter <= miter_max
  //                                    && stat0 > 0.1); niter++){
  //double stat_prev = stat0;
  msh.tag[ithrdfro]++;
  for(int niter = 1; niter <= miter || (miter < 0 && niter < miter_max); niter++){
    stat0 = 0;
    double tloop0 = get_cpu_time();

    qmax_suf = qavg * MAX(10 / (niter * 1.0), 1.0);
    //qmax_suf = 1.0 - (niter - 1) / (double) miter;
    //qmax_suf = 0;
    double stat;
    // 1. Collapse short edges
    t0 = get_cpu_time();
    stat  = collapseShortEdges<MFT,gdim,ideg>(msh, tdim, qmax_suf, &ncoll, ithrdfro, ithrd1, ithrd2, ithrd3);
    stat0 = MAX(stat0,stat);
    t1 = get_cpu_time();
    tcollapse += t1-t0;

    if(DOPRINTS2()){
      writeMesh("v2_collapse_adp" + std::to_string(tdim) + "D"+ std::to_string(niter) + ".meshb",msh);
      msh.met.writeMetricFile("v2_collapse_adp" + std::to_string(tdim) + "D" + std::to_string(niter)+".solb");
      writeBackLinks("v2_collapse_adp_poi2bak" + std::to_string(niter), msh);
    }

    if(DOPRINTS2()){
      getLengthEdges(msh,tdim,-1,ilned,rlned,lenstat);
      CPRINTF1(" - Length qua short = {}\n",lenstat.qua_short);
      CPRINTF1(" -            long  = {}\n",lenstat.qua_long);
      if(DOPRINTS3()) print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");

      getmetquamesh<MFT>(msh,tdim,AsDeg::P1,AsDeg::P1,
                         &iinva,&qmin,&qmax,&qavg,&lquae);
      CPRINTF1(" - Quality min = {:15.7e} \n",qmin);
      CPRINTF1("           max = {:15.7e} \n",qmax);
      CPRINTF1("           avg = {:15.7e} \n",qavg);
      writeField("v2_collapsequa_adp" + std::to_string(tdim) + "D" + std::to_string(niter)+".solb",
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
      t0 = get_cpu_time();
      stat  = swapMesh<MFT,gdim,ideg>(msh, Defaults::swapOptAdapt, &nswap, ithrdfro, ithrd1, ithrd2);
      stat0 = MAX(stat0,stat);
      t1 = get_cpu_time();
      tswap += t1-t0;

      if(DOPRINTS2()){
        writeMesh("v2_swap_adp" + std::to_string(tdim) + "D" + std::to_string(niter)+".meshb",msh);
        msh.met.writeMetricFile("v2_swap_adp" + std::to_string(tdim) + "D" + std::to_string(niter)+".solb");
        writeBackLinks("v2_swap_adp_poi2bak" + std::to_string(niter), msh);
      }
      if(DOPRINTS2()){
        getLengthEdges(msh,tdim,-1,ilned,rlned,lenstat);
        CPRINTF1(" - Length qua short = {}\n",lenstat.qua_short);
        CPRINTF1(" -            long  = {}\n",lenstat.qua_long);
        if(DOPRINTS3()) print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
        getmetquamesh<MFT>(msh,tdim,AsDeg::P1,AsDeg::P1,
                           &iinva,&qmin,&qmax,&qavg,&lquae);
        CPRINTF2(" - Quality min = {:15.7e} \n",qmin);
        CPRINTF2("           max = {:15.7e} \n",qmax);
        CPRINTF2("           avg = {:15.7e} \n",qavg);
        if(DOPRINTS2()) writeField("v2_swapqua_adp" + std::to_string(tdim) + "D" + std::to_string(niter)+".solb",
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

    t0 = get_cpu_time();
    stat  = insertLongEdges<MFT,gdim,ideg>(msh, tdim, &ninser,ithrd1, ithrd2);
    stat0 = MAX(stat0,stat);
    t1 = get_cpu_time();
    tinsert += t1-t0;

    //printf("Debug wait here\n");
    //wait();

    msh.cleanup();
    #ifndef NDEBUG
    check_topo(msh,1);
    #endif

    if(DOPRINTS2()){
      writeMesh("v2_insert_adp"  + std::to_string(tdim) + "D" + std::to_string(niter)+".meshb",msh);
      msh.met.writeMetricFile("v2_insert_adp" + std::to_string(tdim) + "D" + std::to_string(niter)+".solb");
      writeBackLinks("v2_insert_adp_poi2bak" + std::to_string(niter), msh);
    }
    if(DOPRINTS2()){
      getLengthEdges(msh,tdim,-1,ilned,rlned,lenstat);
      CPRINTF1(" - Length qua short = {}\n",lenstat.qua_short);
      CPRINTF1(" -            long  = {}\n",lenstat.qua_long);
      if(DOPRINTS3()) print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
      CPRINTF2("------------------------------------------------------------\n");
      CPRINTF2("- iteration {} insertions stat = {:.2e} time = {:.2e}s \n",niter,stat,t1-t0);
      CPRINTF2("------------------------------------------------------------\n");
    }

    if(msh.param->opt_unif && tdim == msh.get_tdim()){
      // 4. Smoothing (heuristic) -> fast but bad; improve
      t0 = get_cpu_time();
      #ifdef TESTQUAFSIZESHAPE
      double stat = smoothInterior_Ball<MFT>(msh,QuaFun::SizeShape,ithrd1,ithrd2);
      #else
      double stat = smoothInterior_Ball<MFT>(msh,QuaFun::Unit,ithrd1,ithrd2);
      #endif
      stat0 = MAX(stat, stat0);
      t1 = get_cpu_time();
      if(DOPRINTS2()) writeMesh("v2_unif_adp"+ std::to_string(niter)+".meshb",msh);
      if(DOPRINTS2()) msh.met.writeMetricFile("v2_unif_adp"+ std::to_string(niter)+".solb");
      if(DOPRINTS2()) writeBackLinks("v2_unif_adp_poi2bak" + std::to_string(niter), msh);
      if(DOPRINTS2()){
        CPRINTF2("------------------------------------------------------------\n");
        CPRINTF2("- iteration {} unif ball stat = {:.2e} time = {:.2e}s \n",niter,stat,t1-t0);
        CPRINTF2("------------------------------------------------------------\n");
        getmetquamesh<MFT>(msh,tdim,AsDeg::P1,AsDeg::P1,
                           &iinva,&qmin,&qmax,&qavg,&lquae);
        CPRINTF2(" - Quality min = {:15.7e} \n",qmin);
        CPRINTF2("           max = {:15.7e} \n",qmax);
        CPRINTF2("           avg = {:15.7e} \n",qavg);
      }
    }


    if(msh.param->adp_smoo_len){
      t0 = get_cpu_time();
      smoothMeshLength(msh, tdim, ithrd1, ithrd2);
      t1 = get_cpu_time();

      if(DOPRINTS2()){
        writeMesh("v2_smoolen_adp"  + std::to_string(tdim) + "D" + std::to_string(niter)+".meshb",msh);
        msh.met.writeMetricFile("v2_smoolen_adp" + std::to_string(tdim) + "D" + std::to_string(niter)+".solb");
        writeBackLinks("v2_smoolen_adp_poi2bak" + std::to_string(niter), msh);
      }
      if(DOPRINTS2()){
        getLengthEdges(msh,tdim,-1,ilned,rlned,lenstat);
        CPRINTF1(" - Length qua short = {}\n",lenstat.qua_short);
        CPRINTF1(" -            long  = {}\n",lenstat.qua_long);
        if(DOPRINTS3()) print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length");
        CPRINTF2("------------------------------------------------------------\n");
        CPRINTF2("- iteration {} length smoothing stat = {:.2e} time = {:.2e}s \n",niter,stat,t1-t0);
        CPRINTF2("------------------------------------------------------------\n");
      }

      #ifndef NDEBUG
      check_topo(msh,0);
      #endif
    }


    if(msh.idim == tdim){
      t0 = get_cpu_time();
      int noper = reinsertFlat<MFT,gdim,ideg>(msh);
      t1 = get_cpu_time();
      msh.cleanup();
      stat  = noper / (double) msh.nface;
      stat0 = MAX(stat0, stat);
      if(DOPRINTS1()){
        if(DOPRINTS2() && noper >= 0) writeMesh("v2_flat_opt" + std::to_string(tdim) + "D" + std::to_string(niter)+".meshb",msh);
        if(DOPRINTS2() && noper >= 0) msh.met.writeMetricFile("v2_flat_opt" + std::to_string(tdim) + "D" + std::to_string(niter)+".solb");
        CPRINTF1("------------------------------------------------------------\n");
        CPRINTF1("- iteration {} flat collapse noper = {} stat = {:.2e} time = {:.2e}s \n",niter,noper,stat,t1-t0);
        CPRINTF1("------------------------------------------------------------\n");
      }

      #ifndef NDEBUG
      check_topo(msh,0);
      #endif

    }else{
      CPRINTF1("## reinsertFlat disabled in case gdim = {} tdim = {} \n", msh.idim, tdim);
    }


    getLengthEdges(msh,tdim,-1,ilned,rlned,lenstat);
    int ndigit = ceil(log10((double)msh.npoin
                        + (double)msh.nelem
                        + (double)msh.nface
                        + (double)msh.nedge) ) + 1;

    double tloop1 = get_cpu_time();

    std::string fmt =
    "{}-- Adp loop {:3} / {:3} dim {} time {:.2e}s "
    "{:" + std::to_string(ndigit) + "} inser "
    "{:" + std::to_string(ndigit) + "} coll "
    "{:" + std::to_string(ndigit) + "} swap, "
    + "{:7.2f}% unit, op stat = {:.2e} \n";
    // "-- Adp loop {:3} / {:3} inser {} coll {} swap {}, {:.2f}% unit, op stat = {:.2e} \n"
    //CPRINTF1(fmt.c_str(), niter, miter, tloop1 - tloop0,
    //         msh.npoin, ninser, ncoll, nswap,
    //         100*lenstat.prop_unit, stat0);
    if(DOPRINTS1()) fmt::print(LOGFILE__, fmt.c_str(), spaces_string__,
             niter,miter, tdim, tloop1 - tloop0, ninser,ncoll,nswap, 100*lenstat.prop_unit,stat0);

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
      // If we continue, unconstrain the points now
      // We notice that the boundary in 3D, and the whole mesh in 2D
      // improves by unconstraining the points, but not the interior in 3D.
      // So, regardless of dimension, we unconstrain only dim <= 2 points.
      bool uncstr = false;
      if(++irestart < MAX(2,msh.param->adp_opt_niter)){
        for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
          if(msh.isdeadpoint(ipoin)) continue;
          if(msh.getpoitdim(ipoin) > 2) continue;
          msh.poicstr[ipoin] = false;
          uncstr = true;
        }
      }

      CPRINTF1(" - low stat = {:.2e} break or optimize\n",stat0);
      if(uncstr) CPRINTF1(" - unconstrained points\n");

      if(niter >= miter -1) break;
      if(msh.param->opt_niter > 0 &&
        (iopt_niter < msh.param->adp_opt_niter|| msh.param->adp_opt_niter < 0)
         && !msh.param->opt_unif){
        iopt_niter++;
        double tsmo0 = get_cpu_time();
        stat = optimMesh();
        double tsmo1 = get_cpu_time();
        tsmooth += tsmo1 - tsmo0;
        msh.tag[ithrdfro]++;
        if(DOPRINTS2()){
          writeMesh("v2_optim_adp" + std::to_string(tdim) + "D" + std::to_string(iopt_niter), msh);
          msh.met.writeMetricFile("v2_optim_adp" + std::to_string(tdim) + "D" + std::to_string(iopt_niter));
          writeBackLinks("v2_optim_adp_poi2bak" + std::to_string(niter), msh);
        }
        if(stat < minstat){
          CPRINTF1(" - low optim stat {:.2e} break\n",stat);
          break;
        }

        //msh.poicstr.fill(false);
      }else if (!uncstr){
        break;
      }
    }

  }// for niter


  ttotal = get_cpu_time() - ttotal;

  msh.cleanup();

  CPRINTF1("-- Adaptation dim {} end total time = {:.2e}s \n",tdim,ttotal);
  CPRINTF1(" - insertion time = {:.2e}s \n",tinsert);
  CPRINTF1(" -  collapse time = {:.2e}s \n",tcollapse);
  CPRINTF1(" -      swap time = {:.2e}s \n",tswap);
  CPRINTF1(" - smoothing time = {:.2e}s \n",tsmooth);

  if(DOPRINTS1() || DOPRINTS3()){
    std::string fname = "adapt_end_" + std::to_string(tdim);
    writeMesh(fname + ".meshb",msh);
    msh.met.writeMetricFile(fname + ".solb");
    //if(DOPRINTS3()){
    //  getmetquamesh<MFT>(msh,tdim,AsDeg::P1,AsDeg::P1,
    //                     &iinva,&qmin,&qmax,&qavg,&lquae);
    //  writeField("adapt_end_" + std::to_string(tdim) + ".qua.solb",msh,SolTyp::P0Elt,lquae);
    //}
  }

  if(DOPRINTS1()) statMesh(tdim);

  #if 0
  printf("## DEBUG REMOVE THIS\n");
  double voltot = msh.getDomainVolume();
  int npvol = 0;
  for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
    if(msh.isdeadpoint(ipoin)) continue;
    npvol ++;
  }
  printf("-- Point density = {} ; vol = {} np = {}\n",voltot/npvol, voltot, npvol);
  const double pi = 3.141592653589793238462643383279502884;
  //double dens0 = tdim == 2 ? pi / 4 : 2*pi/3;
  double dens0 = tdim == 2 ? pi / 4 : 1.0/sqrt(2);
  printf(" - expected = {} err = {}\n",dens0,abs(dens0 - voltot/npvol));
  wait();
  #endif

#endif
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void MetrisRunner::adaptMesh0<MetricFieldFE        ,2,n>(int tdim);\
template void MetrisRunner::adaptMesh0<MetricFieldFE        ,3,n>(int tdim);\
template void MetrisRunner::adaptMesh0<MetricFieldAnalytical,2,n>(int tdim);\
template void MetrisRunner::adaptMesh0<MetricFieldAnalytical,3,n>(int tdim);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

// =============================================================================================== //
// =============================================================================================== //

template<class MFT,int gdim,int ideg>
void MetrisRunner::adaptMeshQuality0(int tdim){
  #ifdef TESTQUALITYALGO
  if(tdim > gdim) return;

  GETVDEPTH(this->param);

  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*msh_g);

  msh.cleanup();

  CPRINTF1("-- START adaptMeshQuality0 tdim = {}\n",tdim);
  if(DOPRINTS1()){
    writeMesh("debug_adapt_inp", msh);
    msh.met.writeMetricFile("debug_adapt_inp");
  }

  msh.met.setSpace(MetSpace::Exp);
  msh.setBasis(FEBasis::Lagrange);

  #ifdef DIAGNOSIS_QUALALGO
  std::string diagnosisFile = "diagnosisQualAlgo.txt";
  std::fstream foutput;
  foutput.open(diagnosisFile, std::fstream::out);
  METRIS_ASSERT_MSG(foutput.good(), "Error opening file: " + diagnosisFile);

  foutput << std::setw(4) << "Iter"
          << std::setw(30) << "qualError"
          << std::setw(30) << "xavgWorst"
          << std::setw(30) << "yavgWorst"
          << std::setw(30) << "trySmooth"
          << std::setw(30) << "statusSmooth"
          << std::setw(30) << "tryInsert"
          << std::setw(30) << "statusInsert"
          << std::setw(30) << "tryCollapse"
          << std::setw(30) << "statusCollapse"
          << std::setw(30) << "qualOp"
          << std::setw(30) << "xavgOp"
          << std::setw(30) << "yavgOp"
          << std::endl;
  #endif

  // This is the common thread for all routines. Tagged elements are ignored
  const int ithrdfro = 0;
  const int ithrd1 = 1;
  const int ithrd2 = 2;
  const int ithrd3 = 3;
  msh.tag[ithrdfro]++;

  const double alpha = 0.1;
  const double badX = 100;
  const double lengthThreshold = sqrt(2);

  // initial number of elements
  const int nentt0 = msh.nentt(tdim);

  // initial quality array
  bool iinva = false; double qmin = 0, qmax = 0, qavg = 0;
  dblAr1 lquae(nentt0);

  // initial quality computation
  getmetquamesh<MFT, QuaFun::SizeShape>(msh,tdim,AsDeg::P1,AsDeg::P1,&iinva,&qmin,&qmax,&qavg,&lquae);

  std::vector<int> sortedIDs(nentt0);
  std::iota(sortedIDs.begin(), sortedIDs.end(), 0);

  // sort from max (worst) to min (best) quality. The handler is expecting like this
  std::sort(sortedIDs.begin(), sortedIDs.end(),
            [&](int a, int b){ return lquae[a] > lquae[b]; });

  // initialize handler for top X% worst, K, and remainder, R
  BadEntHandler handlerTopX(tdim, badX, alpha);

  const intAr2& ent2poi = msh.ent2poi(tdim);
  const intAr2& ent2ent = msh.ent2ent(tdim);
  intAr2& ent2tag = msh.ent2tag(tdim);

  handlerTopX.setCallbacks(
                            [&](int ientt){ return lquae[ientt]; },
                            [&](int ientt){ return isdeadent(ientt,ent2poi); });

  // builds K and R
  handlerTopX.seedFromSortedIDs(sortedIDs);

  // ============================================= //
  // main loop:
  // - traverse K from worst to best
  // - if successful operation, update K and restart
  // - ends when reaching end of K
  // ============================================= //

  // declare and reserve memory for cavity and cavity work array
  // costly to do at lower level, just declare once here and re-use
  MshCavity cav(100,100,1);
  CavWrkArrs work;

  // Error counting:
  const int mcaverr = CAV_ERR_NERROR;
  intAr1 lcaverr(mcaverr);
  const int minserr = INS2D_ERR_NERROR;
  intAr1 linserr(minserr);

  int ntrySmoothing = 0; int nSuccessSmoothing = 0;
  int ntryInsert = 0; int nSuccessInsert = 0;
  int ntryCollapse = 0; int nSuccessCollapse = 0;

  int iter = 0;
  int smooStreak = 0;
  while (true){

    bool didOperation = false;

    for (auto itK = handlerTopX.K.begin(); itK != handlerTopX.K.end(); itK++){

      const int ientt = itK->ientt;
      // if (itK->genStamp != handlerTopX.currentGen) continue;

      iter++;
      #ifdef DIAGNOSIS_QUALALGO
      const int nentt = msh.nentt(tdim);
      lquae.set_n(nentt);
      getmetquamesh<MFT, QuaFun::SizeShape>(msh,tdim,AsDeg::P1,AsDeg::P1,&iinva,&qmin,&qmax,&qavg,&lquae);

      std::vector<int> enttList(nentt);
      std::iota(enttList.begin(), enttList.end(), 0);

      std::sort(enttList.begin(), enttList.end(),
            [&](int a, int b){ return lquae[a] > lquae[b]; });

      int iworst = 0;
      while(isdeadent(enttList[iworst],ent2poi)) iworst++;

      const int worstEntt = enttList[iworst];
      METRIS_ASSERT(!isdeadent(worstEntt,ent2poi));
      const double worstQual = lquae[worstEntt];

      double xavgWorst = 0.;
      double yavgWorst = 0.;
      for (int ii = 0; ii < tdim + 1; ii++){

        const int ipoin = ent2poi(worstEntt, ii);

        xavgWorst += msh.coord(ipoin,0);
        yavgWorst += msh.coord(ipoin,1);
      }
      xavgWorst /= double(tdim+1);
      yavgWorst /= double(tdim+1);

      foutput << std::setw(4) << iter
              << std::setw(30) << std::setprecision(16) << std::scientific << worstQual
              << std::setw(30) << std::setprecision(16) << std::scientific << xavgWorst
              << std::setw(30) << std::setprecision(16) << std::scientific << yavgWorst;

      #endif



      #ifdef DIAGNOSIS_QUALALGO

      const bool trySmoo = 1;
      bool tryIns = 0;
      bool tryColl = 0;

      bool statusSmoo = 0;
      bool statusIns = 0;
      bool statusColl = 0;

      const double qualEntt = itK->qentt;

      double xavg = 0.;
      double yavg = 0.;
      for (int ii = 0; ii < tdim + 1; ii++){

        const int ipoin = ent2poi(ientt, ii);

        xavg += msh.coord(ipoin,0);
        yavg += msh.coord(ipoin,1);
      }
      xavg /= double(tdim+1);
      yavg /= double(tdim+1);
      #endif

      // we assume that if an operation is successful,
      // the operator takes care of informing the handler:
      // - new/modified entities, with their qualities
      // - deleted entities

      // ----------------------------- //
      // 1. Smoothing
      // ----------------------------- //

      ntrySmoothing++;

      double statSmoothing = smoothElement_Ball<MFT>(msh,ientt,handlerTopX,QuaFun::SizeShape,ithrd1,ithrd2);

      if (statSmoothing > 0){
        smooStreak++;
        nSuccessSmoothing++;
        didOperation = true;
        std::cout << "Successful smoothing in ientt = " << ientt << std::endl;
        handlerTopX.updateK(ientt,ent2ent,ent2tag,msh.tag[ithrd1]+1,ithrd1);
        #ifdef DIAGNOSIS_QUALALGO
        statusSmoo = 1;
        statusIns = 0;
        statusColl = 0;
        // goto LOGSTATUS;
        #else
        break;
        #endif
      }

      // ----------------------------- //
      //  -- End smoothing --
      // ----------------------------- //

      // ----------------------------- //
      // 2. Collapse or Insert
      // ----------------------------- //

      // we'll compute edge lengths in ientt and
      // use it to decide if we do collapse or insertion

      const int nedgl = tdim * (tdim+1)/2;
      const intAr2& ent2poi = msh.ent2poi(tdim);

      const intAr2 lnoed(nedgl, 2,
                         tdim == 1 ? lnoed1[0] :
                         tdim == 2 ? lnoed2[0] :
                         lnoed3[0]);

      double shortestLen = 1e30; int iedMin = -1;
      double longestLen  = 0;    int iedMax = -1;

      dblAr1 edgeLengths(nedgl);

      for (int ied = 0; ied < nedgl; ied++){
        double sz[2];
        double elen = getlenedg_geosz<MFT,gdim,ideg>(msh,ientt,tdim,ied,sz);
        edgeLengths[ied] = elen;

        if (elen > longestLen)  { longestLen  = elen; iedMax = ied; }
        if (elen < shortestLen) { shortestLen = elen; iedMin = ied; }
      }

      if (shortestLen >= lengthThreshold
          #ifdef DIAGNOSIS_QUALALGO
          && !(statusSmoo || statusIns || statusColl)
          #endif
         ){

        // all edges longer than long threshold, insert in longest

        // ----------------------------- //
        //  Insertion in longest edge
        // ----------------------------- //

        ntryInsert++;
        #ifdef DIAGNOSIS_QUALALGO
        tryIns = 1;
        #endif

        INCVDEPTH(msh.param);

        // just for printing info
        const int ip1 = ent2poi(ientt, lnoed(iedMax,0));
        const int ip2 = ent2poi(ientt, lnoed(iedMax,1));
        CPRINTF1(" - enact ins ientt = {} ied = {} edg {} {}\n",ientt,iedMax,ip1,ip2);

        EdgeSeed insertionSeed(msh, cav, tdim, tdim, ientt, iedMax);
        int ierro = insertEdge(msh,insertionSeed,0,false,
                               cav,work,lcaverr,handlerTopX,ithrd1,ithrd2);

        if (ierro <= 0){
          // success
          smooStreak = 0;
          nSuccessInsert++;
          didOperation = true;
          msh.poicstr[cav.ipins] = false;
          std::cout << "Successful insertion in ientt = " << ientt << ", ied = " << iedMax << ", edg (" << ip1 << "," << ip2 << ")" << std::endl;
          handlerTopX.updateK(ientt,ent2ent,ent2tag,msh.tag[ithrd1]+1,ithrd1);
           #ifdef DIAGNOSIS_QUALALGO
          statusSmoo = 0;
          statusIns = 1;
          statusColl = 0;
          // goto LOGSTATUS;
          #else
          break;
          #endif
        }
        else {
          CPRINTF2(" # longest-edge insertion failed ierro = {}\n", ierro);

          // steiner fallback
          // TODO: for now don't do this.
          // TODO: if we want to do this, needs to be modified to inform the handler
          // TODO: about deleted/new entities, similar to insertEdge
          // const bool doSteiner = false;

          // if (doSteiner && tdim <= 2 && insertionSeed.tdimp <= tdim && insertionSeed.tdimp < msh.get_tdim()){

          //   CPRINTF1(" -> try Steiner point insertion\n");
          //   int ierro_Steiner = insertSteiner(msh,insertionSeed,cav,work,lcaverr,ithrd1,ithrd2);

          //   if (ierro_Steiner == -1){
          //     CPRINTF1(" - insertSteiner succeeded, retry edge split\n");

          //     // rebuild the cavity around the original edge and retry insertEdge
          //     int iopen;
          //     intAr1 dum;
          //     shell(msh, insertionSeed.ipedg[0], insertionSeed.ipedg[1],
          //           tdim, ientt, dum, dum, cav.lctet, &iopen);
          //     METRIS_ASSERT(cav.lctet.get_n() > 0);

          //     // retry insertion
          //     ierro = insertEdge(msh,insertionSeed,0,false,
          //                       cav,work,lcaverr,handlerTopX,ithrd1,ithrd2);

          //     if (ierro <= 0){
          //       // success
          //       didOperation = true;
          //       msh.poicstr[cav.ipins] = false;
          //       handlerTopX.updateK(ientt,ent2ent,ent2tag,msh.tag[ithrd1]+1,ithrd1);;
          //       break;
          //     }
          //     else {
          //       CPRINTF2(" # retry after Steiner failed ierro = %d\n", ierro);
          //     }
            // }
            // else {
            //   CPRINTF1(" # insertSteiner failed ierro %d\n", ierro_Steiner);
            // }
          // }
        }

        // ----------------------------- //
        //  End insertion in longest edge
        // ----------------------------- //

      }
      else if (longestLen <= 1./lengthThreshold
               #ifdef DIAGNOSIS_QUALALGO
               && !(statusSmoo || statusIns || statusColl)
               #endif
              ){

        // all edges shorter than short threshold, collapse shortest

        // ----------------------------- //
        //  Collapse shortest edge
        // ----------------------------- //

        ntryCollapse++;
        #ifdef DIAGNOSIS_QUALALGO
        tryColl = 1;
        #endif

        INCVDEPTH(msh.param);

        // just for printing info
        const int ip1 = ent2poi(ientt, lnoed(iedMin,0));
        const int ip2 = ent2poi(ientt, lnoed(iedMin,1));
        CPRINTF1(" - enact collapse ientt = {} ied = {} edg {} {}\n",ientt,iedMin,ip1,ip2);

        int ierro = collapseEdge<MFT>(msh,tdim,ientt,iedMin,0,cav,work,lcaverr,handlerTopX,ithrd1,ithrd2,ithrd3);

        if (ierro == 0){
          smooStreak = 0;
          nSuccessCollapse++;
          didOperation = true;
          std::cout << "Successful collapse in ientt = " << ientt << ", ied = " << iedMin << ", edg (" << ip1 << "," << ip2 << ")" << std::endl;
          handlerTopX.updateK(ientt,ent2ent,ent2tag,msh.tag[ithrd1]+1,ithrd1);
           #ifdef DIAGNOSIS_QUALALGO
          statusSmoo = 0;
          statusIns = 0;
          statusColl = 1;
          // goto LOGSTATUS;
          #else
          break;
          #endif
        }
        else{
          CPRINTF2(" # shortest-edge collapse failed ierro = {}\n", ierro);
        }

        // ----------------------------- //
        //  End collapse shortest edge
        // ----------------------------- //

      }
      else if (log(longestLen) >= log(1./shortestLen)
               #ifdef DIAGNOSIS_QUALALGO
               && !(statusSmoo || statusColl || statusIns)
               #endif
              ){

        // longest edge deviates from one more than shortest edge

        // ----------------------------- //
        //  Insertion in longest edge
        // ----------------------------- //

        ntryInsert++;
        #ifdef DIAGNOSIS_QUALALGO
        tryIns = 1;
        #endif

        INCVDEPTH(msh.param);

        // just for printing info
        const int ip1 = ent2poi(ientt, lnoed(iedMax,0));
        const int ip2 = ent2poi(ientt, lnoed(iedMax,1));
        CPRINTF1(" - enact ins ientt = {} ied = {} edg {} {}\n",ientt,iedMax,ip1,ip2);

        EdgeSeed insertionSeed(msh, cav, tdim, tdim, ientt, iedMax);
        int ierro = insertEdge(msh,insertionSeed,0,false,
                               cav,work,lcaverr,handlerTopX,ithrd1,ithrd2);

        if (ierro <= 0){
          // success
          smooStreak = 0;
          nSuccessInsert++;
          didOperation = true;
          msh.poicstr[cav.ipins] = false;
          std::cout << "Successful insertion in ientt = " << ientt << ", ied = " << iedMax << ", edg (" << ip1 << "," << ip2 << ")" << std::endl;
          handlerTopX.updateK(ientt,ent2ent,ent2tag,msh.tag[ithrd1]+1,ithrd1);
           #ifdef DIAGNOSIS_QUALALGO
          statusSmoo = 0;
          statusIns = 1;
          statusColl = 0;
          // goto LOGSTATUS;
          #else
          break;
          #endif
        }
        else {
          CPRINTF2(" # longest-edge insertion failed ierro = {}\n", ierro);

          // steiner fallback
          // TODO: for now don't do this.
          // TODO: if we want to do this, needs to be modified to inform the handler
          // TODO: about deleted/new entities, similar to insertEdge
          // const bool doSteiner = false;

          // if (doSteiner && tdim <= 2 && insertionSeed.tdimp <= tdim && insertionSeed.tdimp < msh.get_tdim()){

          //   CPRINTF1(" -> try Steiner point insertion\n");
          //   int ierro_Steiner = insertSteiner(msh,insertionSeed,cav,work,lcaverr,ithrd1,ithrd2);

          //   if (ierro_Steiner == -1){
          //     CPRINTF1(" - insertSteiner succeeded, retry edge split\n");

          //     // rebuild the cavity around the original edge and retry insertEdge
          //     int iopen;
          //     intAr1 dum;
          //     shell(msh, insertionSeed.ipedg[0], insertionSeed.ipedg[1],
          //           tdim, ientt, dum, dum, cav.lctet, &iopen);
          //     METRIS_ASSERT(cav.lctet.get_n() > 0);

          //     // retry insertion
          //     ierro = insertEdge(msh,insertionSeed,0,false,
          //                       cav,work,lcaverr,handlerTopX,ithrd1,ithrd2);

          //     if (ierro <= 0){
          //       // success
          //       didOperation = true;
          //       msh.poicstr[cav.ipins] = false;
          //       handlerTopX.updateK(ientt,ent2ent,ent2tag,msh.tag[ithrd1]+1,ithrd1);;
          //       break;
          //     }
          //     else {
          //       CPRINTF2(" # retry after Steiner failed ierro = %d\n", ierro);
          //     }
            // }
            // else {
            //   CPRINTF1(" # insertSteiner failed ierro %d\n", ierro_Steiner);
            // }
          // }
        }

        // ----------------------------- //
        //  End insertion in longest edge
        // ----------------------------- //
      }
      else if (1
              #ifdef DIAGNOSIS_QUALALGO
              && !(statusSmoo || statusIns || statusColl)
              #endif
              ){

        // shortest edge deviates from one more than longest edge

        // ----------------------------- //
        //  Collapse shortest edge
        // ----------------------------- //

        ntryCollapse++;
        #ifdef DIAGNOSIS_QUALALGO
        tryColl = 1;
        #endif

        INCVDEPTH(msh.param);

        // just for printing info
        const int ip1 = ent2poi(ientt, lnoed(iedMin,0));
        const int ip2 = ent2poi(ientt, lnoed(iedMin,1));
        CPRINTF1(" - enact collapse ientt = {} ied = {} edg {} {}\n",ientt,iedMin,ip1,ip2);

        int ierro = collapseEdge<MFT>(msh,tdim,ientt,iedMin,0,cav,work,lcaverr,handlerTopX,ithrd1,ithrd2,ithrd3);

        if (ierro == 0){
          smooStreak = 0;
          nSuccessCollapse++;
          didOperation = true;
          std::cout << "Successful collapse in ientt = " << ientt << ", ied = " << iedMin << ", edg (" << ip1 << "," << ip2 << ")" << std::endl;
          handlerTopX.updateK(ientt,ent2ent,ent2tag,msh.tag[ithrd1]+1,ithrd1);
           #ifdef DIAGNOSIS_QUALALGO
          statusSmoo = 0;
          statusIns = 0;
          statusColl = 1;
          // goto LOGSTATUS;
          #else
          break;
          #endif
        }
        else{
          CPRINTF2(" # shortest-edge collapse failed ierro = {}\n", ierro);
        }

        // ----------------------------- //
        //  End collapse shortest edge
        // ----------------------------- //
      }

      #ifdef DIAGNOSIS_QUALALGO
      // LOGSTATUS:
      foutput << std::setw(30) << trySmoo
              << std::setw(30) << statusSmoo
              << std::setw(30) << tryIns
              << std::setw(30) << statusIns
              << std::setw(30) << tryColl
              << std::setw(30) << statusColl
              << std::setw(30) << std::scientific << qualEntt
              << std::setw(30) << std::scientific << xavg
              << std::setw(30) << std::scientific << yavg
              << std::endl;

      if (statusSmoo || statusIns || statusColl ) break;
      #endif
    }
    if (!didOperation || smooStreak >= 2500 || iter >= 8000 ) break;
  }

  std::cout << "iter = " << iter << std::endl;
  std::cout << "-- END adaptMeshQuality tdim = " << tdim << std::endl;
  std::cout << "-- ntrySmoothing = " << ntrySmoothing << ",  nSuccessSmoothing = " << nSuccessSmoothing << std::endl;
  std::cout << "-- ntryInsert = " << ntryInsert << ",  nSuccessInsert = " << nSuccessInsert << std::endl;
  std::cout << "-- ntryCollapse = " << ntryCollapse << ",  nSuccessCollapse = " << nSuccessCollapse << std::endl;

  #ifdef DIAGNOSIS_QUALALGO
  foutput.close();
  #endif

#endif
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void MetrisRunner::adaptMeshQuality0<MetricFieldFE        ,2,n>(int tdim);\
template void MetrisRunner::adaptMeshQuality0<MetricFieldFE        ,3,n>(int tdim);\
template void MetrisRunner::adaptMeshQuality0<MetricFieldAnalytical,2,n>(int tdim);\
template void MetrisRunner::adaptMeshQuality0<MetricFieldAnalytical,3,n>(int tdim);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

}//end namespace