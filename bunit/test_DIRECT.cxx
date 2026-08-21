//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_DIRECT

#include <boost/test/included/unit_test.hpp> 
#include <math.h>
#include <random>
#include <filesystem>

#include "Optimization/DIRECT.hxx"
#include "types.hxx"
#include "low_geo/misc.hxx"
#include "low_geo/measure.hxx"
#include "io_libmeshb.hxx"
#include "low_eval.hxx"
#include "utils/mprintf.hxx"

#include <boost/timer/progress_display.hpp>

using namespace Metris;

// If idim == 2, a shell
// otherwise a line of triangles
void gen_blob(int idim, int nblob, dblAr2& coord, intAr2& ent2poi);

void dump_init(int idim, std::string mshname, intAr2 &ent2poi, dblAr2 &coord);
void dump_tess(DIBLOB_args& args, intAr2& ent2poi, dblAr2 &coord);

double funsq(double *coop){
  return coop[0]*coop[0] + coop[1]*coop[1] + coop[2]*coop[2];
}

double funsq2(int idim, double *coop, double *coor0){
  if(idim == 2) return geterrl2<2>(coop,coor0);
  else          return geterrl2<3>(coop,coor0);
}

BOOST_AUTO_TEST_CASE(test_DIRECT) 
{

  auto costfun = funsq2;

  std::filesystem::create_directory("tmp");

  for(int idim = 2; idim <= 2; idim++){


    const int nsamp = 50;
    dblAr2 bary(nsamp,idim+1);
    std::uniform_real_distribution<double> unif(0.0,1.0);
    std::default_random_engine rng(0);
    for(int ii = 0; ii < nsamp; ii++){
      double sum = 0;
      do{
        for(int jj = 0; jj < idim+1; jj++){
          bary(ii,jj) = unif(rng);
          sum += bary(ii,jj);
        }
      }while(abs(sum) < 1.0e-16);

      for(int jj = 0; jj < idim+1; jj++){
        bary(ii,jj) /= sum;
      }
    }


    MetrisParameters param;
    param.iverb   = 1;
    param.ivdepth = 0;
    param.outmPrefix = "tmp/";
    dblAr2 coord;
    intAr2 ent2poi;
    double max_errx = -1.0, avg_errx = 0, max_errf = -1.0, avg_errf = 0;
    double avg_neval_tot = 0.0;
    int max_neval_tot = -1;
    unsigned long int navg = 0;

    // A value of 10 seems reasonable. 
    for(int nloc_switch : {-1, 10}){ // 3, 5, 
      printf("-- START tests with nloc_switch = %d \n", nloc_switch);

      double max_errx2 = -1.0, avg_errx2 = 0, max_errf2 = -1.0, avg_errf2 = 0;
      double avg_neval_tot2 = 0.0;
      int max_neval_tot2 = -1;
      unsigned long int navg2 = 0;
      for(int nblob = 1; nblob <= 10; nblob++){
        gen_blob(idim,nblob,coord,ent2poi);
        dump_init(idim,"tmp/init"+std::to_string(idim)+".meshb",ent2poi,coord);

        //boost::timer::progress_display progress(nblob*nsamp);

        double max_errx1 = -1.0, avg_errx1 = 0, max_errf1 = -1.0, avg_errf1 = 0;
        double avg_neval_tot1 = 0.0;
        int max_neval_tot1 = -1;
        unsigned long int navg1 = 0;
        for(int ishell = 0; ishell < nblob; ishell++){

          for(int isamp = 0; isamp < nsamp; isamp++){
            INCVDEPTH((&param));

            //++progress;

            double coor0[3];
            if(idim == 2){
              eval2<2,1>(coord, ent2poi[ishell], FEBasis::Lagrange, DifVar::None, DifVar::None,
                         bary[isamp], coor0, NULL, NULL);
            }else{
              eval3<3,1>(coord, ent2poi[ishell], FEBasis::Lagrange, DifVar::None, DifVar::None,
                         bary[isamp], coor0, NULL, NULL);
            }


            DIBLOB_args args(param);
            args.miter = 300;
            args.ftol  = 1.0e-6;
            args.fscale= 1;
            args.nloc_switch = nloc_switch;

            double xmin[3], fmin;
            int ifmin;
            dblAr2 peval;
            intAr2 leval;
            dblAr1 feval;
            int neval_tot = 0;
            while(true){
              INCVDEPTH((&param));

              try{
                DIBLOB(args, idim, nblob, leval, peval, feval, &ifmin, xmin, &fmin);
              }catch(const MetrisExcept &e){
                args.iflag = -1;
              }

              CPRINTF1(" - iter %d/%d fmin %e neval %d\n",args.niter,args.miter,fmin,leval.get_n());
              if(DOPRINTS3()) dump_tess(args, ent2poi, coord);

              if(args.iflag <= 0){
                break;
              }


              int neval = leval.get_n();
              bool isineval = false;
              for(int ieval = 0; ieval < neval; ieval++){
                INCVDEPTH((&param));
                int ielem = leval(ieval,1);

                double coop[3];
                if(idim == 2){
                  eval2<2,1>(coord, ent2poi[ielem], FEBasis::Lagrange, DifVar::None, DifVar::None,
                             peval[ieval], coop, NULL, NULL);
                }else{
                  eval3<3,1>(coord, ent2poi[ielem], FEBasis::Lagrange, DifVar::None, DifVar::None,
                             peval[ieval], coop, NULL, NULL);
                }
                feval[ieval] = costfun(idim,coop,coor0);
                neval_tot++;

                if(isineval) continue;
                // Check if the sought coordinates are in one of the eval requests.
                // With this convex function, this should always be the case. 
                int ieloc = leval(ieval,0);
                //if(idim == 3){
                //  double meas0 = getmeasentP1<3>(args.ent2pol[ieloc], args.coorl);
                //  if(meas0 < 0){
                //    printf("## NEGATIVE TESS ELEMENT %e \n",meas0);
                //    continue;
                //    break;
                //  }
                //  isineval = isintetP1(args.coorl[args.ent2pol(ieloc,0)], args.coorl[args.ent2pol(ieloc,1)],
                //                       args.coorl[args.ent2pol(ieloc,2)], args.coorl[args.ent2pol(ieloc,3)],
                //                       bary[isamp]);
                //  printf("Debug vol %e looking for %f %f %f %f \n",meas0,bary[isamp][0],bary[isamp][1],bary[isamp][2],bary[isamp][3]);
                //  for(int ii = 0; ii < 4; ii++){
                //    printf(" %d = %d : ",ii,args.ent2pol(ieloc,ii));
                //    dblAr1(idim,args.coorl[args.ent2pol(ieloc,ii)]).print();
                //  }
                //}
              }
              //if(!isineval){
                //printf("## iter %d / %d EVALS DO NOT CONTAIN BARY\n",args.niter,args.miter);
              //  break;
              //}
            }// while true

            if(DOPRINTS1()){
              CPRINTF1(" - neval %d fmin = %e in tet %d at ",neval_tot,fmin,ifmin);
              dblAr1(idim, xmin).print();
              CPRINTF1(" - Optimal solution is in tetra %d at ", ishell);
              dblAr1(idim,bary[isamp]).print();
            }

            max_neval_tot = METRIS_MAX(max_neval_tot, neval_tot);
            avg_neval_tot += (double) neval_tot;

            max_neval_tot1 = METRIS_MAX(max_neval_tot1, neval_tot);
            avg_neval_tot1 += (double) neval_tot;

            max_neval_tot2 = METRIS_MAX(max_neval_tot2, neval_tot);
            avg_neval_tot2 += (double) neval_tot;

            double errf = abs(fmin);
            double errx = idim == 2 ? sqrt(geterrl2<2>(xmin,bary[isamp])) 
                                    : sqrt(geterrl2<3>(xmin,bary[isamp]));
            //BOOST_TEST(ifmin == ishell);
            //BOOST_TEST(errf < args.ftol);
            max_errx = METRIS_MAX(errx, max_errx);
            max_errf = METRIS_MAX(errf, max_errf);
            avg_errx += errx;
            avg_errf += errf;
            navg++;

            max_errx1 = METRIS_MAX(errx, max_errx1);
            max_errf1 = METRIS_MAX(errf, max_errf1);
            avg_errx1 += errx;
            avg_errf1 += errf;
            navg1++;

            max_errx2 = METRIS_MAX(errx, max_errx2);
            max_errf2 = METRIS_MAX(errf, max_errf2);
            avg_errx2 += errx;
            avg_errf2 += errf;
            navg2++;

          }// for isamp

        }// for ishell


        avg_errx1 /= navg1;
        avg_errf1 /= navg1;
        avg_neval_tot1 /= navg1;

        printf("     - END tests nblob %d max x error %e f error %e neval %d\n",nblob,max_errx1,max_errf1,max_neval_tot1);
        printf("     - END tests nblob %d avg x error %e f error %e neval %d\n",nblob,avg_errx1,avg_errf1,(int)avg_neval_tot1);

      }// for nblob = 3
      avg_errx2 /= navg2;
      avg_errf2 /= navg2;
      avg_neval_tot2 /= navg2;

      printf("   - END tests nloc %d max x error %e f error %e neval %d\n",nloc_switch,max_errx2,max_errf2,max_neval_tot2);
      printf("   - END tests nloc %d avg x error %e f error %e neval %d\n",nloc_switch,avg_errx2,avg_errf2,(int)avg_neval_tot2);

    }// for nloc_switch

    avg_errx /= navg;
    avg_errf /= navg;
    avg_neval_tot /= navg;

    printf("-- END tests dim %d max x error %e f error %e neval %d\n",idim,max_errx,max_errf,max_neval_tot);
    printf("-- END tests dim %d avg x error %e f error %e neval %d\n",idim,avg_errx,avg_errf,(int)avg_neval_tot);

  }

}// end boost test case


void gen_blob(int idim, int nblob, dblAr2& coord, intAr2& ent2poi){

  ent2poi.allocate(nblob,idim+1);
  ent2poi.set_n(nblob);

  if(idim == 2){
    // Band of square triangles
    int npoin = nblob + 2; // Each new triangle is made by adding one point. 
    coord.allocate(npoin,3);
    coord.set_n(npoin);
    // Even nodes are (ipoin/2,0)
    // odd nodes are (ipoin/2,1)
    for(int ipoin = 0; ipoin < npoin; ipoin++){
      coord(ipoin,0) = ipoin/2;
      coord(ipoin,1) = ipoin%2;
    }

    for(int ielem = 0; ielem < nblob; ielem++){
      ent2poi(ielem,0) = ielem + 0;
      ent2poi(ielem,1) = ielem + 1;
      ent2poi(ielem,2) = ielem + 2;
    }
  }else{//endif idim == 2


    if(nblob == 1){
      coord.allocate(4,3);
      coord.set_n(4);
    }else{
      coord.allocate(2+nblob,3);
      coord.set_n(2+nblob);
    }

    coord(0,0) = 0;
    coord(0,1) = 0;
    coord(0,2) = -1;

    coord(1,0) = 0;
    coord(1,1) = 0;
    coord(1,2) = 1;

    if(nblob == 1){
      coord(2,0) = 0;
      coord(2,1) = 1;
      coord(2,2) = 0;

      coord(3,0) = 1;
      coord(3,1) = 0;
      coord(3,2) = 0;

      ent2poi(0,0) = 0;
      ent2poi(0,1) = 1;
      ent2poi(0,2) = 2;
      ent2poi(0,3) = 3;

      return;
    }

    const double pi = 3.141592653589793238462643383279502884;
    for(int ipoin = 2; ipoin < nblob+2; ipoin++){
      double angle = 2*pi*(ipoin-2) / (double) nblob;
      coord(ipoin,0) = cos(angle);
      coord(ipoin,1) = sin(angle);
      coord(ipoin,2) = 0;
    }

    for(int ielem = 0; ielem < nblob; ielem++){
      ent2poi(ielem,0) = 0;
      ent2poi(ielem,1) = 1;
      ent2poi(ielem,2) = ielem+2;
      ent2poi(ielem,3) = ielem < nblob - 1 ? ielem+3 : 2;
      double meas0 = getmeasentP1<3>(ent2poi[ielem], coord);
      METRIS_ENFORCE(meas0 > 1.0e-1);
    }
  }// endif idim == 3
}




void dump_init(int idim, std::string mshname, intAr2 &ent2poi, dblAr2 &coord){

  int64_t libIdx;
  libIdx = MetrisOpenMeshFile<GmfWrite>(mshname.c_str(), idim);

  // Write points
  int npoin = coord.get_n();
  intAr1 lpoic(npoin);
  lpoic.fill(0);
  GmfSetKwd(libIdx, GmfVertices, npoin);
  GmfSetBlock(libIdx, GmfVertices, 1, npoin, 0, NULL, NULL,
    GmfDoubleVec, idim, &coord(0,0), &coord[npoin-1][0],
    GmfInt            , &lpoic[0]      , &lpoic[npoin-1]);

  // Write elements
  int nelem = ent2poi.get_n();
  int nnode = idim+1;
  // We need to copy because the DIRECT tess has some aux values in there
  intAr2 ent2po2(nelem,nnode);
  intAr1 ent2ref(nelem);
  ent2ref.fill(1);
  for(int ielem = 0; ielem < nelem; ielem++){
    for(int ii = 0; ii < nnode ; ii++){
      ent2po2(ielem,ii) = ent2poi(ielem,ii)+1;
    }
  }
  int eKwd = idim == 2 ? libmeshb::faceKwds[1] : libmeshb::elemKwds[1];
  GmfSetKwd( libIdx, eKwd, nelem);
  GmfSetBlock(libIdx, eKwd, 1, nelem, 0, NULL, NULL,
    GmfIntVec, nnode, &ent2po2(0,0), &ent2po2(nelem-1,0),
    GmfInt   ,        &ent2ref[0]       , &ent2ref[nelem-1]);

  GmfCloseMesh( libIdx );
}





void dump_tess(DIBLOB_args& args, intAr2& ent2poi, dblAr2 &coord){

  std::string mshname = "tmp/tess"+std::to_string(args.niter)+".meshb";
  int idim = args.ent2pol.get_stride() - 3;
  METRIS_ENFORCE(idim == ent2poi.get_stride() - 1);

  int64_t libIdx;
  libIdx = MetrisOpenMeshFile<GmfWrite>(mshname.c_str(), idim);

  int nelem = args.ent2pol.get_n();
  int npoin = args.coorl.get_n();

  // Write points
  intAr1 lpoic(npoin);
  lpoic.fill(0);
  // Evaluate coords. The args.coorl are storing barycentric coordinates. 
  intAr1 lpoie(npoin);
  lpoie.fill(-1);
  dblAr2 coorp(npoin,idim); // physical coordinates

  int npoiw = 0;
  for(int ielem = 0; ielem < nelem; ielem++){
    int ieglo = args.ent2pol(ielem,idim+2);
    METRIS_ENFORCE(ieglo >= 0 && ieglo < ent2poi.get_n());

    for(int iver = 0; iver < idim + 1; iver++){
      int ipoin = args.ent2pol(ielem,iver);
      if(lpoie[ipoin] >= 0) continue;
      npoiw++;
      lpoie[ipoin] = ieglo;

      double bary[4];
      double sum = 0;
      for(int ii = 0; ii < idim; ii++){
        bary[ii] = args.coorl(ipoin,ii);
        sum += bary[ii];
      }
      bary[idim] = 1 - sum;

      if(idim == 2){
        eval2<2,1>(coord, ent2poi[ieglo], FEBasis::Lagrange, DifVar::None, DifVar::None,
                   bary, coorp[ipoin], NULL, NULL);
      }else{
        eval3<3,1>(coord, ent2poi[ieglo], FEBasis::Lagrange, DifVar::None, DifVar::None,
                   bary, coorp[ipoin], NULL, NULL);
      }
    }
  }
  printf("%d / %d points converted to physical coordinates\n",npoiw,npoin);


  GmfSetKwd(libIdx, GmfVertices, npoin);
  GmfSetBlock(libIdx, GmfVertices, 1, npoin, 0, NULL, NULL,
    GmfDoubleVec, idim, &coorp(0,0), &coorp[npoin-1][0],
    GmfInt            , &lpoic[0]  , &lpoic[npoin-1]);

  // Write elements
  int nnode = idim+1;
  printf("nnode = %d \n",nnode);
  // We need to copy because the DIRECT tess has some aux values in there
  intAr2 ent2po2(nelem,nnode);
  intAr1 ent2ref(nelem);
  for(int ielem = 0; ielem < nelem; ielem++){
    for(int ii = 0; ii < nnode ; ii++){
      ent2po2(ielem,ii) = args.ent2pol(ielem,ii)+1;
    }
    ent2ref[ielem] = args.ent2pol(ielem,idim+1); // write the level
  }
  int eKwd = idim == 2 ? libmeshb::faceKwds[1] : libmeshb::elemKwds[1];
  GmfSetKwd( libIdx, eKwd, nelem);
  GmfSetBlock(libIdx, eKwd, 1, nelem, 0, NULL, NULL,
    GmfIntVec, nnode, &ent2po2(0,0), &ent2po2(nelem-1,0),
    GmfInt   ,        &ent2ref[0]       , &ent2ref[nelem-1]);

  GmfCloseMesh( libIdx );
}

