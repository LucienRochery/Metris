//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include "common_setup.hxx"

#include "../src/low_eval.hxx"
#include "../src/low_eval_d.hxx"
#include "../src/low_evalS.hxx"
#include "../src/low_eval_d_SurrealS.hxx"

#include "../src/utils/CT_loop.hxx"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
#include "../src/Mesh/Mesh.hxx"

#include "../src/low_eval_d.hxx"

#include <random>
#include <boost/hana.hpp> 
namespace hana = boost::hana;
using namespace hana::literals;

namespace utf = boost::unit_test;


namespace Metris{

typedef  MetricFieldAnalytical MFT;

// 1.0e-8 relative error -> 1.0e-6% utf::tolerance()
BOOST_AUTO_TEST_CASE(test_eval3_d) 
{
//METRIS_MAX_DEG


  constexpr int nbase = MAX((int)FEBasis::Bezier, (int)FEBasis::Lagrange) + 1;
  std::string basname[nbase];
  basname[(int)FEBasis::Bezier]   = "Bez";
  basname[(int)FEBasis::Lagrange] = "Lag";

  std::vector<std::string> meshes = {
                                      "../cases/1200_p1.meshb"
                                     ,"../cases/2D/square.p1.10"
                                     ,"../cases/2D/square.circmet.5k.curved.meshb"
                                     ,"../cases/2D/square.circmet.50.curved.meshb"
                                     ,"../cases/1200_p2.meshb"
                                     ,"../cases/curved_p2.meshb"
                                     ,"../cases/2D/square.p2.100k"
                                     #if METRIS_MAX_DEG >= 3
                                     ,"../cases/2D/square.circmet.50.curved.meshb -t 3"
                                     ,"../cases/2D/square.circmet.5k.curved.meshb -t 3"
                                     ,"../cases/1200_p3.meshb"
                                     ,"../cases/curved_p3.meshb"
                                     #if METRIS_MAX_DEG >= 4
                                     ,"../cases/1200_p4.meshb"
                                     ,"../cases/curved_p4.meshb"
                                     #endif
                                     #endif
                                   };


  double tol = 1.0e-11;

  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(0);

  double dum;

  for(auto s : meshes)
  {
    std::cout<<"\n\n\n -- CASE : Mesh "<<s<<std::endl;

    try{

    cargHandler arg("-in " + s + "  -anamet 1");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);

    msh.cleanup();

    for(FEBasis ibasis : {FEBasis::Lagrange, FEBasis::Bezier}){
      msh.setBasis(ibasis);
      printf("\n\n-- Running as FEBasis %s\n",basname[(int)ibasis].c_str());

      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      CT_FOR0_INC(2,3,idim){if(idim == msh.idim){

        constexpr int gdim = idim;
        constexpr int tdim = idim;

        constexpr int njmat = gdim * tdim;
        constexpr int nnode = getnnode(tdim, ideg);

        const int nentt = msh.nentt(tdim);
        const intAr2 &ent2poi = msh.ent2poi(tdim);
    
        constexpr auto evalf = idim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;
        constexpr auto evalf_d0 = idim == 2 ? eval2_d<gdim,ideg,0,gdim> : eval3_d<gdim,ideg,0,gdim>;


        double eval[idim], eva2[idim], jmat[idim*idim], jmat2[idim*idim], deval[idim][idim],djmat[idim][idim*idim];
        double evals[idim+1][idim], jmats[idim+1][idim*idim];
        double dx = 2.0;
        double dfld0[idim][idim];
        double dfld1[idim][idim];
        for(int i = 0; i < idim ;i++)
          for(int j = 0; j < idim ;j++)
            dfld0[i][j] = (i==j);


        for(int i = 0; i < idim ;i++)
          for(int j = 0; j < idim ;j++)
            dfld1[i][j] = unif(rng);


        constexpr int nsamp = nnode;
        double bary[nsamp][idim+1];
        constexpr auto ordent = ORDELT(tdim);
        for(int ii = 0 ;ii < nsamp; ii++)
          for(int jj = 0; jj < tdim + 1; jj++) 
            bary[ii][jj] = ordent[ideg][ii][jj] / ((double) ideg);

        dum = 0.0;
    
    
        printf("\n-- 1. Test evaluation and jmat (not diff)\n\n");
        for(int ientt = 0; ientt < nentt; ientt++){
          if(isdeadent(ientt,ent2poi)) continue;
          constexpr int ivar = 0;

          double epsent = getepsent<gdim>(msh,tdim,ientt);
    
          // Test eval3_direct
          double errL2 = 0;
          for(int isamp = 0; isamp < nsamp; isamp++){
            evalf(msh.coord,ent2poi[ientt],
                  ibasis, DifVar::Bary, DifVar::None, 
                  bary[isamp],eval,jmat,NULL);
            eval_d_direct<gdim,tdim,ideg,ivar>(msh.coord,ent2poi[ientt],
                                              ibasis, DifVar::Bary, DifVar::None, 
                                              bary[isamp],eva2,jmat2,NULL,deval[0],djmat[0],NULL);
            double errloc_eval = geterrl2<gdim>(eval,eva2);
            double errloc_jmat = geterrl2<njmat>(jmat,jmat2);
            BOOST_REQUIRE(errloc_eval <= tol*tol*epsent*epsent);
            BOOST_REQUIRE(errloc_jmat <= tol*tol*epsent*epsent);
            errL2 += errloc_eval + errloc_jmat;
          }
          BOOST_REQUIRE(errL2 < tol*tol*12*nsamp*epsent*epsent);
    
          // Test eval3_SurrealS
          if(ibasis == FEBasis::Bezier){
            errL2 = 0;
            double errL2_2 = 0;
            double errL2_3 = 0;
            for(int isamp = 0; isamp < nsamp; isamp++){
              evalf(msh.coord,ent2poi[ientt],
                    ibasis, DifVar::Bary, DifVar::None, 
                    bary[isamp],eval,jmat,NULL);
              eval_d_SurrealS<idim,idim,ideg,ivar>(msh.coord,ent2poi[ientt],
                                                  ibasis, DifVar::Bary, DifVar::None, 
                                                  bary[isamp],eva2,jmat2,NULL,deval[0],djmat[0],NULL);
              double errloc_eval = geterrl2<gdim>(eval,eva2);
              double errloc_jmat = geterrl2<njmat>(jmat,jmat2);
              BOOST_REQUIRE(errloc_eval <= tol*tol*epsent*epsent);
              BOOST_REQUIRE(errloc_jmat <= tol*tol*epsent*epsent);
              errL2 += errloc_eval + errloc_jmat;


              eval_d_SurrealS_simple<idim,idim,ideg,idim>(msh.coord,ent2poi[ientt],
                                                  ibasis, DifVar::Bary, DifVar::None, 
                                                  bary[isamp],ivar,eva2,jmat2,
                                                  NULL,deval[0],djmat[0],NULL);
              errloc_eval = geterrl2<gdim>(eval,eva2);
              errloc_jmat = geterrl2<njmat>(jmat,jmat2);
              BOOST_REQUIRE(errloc_eval <= tol*tol*epsent*epsent);
              BOOST_REQUIRE(errloc_jmat <= tol*tol*epsent*epsent);
              errL2_2 += errloc_eval + errloc_jmat;
              


              eval_d_SurrealS_bcast<idim,idim,ideg,ivar,1>(msh.coord,ent2poi[ientt],
                                                  ibasis, DifVar::Bary, DifVar::None, 
                                                  bary[isamp],eva2,jmat2,
                                                  NULL,deval[0],djmat[0],NULL);
              errloc_eval = geterrl2<gdim>(eval,eva2);
              errloc_jmat = geterrl2<njmat>(jmat,jmat2);
              BOOST_REQUIRE(errloc_eval <= tol*tol*epsent*epsent);
              BOOST_REQUIRE(errloc_jmat <= tol*tol*epsent*epsent);
              errL2_3 += errloc_eval + errloc_jmat;


            }
            BOOST_REQUIRE(errL2 < tol*tol*12*nsamp*epsent*epsent);
            BOOST_REQUIRE(errL2_2 < tol*tol*12*nsamp*epsent*epsent);
            BOOST_REQUIRE(errL2_3 < tol*tol*12*nsamp*epsent*epsent);
          }
        }// for ientt
    

        printf("\n-- 2. Test eval derivatives\n\n");
        // Test derivatives
        // Linear function: exactly given by diffs
        for(int ientt = 0; ientt < nentt; ientt++){
          if(isdeadent(ientt,ent2poi)) continue;
          double deval_disc[idim][idim], djmat_disc[idim][idim*idim];
      
          CT_FOR0_EXC(0,nnode,ivar){

            int ipoin = ent2poi(ientt,ivar);
            double coor0[gdim];
            for(int ii = 0; ii < gdim; ii++) coor0[ii] = msh.coord(ipoin,ii);
    
            for(int isamp = 0; isamp < nsamp; isamp++){
              eval_d_direct<gdim,tdim,ideg,ivar>(msh.coord,ent2poi[ientt],
                                           ibasis, DifVar::Bary, DifVar::None, 
                                           bary[isamp],evals[0],jmats[0],
                                           NULL,deval[0],djmat[0],NULL);
              for(int ii = 0 ; ii < gdim; ii++){
                for(int jj = 0; jj < idim; jj++) msh.coord(ipoin,jj) = coor0[jj];
                msh.coord(ipoin,ii) = coor0[ii] + dx;
                evalf(msh.coord,ent2poi[ientt],
                      ibasis, DifVar::Bary, DifVar::None, 
                      bary[isamp],evals[ii+1],jmats[ii+1],NULL);
              }
              for(int ii = 0; ii < gdim; ii++) msh.coord(ipoin,ii) = coor0[ii];
              for(int ii = 0; ii < idim; ii++)
                for(int jj = 0; jj < idim; jj++)
                  deval_disc[jj][ii] = (evals[jj+1][ii] - evals[0][ii])/dx;
              
              for(int ii = 0; ii < idim*idim; ii++)
                for(int jj = 0; jj < idim; jj++)
                  djmat_disc[jj][ii] = (jmats[jj+1][ii] - jmats[0][ii])/dx;
              

              double errdeval = sqrt(geterrl2<idim*idim>(deval_disc[0],deval[0]));
              double errdjmat = sqrt(geterrl2<idim*njmat>(djmat_disc[0],djmat[0]));

              BOOST_REQUIRE(errdeval < tol);
              BOOST_REQUIRE(errdjmat < tol);

              eval_d_direct<gdim,tdim,ideg,ivar>(msh.coord,ent2poi[ientt],
                                           ibasis, DifVar::Bary, DifVar::None, 
                                           bary[isamp],evals[0],jmats[0],
                                           NULL,deval[0],djmat[0],NULL,dfld0[0]);
              errdeval = sqrt(geterrl2<idim*idim>(deval_disc[0],deval[0]));
              errdjmat = sqrt(geterrl2<idim*njmat>(djmat_disc[0],djmat[0]));
              BOOST_REQUIRE(errdeval < tol);
              BOOST_REQUIRE(errdjmat < tol);
    

              if(ibasis != FEBasis::Bezier) continue;

              for(const double* dfld : {dfld0[0], (double*)NULL}){

                eval_d_SurrealS<idim,idim,ideg,ivar>(msh.coord,ent2poi[ientt],
                                               ibasis, DifVar::Bary, DifVar::None, 
                                               bary[isamp],evals[0],jmats[0],
                                               NULL,deval[0],djmat[0],NULL, dfld);
                errdeval = sqrt(geterrl2<idim*idim>(deval_disc[0],deval[0]));
                errdjmat = sqrt(geterrl2<idim*njmat>(djmat_disc[0],djmat[0]));
                BOOST_REQUIRE(errdeval < tol);
                BOOST_REQUIRE(errdjmat < tol);


                eval_d_SurrealS_simple<idim,idim,ideg,idim>(msh.coord,ent2poi[ientt],
                                               ibasis, DifVar::Bary, DifVar::None, 
                                               bary[isamp], ivar, evals[0],jmats[0],
                                               NULL,deval[0],djmat[0],NULL, dfld);
                errdeval = sqrt(geterrl2<idim*idim>(deval_disc[0],deval[0]));
                errdjmat = sqrt(geterrl2<idim*njmat>(djmat_disc[0],djmat[0]));
                BOOST_REQUIRE(errdeval < tol);
                BOOST_TEST(errdjmat < tol);


                //eval_d_SurrealS_bcast<idim,idim,ideg,ivar,1>(msh.coord,ent2poi[ientt],
                //                               ibasis, DifVar::Bary, DifVar::None, 
                //                               bary[isamp], evals[0],jmats[0],NULL,
                //                               deval[0],djmat[0],NULL, dfld);
                //errdeval = sqrt(geterrl2<idim*idim>(deval_disc[0],deval[0]));
                //errdjmat = sqrt(geterrl2<idim*njmat>(djmat_disc[0],djmat[0]));
                //BOOST_REQUIRE(errdeval < tol);
                //BOOST_REQUIRE(errdjmat < tol);

              }

                
            }// for isamp
    
          }CT_FOR1(ivar);
        } 
    
    
    
    
    
        //int ntar = msh.idim == 3 ? 1e6 : 1e7;
        int ntar = 1e6;
    
        printf("-- Start benchmarks ideg = %d ilag = %d n run > %8e \n",ideg,ibasis,(double)ntar);
        printf("-1 No   dfld matrix\n");
        double ps1[2],ps2[2],ps3[2],ps4[2],ps5[2];
        for(int ijac = 0; ijac <= 1; ijac++){
          DifVar ideriv = ijac == 0 ? DifVar::None : DifVar::Bary;

          dum = 0;
          int nrep = (int) (ntar / nentt);
          
          double t0 = get_wall_time();
          for(int irep = 0; irep < nrep; irep++)
            for(int ientt = 0; ientt < nentt; ientt++){
              if(isdeadent(ientt,ent2poi)) continue;
              for(int isamp = 0; isamp < nsamp; isamp++){
                evalf(msh.coord,ent2poi[ientt],
                      ibasis, ideriv, DifVar::None, 
                      bary[isamp],eval,jmat,NULL);
                dum += eval[0]+eval[1];
                if(abs(dum) > 1) dum = 1.0 / dum;
              }
            }
          double t1 = get_wall_time();
          for(int irep = 0; irep < nrep; irep++)
            for(int ientt = 0; ientt < nentt; ientt++){
              if(isdeadent(ientt,ent2poi)) continue;
              constexpr int ivar = 0;
              for(int isamp = 0; isamp < nsamp; isamp++){
                eval_d_direct<gdim,tdim,ideg,ivar>(msh.coord,ent2poi[ientt],
                                             ibasis, ideriv, DifVar::None, 
                                             bary[isamp],eval,jmat,NULL,deval[0],djmat[0],NULL);
                dum += eval[0]+eval[1];
                if(abs(dum) > 1) dum = 1.0 / dum;
              }
            }
          double t2 = get_wall_time();
          if(ibasis == FEBasis::Bezier){
            for(int irep = 0; irep < nrep; irep++)
              for(int ientt = 0; ientt < nentt; ientt++){
                if(isdeadent(ientt,ent2poi)) continue;
                constexpr int ivar = 0;
                for(int isamp = 0; isamp < nsamp; isamp++){
                  eval_d_SurrealS<idim,idim,ideg,ivar>(msh.coord,ent2poi[ientt],
                                                 ibasis, ideriv, DifVar::None, 
                                                 bary[isamp],eval,jmat,NULL,deval[0],djmat[0],NULL);
                  dum += eval[0]+eval[1];
                  if(abs(dum) > 1) dum = 1.0 / dum;
                }
              }
          }
          double t3 = get_wall_time();

          if(ibasis == FEBasis::Bezier){
            for(int irep = 0; irep < nrep; irep++)
              for(int ientt = 0; ientt < nentt; ientt++){
                if(isdeadent(ientt,ent2poi)) continue;
                constexpr int ivar = 0;
                for(int isamp = 0; isamp < nsamp; isamp++){
                  eval_d_SurrealS_simple<idim,idim,ideg,idim>(msh.coord,ent2poi[ientt],
                                                 ibasis, ideriv, DifVar::None, 
                                                 bary[isamp],ivar,eval,jmat,NULL,deval[0],djmat[0],NULL);
                  dum += eval[0]+eval[1];
                  if(abs(dum) > 1) dum = 1.0 / dum;
                }
              }
          }

          double t4 = get_wall_time();


          //if(ibasis == FEBasis::Bezier){
          //  for(int irep = 0; irep < nrep; irep++)
          //    for(int ientt = 0; ientt < nentt; ientt++){
          //      if(isdeadent(ientt,ent2poi)) continue;
          //      constexpr int ivar = 0;
          //      for(int isamp = 0; isamp < nsamp; isamp++){
          //        eval_d_SurrealS_bcast<idim,idim,ideg,ivar,1>(msh.coord,ent2poi[ientt],
          //                                       ibasis, ideriv, DifVar::None, 
          //                                       bary[isamp],eval,jmat,NULL,deval[0],djmat[0],NULL);
          //        dum += eval[0]+eval[1];
          //        if(abs(dum) > 1) dum = 1.0 / dum;
          //      }
          //    }
          //}

          //double t5 = get_wall_time();

          ps1[ijac] = (nentt*nrep*nsamp/(1e6 * (t1-t0)));
          ps2[ijac] = (nentt*nrep*nsamp/(1e6 * (t2-t1)));
          ps3[ijac] = (nentt*nrep*nsamp/(1e6 * (t3-t2)));
          ps4[ijac] = (nentt*nrep*nsamp/(1e6 * (t4-t3)));
          //ps5[ijac] = (nentt*nrep*nsamp/(1e6 * (t5-t4)));
          //printf("Debug t0, t1 %23.16e %23.16e diff %23.16e \n", t0, t1, t1-t0);
    
  //    writeMesh(std::string("ccoeff-in.1.") + std::to_string(ideg) + ".mesh",f.msh);
        }

        if(ibasis == FEBasis::Lagrange){
          printf("(%1.0f) eval3/s = %6.1f/%6.1fM/s direct = %6.1f/%6.1fM/s Surreal = x\n",dum,ps1[0],ps1[1],ps2[0],ps2[1]);
        }else{
          //printf("(%1.0f) eval3/s = %6.1f/%6.1fM/s direct = %6.1f/%6.1fM/s Surreal = %6.1f/%6.1fM/s SSimple = %6.1f/%6.1fM/s Sbcast = %6.1f/%6.1fM/s\n",dum,ps1[0],ps1[1],ps2[0],ps2[1],ps3[0],ps3[1],ps4[0],ps4[1],ps5[0],ps5[1]);
          printf("(%1.0f) eval3/s = %6.1f/%6.1fM/s direct = %6.1f/%6.1fM/s Surreal = %6.1f/%6.1fM/s SSimple = %6.1f/%6.1fM/s\n",dum,ps1[0],ps1[1],ps2[0],ps2[1],ps3[0],ps3[1],ps4[0],ps4[1]);
        }


        printf("-2 With dfld matrix\n");
        
        for(int ijac = 0; ijac <= 1; ijac++){
          DifVar ideriv = ijac == 0 ? DifVar::None : DifVar::Bary;

          dum = 0;
          int nrep = (int) (ntar / nentt);
          
          double t0 = get_wall_time();
          for(int irep = 0; irep < nrep; irep++)
            for(int ientt = 0; ientt < nentt; ientt++){
              if(isdeadent(ientt,ent2poi)) continue;
              for(int isamp = 0; isamp < nsamp; isamp++){
                evalf(msh.coord,ent2poi[ientt],
                      ibasis, ideriv, DifVar::None, 
                      bary[isamp],eval,jmat,NULL);
                dum += eval[0]+eval[1];
                if(abs(dum) > 1) dum = 1.0 / dum;
              }
            }
          double t1 = get_wall_time();
          for(int irep = 0; irep < nrep; irep++)
            for(int ientt = 0; ientt < nentt; ientt++){
              if(isdeadent(ientt,ent2poi)) continue;
              constexpr int ivar = 0;
              for(int isamp = 0; isamp < nsamp; isamp++){
                eval_d_direct<gdim,tdim,ideg,ivar>(msh.coord,ent2poi[ientt],
                                             ibasis, ideriv, DifVar::None, 
                                             bary[isamp],eval,jmat,NULL,deval[0],djmat[0],NULL,dfld1[0]);
                dum += eval[0]+eval[1];
                if(abs(dum) > 1) dum = 1.0 / dum;
              }
            }
          double t2 = get_wall_time();
          if(ibasis == FEBasis::Bezier){
            for(int irep = 0; irep < nrep; irep++)
              for(int ientt = 0; ientt < nentt; ientt++){
                if(isdeadent(ientt,ent2poi)) continue;
                constexpr int ivar = 0;
                for(int isamp = 0; isamp < nsamp; isamp++){
                  eval_d_SurrealS<idim,idim,ideg,ivar>(msh.coord,ent2poi[ientt],
                                                 ibasis, ideriv, DifVar::None, 
                                                 bary[isamp],eval,jmat,NULL,deval[0],djmat[0],NULL,dfld1[0]);
                  dum += eval[0]+eval[1];
                  if(abs(dum) > 1) dum = 1.0 / dum;
                }
              }
          }
          double t3 = get_wall_time();
          if(ibasis == FEBasis::Bezier){
            for(int irep = 0; irep < nrep; irep++)
              for(int ientt = 0; ientt < nentt; ientt++){
                if(isdeadent(ientt,ent2poi)) continue;
                constexpr int ivar = 0;
                for(int isamp = 0; isamp < nsamp; isamp++){
                  eval_d_SurrealS_simple<idim,idim,ideg,idim>(msh.coord,ent2poi[ientt],
                                                 ibasis, ideriv, DifVar::None, 
                                                 bary[isamp],ivar,eval,jmat,NULL,deval[0],djmat[0],NULL,dfld1[0]);
                  dum += eval[0]+eval[1];
                  if(abs(dum) > 1) dum = 1.0 / dum;
                }
              }
          }
          double t4 = get_wall_time();
          ps1[ijac] = (nentt*nrep*nsamp/(1e6 * (t1-t0)));
          ps2[ijac] = (nentt*nrep*nsamp/(1e6 * (t2-t1)));
          ps3[ijac] = (nentt*nrep*nsamp/(1e6 * (t3-t2)));
          ps4[ijac] = (nentt*nrep*nsamp/(1e6 * (t4-t3)));
          //printf("Debug t0, t1 %23.16e %23.16e diff %23.16e \n", t0, t1, t1-t0);
    
  //    writeMesh(std::string("ccoeff-in.1.") + std::to_string(ideg) + ".mesh",f.msh);
        }

        if(ibasis == FEBasis::Lagrange){
          printf("(%3.1f) eval3/s = %6.1f/%6.1fM/s direct = %6.1f/%6.1fM/s Surreal = x\n",dum,ps1[0],ps1[1],ps2[0],ps2[1]);
        }else{
          //printf("(%3.1f) eval3/s = %6.1f/%6.1fM/s direct = %6.1f/%6.1fM/s Surreal = %6.1f/%6.1fM/s\n",dum,ps1[0],ps1[1],ps2[0],ps2[1],ps3[0],ps3[1]);
          printf("(%1.0f) eval3/s = %6.1f/%6.1fM/s direct = %6.1f/%6.1fM/s Surreal = %6.1f/%6.1fM/s SSimple = %6.1f/%6.1fM/s\n",dum,ps1[0],ps1[1],ps2[0],ps2[1],ps3[0],ps3[1],ps4[0],ps4[1]);
        }
      }}CT_FOR1(idim);
      }}CT_FOR1(ideg);
    }

    }catch(const MetrisExcept& e){
      printf("## Failed to load case, possibly missing file\n");
    }

  }// end for mesh

}



} // end namespace