//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_eval_d

#include <boost/test/included/unit_test.hpp> 
#include "common_setup.hxx"

#include "../src/low_eval.hxx"
#include "../src/low_eval_d.hxx"
#include "../src/low_evalS.hxx"
#include "../src/low_eval_d_SurrealS.hxx"

#include "../src/utils/CT_loop.hxx"
#include "../src/utils/mprintf.hxx"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
#include "../src/Mesh/Mesh.hxx"

#include "../src/low_eval_d.hxx"

#include <random>


namespace Metris{

typedef  MetricFieldAnalytical MFT;

// 1.0e-8 relative error -> 1.0e-6% utf::tolerance()
BOOST_AUTO_TEST_CASE(test_eval_d) 
{
//METRIS_MAX_DEG


  constexpr int nbase = MAX((int)FEBasis::Bezier, (int)FEBasis::Lagrange) + 1;
  std::string basname[nbase];
  basname[(int)FEBasis::Bezier]   = "Bez";
  basname[(int)FEBasis::Lagrange] = "Lag";

  std::vector<std::string> meshes = {
                                      METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
                                     ,METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
                                     ,METRIS_CASES_DIR "/unit/2D/square/iso.p1.100k"
                                     ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.500"
                                     ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.5k"
                                     ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k -t 2"
                                     ,METRIS_CASES_DIR "/unit/3D/cube/curved.p2.2k"
                                     #if METRIS_MAX_DEG >= 3
                                     ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.500 -t 3"
                                     ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.5k -t 3"
                                     ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k -t 3"
                                     ,METRIS_CASES_DIR "/unit/3D/cube/curved.p2.2k -t 3"
                                     #if METRIS_MAX_DEG >= 4
                                     ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.500 -t 4"
                                     ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.5k -t 4"
                                     ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k -t 4"
                                     ,METRIS_CASES_DIR "/unit/3D/cube/curved.p2.2k -t 4"
                                     #endif
                                     #endif
                                   };


  double tol = 1.0e-11;

  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(0);

  double dum;

  MetrisParameters prtparam;
  prtparam.iverb = 5;
  prtparam.ivdepth = 5;

  for(auto s : meshes)
  {

    try{

    cargHandler arg("-in " + s + "  -anamet 1 -verb 0");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);

    INCVDEPTH((&prtparam));
    CPRINTF1("-- Mesh : %s\n",s.c_str());

    msh.cleanup();

    for(FEBasis ibasis : {FEBasis::Lagrange, FEBasis::Bezier}){
      INCVDEPTH((&prtparam));
      msh.setBasis(ibasis);
      CPRINTF1(" - Running as FEBasis %s\n",basname[(int)ibasis].c_str());

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
    
    
        CPRINTF1(" - 1. Test evaluation and jmat (not diff)\n");
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
    

        CPRINTF1(" - 2. Test eval derivatives\n");
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
    
    
      }}CT_FOR1(idim);
      }}CT_FOR1(ideg);
    }

    }catch(const MetrisExcept& e){
      printf("## Failed to load case, possibly missing file\n");
    }

  }// end for mesh

}



} // end namespace