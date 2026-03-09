//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE bench_eval_d

#include "common_setup.hxx"


namespace Metris{

typedef  MetricFieldAnalytical MFT;

// 1.0e-8 relative error -> 1.0e-6% utf::tolerance()
BOOST_AUTO_TEST_CASE(bench_eval_d) 
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
      CPRINTF1(" - Running as FEBasis {}\n",basname[(int)ibasis].c_str());

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
    
    
    
    
        //int ntar = msh.idim == 3 ? 1e6 : 1e7;
        int ntar = 1e6;
    
        CPRINTF1("-- Start benchmarks ideg = %d ilag = %d n run > %8e \n",ideg,(int)ibasis,(double)ntar);

        CPRINTF1("   - 1: No   dfld matrix\n");
        double ps1[2],ps2[2],ps3[2],ps4[2],ps5[2];
        for(int ijac = 0; ijac <= 1; ijac++){
          DifVar ideriv = ijac == 0 ? DifVar::None : DifVar::Bary;

          dum = 0;
          int nrep = (int) (ntar / nentt);
          
          double t0 = get_cpu_time();
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
          double t1 = get_cpu_time();
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
          double t2 = get_cpu_time();
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
          double t3 = get_cpu_time();

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

          double t4 = get_cpu_time();

          ps1[ijac] = (nentt*nrep*nsamp/(1e6 * (t1-t0)));
          ps2[ijac] = (nentt*nrep*nsamp/(1e6 * (t2-t1)));
          ps3[ijac] = (nentt*nrep*nsamp/(1e6 * (t3-t2)));
          ps4[ijac] = (nentt*nrep*nsamp/(1e6 * (t4-t3)));
        }

        if(ibasis == FEBasis::Lagrange){
          CPRINTF1("    (%1.0f) eval3/s = %6.1f/%6.1fM/s direct = %6.1f/%6.1fM/s Surreal = x\n",dum,ps1[0],ps1[1],ps2[0],ps2[1]);
        }else{
          CPRINTF1("    (%1.0f) eval3/s = %6.1f/%6.1fM/s direct = %6.1f/%6.1fM/s Surreal = %6.1f/%6.1fM/s SSimple = %6.1f/%6.1fM/s\n",dum,ps1[0],ps1[1],ps2[0],ps2[1],ps3[0],ps3[1],ps4[0],ps4[1]);
        }


        CPRINTF1("   - 2: with dfld matrix\n");
        
        for(int ijac = 0; ijac <= 1; ijac++){
          DifVar ideriv = ijac == 0 ? DifVar::None : DifVar::Bary;

          dum = 0;
          int nrep = (int) (ntar / nentt);
          
          double t0 = get_cpu_time();
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
          double t1 = get_cpu_time();
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
          double t2 = get_cpu_time();
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
          double t3 = get_cpu_time();
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
          double t4 = get_cpu_time();
          ps1[ijac] = (nentt*nrep*nsamp/(1e6 * (t1-t0)));
          ps2[ijac] = (nentt*nrep*nsamp/(1e6 * (t2-t1)));
          ps3[ijac] = (nentt*nrep*nsamp/(1e6 * (t3-t2)));
          ps4[ijac] = (nentt*nrep*nsamp/(1e6 * (t4-t3)));
          //printf("Debug t0, t1 %23.16e %23.16e diff %23.16e \n", t0, t1, t1-t0);
    
  //    writeMesh(std::string("ccoeff-in.1.") + std::to_string(ideg) + ".mesh",f.msh);
        }

        if(ibasis == FEBasis::Lagrange){
          CPRINTF1("    (%3.1f) eval3/s = %6.1f/%6.1fM/s direct = %6.1f/%6.1fM/s Surreal = x\n",
                   dum,ps1[0],ps1[1],ps2[0],ps2[1]);
        }else{
          CPRINTF1("    (%3.1f) eval3/s = %6.1f/%6.1fM/s direct = %6.1f/%6.1fM/s Surreal = %6.1f/%6.1fM/s SSimple = %6.1f/%6.1fM/s\n",
                   dum,ps1[0],ps1[1],ps2[0],ps2[1],ps3[0],ps3[1],ps4[0],ps4[1]);
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