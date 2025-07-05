//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include <bunit/common_setup.hxx>

#include <boost/timer/progress_display.hpp>

#include "../src/ho_constants.hxx"
#include "../src/utils/aux_misc.hxx"
#include "../src/quality/low_metqua.hxx"
#include "../SANS/Surreal/SurrealS.h"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
#include "../src/low_geo.hxx"

#include <boost/hana.hpp> 
namespace hana = boost::hana;
using namespace hana::literals;

namespace utf = boost::unit_test;

using namespace Metris;

typedef MetricFieldAnalytical MFT;

// -- Test metqua3_xi_d 
// Constant metric fields should yield reliable derivatives in all cases 
// In non constant metric fields, derivatives only defined for DoFs in back element
// interiors... 

BOOST_AUTO_TEST_CASE(test_eval3) 
{

  // bool is whether straight
  std::vector<std::pair<std::string,bool>> meshes = {
     {METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k",true}
    ,{METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k -t 2",true}
    ,{METRIS_CASES_DIR "/unit/3D/cube/curved.p2.2k",false}
    ,{METRIS_CASES_DIR "/unit/2D/square/circmet.p2.500",false}
    #if METRIS_MAX_DEG >= 3
    ,{METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k -t 3",true}
    ,{METRIS_CASES_DIR "/unit/3D/cube/curved.p2.2k -t 3",false}
    #if METRIS_MAX_DEG >= 4
    ,{METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k -t 4",true}
    ,{METRIS_CASES_DIR "/unit/3D/cube/curved.p2.2k -t 4",false}
    #endif
    #endif
    }; 

    const double tol   = 1.0e-11;

    constexpr int nbase = MAX((int)FEBasis::Bezier, (int)FEBasis::Lagrange) + 1;
    std::string basname[nbase];
    basname[(int)FEBasis::Bezier]   = "Bez";
    basname[(int)FEBasis::Lagrange] = "Lag";


    const int mdx = 256;
    const double dx0   = 1.0e-1;
    const double dx1   = 1.0e-8;
    const double qdx   = 10.0;
    const double minsl = 0.5;
    int ndx = 0;
    for(double dx = dx0; dx > dx1; dx /= qdx) ndx++;
    
    if(ndx > mdx) METRIS_THROW_MSG(SMemExcept(),"Increase mdx")
    double err3dx[nbase][mdx], err6dx[nbase][mdx], logdx[mdx];

    for(auto testcase : meshes)
    { 
      std::string s = testcase.first;
      bool istr8    = testcase.second;

      try{
        cargHandler arg("-in " + s + "  -anamet 1");
        MetrisRunner run(arg.c, arg.v);
        Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);


        std::cout<<"\n\n------------------------------------------------\n";
        std::cout<<"Mesh "<<s<<"\n";
        std::cout<<"------------------------------------------------\n";


        CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
          CT_FOR0_INC(2,3,idim){if(idim == msh.idim){
            constexpr int nvar = idim;
            constexpr int nhess = (idim * (idim + 1)) / 2;

            printf("------- ideg = %d \n",ideg);

            constexpr auto ordent = ORDELT(idim);
            constexpr int nnode = getnnode(idim,ideg);
            constexpr auto evalf  =  idim == 2 ? eval2<idim,ideg> : eval3<idim,ideg>;
            constexpr auto evalP1 =  idim == 2 ? eval2<idim,1> : eval3<idim,1>;
            constexpr auto evalM  =  idim == 2 ? eval2<nhess,1> : eval3<nhess,1>;

            const intAr2 &ent2poi = msh.ent2poi(idim);
            const int nentt = msh.nentt(idim);

            const int nsamp = nnode;
            double bar0[nsamp][idim+1];
            for(int ii = 0 ; ii < nsamp ; ii++){
              for(int jj = 0; jj < idim + 1; jj++) 
                bar0[ii][jj] = ordent[ideg][ii][jj] / ((double) ideg);
            }

            double coob[idim],cool[idim],metl[nhess],metb[nhess];
            double jmat3b[idim*idim],jmat3l[idim*idim],jmat6l[idim*nhess],jmat6b[idim*nhess];
            double jmat3[idim*idim], jmat3_disc[idim][idim];
            double jmat6[nhess*idim], jmat6_disc[idim][nhess];
            double eva3s[1+idim][idim], eva6s[1+idim][nhess];
            double jmat3s[1+idim][idim*idim], jmat6s[1+idim][idim*nhess];
            double hmat3[nhess*idim], hmat6[nhess*nhess];
            double hmat3_disc[nhess][idim];
            double hmat6_disc[nhess][nhess];

        // START body

        // Test values
        // Linear function: exactly given by diffs

            if(istr8){
          // NB: although meshes are straight, this doesn't mean metrics are "straight"
          // i.e. the metric field itself is not P1 due to how it is constructed at the 
          // edges by averaging with neighbouring elements. 
              printf("-- Start value tests.\n");
              for(int ielem = 0; ielem < msh.nelem; ielem++){
                for(int isamp = 0; isamp < nsamp; isamp++){
                  evalf(msh.coord,ent2poi[ielem],FEBasis::Bezier  ,DifVar::Bary,DifVar::None,bar0[isamp],coob,jmat3b,NULL);

                  evalf(msh.coord,ent2poi[ielem],FEBasis::Lagrange,DifVar::Bary,DifVar::None,bar0[isamp],cool,jmat3l,NULL);

                  double errc = sqrt(geterrl2<3>(coob,cool));
                  double nrmc = sqrt(getnrml2<3>(coob));
                  double errj3 = sqrt(geterrl2<9>(jmat3b,jmat3l));
                  double nrmj3 = sqrt(getnrml2<9>(jmat3b));

                  BOOST_REQUIRE_MESSAGE(errc < tol*nrmc,
                    "Large coord eval error:"<<errc/nrmc<<"\n");
                  BOOST_REQUIRE_MESSAGE(errj3 < tol*nrmj3,
                    "Large coord jmat error:"<<errj3/nrmj3<<"\n");
                }
              }   
              printf("-- Value tests passed.\n");

            }


            printf("-- Start diff tests.\n");
            for(int ielem = 0; ielem < msh.nelem; ielem++){
              for(int isamp = 0; isamp < nsamp; isamp++){

                for(auto ibasis : {FEBasis::Lagrange, FEBasis::Bezier}){
                  evalf(msh.coord,ent2poi[ielem],ibasis,DifVar::Bary,DifVar::None,
                    bar0[isamp],eva3s[0],jmat3,NULL);
                  evalM(msh.met.rfld,ent2poi[ielem],ibasis,DifVar::Bary,DifVar::None,
                    bar0[isamp],eva6s[0],jmat6,NULL);

                  int idx = 0;
                  for(double dx = dx0; dx > dx1; dx /= qdx){
                    for(int ii = 0; ii < idim; ii++){
                      double bary[idim+1];
                      for(int jj = 0; jj < idim + 1; jj++) bary[jj] = bar0[isamp][jj];
                        bary[ii + 1] += dx;
                      bary[0] -= dx;
                      evalf(msh.coord,ent2poi[ielem],ibasis,DifVar::None,DifVar::None,
                        bary,eva3s[ii+1],NULL,NULL);
                      evalM(msh.met.rfld  ,ent2poi[ielem],ibasis,DifVar::None,DifVar::None,
                        bary,eva6s[ii+1],NULL,NULL);
                    }
                    for(int i = 0; i < idim; i++){
                      for(int jj = 0; jj < idim; jj++)
                        jmat3_disc[jj][i] = (eva3s[jj+1][i] - eva3s[0][i])/dx;
                    }
                    for(int i = 0; i < nhess; i++){
                      for(int jj = 0; jj < idim; jj++)
                        jmat6_disc[jj][i] = (eva6s[jj+1][i] - eva6s[0][i])/dx;
                    }
                    double err3 = sqrt(geterrl2<idim*idim >(jmat3_disc[0],jmat3));
                    double err6 = sqrt(geterrl2<idim*nhess>(jmat6_disc[0],jmat6));
                    err3dx[(int)ibasis][idx] = MAX(err3,1.0e-32);
                    err6dx[(int)ibasis][idx] = MAX(err6,1.0e-32);
                    logdx[idx]  = log(dx);
                    idx++;
                  }
                }
                if(!istr8){
                  double minerr3[2], minerr6[2];
                  for(auto ibasis : {FEBasis::Lagrange, FEBasis::Bezier}){
                    double minerr3 = 1.0e32;
                    double minerr6 = 1.0e32;
                    for(int idx = 0; idx < ndx; idx++){
                      minerr3 = MIN(minerr3, err3dx[(int)ibasis][idx]);
                      minerr6 = MIN(minerr6, err6dx[(int)ibasis][idx]);
                      err3dx[(int)ibasis][idx] = log(err3dx[(int)ibasis][idx]);
                      err6dx[(int)ibasis][idx] = log(err6dx[(int)ibasis][idx]);
                    }

                    if(minerr3 >= tol){
                      double slop3 = linearRegression(ndx,logdx,err3dx[(int)ibasis]);
                      BOOST_CHECK_MESSAGE(minsl < slop3,
                      " Coord slope ("<<basname[(int)ibasis]<<") "<<slop3<<" under minimum "<<minsl<<"\n");
                    }
                    if(minerr6 >= tol){
                      double slop6 = linearRegression(ndx,logdx,err6dx[(int)ibasis]);
                      BOOST_CHECK_MESSAGE(minsl < slop6,
                      " Metric slope ("<<basname[(int)ibasis]<<") "<<slop6<<" under minimum "<<minsl<<"\n");
                    }
                  }
                  for(auto ibasis : {FEBasis::Lagrange, FEBasis::Bezier}){
                    for(int idx = 0; idx < ndx; idx++){
                      err3dx[(int)ibasis][idx] = log(err3dx[(int)ibasis][idx]);
                      err6dx[(int)ibasis][idx] = log(err6dx[(int)ibasis][idx]);
                    }
                  }
                }else{
                  BOOST_CHECK_MESSAGE(err3dx[(int)FEBasis::Bezier][0] < tol,
                    "Large err3 "<< err3dx[(int)FEBasis::Bezier][0]<<" Bézier ideg = "<<ideg<<"\n");
                  BOOST_CHECK_MESSAGE(err3dx[(int)FEBasis::Lagrange][0] < tol,
                    "Large err3 "<< err3dx[(int)FEBasis::Lagrange][0]<<" Lagrange ideg = "<<ideg<<"\n");
                  BOOST_CHECK_MESSAGE(err6dx[(int)FEBasis::Bezier][0] < tol,
                    "Large err6 "<< err6dx[(int)FEBasis::Bezier][0]<<" Bézier ideg = "<<ideg<<"\n");
                  BOOST_CHECK_MESSAGE(err6dx[(int)FEBasis::Lagrange][0] < tol,
                    "Large err6 "<< err6dx[(int)FEBasis::Lagrange][0]<<" Lagrange ideg = "<<ideg<<"\n");
                } 



            /* Test second derivatives */
                for(auto ibasis : {FEBasis::Lagrange, FEBasis::Bezier}){
                  evalf(msh.coord,ent2poi[ielem],ibasis,DifVar::Bary,DifVar::Bary,
                    bar0[isamp],eva3s[0],jmat3s[0],hmat3);
                  evalM(msh.met.rfld  ,ent2poi[ielem],ibasis,DifVar::Bary,DifVar::Bary,
                    bar0[isamp],eva6s[0],jmat6s[0],hmat6);
                  int idx = 0;
                  for(double dx = dx0; dx > dx1; dx /= qdx){
                    for(int ii = 0; ii < idim; ii++){
                      double bary[idim+1];
                      for(int jj = 0; jj < idim + 1; jj++) bary[jj] = bar0[isamp][jj];
                        bary[ii + 1] += dx;
                      bary[0]      -= dx;
                      evalf(msh.coord,ent2poi[ielem],ibasis,DifVar::Bary,DifVar::None,
                        bary,eva3s[ii+1],jmat3s[ii+1],NULL);
                      evalM(msh.met.rfld  ,ent2poi[ielem],ibasis,DifVar::Bary,DifVar::None,
                        bary,eva6s[ii+1],jmat6s[ii+1],NULL);
                    }
                    for(int idif1 = 0; idif1 < idim; idif1++){
                      for(int idif2 = 0; idif2 < idim; idif2++){
                        for(int k = 0; k < idim; k++){
                          hmat3_disc[sym2idx(idif1,idif2)][k] 
                          = (jmat3s[1+idif1][idim*idif2+k] - jmat3s[0][idim*idif2+k])/dx;
                        }
                      }
                    }
                    for(int idif1 = 0; idif1 < idim; idif1++){
                      for(int idif2 = 0; idif2 < idim; idif2++){
                        for(int k = 0; k < nhess; k++){
                          hmat6_disc[sym2idx(idif1,idif2)][k] 
                          = (jmat6s[1+idif1][nhess*idif2+k] - jmat6s[0][nhess*idif2+k])/dx;
                        }
                      }
                    }
                    double err3 = sqrt(geterrl2<idim*nhess  >(hmat3_disc[0],hmat3));
                    double err6 = sqrt(geterrl2<nhess*nhess >(hmat6_disc[0],hmat6));
                    err3dx[(int)ibasis][idx] = MAX(err3,1.0e-32);
                    err6dx[(int)ibasis][idx] = MAX(err6,1.0e-32);
                    logdx[idx]  = log(dx);
                    idx++;

                    if(istr8){
                      double nrm3 = getnrml2<idim*nhess>(hmat3);
                      BOOST_CHECK_MESSAGE(nrm3 < tol, "hmat3 norm "<<nrm3<<" on straight mesh "<<s<<" ibasis = " << (int)ibasis);
                      if(nrm3>=tol){
                        printf("Debug hmat3\n");
                        dblAr2(nhess,idim,hmat3).print();
                        printf("Debug jmats0 jmats1:\n");
                        dblAr2(idim,idim,jmat3s[0]).print();
                        dblAr2(idim,idim,jmat3s[1]).print();
                        double dbg[idim*idim];
                        evalP1(msh.coord,ent2poi[ielem],ibasis,DifVar::Bary,DifVar::None,
                         bar0[isamp],eva3s[1],dbg,NULL);
                        printf("Eval as P1\n");
                        dblAr2(idim,idim,dbg).print();
                        evalf(msh.coord,ent2poi[ielem],ibasis,DifVar::Bary,DifVar::Bary,
                          bar0[isamp],eva3s[0],dbg,hmat3);
                        printf("Eval as P3 with dif2: \n");
                        dblAr2(idim,idim,dbg).print();
                        printf("eval = "); dblAr1(idim,eva3s[0]).print();
                        evalf(msh.coord,ent2poi[ielem],ibasis,DifVar::Bary,DifVar::None,
                          bar0[isamp],eva3s[0],dbg,hmat3);
                        printf("Eval as P3 no dif2: \n");
                        dblAr2(idim,idim,dbg).print();
                        printf("eval = "); dblAr1(idim,eva3s[0]).print();
                        wait();

                      }
                    }
                  }
                }
                if(!istr8){
                  double minerr3[2], minerr6[2];
                  for(auto ibasis : {FEBasis::Lagrange, FEBasis::Bezier}){
                    double minerr3 = 1.0e32;
                    double minerr6 = 1.0e32;
                    for(int idx = 0; idx < ndx; idx++){
                      minerr3 = MIN(minerr3, err3dx[(int)ibasis][idx]);
                      minerr6 = MIN(minerr6, err6dx[(int)ibasis][idx]);
                      err3dx[(int)ibasis][idx] = log(err3dx[(int)ibasis][idx]);
                      err6dx[(int)ibasis][idx] = log(err6dx[(int)ibasis][idx]);
                    }

                    if(minerr3 >= tol){
                      double slop3 = linearRegression(ndx,logdx,err3dx[(int)ibasis]);
                      BOOST_CHECK_MESSAGE(minsl < slop3,
                      " D2 coord slope ("<<basname[(int)ibasis]<<") "<<slop3<<" under minimum "<<minsl<<"\n");
                    }
                    if(minerr6 >= tol){
                      double slop6 = linearRegression(ndx,logdx,err6dx[(int)ibasis]);
                      BOOST_CHECK_MESSAGE(minsl < slop6,
                      " D2 metric slope ("<<basname[(int)ibasis]<<") "<<slop6<<" under minimum "<<minsl<<"\n");
                    }
                  }
                }


              }
            }   
            printf("-- Diff tests passed.\n");
      // END body
        }}CT_FOR1(idim);
      }}CT_FOR1(ideg);

    }catch(const MetrisExcept& e){
      printf("-- Degree not compiled for or mesh does not exist\n");
      continue;
    }

  }// fortest cases

}// end boost test case











