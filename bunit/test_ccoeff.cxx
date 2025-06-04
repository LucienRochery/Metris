//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include "common_setup.hxx"

#include "../src/Metris.h"
//#include "../src/codegen_ccoef.hxx"
//#include "../src/low_ccoef.hxx"
//#include "../src/msh_lag2bez.hxx"
//#include "../src/linalg/det.hxx"


#include <boost/hana.hpp> 
#include <boost/test/included/unit_test.hpp> 
#include <cmath>
namespace hana = boost::hana;
using namespace hana::literals;

namespace utf = boost::unit_test;

using namespace Metris;

typedef MetricFieldAnalytical MFT;

// 1.0e-8 relative error -> 1.0e-6% utf::tolerance()
BOOST_AUTO_TEST_CASE(ccoef, * utf::tolerance(double(1.0e-6)) ) 
{ 

  // bool is whether straight
  std::vector<std::pair<std::string,bool>> meshes = {
  {"../cases/2D/square.p1.100.meshb",true},
  {"../cases/2D/square.circmet.5k.curved",false},
  {"../cases/1200_p2.meshb",true},
  #if METRIS_MAX_DEG >= 3
  {"../cases/1200_p3.meshb",true},
  #endif
  #if METRIS_MAX_DEG >= 4
  {"../cases/1200_p4.meshb",true},
  #endif
  {"../cases/curved_p2.meshb",false}
  #if METRIS_MAX_DEG >= 3
  ,{"../cases/curved_p3.meshb",false}
  #endif
  #if METRIS_MAX_DEG >= 4
  ,{"../cases/curved_p4.meshb",false}
  #endif
  };


  double tol = 1.0e-12;
  for(auto testcase : meshes)
  {
    std::string s = testcase.first;
    bool istr8    = testcase.second;

    cargHandler arg("-in " + s + "  -anamet 1" 
    #ifdef NDEBUG 
      +" -verb 0"
    #endif
      );
    
    MeshTestSetup<MFT> f(arg.c, arg.v);
    Mesh<MFT> &msh = *(f.msh);
    msh.cleanup();

    std::cout<<"\n\n-- Mesh "<<s<<" dim "<<msh.idim<<" deg "<<msh.curdeg<<"\n";

    int ps;

    CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
    CT_FOR0_INC(2,3,idim){if(idim == msh.idim){

      msh.setBasis(FEBasis::Bezier);

      constexpr auto ordent = ORDELT(idim);
      constexpr int jdeg = idim*(ideg-1);
      constexpr int nnodj = idim == 2 ? getnnod2(jdeg) : getnnod3(jdeg);
      const int vol0 = ifact<idim>();
      const intAr2 &ent2poi = msh.ent2poi(idim);
      const int nentt = msh.nentt(idim);
      double ccoef[nnodj], ccoef2[nnodj];

      constexpr auto ccoef_genbez = idim == 2 ? ccoef_genbez2<ideg> : ccoef_genbez3<ideg>;
      constexpr auto evalf =  idim == 2 ? eval2<idim,ideg> : eval3<idim,ideg>;
      constexpr auto evalj =  idim == 2 ? eval2<1   ,jdeg> : eval3<1   ,jdeg>;

      #ifdef NDEBUG
        printf("-- Release target: benchmarks\n");
        std::cout<<"Mesh "<<s<<"\n";
        double t0,t1,dum[3]={0.0};
        double jmin,jmax;

        printf("--- Codegen \n");
        jmin = 1.0e30;
        jmax =-1.0e30;
        t0 = get_wall_time();
        for(int ielem = 0; ielem < nentt; ielem++){
          ccoef_genbez(ent2poi,msh.coord,ielem,ccoef);
          dum[0] += ccoef[0];
          double vol = getmeasentP1<idim>(ent2poi[ielem],msh.coord)*vol0;
          for(int i = 0; i < nnodj; i++){
            jmin = MIN(ccoef[i] / vol, jmin);
            jmax = MAX(ccoef[i] / vol, jmax);
          }
        }
        t1 = get_wall_time();
        ps = (int)(nentt/(t1-t0));
        ps /= 1000;
        printf(" %2.0f P%d Full elt coefs %dk/s \n",dum[0],ideg,ps);
        printf("   %15.8e < J_K < %15.8e \n",jmin,jmax);
        if(istr8){
          BOOST_TEST(jmin == 1.0);
          BOOST_TEST(jmax == 1.0);
        }
      
        printf("--- Manual compute at nodes \n");
        jmin = 1.0e30;
        jmax =-1.0e30;
        t0 = get_wall_time();
        for(int ielem = 0; ielem < nentt; ielem++){
          ccoef_eval<idim,idim,ideg>(msh.getBasis(),ent2poi,msh.coord,ielem,NULL,ccoef);
          dum[0] += ccoef[0];
          double vol = getmeasentP1<idim>(ent2poi[ielem],msh.coord)*vol0;
          for(int i = 0; i < nnodj; i++){
            jmin = MIN(ccoef[i] / vol, jmin);
            jmax = MAX(ccoef[i] / vol, jmax);
          }
        }
        t1 = get_wall_time();
        ps = (int)(nentt/(t1-t0));
        ps /= 1000;
        printf(" %2.0f P%d Full elt coefs %dk/s \n",dum[0],ideg,ps);
        printf("   %15.8e < J_K < %15.8e \n",jmin,jmax);
        if(istr8){
          BOOST_TEST(jmin == 1.0);
          BOOST_TEST(jmax == 1.0);
        }
      #endif



      bool first = true;
      double err2 = 0.0;
      double erri = -1.0e30;
      double err2_ev1 = 0.0;
      double erri_ev1 = -1.0e30;
      double err2_ev2 = 0.0;
      double erri_ev2 = -1.0e30;
      int lcoef[nnodj];
      for(int i=0; i<nnodj; i++){
        lcoef[i] = i;
      }
  
      bool firstnan = true;
      for(int ielem = 0; ielem < nentt; ielem++){
        ccoef_genbez(ent2poi,msh.coord,ielem,ccoef);
        double vol = vol0 * getmeasentP1<idim>(ent2poi[ielem],msh.coord);
        //double vol = det3_vdif(msh.coord[ent2poi(ielem,1)],msh.coord[ent2poi(ielem,0)]
        //                      ,msh.coord[ent2poi(ielem,2)],msh.coord[ent2poi(ielem,0)]
        //                      ,msh.coord[ent2poi(ielem,3)],msh.coord[ent2poi(ielem,0)]); 
        ccoef_eval<idim,idim,ideg>(msh.getBasis(),ent2poi,msh.coord,ielem,NULL,ccoef2);
        err2 += geterrl2<nnodj>(ccoef2,ccoef)/vol/vol;
        for(int i = 0; i < nnodj; i++){
          //int idx1 = ordtet.s[jdeg][i][0];
          //int idx2 = ordtet.s[jdeg][i][1];
          //int idx3 = ordtet.s[jdeg][i][2];
          //int idx4 = ordtet.s[jdeg][i][3];
          erri = erri > abs(ccoef[i]-ccoef2[i])/vol ? erri : abs(ccoef[i]-ccoef2[i])/vol;
          //if(first && abs(ccoef[i]-ccoef2[i]) > tol)printf("dbg %d (%d%d%d%d) ccoef = %20.15e ccoef2 = %20.15e err = %20.15e addr1 %p \n",i,
          //  idx1,idx2,idx3,idx4,ccoef[i],ccoef2[i],abs(ccoef[i]-ccoef2[i]),(void *)&ccoef[i]);
        }
        if( (std::isnan(err2) || std::isnan(erri)) && firstnan){
          firstnan = false;
          printf("Error became NaN\n");
          printf("Printing ccoef1:\n");
          for(int i = 0; i < nnodj ;i++){
            //int idx1 = ordtet.s[jdeg][i][0];
            //int idx2 = ordtet.s[jdeg][i][1];
            //int idx3 = ordtet.s[jdeg][i][2];
            //int idx4 = ordtet.s[jdeg][i][3];
            //printf("NaN %d (%d%d%d%d) ccoef = %20.15e ccoef2 = %20.15e err = %20.15e addr1 %p \n",i,
            //idx1,idx2,idx3,idx4,ccoef[i],ccoef2[i],abs(ccoef[i]-ccoef2[i]),(void *)&ccoef[i]);
            printf("NaN %d (",i);
            for(int jj = 0; jj < idim + 1; jj++) printf("%d",ordent[jdeg][i][jj]);
            printf(") ccoef = %20.15e ccoef2 = %20.15e err = %20.15e addr1 %p \n",
                  ccoef[i],ccoef2[i],abs(ccoef[i]-ccoef2[i]),(void *)&ccoef[i]);
          }
        }
        /* Compare Jacobian determinant to evaluations using ccoef */

        double jmat[idim*idim],dum[idim],bary[idim+1];
        for(int irnk = 0; irnk < nnodj; irnk++){
          for(int ii = 0; ii < idim + 1; ii++) 
            bary[ii] = ordent[jdeg][irnk][ii]/(1.0*jdeg);
          //bary[0] = ordtet.s[jdeg][irnk][0]/(1.0*jdeg);
          //bary[1] = ordtet.s[jdeg][irnk][1]/(1.0*jdeg);
          //bary[2] = ordtet.s[jdeg][irnk][2]/(1.0*jdeg);
          //bary[3] = ordtet.s[jdeg][irnk][3]/(1.0*jdeg);
          //int idx1 = ordtet.s[jdeg][irnk][0];
          //int idx2 = ordtet.s[jdeg][irnk][1];
          //int idx3 = ordtet.s[jdeg][irnk][2];
          //int idx4 = ordtet.s[jdeg][irnk][3];
          //eval3<3,ideg,ilag,1,0>(msh.coord,ent2poi[ielem],bary,dum,jmat,NULL);
          evalf(msh.coord,ent2poi[ielem],msh.getBasis(),DifVar::Bary,DifVar::None,
                bary,dum,jmat,NULL);
          double det = detmat<idim>(jmat);
          double detc1; 
          //eval3<1,jdeg,ilag,0,0>(dblAr2(nnodj,1,ccoef),lcoef,bary,&detc1,NULL,NULL);
          dblAr2 ccoef_(nnodj,1,ccoef);
          evalj(ccoef_,lcoef,FEBasis::Bezier,DifVar::None,DifVar::None,
                bary,&detc1,NULL,NULL);
          err2_ev1 += geterrl2<1>(&detc1,&det)/vol/vol;
          erri_ev1 = erri_ev1 > abs(det-detc1)/vol ? erri_ev1 : abs(det-detc1)/vol;


          double detc2; 
          //eval3<1,jdeg,ilag,0,0>(dblAr2(nnodj,1,ccoef2),lcoef,bary,&detc2,NULL,NULL);
          dblAr2 ccoef2_(nnodj,1,ccoef2);
          evalj(ccoef2_,lcoef,FEBasis::Bezier,DifVar::None,DifVar::None,
                bary,&detc2,NULL,NULL);
          err2_ev2 += geterrl2<1>(&detc2,&det)/vol/vol;
          erri_ev2 = erri_ev2 > abs(det-detc2)/vol ? erri_ev2 : abs(det-detc2)/vol;
          //if(first && (abs(detc1-det) > tol || abs(detc2-det) > tol) )printf("(eval) %d (%d%d%d%d) detc2 = %12.5e err1 = %12.5e err2 = %12.5e bary = %f %f %f %f \n",irnk,
          //  idx1,idx2,idx3,idx4,detc2,abs(detc1-det),abs(detc2-det),bary[0],bary[1],bary[2],bary[3]);
          if(detc1>1.0e10){
            printf("## VERY LARGE VALUE %f \n",detc1);
            printf("ccoef: \n");
            dblAr1(nnodj,ccoef).print(nnodj);
            printf("\n");
            printf("Deg = %d nnode jac %d \n",ideg,nnodj);
          }
        }

        //if(first) wait();
        first = false;
      }



      msh.setBasis(FEBasis::Lagrange);
      /* ccoef_genbez only for Bézier but ccoef_eval3 can be tried in Lagrange */

      double erri_ev2_lag = -1.0e30;
      double err2_ev2_lag = 0.0;
      first = true;
      for(int ielem = 0; ielem < nentt; ielem++){
        double vol = vol0 * getmeasentP1<idim>(ent2poi[ielem],msh.coord);
        //double vol = det3_vdif(msh.coord[ent2poi(ielem,1)],msh.coord[ent2poi(ielem,0)]
        //                      ,msh.coord[ent2poi(ielem,2)],msh.coord[ent2poi(ielem,0)]
        //                      ,msh.coord[ent2poi(ielem,3)],msh.coord[ent2poi(ielem,0)]); 
        ccoef_eval<idim,idim,ideg>(msh.getBasis(),ent2poi,msh.coord,ielem,NULL,ccoef2);
        /* Compare Jacobian determinant to evaluations using ccoef */
        double jmat[idim*idim],dum[idim],bary[idim+1];
        for(int irnk = 0; irnk < nnodj; irnk++){
          for(int ii = 0; ii < idim + 1; ii++) 
            bary[ii] = ordent[jdeg][irnk][ii]/(1.0*jdeg);
          //bary[0] = ordtet.s[jdeg][irnk][0]/(1.0*jdeg);
          //bary[1] = ordtet.s[jdeg][irnk][1]/(1.0*jdeg);
          //bary[2] = ordtet.s[jdeg][irnk][2]/(1.0*jdeg);
          //bary[3] = ordtet.s[jdeg][irnk][3]/(1.0*jdeg);
          //int idx1 = ordtet.s[jdeg][irnk][0];
          //int idx2 = ordtet.s[jdeg][irnk][1];
          //int idx3 = ordtet.s[jdeg][irnk][2];
          //int idx4 = ordtet.s[jdeg][irnk][3];
          //eval3<3,ideg,1,1,0>(msh.coord,ent2poi[ielem],bary,dum,jmat,NULL);
          evalf(msh.coord,ent2poi[ielem],msh.getBasis(),DifVar::Bary,DifVar::None,
                bary,dum,jmat,NULL);
          double det = detmat<idim>(jmat);
          double detc2; 
          //eval3<1,jdeg,0,0,0>(dblAr2(nnodj,1,ccoef2),lcoef,bary,&detc2,NULL,NULL);
          dblAr2 ccoef2_(nnodj,1,ccoef2);
          evalj(ccoef2_,lcoef,FEBasis::Bezier,DifVar::None,DifVar::None,
                bary,&detc2,NULL,NULL);
          err2_ev2_lag += geterrl2<1>(&detc2,&det)/vol/vol;
          erri_ev2_lag = erri_ev2_lag > abs(det-detc2)/vol ? erri_ev2_lag : abs(det-detc2)/vol;
          //if(first && abs(detc2-det) > tol )printf("(lagval) %d (%d%d%d%d) detc2 = %12.5e err = %12.5e \n",irnk,
          //                          idx1,idx2,idx3,idx4,detc2,abs(detc2-det));
        }

        //if(first) wait();
        first = false;
      }




      //printf("-- Results %s\n",s.c_str());
      err2 /= nnodj*nentt;
      err2 = sqrt(err2);
      err2_ev1 /= nnodj*nentt;
      err2_ev1 = sqrt(err2_ev1);
      err2_ev2 /= nnodj*nentt;
      err2_ev2 = sqrt(err2_ev2);
      BOOST_TEST(erri <= tol);
      BOOST_TEST(err2 <= tol);
      BOOST_TEST(erri_ev1 <= tol);
      BOOST_TEST(err2_ev1 <= tol);
      BOOST_TEST(erri_ev2 <= tol);
      BOOST_TEST(err2_ev2 <= tol);
      //printf("   %15.8e < J_K < %15.8e \n",jmin,jmax);
      //printf("-- Bézier\n");
      //printf("Debug erri = %20.15e err2 = %20.15e\n",erri,err2);
      //printf("Dbg ev1 erri = %20.15e err2 = %20.15e\n",erri_ev1,err2_ev1);
      //printf("Dbg ev2 erri = %20.15e err2 = %20.15e\n",erri_ev2,err2_ev2);
      //printf("-- Lagrange\n");
      //printf("Dbg ev2 erri = %20.15e err2 = %20.15e\n",erri_ev2_lag,err2_ev2_lag);

      //printf("\n\n");
    }}CT_FOR1(idim);
    }}CT_FOR1(ideg);
  }

}