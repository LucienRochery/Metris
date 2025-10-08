//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#define BOOST_TEST_MODULE test_metqua_d2

#include <boost/test/included/unit_test.hpp>
#include <random>

#include "gen_bary.hxx"

#include "utils/CT_loop.hxx"
#include "utils/mprintf.hxx"
#include "quality/low_metqua.hxx"
#include "quality/low_metqua_d.hxx"
#include "metris_options.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "Mesh/Mesh.hxx"

#include "quality/low_metqua.hxx"
#include "quality/quafun_tradet.hxx"
#include "quality/quafun.hxx"

#include "SANS/Surreal/SurrealS.h"

namespace Metris{

typedef MetricFieldAnalytical MFT;
typedef double ftype;
typedef std::pair<AsDeg,AsDeg> AsDegPair;

BOOST_AUTO_TEST_CASE(test_metqua_d2)
{
  std::cout << "================================" << std::endl;
  std::cout << "Testing with quafun_distortion" << std::endl;
  std::cout << "================================" << std::endl;

  std::vector<std::string> meshes =
  {METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
  ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
  ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.5k"
  };


  for(auto testcase : meshes){
    std::string mesh_name = testcase;
    std::cout<<"-- Test case mesh = "<<mesh_name<<std::endl;
    cargHandler arg("-in "+mesh_name+" -verb 0 -anamet 1 -sclmet 1");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
    msh.met.setSpace(MetSpace::Log);

    dblAr2 bary[4];
    const int nsamp = 20;
    genBary(nsamp, 2, bary[2]);
    genBary(nsamp, 3, bary[3]);

    GETVDEPTH(msh.param);

    CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
      constexpr int nhess = (gdim*(gdim+1))/2;
      // Use Boost PP to loop between options in the future
      constexpr auto iquaf = QuaFun::Distortion;


      CT_FOR0_INC(2,gdim,tdim){if(tdim <= msh.get_tdim()){
        INCVDEPTH(msh.param);
        MPRINTF("-- Tests tdim = {} gdim = {}\n",tdim,gdim);
        constexpr auto quafun_xi   = get_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
        constexpr auto d_quafun_xi = get_d_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
        const intAr2 &ent2poi = msh.ent2poi(tdim);
        const int nnode = getnnode(tdim, ideg);

        {//INCVDEPTH
        INCVDEPTH(msh.param);
        MPRINTF("-- Test 1: quafun_xi == d_quafun_xi\n");
        for(AsDegPair asdeg : {AsDegPair({AsDeg::P1, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::Pk})}){
          std::string asdmsh_str = asdeg.first  == AsDeg::P1 ? "P1" : "Pk";
          std::string asdmet_str = asdeg.second == AsDeg::P1 ? "P1" : "Pk";
          double err_min = 1.0e30;
          double err_max = -1.0;
          double err_avg = 0;
          double ntest = 0;
          for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
            for(int isamp = 0; isamp < nsamp; isamp++){
              double qua = quafun_xi(msh, asdeg.first, asdeg.second, ent2poi[ientt],
                                     bary[tdim][isamp], NULL);
              double dquael[gdim], hquael[nhess];
              for(int ivar = 0; ivar < nnode; ivar++){
                double qua_d = d_quafun_xi(msh, asdeg.first, asdeg.second, ent2poi[ientt],
                                           bary[tdim][isamp], NULL,
                                           ivar, msh.getBasis(), DifVar::None,
                                           dquael, hquael);
                double err = abs(qua - qua_d);
                err_min = MIN(err_min, err);
                err_max = MAX(err_max, err);
                err_avg += err;
                ntest += 1;
              }
            }// for isamp
          }// for ientt
          if(ntest == 0) continue;
          err_avg /= ntest;
          MPRINTF(" - DONE using asdmsh = {}, asdmet = {}:"
                  " got err min = {}, avg = {}, max = {}\n",
                  asdmsh_str.c_str(),asdmet_str.c_str(),
                  err_min,err_avg,err_max);
          BOOST_TEST(err_max < 1.0e-12);
        }// for AsDegPair



        MPRINTF("-- Test 2.1: d_quafun_tradet_SurrealS derivatives (identity)\n");
        {// using
        constexpr int nvar = gdim;
        auto d_quafun = d_quafun_tradet_SurrealS<MFT,gdim,tdim,gdim>;
        auto   quafun =   quafun_tradet         <MFT,gdim,tdim,ftype>;
        SANS::DLA::MatrixS<gdim,nvar,double> dpoint = SANS::DLA::Identity();
        //for(int ii = 0; ii < gdim; ii++){
        //  for(int jj = 0; jj < gdim; jj++){
        //    dpoint(ii,jj) = ii == jj;
        //  }
        //}
        for(AsDegPair asdeg : {AsDegPair({AsDeg::P1, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::Pk})}){
          double ntest = 0;
          for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
            for(int isamp = 0; isamp < nsamp; isamp++){
              SANS::SurrealS<nvar, ftype> tra, det;
              ftype tra0, det0;

              quafun(msh, asdeg.first, asdeg.second,
                     ent2poi[ientt], bary[tdim][isamp],
                     NULL,
                     &tra0, &det0);
              for(int ivar = 0; ivar < nnode; ivar++){
                d_quafun(msh, asdeg.first, asdeg.second,
                         ent2poi[ientt], bary[tdim][isamp],
                         ivar, msh.getBasis(), DifVar::None,
                         NULL,
                         tra, det, &dpoint);
                ntest += 1;
              }
            }// for isamp
          }// for ientt
        }// for AsDegPair
        }// using

        }// INCVDEPTH

      }}CT_FOR1(tdim);
    }}CT_FOR1(ideg);
    }}CT_FOR1(gdim);


  }// for testcase


}// BOOST_AUTO_TEST_CASE

BOOST_AUTO_TEST_CASE(test_quafun_sizeshape)
{

  std::cout << "================================" << std::endl;
  std::cout << "Testing with quafun_sizeshape" << std::endl;
  std::cout << "================================" << std::endl;

  std::vector<std::string> meshes =
  {METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
  // ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
  // ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.5k"
  };


  for(auto testcase : meshes){
    std::string mesh_name = testcase;
    std::cout<<"-- Test case mesh = "<<mesh_name<<std::endl;
    cargHandler arg("-in "+mesh_name+" -verb 0 -anamet 1 -sclmet 1");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
    msh.met.setSpace(MetSpace::Log);

    dblAr2 bary[4];
    const int nsamp = 20;
    genBary(nsamp, 2, bary[2]);
    genBary(nsamp, 3, bary[3]);

    GETVDEPTH(msh.param);

    CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
      constexpr int nhess = (gdim*(gdim+1))/2;
      // Use Boost PP to loop between options in the future
      constexpr auto iquaf = QuaFun::SizeShape;


      CT_FOR0_INC(2,gdim,tdim){if(tdim <= msh.get_tdim()){
        INCVDEPTH(msh.param);
        MPRINTF("-- Tests tdim = {} gdim = {}\n",tdim,gdim);
        constexpr auto quafun_xi   = get_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
        constexpr auto d_quafun_xi = get_d_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
        const intAr2 &ent2poi = msh.ent2poi(tdim);
        const int nnode = getnnode(tdim, ideg);

        {//INCVDEPTH
        INCVDEPTH(msh.param);
        MPRINTF("-- Test 1: quafun_xi == d_quafun_xi\n");
        for(AsDegPair asdeg : {AsDegPair({AsDeg::P1, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::Pk})}){
          std::string asdmsh_str = asdeg.first  == AsDeg::P1 ? "P1" : "Pk";
          std::string asdmet_str = asdeg.second == AsDeg::P1 ? "P1" : "Pk";
          double err_min = 1.0e30;
          double err_max = -1.0;
          double err_avg = 0;
          double ntest = 0;
          for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
            for(int isamp = 0; isamp < nsamp; isamp++){
              double qua = quafun_xi(msh, asdeg.first, asdeg.second, ent2poi[ientt],
                                     bary[tdim][isamp], NULL);
              double dquael[gdim], hquael[nhess];
              for(int ivar = 0; ivar < nnode; ivar++){
                double qua_d = d_quafun_xi(msh, asdeg.first, asdeg.second, ent2poi[ientt],
                                           bary[tdim][isamp], NULL,
                                           ivar, msh.getBasis(), DifVar::None,
                                           dquael, hquael);
                double err = abs(qua - qua_d);
                err_min = MIN(err_min, err);
                err_max = MAX(err_max, err);
                err_avg += err;
                ntest += 1;
              }
            }// for isamp
          }// for ientt
          if(ntest == 0) continue;
          err_avg /= ntest;
          MPRINTF(" - DONE using asdmsh = {}, asdmet = {}:"
                  " got err min = {}, avg = {}, max = {}\n",
                  asdmsh_str.c_str(),asdmet_str.c_str(),
                  err_min,err_avg,err_max);
          BOOST_TEST(err_max < 1.0e-12);
        }// for AsDegPair

        MPRINTF("-- Test 2: Check quality value for knwon P1 cases (unit regular, scaled regular, anisotropic example)\n");
        if (gdim == 2 && tdim == 2 && ideg == 1 && msh.nentt(tdim) > 0) {
          const int ientt = 0;
          intAr2& ent2poi = msh.ent2poi(2);

          // save old coords to restore later
          dblAr2 oldCoord(3,2);
          for (int ipoi = 0; ipoi < 3; ipoi++) {
            oldCoord(ipoi,0) = msh.coord(ent2poi(ientt,ipoi),0);
            oldCoord(ipoi,1) = msh.coord(ent2poi(ientt,ipoi),1);
          }

          // ===================
          // Test unit regular
          // ===================

          // overwrite with unit equilateral triangle
          msh.coord(ent2poi(ientt,0),0) = 0.; msh.coord(ent2poi(ientt,0),1) = 0.;
          msh.coord(ent2poi(ientt,1),0) = 1.; msh.coord(ent2poi(ientt,1),1) = 0.;
          msh.coord(ent2poi(ientt,2),0) = 0.5; msh.coord(ent2poi(ientt,2),1) = sqrt(3.)/2.;

          // barycenter
          const double btri[3] = { 1./3., 1./3., 1./3. };

          // use identity metric for this test
          const double metI[3] = { 1., 0., 1.};

          double trueQuality = 1.;

          for (int powerSign : {+1, -1}) {
            msh.param->opt_power = powerSign;

            // value-only
            double qv = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

            // returns value as well, and additionally compute derivatives
            double dquael_loc[gdim], hquael_loc[(gdim*(gdim+1))/2];
            double qv_d = d_quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri,
                                    metI, /*ivar=*/0, msh.getBasis(), DifVar::None,
                                    dquael_loc, hquael_loc);

            BOOST_CHECK_SMALL(qv - trueQuality, 1e-13);
            BOOST_CHECK_SMALL(qv_d - trueQuality, 1e-13);
          }

          // ===================
          // Test scaled regular
          // ===================

          const double scale = 8.;

          msh.coord(ent2poi(ientt,0),0) *= scale; msh.coord(ent2poi(ientt,0),1) *= scale;
          msh.coord(ent2poi(ientt,1),0) *= scale; msh.coord(ent2poi(ientt,1),1) *= scale;
          msh.coord(ent2poi(ientt,2),0) *= scale; msh.coord(ent2poi(ientt,2),1) *= scale;

          trueQuality = 1./2. * (pow(scale,4.) + pow(scale,-4.));

          for (int powerSign : {+1, -1}) {
            msh.param->opt_power = powerSign;

            // value-only
            double qv = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

            // returns value as well, and additionally compute derivatives
            double dquael_loc[gdim], hquael_loc[(gdim*(gdim+1))/2];
            double qv_d = d_quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri,
                                    metI, /*ivar=*/0, msh.getBasis(), DifVar::None,
                                    dquael_loc, hquael_loc);

            if (powerSign == 1){
              BOOST_CHECK_SMALL(qv - trueQuality, 1e-12);
              BOOST_CHECK_SMALL(qv_d - trueQuality, 1e-12);
            }
            else{
              BOOST_CHECK_SMALL(qv - 1./trueQuality, 1e-12);
              BOOST_CHECK_SMALL(qv_d - 1./trueQuality, 1e-12);
            }
          }

          // ===================
          // Test anisotropic
          // ===================

          // vertices: (0,0) - (1,0) - (0.5,0.1)

          msh.coord(ent2poi(ientt,0),0) = 0.; msh.coord(ent2poi(ientt,0),1) = 0.;
          msh.coord(ent2poi(ientt,1),0) = 1.; msh.coord(ent2poi(ientt,1),1) = 0.;
          msh.coord(ent2poi(ientt,2),0) = 0.5; msh.coord(ent2poi(ientt,2),1) = 0.1;

          const dblAr1 e1({msh.coord(ent2poi(ientt,1),0) - msh.coord(ent2poi(ientt,0),0), msh.coord(ent2poi(ientt,1),1) - msh.coord(ent2poi(ientt,0),1)});
          const dblAr1 e2({msh.coord(ent2poi(ientt,2),0) - msh.coord(ent2poi(ientt,0),0), msh.coord(ent2poi(ientt,2),1) - msh.coord(ent2poi(ientt,0),1)});

          double tra = 4./3. * (e1[0]*e1[0] + e1[1]*e1[1] +  e2[0]*e2[0] + e2[1]*e2[1] - e1[0]*e2[0] - e1[1]*e2[1]);
          double det = 4./3. * (e1[0]*e2[1] - e1[1]*e2[0])*(e1[0]*e2[1] - e1[1]*e2[0]);

          trueQuality = 1./(2*2*2) * tra * tra * (1. + 1./(det*det));

          for (int powerSign : {+1, -1}) {
            msh.param->opt_power = powerSign;

            // value-only
            double qv = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

            // returns value as well, and additionally compute derivatives
            double dquael_loc[gdim], hquael_loc[(gdim*(gdim+1))/2];
            double qv_d = d_quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri,
                                    metI, /*ivar=*/0, msh.getBasis(), DifVar::None,
                                    dquael_loc, hquael_loc);

            if (powerSign == 1){
              BOOST_CHECK_SMALL(qv - trueQuality, 1e-12);
              BOOST_CHECK_SMALL(qv_d - trueQuality, 1e-12);
            }
            else{
              BOOST_CHECK_SMALL(qv - 1./trueQuality, 1e-12);
              BOOST_CHECK_SMALL(qv_d - 1./trueQuality, 1e-12);
            }
          }

          // Restore original coords
          for (int ipoin = 0; ipoin < 3; ipoin++) {
            msh.coord(ent2poi(ientt,ipoin),0) = oldCoord(ipoin,0);
            msh.coord(ent2poi(ientt,ipoin),1) = oldCoord(ipoin,1);
          }


        }
        MPRINTF("-- Test 3: Check quality derivatives using 1st order FD\n");
        if (gdim == 2 && tdim == 2 && ideg == 1 && msh.nentt(tdim) > 0) {

          // use identity metric for this test
          const double metI[3] = {1., 0., 1.};
          const double btri[3]  = {1./3., 1./3., 1./3.};

          const int nnode = getnnode(tdim, ideg);

          // check for several elements (not entire mesh)
          const int ncheck_ent = std::min(10, msh.nentt(tdim));

          // step sizes
          const int nh = 7;
          const double hvec[7] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};

          // to save max error
          double err_max_by_h[7] = {0,0,0,0,0,0,0};

          intAr2& ent2poi = msh.ent2poi(2);

          for (int ientt = 0; ientt < ncheck_ent; ientt++) {

            // Compute analytic value and gradient
            double qua0;
            dblAr1 grad_analytic(6); // 3 nodes with 2 components each

            // we will need the coords backup to restore after each perturbation
            dblAr2 coord_backup(3,2);
            for (int ipoin = 0; ipoin < 3; ipoin++) {
              for (int ii = 0; ii < gdim; ii++)
                coord_backup(ipoin,ii) = msh.coord(ent2poi(ientt,ipoin), ii);
            }

            // get base value once
            qua0 = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

            // fill analytic gradient per node
            for (int ivar = 0; ivar < nnode; ivar++) {
              double dqa[gdim], hqa[(gdim*(gdim+1))/2];
              double qua_d = d_quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI,
                                        ivar, msh.getBasis(), DifVar::None,
                                        dqa, hqa);
              // consistency : qua_d == qua0
              BOOST_CHECK_SMALL(qua_d - qua0, 1e-12);
              for (int ii = 0; ii < gdim; ii++)
                grad_analytic[ivar*gdim + ii] = dqa[ii];
            }

            // now finite differences for each coordinate
            for (int ivar = 0; ivar < nnode; ivar++) {
              for (int ii = 0; ii < gdim; ii++) {

                double prev_err = std::numeric_limits<double>::quiet_NaN();

                for (int ih = 0; ih < nh; ih++) {
                  const double h = hvec[ih];

                  // perturb x(ivar,k) by +h
                  msh.coord(ent2poi(ientt,ivar), ii) = coord_backup(ivar,ii) + h;

                  // recompute quality value
                  double quah = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

                  // forward finite difference
                  const double fd = (quah - qua0) / h;

                  // error against analytic
                  const double an = grad_analytic[ivar*gdim + ii];

                  const double err = std::abs(fd - an);
                  err_max_by_h[ih] = std::max(err_max_by_h[ih], err);

                  // crude rate check: first order
                  if (ih > 0 && std::isfinite(prev_err) && prev_err > 0) {
                    const double h_prev = hvec[ih-1];
                    const double order = std::log(prev_err/err) / std::log(h_prev/h);
                  }

                  // restore the perturbed coordinate
                  msh.coord(ent2poi(ientt,ivar), ii) = coord_backup(ivar,ii);
                  prev_err = err;
                } // hvec
              } // ii
            } // ivar
          } // ientt

          // print max errors per h
          MPRINTF(" - Gradient 1st order FD max errors by h:\n");
          for (int ih = 0; ih < nh; ih++) MPRINTF("  - h = {:.1e}, max err = {:.1e}\n",hvec[ih],err_max_by_h[ih]);

        }

        MPRINTF("-- Test 4: Check quality derivatives using 2nd order FD\n");
        if (gdim == 2 && tdim == 2 && ideg == 1 && msh.nentt(tdim) > 0) {

          // use identity metric for this test
          const double metI[3] = {1., 0., 1.};
          const double btri[3]  = {1./3., 1./3., 1./3.};

          const int nnode = getnnode(tdim, ideg);

          // check for several elements (not entire mesh)
          const int ncheck_ent = std::min(10, msh.nentt(tdim));

          // step sizes
          const int nh = 7;
          const double hvec[7] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};

          // to save max error
          double err_max_by_h[7] = {0,0,0,0,0,0,0};

          intAr2& ent2poi = msh.ent2poi(2);

          for (int ientt = 0; ientt < ncheck_ent; ientt++) {

            // Compute analytic value and gradient
            double qua0;
            dblAr1 grad_analytic(6); // 3 nodes with 2 components each

            // we will need the coords backup to restore after each perturbation
            dblAr2 coord_backup(3,2);
            for (int ipoin = 0; ipoin < 3; ipoin++) {
              for (int ii = 0; ii < gdim; ii++)
                coord_backup(ipoin,ii) = msh.coord(ent2poi(ientt,ipoin), ii);
            }

            // get base value once
            qua0 = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

            // fill analytic gradient per node
            for (int ivar = 0; ivar < nnode; ivar++) {
              double dqa[gdim], hqa[(gdim*(gdim+1))/2];
              double qua_d = d_quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI,
                                        ivar, msh.getBasis(), DifVar::None,
                                        dqa, hqa);
              // consistency : qua_d == qua0
              BOOST_CHECK_SMALL(qua_d - qua0, 1e-12);
              for (int ii = 0; ii < gdim; ii++)
                grad_analytic[ivar*gdim + ii] = dqa[ii];
            }

            // now finite differences for each coordinate
            for (int ivar = 0; ivar < nnode; ivar++) {
              for (int ii = 0; ii < gdim; ii++) {

                double prev_err = std::numeric_limits<double>::quiet_NaN();

                for (int ih = 0; ih < nh; ih++) {
                  const double h = hvec[ih];

                  // perturb x(ivar,k) by +h
                  msh.coord(ent2poi(ientt,ivar), ii) = coord_backup(ivar,ii) + h;
                  // compute quality value
                  double quahplus = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);
                  // perturb x(ivar,k) by -h
                  msh.coord(ent2poi(ientt,ivar), ii) = coord_backup(ivar,ii) - h;
                  // compute quality value
                  double quahminus = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

                  // forward finite difference
                  const double fd = (quahplus - quahminus) / (2.*h);

                  // error against analytic
                  const double an = grad_analytic[ivar*gdim + ii];

                  const double err = std::abs(fd - an);
                  err_max_by_h[ih] = std::max(err_max_by_h[ih], err);

                  // crude rate check: first order
                  if (ih > 0 && std::isfinite(prev_err) && prev_err > 0) {
                    const double h_prev = hvec[ih-1];
                    const double order = std::log(prev_err/err) / std::log(h_prev/h);
                  }

                  // restore the perturbed coordinate
                  msh.coord(ent2poi(ientt,ivar), ii) = coord_backup(ivar,ii);
                  prev_err = err;
                } // hvec
              } // ii
            } // ivar
          } // ientt

          // print max errors per h
          MPRINTF(" - Gradient 2nd order FD max errors by h:\n");
          for (int ih = 0; ih < nh; ih++) MPRINTF("  - h = {:.1e}, max err = {:.1e}\n",hvec[ih],err_max_by_h[ih]);

        }
        }// INCVDEPTH

      }}CT_FOR1(tdim);
    }}CT_FOR1(ideg);
    }}CT_FOR1(gdim);

  }// for testcase

}// BOOST_AUTO_TEST_CASE

}// namespace