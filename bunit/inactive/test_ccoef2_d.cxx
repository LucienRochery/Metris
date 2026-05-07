//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test

//#include <src/codegen_ccoef_d.hxx>

#include <src/low_ccoef_d.hxx>
#include <src/low_geo/ccoef.hxx>

#include <boost/test/included/unit_test.hpp>
#include <bunit/common_setup.hxx>

#include <boost/timer/progress_display.hpp>

#include "utils/aux_misc.hxx"
#include "quality/low_metqua.hxx"
#include "quality/low_metqua_d.hxx"

#include "..libs/SANS/Surreal/SurrealS.h"
#include "..libs/SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"


#include <boost/hana.hpp>
namespace hana = boost::hana;
using namespace hana::literals;

namespace utf = boost::unit_test;

namespace Metris{


typedef MetricFieldFE MFT;


BOOST_AUTO_TEST_CASE(test_ccoef2_d)
{
  std::vector<std::string> meshes = {
   "../cases/2D/square.circmet.50.curved.meshb",
   "../cases/2D/square.circmet.5k.curved.meshb"
  };


  for(auto s : meshes)
  {

    std::cout<<"Mesh "<<s<<std::endl;

    cargHandler arg("-i "+s);
    MeshTestSetup<MFT> f(arg.c,arg.v);
    Mesh<MFT> &msh = *(f.msh);

    msh.setBasis(FEBasis::Bezier);



    CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
    CT_FOR0_INC(2,3,idim){if(idim == msh.idim){


      constexpr int jdeg = idim*(ideg-1);
      constexpr int nnodj = idim == 2 ? getnnod2(jdeg) : getnnod3(jdeg);
      const int nnode = msh.nnode(idim);

      //constexpr auto d_ccoef_genbez = idim == 2 ? d_ccoef_genbez2<ideg> : d_ccoef_genbez3<ideg>;

      const int nentt = msh.nentt(idim);
      const intAr2& ent2poi = msh.ent2poi(idim);

      const double h = 1;
      const double tol_relative = 1e-10;
      const double tol0 = 1.0e-14;
      int num_total = 0;
      int num_pass = 0;

      for (int ientt = 0; ientt < nentt; ientt++){

        if(isdeadent(ientt,ent2poi)) continue;

        double ccoef[nnodj], dum_ccoef[nnodj];
        dblAr2 d_ccoef_x(nnodj,nnode), d_ccoef_y(nnodj,nnode);

        //d_ccoef_genbez(ent2poi, msh.coord, ientt, 0, d_ccoef_x);
        //d_ccoef_genbez(ent2poi, msh.coord, ientt, 1, d_ccoef_y);

        getccoef_dcoord<idim,ideg>(msh,ientt,0,ccoef,d_ccoef_x);
        getccoef_dcoord<idim,ideg>(msh,ientt,1,ccoef,d_ccoef_y);

        for (int i = 0; i < nnodj; i++) //  N_i
        {
          double ccoef_fd_x[nnodj], ccoef_fd_y[nnodj];
          dblAr2 d_ccoef_fd_x(nnodj, nnodj), d_ccoef_fd_y(nnodj, nnodj);

          for (int j = 0; j < nnode; j++) //   P_j  -> compute D_{P_j} (N_i)
          {
            int idx = ent2poi(ientt,j);

            msh.coord(idx,0) += h;
            getccoef<idim,idim,ideg>(msh, ientt, NULL, ccoef_fd_x);
            msh.coord(idx,0) -= h;

            msh.coord(idx,1) += h;
            getccoef<idim,idim,ideg>(msh, ientt, NULL, ccoef_fd_y);
            msh.coord(idx,1) -= h;

            d_ccoef_fd_x(i,j) = (ccoef_fd_x[i] - ccoef[i]) / h;
            d_ccoef_fd_y(i,j) = (ccoef_fd_y[i] - ccoef[i]) / h;

            num_total++;
            if (std::abs(d_ccoef_fd_x(i,j) - d_ccoef_x(i,j)) < tol_relative
              &&std::abs(d_ccoef_fd_y(i,j) - d_ccoef_y(i,j)) < tol_relative) {
                num_pass++;
            }else{
              std::cout << "d_ccoef_fd_x[" << i << "][" << j << "] = "
              << std::setw(10) << d_ccoef_fd_x(i,j)
              << "     d_ccoef_fd_y[" << i << "][" << j << "] = "
              << std::setw(10) << d_ccoef_fd_y(i,j)
              << "     d_ccoef_x[" << i << "][" << j << "] = "
              << std::setw(10) << d_ccoef_x(i,j)
              << "     d_ccoef_y[" << i << "][" << j << "] = "
              << std::setw(10) << d_ccoef_y(i,j) << std::endl;
              wait();
            }

            if(abs(d_ccoef_x(i,j)) < tol0){
              BOOST_TEST(abs(d_ccoef_fd_x(i,j)) < tol0 );
            }else{
              BOOST_TEST(abs(d_ccoef_fd_x(i,j) - d_ccoef_x(i,j)) < tol_relative * abs(d_ccoef_x(i,j)));
            }

            if(abs(d_ccoef_y(i,j)) < tol0){
              BOOST_TEST(abs(d_ccoef_fd_y(i,j)) < tol0 );
            }else{
              BOOST_TEST(abs(d_ccoef_fd_y(i,j) - d_ccoef_y(i,j)) < tol_relative * abs(d_ccoef_y(i,j)));
            }

            //BOOST_CHECK_CLOSE(d_ccoef_fd_x[i][j], d_ccoef_x[i][j], tol_relative);
            //BOOST_CHECK_CLOSE(d_ccoef_fd_y[i][j], d_ccoef_y[i][j], tol_relative);

          }
        }
      }
      double pct_pass = (static_cast<double>(num_pass) / num_total) * 100.0;
      std::cout << "Percentage of tests passed: " << pct_pass << "%" << std::endl;

    }}CT_FOR1(idim);
    }}CT_FOR1(ideg);
  }
}
}


