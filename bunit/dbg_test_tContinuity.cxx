//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_tContinuity

#include "common_setup.hxx"


using namespace Metris;

typedef MetricFieldAnalytical MFT;


BOOST_AUTO_TEST_CASE(test_tContinuity) 
{
  std::vector<std::string> meshes = {
    METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
   ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
  };

  for(auto mesh : meshes){
    cargHandler arg("-in "+mesh+" -verb 0 -anamet 1 -sclmet 1");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);

    fmt::print("-- Testing {}D\n",msh.idim);

    const int tdim = msh.get_tdim();
    const int nedgl = (tdim*(tdim+1))/2;
    const intAr2 lnoed(nedgl,2,tdim == 2 ? lnoed2[0] : lnoed3[0]);
    const int nnmet = (msh.idim*(msh.idim+1))/2;
    const intAr2& ent2poi = msh.ent2poi(tdim);
    const int nentt = msh.nentt(tdim);

    HshTab_I2I ledge;
    ledge.reserve(nentt);
    for(int ientt = 0; ientt < nentt; ientt++){
      if(isdeadent(ientt,ent2poi)) continue;
      for(int ied = 0; ied < nedgl; ied++){
        int ip1 = ent2poi(ientt, lnoed(ied,0));
        int ip2 = ent2poi(ientt, lnoed(ied,1));
        auto key = stup2(ip1,ip2);
        
        // Check edge already seen
        if(ledge.find(key) != ledge.end()) continue;
        
        ledge[key] = ientt;
      }
    }

    double max_err_mdu    = -1;
    double max_err_nrmmdu = -1;

    intAr1 lsedg(1);
    intAr1 lsfac(2);
    intAr1 lstet(100);
    double bary[4]; // if considering HO elements, build bary to match the edge
    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0/(tdim + 1);
    dblAr2 lsmet(100, nnmet);
    double eigval[3], eigvec[9];
    for(auto iedge : ledge){
      int ientt = iedge.second;
      auto key = iedge.first;
      int ip1 = std::get<0>(key);
      int ip2 = std::get<1>(key);
      int iopen;
      shell(msh, ip1, ip2, tdim, ientt,
            lsedg, lsfac, lstet, &iopen);
      intAr1 &lsent = tdim == 2 ? lsfac : lstet;
      if(lsent.get_n() == 1) continue;

      lsmet.set_n(lsent.get_n());

      for(int ishell = 0; ishell < lsent.get_n(); ishell++){
        int ieshl = lsent[ishell];
        int ierro;
        MSH_DIM_DEG0(msh){
          CT_FOR0_INC(2,gdim,tdim_c){if(tdim == tdim_c)
            ierro = getintmetxi<gdim,tdim_c,ideg>(msh.coord, ent2poi[ieshl], msh.getBasis(), bary, lsmet[ishell]);
            geteigsym<gdim>(lsmet[ishell], eigval, eigvec);
            for(int ii = 0; ii < gdim; ii++) eigval[ii] = sqrt(eigval[ii]);
            eig2met<gdim>(eigval,eigvec,lsmet[ishell]);
          }CT_FOR1(tdim_c);
        }MSH_DIM_DEG1();
        METRIS_ENFORCE(ierro == 0);
      }

      // Take tangent vector (simply X1 - X2) and check t-continuity
      double du[3];
      for(int ii = 0; ii < msh.idim; ii++) du[ii] = msh.coord(ip1,ii) - msh.coord(ip2,ii);

      double mduref[3];
      if(msh.idim == 2) symXvec<2>(lsmet[0], du, mduref);
      else              symXvec<3>(lsmet[0], du, mduref);

      double nrmmdu_ref = msh.idim == 2 ? sqrt(getnrml2<2>(mduref)) 
                                        : sqrt(getnrml2<3>(mduref));
      //fmt::print("met0 = {} mduref = {}\n",dblAr1(msh.idim,du),dblAr1(nnmet,lsmet[0]),dblAr1(msh.idim,mduref));

      double max_err_mdu0   = -1;
      double max_err_nrmmdu0= -1;
      for(int ishell = 1; ishell < lsent.get_n(); ishell++){
        double mdutest[3];
        if(msh.idim == 2) symXvec<2>(lsmet[ishell], du, mdutest);
        else              symXvec<3>(lsmet[ishell], du, mdutest);

        //fmt::print("met[{}] = {} mdutest = {}\n",ishell,dblAr1(nnmet,lsmet[ishell]),dblAr1(msh.idim,mdutest));

        double err = msh.idim == 2 ? geterrl2<2>(mduref, mdutest) 
                                   : geterrl2<3>(mduref, mdutest);

        double nrmmdu_test = msh.idim == 2 ? sqrt(getnrml2<2>(mdutest)) 
                                           : sqrt(getnrml2<3>(mdutest));

        max_err_mdu0   = MAX(max_err_mdu0,   sqrt(err));
        max_err_nrmmdu0= MAX(max_err_nrmmdu0, abs(nrmmdu_test - nrmmdu_ref));

        BOOST_CHECK_CLOSE(nrmmdu_test, nrmmdu_ref, 1.0e-12);
      }
      //wait();
      
      max_err_mdu    = MAX(max_err_mdu,   max_err_mdu0);
      max_err_nrmmdu = MAX(max_err_nrmmdu, max_err_nrmmdu0);

    }

    fmt::print(" - Over all elements, max error on Mdu = {:.2}, on ||Mdu|| = {:.2}\n",
             max_err_mdu, max_err_nrmmdu);
  }

}// end boost test case
