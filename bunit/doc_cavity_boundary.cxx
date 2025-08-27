//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_doc_cavity_boundary

#include "common_setup.hxx"
#include "doc_cavity_boundary.hxx"
#include "doc_ball.hxx"

using namespace Metris;
typedef MetricFieldAnalytical MFT;

BOOST_AUTO_TEST_CASE(test_doc_cavity_boundary)
{

  const int ithread = 0;

  // bool is whether straight
  std::vector<std::string> meshes = {
   METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
  //,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k" // if 3D is implemented
  };

  for(std::string mesh : meshes){

    // dummy analytical metric
    cargHandler arg("-in " + mesh + " -anamet 1 -verb 0");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
    msh.cleanup(); // to avoid dead elements in case the file has some

    // It's best to declare these objects outside loops as they
    // allocate memory. They can be reset using .set_n(0); this should
    // ideally be done inside the routines that use them.
    intAr1 lcavel(100);
    intAr2 lcavbdry(100,2);
    for(int iface = 0; iface < msh.nface; iface++){
      for(int iver = 0; iver < 3; iver++){

        int ipoin = msh.fac2poi(iface,iver);

        // gather ball of ipoin as ad-hoc cavity

        int iele0;
        if (msh.poi2ent(ipoin,1) == 1) iele0 = msh.edg2fac[msh.poi2ent(ipoin,0)];
        else if (msh.poi2ent(ipoin,1) == 2) iele0 = msh.poi2ent(ipoin,0);
        else METRIS_THROW_MSG(TODOExcept(), "Implement ball for 3D");

        intAr1 ipoinBall = doc_ball(ipoin,msh.fac2poi,msh.fac2fac,ithread,msh.tag,msh.fac2tag,iele0);

        // Get the cavity boundary
        intAr2 ballBoundary = doc_cavity_boundary(2,ipoinBall,msh.fac2poi,msh.fac2fac,ithread,msh.tag,msh.fac2tag);

        // Check it has as many elements as the cavity (only true in case of the ball and for interior ipoin)
        if (msh.poi2bpo[ipoin] < 0) BOOST_TEST(ipoinBall.get_n() == ballBoundary.get_n());

        // TODO
        // // Check boundary faces are only contained in one cavity
        // // element.
        // // Use getedgfac and getfactet
        // for(int icavbdry = 0; icavbdry < lcavbdry.get_n(); icavbdry++){
        //   int ncav_has_bdry = 0;
        //   for(int icavel : lcavel){
        //     METRIS_THROW_MSG(TODOExcept(), "Check if icavel has the facet");
        //   }
        // }

      }
    }
    std::cout << "msh.npoin = " << msh.npoin << std::endl;
  }
}