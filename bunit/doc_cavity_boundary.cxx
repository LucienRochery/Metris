//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE doc_cavity_boundary

#include "common_setup.hxx"

using namespace Metris;


BOOST_AUTO_TEST_CASE(doc_cavity_boundary) 
{ 

  // bool is whether straight
  std::vector<std::string> meshes = {
   METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
  //,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k" // if 3D is implemented
  };

  for(std::string mesh : meshes){

    // dummy analytical metric
    cargHandler arg("-in " + mesh + " -anamet 1 -verb 0");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *(run.msh_g);
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
        METRIS_THROW_MSG(TODOExcept(), "Implement ball or call ball() in low_topo.hxx");

        // Get the cavity boundary
        METRIS_THROW_MSG(TODOExcept(), "Get cavity boundary");

        // Check it has as many elements as the cavity (only true in case of the ball)

        // Check boundary faces are only contained in one cavity 
        // element. 
        // Use getedgfac and getfactet 
        for(int icavbdry = 0; icavbdry < lcavbdry.get_n(); icavbdry++){
          int ncav_has_bdry = 0;
          for(int icavel : lcavel){
            METRIS_THROW_MSG(TODOExcept(), "Check if icavel has the facet");
          }
        }

      }
    }
  
  }

}