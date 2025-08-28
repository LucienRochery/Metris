//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_doc_ball

#include "common_setup.hxx"
#include "doc_ball.hxx"

using namespace Metris;
typedef MetricFieldAnalytical MFT;

BOOST_AUTO_TEST_CASE(test_doc_ball)
{
  // bool is whether straight
  std::vector<std::string> meshes = {
   METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
  //,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k" // if 3D is implemented
  };

  const int ithread = 0;
  for(std::string mesh : meshes){

    // dummy analytical metric
    cargHandler arg("-in " + mesh + " -anamet 1 -verb 0");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
    msh.cleanup(); // to avoid dead elements in case the file has some

    // Placeholders to use with Metris::ball

    // Initialize ball of faces (triangles) array
    intAr1 lbfac(10);
    // Pass the other two empty so ball does not fill them
    intAr1 lbedg;
    intAr1 lbtet;
    int iopen;
    bool append = false;



    int failedPoin = 0;
    for (int ipoin = 0; ipoin < msh.npoin; ipoin++){

      ball(msh,ipoin,lbedg,lbfac,lbtet,&iopen,append,ithread);

      // Now obtain ball of elements using doc_ball

      // seed for doc_ball
      int iele0;
      if (msh.poi2ent(ipoin,1) == 1) iele0 = msh.edg2fac[msh.poi2ent(ipoin,0)];
      else if (msh.poi2ent(ipoin,1) == 2) iele0 = msh.poi2ent(ipoin,0);
      else METRIS_THROW_MSG(TODOExcept(), "Implement ball for 3D");
// Placeholder to use with doc_ball
    intAr1 ball_doc(100);
      doc_ball(ipoin,msh.fac2poi,msh.fac2fac,ithread,msh.tag,msh.fac2tag,iele0,ball_doc);

      // check that balls match
      for (const int ele1 : ball_doc){
        bool found = false;
        for (const int ele2 : lbfac){
          if (ele1 == ele2){
            found = true;
            break;
          }
        }
        if (!found) failedPoin++;
      }

    }
    BOOST_TEST(failedPoin == 0);
  }
}

