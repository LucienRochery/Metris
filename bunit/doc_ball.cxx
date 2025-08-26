//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_doc_ball

#include "common_setup.hxx"

namespace Metris {

intAr1 doc_ball(const int ipoin, const intAr2& fac2poi, const intAr2& fac2fac,
                const int ithread, int* tag, intAr2r fac2tag,
                const int iele0) {

  // update tag
  const int currentTag = tag[ithread] + 1;
  tag[ithread] = currentTag;

  // initialize array for the ball and add iele0
  intAr1 ball;
  ball.stack(iele0);
  fac2tag(ithread,iele0) = currentTag;

  // little helper to identify if a given element has a given point
  auto eleHasPoi = [&](const int iele, const int ipoin){
    for (int j = 0; j < fac2poi.size2(); j++){
      if (fac2poi(iele,j) == ipoin) return true;
    }
    return false;
  };

  // loop visiting neighbors of elements in ball
  for (int i = 0; i < ball.get_n(); i++){

    // fetch element from ball
    const int ele = ball[i];

    // traverse its neighbours
    for (int j = 0; j < fac2fac.size2(); j++){

      const int neighbourElem = fac2fac(ele,j);

      // not really an element index;
      // e.g. ele might be a boundary element and might not contain valid elem index in all the fac2fac stride
      if (neighbourElem < 0) continue;

      // neighbour already visited
      if (fac2tag(ithread,neighbourElem) == currentTag) continue;

      // add neighbour to the visited elements
      fac2tag(ithread,neighbourElem) = currentTag;

      // neighbour contains ipoin, add it to ball
      if (eleHasPoi(neighbourElem,ipoin)) ball.stack(neighbourElem);
    }
  }

  return ball;
}

}

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

    int failedPoin = 0;
    for (int ipoin = 0; ipoin < msh.npoin; ipoin++){

      // First obtain ball of elements with Metris::ball

      // Initialize ball of faces (triangles) array
      intAr1 lbfac(10);
      // Pass the other two empty so ball does not fill them
      intAr1 lbedg;
      intAr1 lbtet;
      int iopen;
      bool append = false;
      ball(msh,ipoin,lbedg,lbfac,lbtet,&iopen,append,ithread);

      // Now obtain ball of elements using doc_ball

      // cheating a bit, use lbfac to obtain seed
      int iele0 = lbfac[0];

      intAr1 ipoinBall = doc_ball(ipoin,msh.fac2poi,msh.fac2fac,ithread,msh.tag,msh.fac2tag,iele0);

      // check that balls match
      for (const int ele1 : ipoinBall){
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

