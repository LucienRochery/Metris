//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

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
    for (int jj = 0; jj < fac2poi.size2(); jj++){
      if (fac2poi(iele,jj) == ipoin) return true;
    }
    return false;
  };

  // loop visiting neighbors of elements in ball
  for (int ii = 0; ii < ball.get_n(); ii++){

    // fetch element from ball
    const int ele = ball[ii];

    // traverse its neighbours
    for (int jj = 0; jj < fac2fac.size2(); jj++){

      const int neighbourElem = fac2fac(ele,jj);

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

} // namespace Metris
