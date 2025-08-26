//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "types_arrays.hxx"

intAr1 doc_ball(const int ipoin, const intAr2& fac2poi, const intAr2& fac2fac,
                const int ithread, intAr1 tag, intAr2r fac2tag,
                const int iele0) {

  const int currentTag = tag[ithread] + 1;
  tag[ithread] = currentTag;

  // initialize array for the ball and add iele0
  intAr1 ball;
  ball.stack(iele0);
  fac2tag(iele0) = currentTag;

  // little helper to identify if a given element has a given point
  auto eleHasPoi = [&](const int iele, const int ipoi){
    for (int j = 0; j < fac2poi.size2(); j++){
      if (fac2poi(iele,j) == ipoi) return true;
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
      if (fac2tag(neighbourElem) == currentTag) continue;

      // add neighbour to the visited elements
      fac2tag(neighbourElem) = currentTag;

      // neighbour contains ipoi, add it to ball
      if (eleHasPoi(neighbourElem,ipoi)) ball.stack(neighbourElem);
    }
  }

  return ball;
}
