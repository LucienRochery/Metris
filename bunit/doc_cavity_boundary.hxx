//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "common_setup.hxx"

namespace Metris {

void doc_cavity_boundary(const int tdim, const intAr1& lcavel,
                         const intAr2& fac2poi, const intAr2& fac2fac,
                         const int ithread, int* tag, intAr2r& fac2tag,
                         intAr2& cavBoundary){
  // update tag
  const int currentTag = tag[ithread] + 1;
  tag[ithread] = currentTag;

  int nBndFacets = 0;

  // color cavity elements
  for (const int ele : lcavel) fac2tag(ithread,ele) = currentTag;

  // loop over cavity elements and check if neighbors are colored
  for (const int ele : lcavel){
    for (int iface = 0; iface < tdim+1; iface++){

      const int neighborEle = fac2fac(ele,iface);

      // neighbor not colored or invalid index, we have iface is a facet on cavity boundary
      if (neighborEle < 0 || fac2tag(ithread,neighborEle) != currentTag){
        nBndFacets++;
        cavBoundary.set_n(nBndFacets);
        // collect points on facet, use helper lnoed2
        for (int ipoin = 0; ipoin < tdim; ipoin++){
          cavBoundary(nBndFacets-1,ipoin) = fac2poi(ele, lnoed2[iface][ipoin]);
        }
      }
    }
  }
}

} // namespace Metris
