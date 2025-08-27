//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "common_setup.hxx"

namespace Metris {

intAr2 doc_cavity_boundary(const int tdim, const intAr1& lcavel,
                           const intAr2& fac2poi, const intAr2& fac2fac,
                           const int ithread, int* tag, intAr2r fac2tag){
  // update tag
  const int currentTag = tag[ithread] + 1;
  tag[ithread] = currentTag;

  // initialize cavity boundary array
  intAr2 cavBoundary(10,tdim);
  int nFacets = 0;

  // color cavity elements
  for (const int ele : lcavel) fac2tag(ithread,ele) = currentTag;

  // loop over cavity elements and check if neighbors are colored
  for (const int ele : lcavel){
    for (int iface = 0; iface < tdim+1; iface++){

      const int neighborEle = fac2fac(ele,iface);

      // neighbor not colored or invalid index, we have iface is a facet on cavity boundary
      if (neighborEle < 0 || fac2tag(ithread,neighborEle) != currentTag){

        // collect points on facet, use helper lnoed2
        for (int ipoin = 0; ipoin < tdim; ipoin++){
          cavBoundary(nFacets,ipoin) = fac2poi(ele, lnoed2[iface][ipoin]);
        }
        // done with current facet, increment facet number
        nFacets++;
      }
    }
  }

  // resize used slots in cavBoundary to the number of boundary facets
  cavBoundary.set_n(nFacets);

  return cavBoundary;
}

} // namespace Metris
