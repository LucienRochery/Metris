//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "common_setup.hxx"

namespace Metris {

void doc_structured_mesh(const double xmin, const double xmax,
                         const double ymin, const double ymax,
                         const int nx, const int ny,
                         intAr2& fac2poi,
                         intAr1& fac2ref,
                         intAr2& edg2poi,
                         intAr1& edg2ref,
                         dblAr2& coord,
                         intAr2& fac2fac,
                         intAr2& edg2edg,
                         intAr1& geonodes) {

  // points are numbered from the SW cornern and following the x direction
  // helper to get point ID given x and y indices
  auto get_ipoin = [&](const int ix, const int jy) -> int { return jy * nx + ix; };

  // we will subdivide domain into rectangles and then split them SW to NE
  // helper to get faces ID given rectangle indices. face at SE comes first
  auto get_iface = [&](const int ii, const int jj, const bool isSE) -> int {

      // return -1 if out of domain
      if (ii < 0 || ii >= nx-1 || jj < 0 || jj >= ny-1) return -1;

      // define rectangle ID
      const int irectangle = jj * (nx-1) + ii;

      // compute face ID
      if (isSE) return 2 * irectangle;
      else return 2 * irectangle + 1;
  };

  // geometric nodes (corners)
  geonodes.set_n(4);

  geonodes[0] = get_ipoin(0,0); // SW
  geonodes[1] = get_ipoin(nx-1,0); // SE
  geonodes[2] = get_ipoin(nx-1,ny-1); // NE
  geonodes[3] = get_ipoin(0,ny-1); // NW

  // fill coordinates array
  const int npoin = nx*ny;
  coord.allocate(npoin,2);
  coord.set_n(npoin);

  double dx = (xmax-xmin)/(nx-1);
  double dy = (ymax-ymin)/(ny-1);

  for (int ix = 0; ix < nx; ix++)
    for (int jy = 0; jy < ny; jy++) {
      int ipoin = get_ipoin(ix,jy);
      coord(ipoin,0) = xmin + dx*ix;
      coord(ipoin,1) = ymin + dy*jy;
    }

  // subdivide mesh in (nx-1) x (ny-1) rectangles and loop over them
  // to split them and fill fac2poi, edg2poi, fac2fac, edg2edg
  const int nfaces = (nx-1)*(ny-1)*2;
  fac2poi.allocate(nfaces,3);
  fac2poi.set_n(nfaces);
  fac2fac.allocate(nfaces,3);
  fac2fac.set_n(nfaces);
  for (int ii = 0; ii < nx-1; ii++)
    for (int jj = 0; jj < ny-1; jj++) {

      const int ifaceSE = get_iface(ii,jj,true);
      const int ifaceNW = get_iface(ii,jj,false);

      // obtain points in the current rectangle
      const int poiSW = get_ipoin(ii,jj);
      const int poiSE = get_ipoin(ii+1,jj);
      const int poiNE = get_ipoin(ii+1,jj+1);
      const int poiNW = get_ipoin(ii,jj+1);

      // fill fac2poi
      // SE face
      fac2poi(ifaceSE,0) = poiSW;
      fac2poi(ifaceSE,1) = poiSE;
      fac2poi(ifaceSE,2) = poiNE;

      // NW face
      fac2poi(ifaceNW,0) = poiSW;
      fac2poi(ifaceNW,1) = poiNE;
      fac2poi(ifaceNW,2) = poiNW;

      // fill neighbours, i-th entry is neighbor across edge opposite to vertex i

      // SE face
      fac2fac(ifaceSE,0) = get_iface(ii+1,jj,false);
      fac2fac(ifaceSE,1) = ifaceNW;
      fac2fac(ifaceSE,2) = get_iface(ii,jj-1,false);

      // NW face
      fac2fac(ifaceNW,0) = get_iface(ii,jj+1,true);
      fac2fac(ifaceNW,1) = get_iface(ii-1,jj,true);
      fac2fac(ifaceNW,2) = ifaceSE;
    }

  // edge connectivity is particularly simple
  const int nedge = 2*(nx-1) + 2*(ny-1);
  edg2poi.allocate(nedge,2);
  edg2poi.set_n(nedge);
  edg2edg.allocate(nedge,2);
  edg2edg.set_n(nedge);

  edg2ref.set_n(nedge);
  int bndRef = 0;

  int iedge = 0;

  // bottom boundary
  int jy = 0;
  for (int ix = 0; ix < nx-1; ix++, iedge++) {

    const int poi0 = get_ipoin(ix,jy);
    const int poi1 = get_ipoin(ix+1,jy);

    edg2poi(iedge,0) = poi0;
    edg2poi(iedge,1) = poi1;

    edg2ref[iedge] = bndRef;
  }
  bndRef++;

  // right boundary
  int ix = nx-1;
  for (int jy = 0; jy < ny-1; jy++, iedge++) {

    const int poi0 = get_ipoin(ix,jy);
    const int poi1 = get_ipoin(ix,jy+1);

    edg2poi(iedge,0) = poi0;
    edg2poi(iedge,1) = poi1;

    edg2ref[iedge] = bndRef;
  }
  bndRef++;

  // top boundary
  jy = ny-1;
  for (int ix = nx-1; ix > 0; ix--, iedge++) {

    const int poi0 = get_ipoin(ix,jy);
    const int poi1 = get_ipoin(ix-1,jy);

    edg2poi(iedge,0) = poi0;
    edg2poi(iedge,1) = poi1;

    edg2ref[iedge] = bndRef;
  }
  bndRef++;

  // left boundary
  ix = 0;
  for (int jy = ny-1; jy > 0; jy--, iedge++) {

    const int poi0 = get_ipoin(ix,jy);
    const int poi1 = get_ipoin(ix,jy-1);

    edg2poi(iedge,0) = poi0;
    edg2poi(iedge,1) = poi1;

    edg2ref[iedge] = bndRef;
  }

  // fill edg2edg connectivity
  for (int iedg = 0; iedg < nedge; iedg++) {

    edg2edg(iedg,0) = iedg+1;
    edg2edg(iedg,1) = iedg-1;
  }

  // correct ends
  edg2edg(0,1) = nedge-1;
  edg2edg(nedge-1,0) = 0;

  // fill fac2ref
  fac2ref.set_n(nfaces);
  const int faceRef = 0;
  for (int iface = 0; iface < nfaces; iface++) fac2ref[iface] = faceRef;
}

} // namespace Metris
