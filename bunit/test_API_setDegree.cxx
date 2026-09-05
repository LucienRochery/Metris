//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_API_setDegree

#include "common_setup.hxx"

using namespace Metris;

BOOST_AUTO_TEST_CASE(quadratic_to_linear_preserves_connectivity)
{
  constexpr int npoin = 9;
  constexpr int nedge = 4;
  constexpr int nface = 2;

  MetrisAPI api(2, 2, 4, 0, 0, npoin, false, nedge, nface, 0,
                FEBasis::Lagrange, FEBasis::Lagrange, MetSpace::Exp);

  for(int ipoin = 0; ipoin < npoin; ipoin++){
    const double coord[2] = {static_cast<double>(ipoin),
                             -static_cast<double>(ipoin)};
    api.setCoord(ipoin, coord);
  }

  for(int icorn = 0; icorn < 4; icorn++) api.setCorner(icorn, icorn);

  const int edges[nedge][3] = {
    {0, 1, 4},
    {1, 2, 5},
    {2, 3, 6},
    {3, 0, 7}
  };
  for(int iedge = 0; iedge < nedge; iedge++){
    api.setElement(1, iedge, edges[iedge], 1);
  }

  const int faces[nface][6] = {
    {0, 1, 2, 5, 8, 4},
    {0, 2, 3, 6, 7, 8}
  };
  for(int iface = 0; iface < nface; iface++){
    api.setElement(2, iface, faces[iface], 1);
  }

  api.setDegree(1);

  int idim, ideg, ncorn, npoin_out, nedge_out, nface_out, nelem;
  bool imet;
  FEBasis mshbasis, metbasis;
  api.getConstants(&idim, &ideg, &ncorn, &npoin_out, &imet,
                   &nedge_out, &nface_out, &nelem, &mshbasis, &metbasis);

  BOOST_TEST(idim == 2);
  BOOST_TEST(ideg == 1);
  BOOST_TEST(npoin_out == 4);
  BOOST_TEST(nedge_out == nedge);
  BOOST_TEST(nface_out == nface);

  for(int iedge = 0; iedge < nedge; iedge++){
    int nodes[2], iref;
    api.getElement(1, iedge, nodes, &iref);
    BOOST_TEST(nodes[0] == edges[iedge][0]);
    BOOST_TEST(nodes[1] == edges[iedge][1]);
    BOOST_TEST(iref == 1);
  }

  for(int iface = 0; iface < nface; iface++){
    int nodes[3], iref;
    api.getElement(2, iface, nodes, &iref);
    for(int inode = 0; inode < 3; inode++){
      BOOST_TEST(nodes[inode] == faces[iface][inode]);
    }
    BOOST_TEST(iref == 1);
  }

  for(int ipoin = 0; ipoin < npoin_out; ipoin++){
    double coord[2];
    api.getCoord(ipoin, coord);
    BOOST_TEST(coord[0] == static_cast<double>(ipoin));
    BOOST_TEST(coord[1] == -static_cast<double>(ipoin));
  }

  for(int icorn = 0; icorn < ncorn; icorn++){
    int ipoin;
    api.getCorner(icorn, &ipoin);
    BOOST_TEST(ipoin == icorn);
  }
}
