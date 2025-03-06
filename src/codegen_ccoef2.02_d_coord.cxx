//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "codegen_ccoef_d.hxx"
#include "types.hxx"


namespace Metris{


template<int ideg> void d_ccoef_genbez2([[maybe_unused]]const intAr2&__restrict__ fac2poi, [[maybe_unused]] const dblAr2&__restrict__ coord,
                                        [[maybe_unused]] int ielem, [[maybe_unused]] int icoor,
                                        [[maybe_unused]] dblAr2&__restrict__ d_ccoef){}
template<> void d_ccoef_genbez2<2>(const intAr2&__restrict__ fac2poi, const dblAr2&__restrict__ coord,
                                   int ielem, int icoor,
                                   dblAr2&__restrict__ d_ccoef){
  METRIS_ASSERT(icoor == 0 || icoor == 1);
  const int sg = 1 - 2*(icoor%2);

  d_ccoef(0,0) = sg*4*(coord(fac2poi(ielem,5),1-icoor) - coord(fac2poi(ielem,4),1-icoor));

  d_ccoef(0,4) = sg*4*(coord(fac2poi(ielem,0),1-icoor) - coord(fac2poi(ielem,5),1-icoor));

  d_ccoef(0,5) = sg*4*(coord(fac2poi(ielem,4),1-icoor) - coord(fac2poi(ielem,0),1-icoor));

  d_ccoef(1,1) = sg*4*(coord(fac2poi(ielem,3),1-icoor) - coord(fac2poi(ielem,5),1-icoor));

  d_ccoef(1,3) = sg*4*(coord(fac2poi(ielem,5),1-icoor) - coord(fac2poi(ielem,1),1-icoor));

  d_ccoef(1,5) = sg*4*(coord(fac2poi(ielem,1),1-icoor) - coord(fac2poi(ielem,3),1-icoor));

  d_ccoef(2,2) = sg*4*(coord(fac2poi(ielem,4),1-icoor) - coord(fac2poi(ielem,3),1-icoor));

  d_ccoef(2,3) = sg*4*(coord(fac2poi(ielem,2),1-icoor) - coord(fac2poi(ielem,4),1-icoor));

  d_ccoef(2,4) = sg*4*(coord(fac2poi(ielem,3),1-icoor) - coord(fac2poi(ielem,2),1-icoor));

  d_ccoef(3,1) = sg*2*(coord(fac2poi(ielem,2),1-icoor) - coord(fac2poi(ielem,4),1-icoor));

  d_ccoef(3,2) = sg*2*(coord(fac2poi(ielem,5),1-icoor) - coord(fac2poi(ielem,1),1-icoor));

  d_ccoef(3,3) = sg*2*(coord(fac2poi(ielem,4),1-icoor) - coord(fac2poi(ielem,5),1-icoor));

  d_ccoef(3,4) = sg*2*(coord(fac2poi(ielem,1),1-icoor) - coord(fac2poi(ielem,3),1-icoor));

  d_ccoef(3,5) = sg*2*(coord(fac2poi(ielem,3),1-icoor) - coord(fac2poi(ielem,2),1-icoor));

  d_ccoef(4,0) = sg*2*(coord(fac2poi(ielem,3),1-icoor) - coord(fac2poi(ielem,2),1-icoor));

  d_ccoef(4,2) = sg*2*(coord(fac2poi(ielem,0),1-icoor) - coord(fac2poi(ielem,5),1-icoor));

  d_ccoef(4,3) = sg*2*(coord(fac2poi(ielem,4),1-icoor) - coord(fac2poi(ielem,0),1-icoor));

  d_ccoef(4,4) = sg*2*(coord(fac2poi(ielem,5),1-icoor) - coord(fac2poi(ielem,3),1-icoor));

  d_ccoef(4,5) = sg*2*(coord(fac2poi(ielem,2),1-icoor) - coord(fac2poi(ielem,4),1-icoor));

  d_ccoef(5,0) = sg*2*(coord(fac2poi(ielem,1),1-icoor) - coord(fac2poi(ielem,3),1-icoor));

  d_ccoef(5,1) = sg*2*(coord(fac2poi(ielem,4),1-icoor) - coord(fac2poi(ielem,0),1-icoor));

  d_ccoef(5,3) = sg*2*(coord(fac2poi(ielem,0),1-icoor) - coord(fac2poi(ielem,5),1-icoor));

  d_ccoef(5,4) = sg*2*(coord(fac2poi(ielem,5),1-icoor) - coord(fac2poi(ielem,1),1-icoor));

  d_ccoef(5,5) = sg*2*(coord(fac2poi(ielem,3),1-icoor) - coord(fac2poi(ielem,4),1-icoor));
}


} // End namespace
