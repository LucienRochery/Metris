//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_DEFAULTS__
#define __METRIS_DEFAULTS__



namespace Metris{

namespace Defaults{
  // Tolerances
  // - jtol = normalized (by P1 measure) Jacobian control coefficient minimum
  // - vtol: e.g. 2D. Volume is c0 det (P2-P1, P3-P1)
  // and det(u,v) <= ||u|| ||v|| 
  // For symmetry, we consider not the factor l_1 l_2 but (l_1l_2l_3)^{2/3}
  // with l_i the edges. In 3D: (l_1l_2l_3l_4l_5l_6)^(1/2) 
  // Hence an element passes the volume test iff: 
  // det(li lj) >= vtol * (pi_i ||l_i||)^(2/(idim+1))   (idim / nnmet)
  // The 1/idim! factor is not included. 
  const double jtol = 1.0e-6;
  const double vtol = 1.0e-9;
  // Absolute tolerance on edge length 
  const double ltol = 1.0e-12;

  const double geo_lentolfac = 1.01;
  const double geo_abstoledg = 1.0e-12;

  // See MetrisParameters 
  const int qopt_pnorm = 2;
  const int qopt_power = -1;
  const int qopt_niter = 5;
  const int qopt_smoo_niter = 100;

  const int qopt_swap_niter = 100;
  const int qopt_swap_pnorm = qopt_pnorm;
  const double qopt_swap_thres = 0.0;

  // What factor to multiply sizes by when realloc necessary.
  const double mem_growfac = 1.5;
}// namespace Defaults

}// namespace Metris

#endif