
//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_QUAFUN_SIZESHAPE__
#define __METRIS_QUAFUN_SIZESHAPE__
/*
  Pointwise SizeShape objective
    psi_SS,p(A) = (q_SS(A) - 1)^p,
  where q_SS is the scale-aware distortion measure with ideal value one.
*/

#include "../Mesh/MeshFwd.hxx"

namespace Metris{

enum class AsDeg;
enum class FEBasis;
enum class DifVar;

/* ---- Complete SizeShape objective from tr / det ---- */
// Pointwise
template <class MetricFieldType, int gdim, int tdim,
          typename ftype = double>
ftype quafun_sizeshape(Mesh<MetricFieldType> &msh,
                        AsDeg asdmsh, AsDeg asdmet,
                        const int*__restrict__ ent2poi,
                        const double*__restrict__ bary,
                        const double*__restrict__ met);

// Differentiated w.r.t. ielem's ivar-th control point/node.
template <class MetricFieldType, int gdim, int tdim,
           typename ftype = double>
ftype d_quafun_sizeshape(Mesh<MetricFieldType> &msh,
                          AsDeg asdmsh, AsDeg asdmet,
                          const int* ent2poi,
                          const double*__restrict__ bary,
                          const double*__restrict__ met_,
                          int ivar,
                          FEBasis dofbas,
                          DifVar idifmet,
                          ftype*__restrict__ dquael,
                          ftype*__restrict__ hquael);


} // End namespace

#endif // __METRIS_QUAFUN_SIZESHAPE__
