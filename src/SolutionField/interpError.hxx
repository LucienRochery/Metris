//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_INTERPERROR__
#define __METRIS_INTERPERROR__

#include "SolutionField.hxx"

#include <initializer_list>

namespace Metris{

// Evaluate interpolation error over an element. (tdim = msh.get_tdim())
// ideg is mesh degree maximum 2
// pdeg is interpolation degree maximum 2 
// pnorm is interp err pnorm maximum 2 
// Derivatives derr can be empty {} or up to two double* sized resp idim, (idim*(idim+1))/2
template<int idim, int ideg, int pdeg, int pnorm, bool iexact = false>
double interpErr(const SolutionFieldAnalytical &sol, int ielem,
                 int idiff = -1, std::initializer_list<double*> derr = {});

template<int idim, int ideg, int pdeg, int pnorm, bool iexact>
double interpErrBall(const SolutionFieldAnalytical &sol, 
                     const intAr1& lball, const intAr1& lnode,
                     std::initializer_list<double*> derr);

template<int idim, int ideg, int pdeg, int pnorm, bool iexact = false>
double interpErrGlo(const SolutionFieldAnalytical &sol);


} //namespace
#endif