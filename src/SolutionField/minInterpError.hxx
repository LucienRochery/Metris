//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MININTERPERROR__
#define __METRIS_MININTERPERROR__

#include "SolutionField.hxx"

#include <initializer_list>

namespace Metris{
template <class MFT> class Mesh;

template <class MFT>
void minimizeInterpErrglo(Mesh<MFT> &msh, const SolutionFieldAnalytical &sol, 
                          int pdeg, int pnorm, int ithrd1 = 0, int ithrd2 = 1);



// Aux
template<class MFT, int idim, int ideg, int pdeg>
void minimizeInterpErrglo0(Mesh<MFT> &msh, const SolutionFieldAnalytical &sol, 
                           int pnorm, int ithrd1, int ithrd2);


// For nlopt
template<class MFT>
struct interpErrglo_constraint_data{
  interpErrglo_constraint_data(Mesh<MFT> *msh, int pnorm, int ithrd1, int ithrd2){
    this->msh = msh;
    this->pnorm = pnorm;
    this->ithrd1 = ithrd1;
    this->ithrd2 = ithrd2;
  }
  Mesh<MFT>* msh;
  const intAr1 *ldof;
  int pnorm;
  int ithrd1; 
  int ithrd2;
};






template<class MFT, int idim, int ideg, int pdeg>
int minimizeInterpErrloc(Mesh<MFT> &msh, const SolutionFieldAnalytical &sol, 
                         int pnorm, intAr1& lball, intAr1& lnode,
                         double *errLp0, double *errLp1);
} //namespace
#endif