//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

/*
Low level routine for "direct" P1 ball smoothing.
From each (facet, metric) pair, generate remaining vertex to be unit. Then average over ball.
Simplest possible approach.
*/

#ifndef __METRIS_LOW_SMOOBALLDIFF__
#define __METRIS_LOW_SMOOBALLDIFF__

#include "../Mesh/MeshFwd.hxx"
#include "../quality/quafun.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../types.hxx"

namespace Metris{

// same as smooballdirect but gradient descent (interior point)
template<class MetricFieldType, int idim, int ideg>
int smooballdiff(Mesh<MetricFieldType>& msh, int ipoin,
                 const intAr1 &lball,
                 double*__restrict__ qavg0, double*__restrict__ qmax0,
                 double*__restrict__ qavg1, double*__restrict__ qmax1,
                 QuaFun iquaf = QuaFun::Distortion);

// same as the above but for when the point is on a boundary
template<class MetricFieldType, int idim, int ideg>
int smooballdiff_boundary(Mesh<MetricFieldType>& msh, int ipoin, const int cadDim,
                          const intAr1 &lball,
                          double*__restrict__ qavg0, double*__restrict__ qmax0,
                          double*__restrict__ qavg1, double*__restrict__ qmax1,
                          QuaFun iquaf = QuaFun::Distortion);



template<class MFT, int gdim, int ideg>
double smooballdiff_fun(unsigned int nvar, const double *x,
                        double *grad, void *f_data);

template<class MFT>
struct smooballdiff_fun_data{
  smooballdiff_fun_data(Mesh<MFT> &msh_, const intAr1 &lball_, int ipoin_,
                        QuaFun iquaf_, double *xopt_) :
  msh(&msh_), lball(&lball_), ipoin(ipoin_), iquaf(iquaf_), qnrm0(0), qmax0(-1.0e30),
  iqset(false), xopt(xopt_), fopt(1.0e30) {}

  Mesh<MFT> *msh;
  const intAr1 *lball;
  int ipoin;
  QuaFun iquaf;

  double qnrm0, qmax0;
  bool iqset;

  // store best valid iterate
  double *xopt;
  double fopt;
};

// inorm <= infi norm , p > 0 L^p norm (over ball)
template<class MFT, int idim, int ideg>
int smooballdiff_luksan(Mesh<MFT>& msh, int ipoin,
                        const intAr1 &lball,
                        double*__restrict__ qnrm0, double*__restrict__ qmax0,
                        double*__restrict__ qnrm1, double*__restrict__ qmax1,
                        dblAr1 &work,
                        QuaFun iquaf);


// same as smooballdiff but smooths cavity instead of ball: the difference is that we need to put together the reconnected elements on the fly
template<class MetricFieldType, int idim, int ideg>
int smoocavdiff(Mesh<MetricFieldType>& msh, MshCavity& cav,
                   double& quaCav1, double& quaMaxCav1,
                   double& targetWeightCav1,
                   QuaFun iquaf, int ithread);

// Raw completed-P2 released-spoke objective and derivatives used by cavity
// smoothing. The surrounding cavity-growth workspace must have reserved its
// temporary cone element as the last full-dimensional entity. Derivatives are
// with respect to the physical insertion-point coordinate and include every
// affine midpoint dependence.
template<class MetricFieldType, int idim>
bool evaluate_released_p2_cavity_objective(
    Mesh<MetricFieldType>& msh,
    MshCavity& cav,
    QuaFun iquaf,
    int ithread,
    double& numerator,
    double* gradient,
    double* hessian = nullptr);

// boundary variant of smoocavdiff
template<class MetricFieldType, int idim, int ideg>
int smoocavdiff_boundary(Mesh<MetricFieldType>& msh, MshCavity& cav,
                         const int cadDim,
                         double& quaCav1, double& quaMaxCav1,
                         double& targetWeightCav1,
                         QuaFun iquaf, int ithread);
} // end namespace
#endif
