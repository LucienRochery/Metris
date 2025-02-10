//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php
#ifndef __METRIS_LOW_LENEDG__
#define __METRIS_LOW_LENEDG__


namespace Metris{

template <class MFT> class MeshMetric;

// -----------------------------------------------------------------------------
template<int gdim>
double getlenedg(const double x1[], const double x2[], const double metl[]);
// -----------------------------------------------------------------------------
template<int gdim>
double getlenedg(const double dx[], const double metl[]);
// -----------------------------------------------------------------------------
template<int gdim>
double getlenedgsq(const double x1[], const double x2[], const double metl[]);
// -----------------------------------------------------------------------------
template<int gdim>
double getlenedgsq(const double dx[],  const double metl[]);
// -----------------------------------------------------------------------------
// The metric is in log format. Relative tolerance 
template<int gdim>
double getlenedg_log(const double dx[], const double metl[], int miter = 100, double tol = 1.0e-6);



template<class MetricFieldType, int gdim, int ideg>
double getlenedg_geosz(MeshMetric<MetricFieldType> &msh,
                       int ientt, int tdimn, int iedg);
template<class MetricFieldType, int gdim, int ideg>
double getlenedg_geosz(MeshMetric<MetricFieldType> &msh,
                       int ientt, int tdimn, int iedg, double *sz);
template<class MetricFieldType, int gdim, int ideg>
double getlenedg_geosz(MeshMetric<MetricFieldType> &msh,
                       int *edg2pol, double *sz);

// This one assumes lpoi of size ideg + 1
//template<int gdim, int ideg>
//double getlenedg_quad(const int lpoi[], const dblAr2 &coord, const dblAr2 &met, int nquad);

// This one is more recent. Since metric refactoring, using gdim etc. 
template<class MetricFieldType, int gdim, int ideg>
double getlenedg_quad(MeshMetric<MetricFieldType> &msh,
                      int ientt, int tdimn, int iedg, int nquad);


}// end namespace
#endif