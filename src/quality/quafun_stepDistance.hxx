
//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_QUAFUN_STEPDISTANCE__
#define __METRIS_QUAFUN_STEPDISTANCE__

#include "../Mesh/MeshFwd.hxx"

#include <limits>

namespace Metris{

enum class AsDeg;
enum class FEBasis;
enum class DifVar;

// Finite sentinel returned for a ShapeVolume trial whose metric-space matrix
// is too singular to evaluate robustly.  Integrated quality routines detect
// it and return immediately, so cavity logic rejects the operation without
// evaluating geometric weights or derivatives on the bad element.
inline constexpr double step_distance_shape_volume_rejection_quality = 1.e100;

// Running value of the mesh-wide CavityTargetAverage objective. Despite the
// historical name, this variant is the arithmetic mean of the unweighted
// reference-space elemental StepDistance integrals. A local operation is
// accepted only when the resulting mesh-wide mean strictly improves.
struct StepDistanceObjectiveState{
  double numerator = 0.;
  int element_count = 0;
  // Compatibility mirror used by cavity/smoothing plumbing. Each element now
  // has unit objective weight, so target_weight == element_count.
  double target_weight = 0.;
  double best_objective = std::numeric_limits<double>::infinity();
  double* best_objective_storage = nullptr;
  double cavity_global_relative_tolerance = 0.;
  double cavity_global_gain_fraction = 0.;

  double value() const;
  double region_value(double region_numerator,
                      int region_element_count,
                      double region_unit_weight) const;
  double replaced_value(double old_region_numerator,
                        int old_region_element_count,
                        double old_region_unit_weight,
                        double new_region_numerator,
                        int new_region_element_count,
                        double new_region_unit_weight) const;
  bool accepts_replacement(double old_region_numerator,
                           int old_region_element_count,
                           double old_region_unit_weight,
                           double new_region_numerator,
                           int new_region_element_count,
                           double new_region_unit_weight,
                           double relative_improvement = 1.e-7) const;
  void replace(double old_region_numerator,
               int old_region_element_count,
               double old_region_unit_weight,
               double new_region_numerator,
               int new_region_element_count,
               double new_region_unit_weight);
};

bool objective_strictly_improves(double value_new,
                                 double value_old,
                                 double relative_improvement = 1.e-7);

// Compatibility entry point for existing cavity/smoothing call sites. Local
// values, best-so-far state, and the former tolerance budget are deliberately
// ignored: acceptance is strict mesh-wide improvement only.
bool cavity_target_average_global_filter_accepts(
    double local_old,
    double local_new,
    double global_current,
    double global_trial,
    double global_best,
    double old_region_unit_weight,
    double new_region_unit_weight,
    double global_unit_weight,
    double global_relative_tolerance,
    double global_gain_fraction);

// Convert accumulated elemental StepDistance data to a regional objective.
// For cavity_target_average, the objective is the arithmetic mean over the
// supplied number of elements.
double step_distance_region_objective(double elemental_sum,
                                      double element_count,
                                      bool cavity_target_average);

// Convert one elemental value into its additive regional contribution.
double step_distance_region_contribution(double element_value,
                                         double unit_weight,
                                         bool cavity_target_average);

// Evaluate the CavityTargetAverage objective after replacing one triangulated
// subregion while leaving the rest of the mesh unchanged.
double step_distance_replaced_region_objective(
    double global_elemental_sum,
    double old_region_elemental_sum,
    double new_region_elemental_sum,
    double global_element_count,
    double old_region_element_count,
    double new_region_element_count);

// Pointwise
template <class MetricFieldType, int gdim, int tdim,
          typename ftype = double>
ftype quafun_stepDistance(Mesh<MetricFieldType> &msh,
                        AsDeg asdmsh, AsDeg asdmet,
                        const int*__restrict__ ent2poi,
                        const double*__restrict__ bary,
                        const double*__restrict__ met);

// Differentiated w.r.t. ielem's ivar-th control point/node.
template <class MetricFieldType, int gdim, int tdim,
           typename ftype = double>
ftype d_quafun_stepDistance(Mesh<MetricFieldType> &msh,
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

#endif // __METRIS_QUAFUN_STEPDISTANCE__
