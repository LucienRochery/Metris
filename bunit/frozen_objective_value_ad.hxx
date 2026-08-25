// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#ifndef __METRIS_BUNIT_FROZEN_OBJECTIVE_VALUE_AD_HXX__
#define __METRIS_BUNIT_FROZEN_OBJECTIVE_VALUE_AD_HXX__

#include "linalg/eigen.hxx"
#include "libs/SANS/Surreal/SurrealS.h"

#include <cmath>

namespace Metris
{

namespace FrozenObjectiveValueAD
{

template<int nvar>
double primal_value(const SANS::SurrealS<nvar,double> &value)
{
  return value.value();
}

template<int dimension, typename T>
T full_determinant(const T *matrix)
{
  static_assert(dimension == 2 || dimension == 3);
  if constexpr(dimension == 2){
    return matrix[0]*matrix[3] - matrix[1]*matrix[2];
  }else{
    return matrix[0]*(matrix[4]*matrix[8] - matrix[5]*matrix[7])
         - matrix[1]*(matrix[3]*matrix[8] - matrix[5]*matrix[6])
         + matrix[2]*(matrix[3]*matrix[7] - matrix[4]*matrix[6]);
  }
}

template<int dimension>
int packed_symmetric_index(int row, int column)
{
  static_assert(dimension == 2 || dimension == 3);
  if(row > column){
    const int temporary = row;
    row = column;
    column = temporary;
  }
  if constexpr(dimension == 2){
    if(row == 0 && column == 0) return 0;
    if(row == 0 && column == 1) return 1;
    return 2;
  }else{
    if(row == 0 && column == 0) return 0;
    if(row == 0 && column == 1) return 1;
    if(row == 1 && column == 1) return 2;
    if(row == 0 && column == 2) return 3;
    if(row == 1 && column == 2) return 4;
    return 5;
  }
}

template<int dimension, typename T>
void form_metric_gram_matrix(
    const T *regular_jacobian_transpose,
    const double *metric,
    T *gram_matrix)
{
  for(int row = 0; row < dimension; row++){
    for(int column = 0; column < dimension; column++){
      T entry = T(0.0);
      for(int physical_row = 0;
          physical_row < dimension; physical_row++){
        for(int physical_column = 0;
            physical_column < dimension; physical_column++){
          entry += regular_jacobian_transpose[
                       row*dimension + physical_row]
                 * T(metric[packed_symmetric_index<dimension>(
                       physical_row,physical_column)])
                 * regular_jacobian_transpose[
                       column*dimension + physical_column];
        }
      }
      gram_matrix[row*dimension + column] = entry;
    }
  }
}

template<int dimension>
void seed_regular_jacobian(
    const double *regular_jacobian_transpose,
    const double *regular_basis_gradient,
    SANS::SurrealS<dimension,double> *differentiated_jacobian)
{
  // For an active physical point X_q and regular-reference basis gradient g,
  // d(J^T)_{ia}/d(X_q)_k = g_i delta_{ak}. Metric and quadrature theta do not
  // receive derivative seeds, which is exactly the frozen-sample contract.
  for(int row = 0; row < dimension; row++){
    for(int component = 0; component < dimension; component++){
      const int entry = row*dimension + component;
      differentiated_jacobian[entry].value()
          = regular_jacobian_transpose[entry];
      for(int derivative = 0; derivative < dimension; derivative++){
        differentiated_jacobian[entry].deriv(derivative)
            = component == derivative
            ? regular_basis_gradient[row]
            : 0.0;
      }
    }
  }
}

template<int dimension>
SANS::SurrealS<dimension,double> sizeshape_pointwise_value(
    const double *regular_jacobian_transpose,
    const double *metric,
    const double *regular_basis_gradient,
    const double objective_power)
{
  using S = SANS::SurrealS<dimension,double>;
  S jacobian[dimension*dimension];
  S gram_matrix[dimension*dimension];
  seed_regular_jacobian<dimension>(
      regular_jacobian_transpose,regular_basis_gradient,jacobian);
  form_metric_gram_matrix<dimension,S>(jacobian,metric,gram_matrix);

  S trace = S(0.0);
  for(int diagonal = 0; diagonal < dimension; diagonal++){
    trace += gram_matrix[diagonal*dimension + diagonal];
  }
  const S determinant = full_determinant<dimension,S>(gram_matrix);
  const double dimension_power = dimension == 2 ? 4.0 : 27.0;
  S trace_power;
  if constexpr(dimension == 2){
    trace_power = trace*trace;
  }else{
    trace_power = trace*trace*trace;
  }
  const S size_shape_error
      = trace_power*(S(1.0) + S(1.0)/(determinant*determinant))
       /S(2.0*dimension_power) - S(1.0);
  if(objective_power == 1.0) return size_shape_error;
  return pow(size_shape_error,objective_power);
}

template<int dimension, typename T>
void invert_full_matrix(const T *matrix, T *inverse)
{
  static_assert(dimension == 2 || dimension == 3);
  const T determinant = full_determinant<dimension,T>(matrix);
  if constexpr(dimension == 2){
    inverse[0] =  matrix[3]/determinant;
    inverse[1] = -matrix[1]/determinant;
    inverse[2] = -matrix[2]/determinant;
    inverse[3] =  matrix[0]/determinant;
  }else{
    inverse[0] =  (matrix[4]*matrix[8] - matrix[5]*matrix[7])
                /determinant;
    inverse[1] = -(matrix[1]*matrix[8] - matrix[2]*matrix[7])
                /determinant;
    inverse[2] =  (matrix[1]*matrix[5] - matrix[2]*matrix[4])
                /determinant;
    inverse[3] = -(matrix[3]*matrix[8] - matrix[5]*matrix[6])
                /determinant;
    inverse[4] =  (matrix[0]*matrix[8] - matrix[2]*matrix[6])
                /determinant;
    inverse[5] = -(matrix[0]*matrix[5] - matrix[2]*matrix[3])
                /determinant;
    inverse[6] =  (matrix[3]*matrix[7] - matrix[4]*matrix[6])
                /determinant;
    inverse[7] = -(matrix[0]*matrix[7] - matrix[1]*matrix[6])
                /determinant;
    inverse[8] =  (matrix[0]*matrix[4] - matrix[1]*matrix[3])
                /determinant;
  }
}

template<int dimension>
void sizeshape_pointwise_gradient(
    const double *regular_jacobian_transpose,
    const double *metric,
    const double *regular_basis_gradient,
    const double objective_power,
    SANS::SurrealS<dimension,double> *pointwise_gradient)
{
  using S = SANS::SurrealS<dimension,double>;
  S jacobian[dimension*dimension];
  S gram_matrix[dimension*dimension];
  S inverse_gram_matrix[dimension*dimension];
  seed_regular_jacobian<dimension>(
      regular_jacobian_transpose,regular_basis_gradient,jacobian);
  form_metric_gram_matrix<dimension,S>(jacobian,metric,gram_matrix);
  invert_full_matrix<dimension,S>(gram_matrix,inverse_gram_matrix);

  S trace = S(0.0);
  for(int diagonal = 0; diagonal < dimension; diagonal++){
    trace += gram_matrix[diagonal*dimension + diagonal];
  }
  const S determinant = full_determinant<dimension,S>(gram_matrix);

  S trace_power_dimension_minus_one;
  S trace_power_dimension;
  if constexpr(dimension == 2){
    trace_power_dimension_minus_one = trace;
    trace_power_dimension = trace*trace;
  }else{
    trace_power_dimension_minus_one = trace*trace;
    trace_power_dimension = trace*trace*trace;
  }

  S trace_derivative[dimension];
  S inverse_gram_times_basis_gradient[dimension];
  for(int row = 0; row < dimension; row++){
    inverse_gram_times_basis_gradient[row] = S(0.0);
    for(int column = 0; column < dimension; column++){
      inverse_gram_times_basis_gradient[row]
          += inverse_gram_matrix[row*dimension + column]
           * S(regular_basis_gradient[column]);
    }
  }

  S metric_jacobian_inverse_gram_gradient[dimension];
  for(int component = 0; component < dimension; component++){
    trace_derivative[component] = S(0.0);
    metric_jacobian_inverse_gram_gradient[component] = S(0.0);
    for(int physical_component = 0;
        physical_component < dimension; physical_component++){
      S jacobian_times_inverse_gram_gradient = S(0.0);
      for(int regular_component = 0;
          regular_component < dimension; regular_component++){
        trace_derivative[component]
            += S(2.0*regular_basis_gradient[regular_component]
                 * metric[packed_symmetric_index<dimension>(
                       component,physical_component)])
             * jacobian[regular_component*dimension + physical_component];
        jacobian_times_inverse_gram_gradient
            += jacobian[regular_component*dimension + physical_component]
             * inverse_gram_times_basis_gradient[regular_component];
      }
      metric_jacobian_inverse_gram_gradient[component]
          += S(metric[packed_symmetric_index<dimension>(
                   component,physical_component)])
           * jacobian_times_inverse_gram_gradient;
    }
  }

  const double dimension_power = dimension == 2 ? 4.0 : 27.0;
  const S size_shape_error
      = trace_power_dimension
       *(S(1.0) + S(1.0)/(determinant*determinant))
       /S(2.0*dimension_power) - S(1.0);
  const S objective_gradient_scale
      = objective_power == 1.0
      ? S(1.0)
      : S(objective_power)
       *pow(size_shape_error,objective_power - 1.0);

  for(int component = 0; component < dimension; component++){
    const S determinant_derivative
        = S(2.0)*determinant
         *metric_jacobian_inverse_gram_gradient[component];
    const S size_shape_gradient
        = (S(dimension)*trace_power_dimension_minus_one
             *trace_derivative[component]
             *(S(1.0) + S(1.0)/(determinant*determinant))
           - S(2.0)*trace_power_dimension*determinant_derivative
             /(determinant*determinant*determinant))
         /S(2.0*dimension_power);
    pointwise_gradient[component]
        = objective_gradient_scale*size_shape_gradient;
  }
}

template<int dimension, typename T>
void pack_symmetric_matrix(const T *matrix, T *packed_matrix)
{
  packed_matrix[0] = matrix[0];
  packed_matrix[1] = matrix[1];
  packed_matrix[2] = matrix[dimension + 1];
  if constexpr(dimension == 3){
    packed_matrix[3] = matrix[2];
    packed_matrix[4] = matrix[5];
    packed_matrix[5] = matrix[8];
  }
}

template<int dimension>
SANS::SurrealS<dimension,double> step_distance_pointwise_value(
    const double *regular_jacobian_transpose,
    const double *metric,
    const double *regular_basis_gradient,
    const double objective_power,
    const double regularization,
    const bool shape_volume)
{
  using S = SANS::SurrealS<dimension,double>;
  constexpr int packed_count = dimension*(dimension + 1)/2;
  S jacobian[dimension*dimension];
  S gram_matrix[dimension*dimension];
  S packed_gram_matrix[packed_count];
  seed_regular_jacobian<dimension>(
      regular_jacobian_transpose,regular_basis_gradient,jacobian);
  form_metric_gram_matrix<dimension,S>(jacobian,metric,gram_matrix);
  pack_symmetric_matrix<dimension,S>(gram_matrix,packed_gram_matrix);

  S eigenvalues[dimension];
  S eigenvectors[dimension*dimension];
  geteigsym<dimension,S>(
      packed_gram_matrix,eigenvalues,eigenvectors);

  S squared_distance = S(0.0);
  if(shape_volume){
    S determinant = S(1.0);
    S mean_logarithm = S(0.0);
    S logarithms[dimension];
    for(int eigenvalue = 0; eigenvalue < dimension; eigenvalue++){
      determinant *= eigenvalues[eigenvalue];
      logarithms[eigenvalue] = log(eigenvalues[eigenvalue]);
      mean_logarithm += logarithms[eigenvalue]/S(dimension);
    }
    for(int eigenvalue = 0; eigenvalue < dimension; eigenvalue++){
      const S centered_logarithm
          = logarithms[eigenvalue] - mean_logarithm;
      squared_distance += centered_logarithm*centered_logarithm;
    }
    const S volume_coordinate
        = determinant - S(1.0)/determinant;
    squared_distance += volume_coordinate*volume_coordinate
                       /S(4.0*dimension);
  }else{
    for(int eigenvalue = 0; eigenvalue < dimension; eigenvalue++){
      const S logarithm = log(eigenvalues[eigenvalue]);
      squared_distance += logarithm*logarithm;
    }
  }

  return pow(squared_distance + S(regularization*regularization),
             objective_power/2.0)
       - S(std::pow(regularization,objective_power));
}

template<int dimension>
SANS::SurrealS<dimension,double> metric_volume_barrier_value(
    const double *regular_jacobian_transpose,
    const double *metric,
    const double *regular_basis_gradient,
    const double rho0,
    const double beta)
{
  using S = SANS::SurrealS<dimension,double>;
  if(beta <= 0.0 || rho0 <= 0.0) return S(0.0);

  S jacobian[dimension*dimension];
  S gram_matrix[dimension*dimension];
  seed_regular_jacobian<dimension>(
      regular_jacobian_transpose,regular_basis_gradient,jacobian);
  form_metric_gram_matrix<dimension,S>(jacobian,metric,gram_matrix);
  const S rho = sqrt(full_determinant<dimension,S>(gram_matrix));
  if(primal_value(rho) >= rho0) return S(0.0);

  const S logarithmic_collapse = log(S(rho0)/rho);
  const S logarithmic_collapse_squared
      = logarithmic_collapse*logarithmic_collapse;
  return S(beta)*logarithmic_collapse_squared
                *logarithmic_collapse_squared;
}

} // namespace FrozenObjectiveValueAD

} // namespace Metris

#endif
