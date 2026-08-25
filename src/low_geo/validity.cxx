// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#include "validity.hxx"

#include "ccoef.hxx"
#include "measure.hxx"

#include "../Mesh/MeshBase.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../ho_constants.hxx"
#include "../linalg/det.hxx"
#include "../low_eval.hxx"
#include "../utils/fact_pow.hxx"

#include <cmath>

namespace Metris {

template<int gdim, int ideg>
ElementValidityResult classify_element_validity(const MeshBase &msh,
                                                 int element)
{
  static_assert(gdim == 2 || gdim == 3);
  static_assert(ideg >= 1);

  ElementValidityResult result;
  if(msh.param == nullptr) return result;

  constexpr int jacobian_degree = gdim*(ideg - 1);
  constexpr int coefficient_count = getnnode(gdim,jacobian_degree);
  const intAr2 &element_to_point = msh.ent2poi(gdim);

  bool p1_orientation_failed = false;
  const double p1_measure = getmeasentP1<gdim,gdim>(
      msh,element,nullptr,&p1_orientation_failed);
  const double normalization
      = std::abs(p1_measure)*(double)ifact<gdim>();
  if(!(normalization > 0.0) || !std::isfinite(normalization)){
    return result;
  }

  double coefficients[coefficient_count];
  getccoef<gdim,gdim,ideg>(msh,element,nullptr,coefficients);

  bool finite_coefficients = true;
  double lower_bound = 0.0;
  int lower_bound_index = -1;
  for(int coefficient = 0; coefficient < coefficient_count; coefficient++){
    const double normalized_coefficient
        = coefficients[coefficient]/normalization;
    if(!std::isfinite(normalized_coefficient)){
      finite_coefficients = false;
      continue;
    }
    if(lower_bound_index < 0 || normalized_coefficient < lower_bound){
      lower_bound = normalized_coefficient;
      lower_bound_index = coefficient;
    }
  }
  if(finite_coefficients && lower_bound_index >= 0){
    result.normalized_lower_bound = lower_bound;
    result.lower_bound_coefficient_index = lower_bound_index;
  }

  const double threshold = msh.param->jtol;
  if(!p1_orientation_failed && finite_coefficients
     && lower_bound_index >= 0 && lower_bound >= threshold){
    result.status = ElementValidityStatus::Certified;
    return result;
  }

  constexpr auto evaluate_geometry
      = gdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;
  constexpr auto ordering = ORDELT(gdim);
  double barycentric_coordinates[gdim + 1];
  double physical_coordinates[gdim];
  double jacobian[gdim*gdim];

  for(int sample = 0; sample < coefficient_count; sample++){
    if constexpr(jacobian_degree == 0){
      for(int barycentric = 0; barycentric < gdim + 1; barycentric++){
        barycentric_coordinates[barycentric] = 1.0/(gdim + 1.0);
      }
    }else{
      for(int barycentric = 0; barycentric < gdim + 1; barycentric++){
        barycentric_coordinates[barycentric]
            = ordering[jacobian_degree][sample][barycentric]
             /(double)jacobian_degree;
      }
    }

    evaluate_geometry(
        msh.coord,element_to_point[element],msh.getBasis(),
        DifVar::Bary,DifVar::None,barycentric_coordinates,
        physical_coordinates,jacobian,nullptr);
    const double normalized_jacobian
        = detmat<gdim>(jacobian)/normalization;
    if(std::isfinite(normalized_jacobian)
       && normalized_jacobian < threshold){
      result.status = ElementValidityStatus::Invalid;
      result.normalized_witness = normalized_jacobian;
      result.witness_sample_index = sample;
      return result;
    }
  }

  return result;
}

#define BOOST_PP_LOCAL_MACRO(n)                                                \
template ElementValidityResult classify_element_validity<2,n>(                \
    const MeshBase &msh, int element);                                         \
template ElementValidityResult classify_element_validity<3,n>(                \
    const MeshBase &msh, int element);
#define BOOST_PP_LOCAL_LIMITS (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

} // namespace Metris
