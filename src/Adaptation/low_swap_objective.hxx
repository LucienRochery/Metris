//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1

#ifndef __METRIS_LOW_SWAP_OBJECTIVE_HXX__
#define __METRIS_LOW_SWAP_OBJECTIVE_HXX__

#include "../Mesh/Mesh.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../quality/low_metqua.hxx"

#include <cmath>

namespace Metris{

#ifdef TESTQUALITYALGO

struct CompletedSwapObjectiveReplacement{
  double old_numerator = 0.;
  double new_numerator = 0.;
  int old_element_count = 0;
  int new_element_count = 0;
  bool evaluated = false;
  bool accepted = false;
};

template<class MFT, int gdim, int tdim, QuaFun iquaf>
double completed_swap_element_objective(Mesh<MFT>& msh, int element){
  return metqua<MFT,gdim,tdim,iquaf>(
      msh,AsDeg::Pk,AsDeg::P1,element,1.);
}

// Install a final acceptance decision for a high-order swap. The callback is
// reached only after common cavity reconstruction has created every high-order
// coefficient, applied CAD ownership, and certified the completed elements.
template<class MFT, int gdim, int tdim, QuaFun iquaf>
void configure_completed_swap_acceptance(
    Mesh<MFT>& msh,
    const intAr1& old_region,
    StepDistanceObjectiveState* global_objective,
    CavOprOpt& opts,
    CompletedSwapObjectiveReplacement& replacement){
  replacement = {};
  const intAr2& element_to_point = msh.ent2poi(tdim);
  for(const int element : old_region){
    if(isdeadent(element,element_to_point)) continue;
    const double value = completed_swap_element_objective<
        MFT,gdim,tdim,iquaf>(msh,element);
    METRIS_ENFORCE(std::isfinite(value));
    replacement.old_numerator += value;
    replacement.old_element_count++;
  }
  METRIS_ENFORCE(replacement.old_element_count > 0);

  opts.accept_completed_elements =
      [&msh,global_objective,&replacement]
      (int candidate_tdim, int first_new, int end_new){
        replacement.evaluated = true;
        replacement.new_numerator = 0.;
        replacement.new_element_count = 0;
        if(candidate_tdim != tdim){
          replacement.accepted = false;
          return false;
        }

        const intAr2& candidate_to_point = msh.ent2poi(tdim);
        for(int element = first_new; element < end_new; element++){
          if(isdeadent(element,candidate_to_point)) continue;
          const double value = completed_swap_element_objective<
              MFT,gdim,tdim,iquaf>(msh,element);
          if(!std::isfinite(value)){
            replacement.accepted = false;
            return false;
          }
          replacement.new_numerator += value;
          replacement.new_element_count++;
        }
        if(replacement.new_element_count == 0){
          replacement.accepted = false;
          return false;
        }

        if constexpr(iquaf == QuaFun::StepDistance){
          if(msh.param->step_distance_cavity_target_average){
            replacement.accepted = global_objective != nullptr
                ? global_objective->accepts_replacement(
                      replacement.old_numerator,
                      replacement.old_element_count,
                      replacement.old_element_count,
                      replacement.new_numerator,
                      replacement.new_element_count,
                      replacement.new_element_count)
                : objective_strictly_improves(
                      replacement.new_numerator
                          /replacement.new_element_count,
                      replacement.old_numerator
                          /replacement.old_element_count);
            return replacement.accepted;
          }
        }

        replacement.accepted = objective_strictly_improves(
            replacement.new_numerator,replacement.old_numerator);
        return replacement.accepted;
      };
}

template<QuaFun iquaf>
inline void commit_completed_swap_objective(
    const CompletedSwapObjectiveReplacement& replacement,
    StepDistanceObjectiveState* global_objective){
  METRIS_ENFORCE(replacement.evaluated && replacement.accepted);
  if constexpr(iquaf == QuaFun::StepDistance){
    if(global_objective == nullptr) return;
    global_objective->replace(
        replacement.old_numerator,replacement.old_element_count,
        replacement.old_element_count,
        replacement.new_numerator,replacement.new_element_count,
        replacement.new_element_count);
  }else{
    METRIS_ASSERT(global_objective == nullptr);
  }
}

#endif

} // namespace Metris

#endif
