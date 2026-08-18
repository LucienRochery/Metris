//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_METQUA__
#define __METRIS_LOW_METQUA__
/*
	Functions relating to Q_M(K) = 1/n tr(J_K^T J_0^{-T} M J_0^{-1} J_K) / det(J_K^T M J_K)^(1/n)
*/


#include "quafun.hxx"

#include "../aux_exceptions.hxx"
#include "../Mesh/MeshFwd.hxx"

#include <functional>


namespace Metris{

enum class AsDeg;

/* ----------------------------------------------
 Main functions
   ---------------------------------------------- */
// Integrated quality
template <class MetricFieldType, int gdim, int tdim,
          QuaFun iquaf = QuaFun::Distortion,
          typename ftype = double>
ftype metqua(Mesh<MetricFieldType> &msh,
             AsDeg asdmsh, AsDeg asdmet,
             int ielem, ftype difto = 1);

// Unit aggregation weight for the historical CavityTargetAverage plumbing.
// The restored formulation is averaged by element count, not target volume.
template <class MetricFieldType, int gdim, int tdim>
double step_distance_element_target_weight(Mesh<MetricFieldType> &msh,
                                           AsDeg asdmet,
                                           int ielem);

// Build the exact mesh-wide totals used to assess a local replacement.
template <class MetricFieldType, int gdim, int tdim>
StepDistanceObjectiveState step_distance_global_objective_state(
    Mesh<MetricFieldType> &msh,
    AsDeg asdmsh,
    AsDeg asdmet);



/* ----------------------------------------------
 Convenience: return f pointers selecting depending on iquaf
   ---------------------------------------------- */
template<class MFT, int gdim, int tdim, class ftype=double>
std::function<ftype(Mesh<MFT>&,AsDeg,AsDeg,int,ftype)>
get_quafun(QuaFun iquaf){
  if(iquaf == QuaFun::Distortion){
    return metqua<MFT,gdim,tdim,QuaFun::Distortion,ftype>;
  }else if(iquaf == QuaFun::Unit){
    return metqua<MFT,gdim,tdim,QuaFun::Unit,ftype>;
  }else if(iquaf == QuaFun::SizeShape){
    return metqua<MFT,gdim,tdim,QuaFun::SizeShape,ftype>;
  }else if(iquaf == QuaFun::StepDistance){
    return metqua<MFT,gdim,tdim,QuaFun::StepDistance,ftype>;
  }
  else{
    METRIS_THROW_MSG("TODO: cf quafun_")
  }
}



} // End namespace

#endif
