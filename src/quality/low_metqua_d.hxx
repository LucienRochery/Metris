//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_METQUA_D__
#define __METRIS_LOW_METQUA_D__


#include "quafun.hxx"
#include "../Mesh/MeshFwd.hxx"

#include <functional>
#include <type_traits>


namespace Metris{

enum class FEBasis;
enum class DifVar;
enum class AsDeg;

// Integrated quality
template <class MetricFieldType, int gdim, int tdim,
          QuaFun iquaf = QuaFun::Distortion, typename ftype = double>
ftype d_metqua(Mesh<MetricFieldType> &msh,
               AsDeg asdmsh, AsDeg asdmet,
               int ielem, int ivar,
               FEBasis dofbas, DifVar idifmet,
               ftype*__restrict__ dquael, ftype*__restrict__ hquael,
               double difto = 1);

// All derivatives
template <class MetricFieldType, int gdim, int mshdeg, AsDeg asdmet, typename ftype = double>
ftype D_metqua(Mesh<MetricFieldType> &msh, int ielem,
               FEBasis dofbas, DifVar idifmet,
               ftype*__restrict__ dquael, ftype*__restrict__ hquael,
               double difto = 1);

template <class MetricFieldType, int gdim, int mshdeg, AsDeg asdmet, typename ftype = double>
ftype D_metqua(Mesh<MetricFieldType> &msh, const int *ent2poi,
               FEBasis dofbas, DifVar idifmet,
               ftype*__restrict__ dquael, ftype*__restrict__ hquael,
               double difto = 1);

// Differentiated wrt all nodes and coordinates
template <class MFT, int gdim, int mshdeg, AsDeg asdmet, typename ftype>
ftype D_quafun_distortion(Mesh<MFT> &msh,
                  const int* ent2poi,
                  const double*__restrict__ bary,
                  FEBasis dofbas,
                  DifVar idifmet,
                  ftype*__restrict__ dquael,
                  ftype*__restrict__ hquael);







template<class MFT, int gdim, int tdim, class ftype=double>
std::function<ftype(Mesh<MFT>&,AsDeg,AsDeg,int,int,FEBasis,DifVar,
                    ftype*__restrict__,ftype*__restrict__,double)>
get_d_quafun(QuaFun iquaf){
  if(iquaf == QuaFun::Distortion){
    return d_metqua<MFT,gdim,tdim,QuaFun::Distortion,ftype>;
  }else if(iquaf == QuaFun::Unit){
    return d_metqua<MFT,gdim,tdim,QuaFun::Unit,ftype>;
  }else if(iquaf == QuaFun::SizeShape){
    return d_metqua<MFT,gdim,tdim,QuaFun::SizeShape,ftype>;
  }else if(iquaf == QuaFun::StepDistance){
    return d_metqua<MFT,gdim,tdim,QuaFun::StepDistance,ftype>;
  }else if(iquaf == QuaFun::Connection){
    if constexpr(gdim == tdim && std::is_same_v<ftype,double>)
      return d_metqua_connection<MFT,gdim>;
    else METRIS_THROW_MSG("Connection smoothing requires full-dimensional double precision");
  }
  else{
    METRIS_THROW_MSG("TODO: cf quafun_")
  }
}

} // End namespace

#endif
