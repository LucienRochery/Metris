//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#if 0
#include "low_metqua.hxx"
#include "low_metqua_d.hxx"

#include "quafun.hxx"

#include "../Mesh/Mesh.hxx"
#include "../metris_constants.hxx"
#include <functional>

namespace Metris{



template<class MFT, int gdim, int tdim, int ideg, 
          AsDeg asdmsh, AsDeg asdmet, class ftype>
typename QuaFunList<MFT,gdim,tdim,ideg,asdmsh,asdmet,ftype>::quafun_t 
QuaFunList<MFT,gdim,tdim,ideg,asdmsh,asdmet,ftype>
  ::quafun(QuaFun iquaf){
  if(iquaf == QuaFun::Distortion){
    return metqua<MFT,gdim,tdim,ideg,asdmsh,QuaFun::Distortion,ftype>;
  }else{
    METRIS_THROW_MSG(TODOExcept(), "Implement "<<(int)iquaf<<" in QuaFunList");
  }
}

template<class MFT, int gdim, int tdim, int ideg, 
          AsDeg asdmsh, AsDeg asdmet, class ftype>
typename QuaFunList<MFT,gdim,tdim,ideg,asdmsh,asdmet,ftype>::d_quafun_t 
QuaFunList<MFT,gdim,tdim,ideg,asdmsh,asdmet,ftype>::d_quafun(QuaFun iquaf){
  if(iquaf == QuaFun::Distortion){
    return d_metqua<MFT,gdim,ideg,asdmet,ftype>;
  }else{
    METRIS_THROW_MSG(TODOExcept(), "Implement "<<(int)iquaf<<" in QuaFunList");
  }
}

template<class MFT, int gdim, int tdim, int ideg, 
          AsDeg asdmsh, AsDeg asdmet, class ftype>
template<QuaFun iquaf>
constexpr 
typename QuaFunList<MFT,gdim,tdim,ideg,asdmsh,asdmet,ftype>::quafun_t 
QuaFunList<MFT,gdim,tdim,ideg,asdmsh,asdmet,ftype>::quafun(){
  return metqua<MFT,gdim,tdim,ideg,asdmsh,iquaf,ftype>;
}



template<class MFT, int gdim, int tdim, int ideg, 
          AsDeg asdmsh, AsDeg asdmet, class ftype>
template<QuaFun iquaf>
constexpr 
typename QuaFunList<MFT,gdim,tdim,ideg,asdmsh,asdmet,ftype>::quafun_xi_t 
QuaFunList<MFT,gdim,tdim,ideg,asdmsh,asdmet,ftype>::quafun_xi(){
  return quafun_distortion<MFT,gdim,tdim,ideg,asdmet,iquaf,ftype>;
}


}

#endif