//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_QUAFUN__
#define __METRIS_QUAFUN__


#include "quafun_distortion.hxx"
#include "quafun_unit.hxx"

#include "../metris_constants.hxx"
#include "../aux_exceptions.hxx"

//#include <functional>

namespace Metris{

  enum class QuaFun{Distortion,Unit};


  template<class MFT, int gdim, int tdim, 
           QuaFun iquaf, class ftype=double>
  constexpr auto get_quafun_xi(){
    if constexpr(iquaf == QuaFun::Distortion){
      return quafun_distortion<MFT,gdim,tdim,ftype>;
    }else if(iquaf == QuaFun::Unit){
      return quafun_unit<MFT,gdim,tdim,ftype>;
    }else{
      METRIS_THROW_MSG(TODOExcept(), "Implement QuaFun in get_quafun_xi");
    }
    METRIS_THROW_MSG(TODOExcept(), "Implement QuaFun in get_quafun_xi");
    return quafun_distortion<MFT,gdim,tdim,ftype>;
  }

  template<class MFT, int gdim, int tdim, 
           QuaFun iquaf, class ftype=double>
  constexpr auto get_d_quafun_xi(){
    if constexpr(iquaf == QuaFun::Distortion){
      return d_quafun_distortion<MFT,gdim,ftype>;
    }else if(iquaf == QuaFun::Unit){
      return d_quafun_unit<MFT,gdim,ftype>;
    }else{
      METRIS_THROW_MSG(TODOExcept(), "Implement QuaFun in get_quafun_xi");
    }
    METRIS_THROW_MSG(TODOExcept(), "Implement QuaFun in get_quafun_xi");
    return d_quafun_distortion<MFT,gdim,ftype>;
  }


  //template<class MFT> class Mesh;

  //template<class MFT, int gdim, int tdim, int ideg, 
  //          AsDeg asdmsh, AsDeg asdmet, class ftype=double>
  //struct quafun_t{

  //}

  #if 0
  template<class MFT, int gdim, int tdim, int ideg, 
            AsDeg asdmsh, AsDeg asdmet, class ftype=double>
  struct QuaFunList{

    using   quafun_t = std::function<ftype(Mesh<MFT>&,int,int,int,ftype)>;
    using d_quafun_t =
      std::function<ftype(Mesh<MFT>&,int,int,int,FEBasis,DifVar,
                    ftype*__restrict__,ftype*__restrict__,int,double)>;

    using   quafun_xi_t = std::function<ftype(Mesh<MFT>&,
                               const int*__restrict__,const double*__restrict__,
                               int, double*__restrict__)>;
#if 0 

    quafun_t quafun(QuaFun iquaf){
      if(iquaf == QuaFun::Distortion){
        return metqua<MFT,gdim,tdim,ideg,asdmsh,QuaFun::Distortion,ftype>;
      }else{
        METRIS_THROW_MSG(TODOExcept(), "Implement "<<(int)iquaf<<" in QuaFunList");
      }
    }

    d_quafun_t d_quafun(QuaFun iquaf){
      if(iquaf == QuaFun::Distortion){
        return d_metqua<MFT,gdim,ftype>;
      }else{
        METRIS_THROW_MSG(TODOExcept(), "Implement "<<(int)iquaf<<" in QuaFunList");
      }
    }

    template<QuaFun iquaf>
    constexpr quafun_t quafun(){
      return metqua<MFT,gdim,tdim,ideg,asdmsh,iquaf,ftype>;
    }

    template<QuaFun iquaf>
    constexpr quafun_xi_t quafun_xi(){
      return quafun_distortion<MFT,gdim,tdim,iquaf,ftype>;
    }

#endif

    quafun_t quafun(QuaFun iquaf);

    d_quafun_t d_quafun(QuaFun iquaf);

    template<QuaFun iquaf>
    constexpr quafun_t quafun();

    template<QuaFun iquaf>
    constexpr quafun_xi_t quafun_xi();

  };

  #endif
}

#endif