//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_metqua.hxx"

#include "../ho_constants.hxx"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "../aux_topo.hxx"
#include "../Mesh/Mesh.hxx"

#include "../aux_pp_inc.hxx"

namespace Metris{



template <class MFT, int gdim, int tdim, QuaFun iquaf, typename ftype>
ftype metqua(Mesh<MFT> &msh, AsDeg asdmsh, AsDeg asdmet,
             int ielem, int power, int pnorm, ftype difto){
	static_assert(tdim==2 || tdim==3);
  METRIS_ASSERT(!isdeadent(ielem,tdim == 2 ? msh.fac2poi : msh.tet2poi));
  int* ent2pol = tdim == 2 ? msh.fac2poi[ielem] : msh.tet2poi[ielem];
  return metqua0<MFT,gdim,tdim,iquaf,ftype>(msh,asdmsh,asdmet,
                                            ent2pol,power,pnorm,difto);
}

template <class MFT, int gdim, int tdim, 
          QuaFun iquaf, typename ftype>
ftype metqua0(Mesh<MFT> &msh, AsDeg asdmsh, AsDeg asdmet, 
              const int* ent2pol, int power, int pnorm, ftype difto){
  static_assert(gdim==2 || gdim==3);

  double bary[tdim+1];

  ftype qutet; 

  //// Performance impact should be zero
  //constexpr auto quafun_xi
  //  = QuaFunList<MFT,gdim,tdim,ideg,AsDeg::Pk,AsDeg::Pk>{}.
  //    template quafun_xi<iquaf>();

  constexpr auto quafun_xi 
    = get_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();

  const int ideg = msh.curdeg;
  if(asdmet == AsDeg::Pk && ideg > 1){
    qutet = 0.0;

    const int idegj = SMOO_DEGJ(ideg);
    const int nnodj = tdim == 2 ? facnpps[idegj] : tetnpps[idegj];

    constexpr auto ordelt = ORDELT(tdim);

    for(int iquad = 0; iquad < tdim+1; iquad++){
      for(int ii = 0; ii < tdim + 1; ii++){
        bary[ii] = ordelt[idegj][iquad][ii]/((double) (idegj));
      }
      //ftype qua0 = quafun_distortion<MFT,gdim,tdim,ideg,asdmet,ftype>(msh,ent2pol,bary,power);
      ftype qua0 = quafun_xi(msh,asdmet,asdmsh,ent2pol,bary,power,NULL);
      qutet += pow(abs(qua0 - difto),pnorm) / nnodj;
    }
    for(int iquad = tdim+1; iquad < nnodj; iquad++){
      for(int ii = 0; ii < tdim + 1; ii++){
        bary[ii] = ordelt[idegj][iquad][ii]/((double) (idegj));
      }
      //ftype qua0 = quafun_distortion<MFT,gdim,tdim,ideg,asdmet,ftype>(msh,ent2pol,bary,power);
      ftype qua0 = quafun_xi(msh,asdmet,asdmsh,ent2pol,bary,power,NULL);
      qutet += pow(abs(qua0 - difto),pnorm) / nnodj;
    }
  }else{
    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim  + 1);
    //qutet = quafun_distortion<MFT,gdim,tdim,1,asdmet,ftype>(msh,ent2pol,bary,power);
    #ifndef NDEBUG
      try{
    #endif
    qutet = quafun_xi(msh,asdmet,asdmsh,ent2pol,bary,power,NULL);
    #ifndef NDEBUG
      }catch(const MetrisExcept &e){
        printf("## metqua ent2pol \n");
        intAr1(facnpps[ideg], ent2pol).print();
        throw(e);
      }
    #endif
    qutet = pow(abs(qutet - difto),pnorm);
  }
  return qutet;
}

#define EXPAND_TEMPLATE(r,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ),\
                              BOOST_PP_SEQ_ELEM(2, SEQ))
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical) 
#define QUAFUN_SEQ (QuaFun::Distortion)(QuaFun::Unit)


#define INSTANTIATE(MFT_VAL,QUAFUN,FTYPE)\
template FTYPE metqua< MFT_VAL , 2, 2, QUAFUN, FTYPE>\
            (Mesh<MFT_VAL> &msh, AsDeg, AsDeg, \
             int ielem, int power, int pnorm, FTYPE difto);\
template FTYPE metqua< MFT_VAL , 3, 2, QUAFUN, FTYPE>\
            (Mesh<MFT_VAL> &msh, AsDeg, AsDeg, \
             int ielem, int power, int pnorm, FTYPE difto);\
template FTYPE metqua< MFT_VAL , 3, 3, QUAFUN, FTYPE>\
            (Mesh<MFT_VAL> &msh, AsDeg, AsDeg, \
             int ielem, int power, int pnorm, FTYPE difto);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,(MFT_SEQ)(QUAFUN_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE

#define INSTANTIATE(MFT_VAL,QUAFUN,FTYPE)\
template FTYPE metqua0< MFT_VAL , 2, 2, QUAFUN,FTYPE>\
   (Mesh<MFT_VAL> &msh, AsDeg, AsDeg,\
    const int* ent2pol, int power, int pnorm, FTYPE difto);\
template FTYPE metqua0< MFT_VAL , 3, 2, QUAFUN,FTYPE>\
   (Mesh<MFT_VAL> &msh, AsDeg, AsDeg,\
    const int* ent2pol, int power, int pnorm, FTYPE difto);\
template FTYPE metqua0< MFT_VAL , 3, 3, QUAFUN,FTYPE>\
   (Mesh<MFT_VAL> &msh, AsDeg, AsDeg,\
    const int* ent2pol, int power, int pnorm, FTYPE difto);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,(MFT_SEQ)(QUAFUN_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE



#undef EXPAND_TEMPLATE
#undef MFT_SEQ // note these two could go into headers 

} // End namespace
