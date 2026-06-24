//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../quality/msh_metqua.hxx"
#include "../Mesh/Mesh.hxx"
#include "../utils/CT_loop.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/aux_pp_inc.hxx"


namespace Metris{


// Version with default "conformity error" settings that should be called except
// perhaps in specific cases
template <class MFT, QuaFun iquaf>
double getmetquamesh(Mesh<MFT> &msh, int tdim, AsDeg asdmsh, AsDeg asdmet,
                     bool *iinva, double *qmin, double *qmax, double *qavg,
                     dblAr1 *lquae){
  INCVDEPTH(msh.param);

  msh.met.setSpace(MetSpace::Exp);
  *iinva = false;
  int nentt = msh.nentt(tdim);

  double qtot = 0;
  *qmin = 1.0e30;
  *qmax = 0;

  int ncnt = 0;

  if(lquae != NULL){
    lquae->allocate(nentt);
    lquae->set_n(nentt);
  }

  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
  CT_FOR0_INC(2,3,tdim_c){if(tdim == tdim_c){

    const intAr2 &ent2poi = msh.ent2poi(tdim);

    for(int ientt = 0; ientt < nentt; ientt++){
      INCVDEPTH(msh.param);

      if(isdeadent(ientt,ent2poi)) continue;
      ncnt ++;
      double quent = 0;
      //try{
        if(msh.idim == 2){
          METRIS_ASSERT(tdim_c == 2);
          if constexpr(tdim_c == 2){ // avoid bad instantiations
            quent = metqua<MFT,2,tdim_c,iquaf>(msh,asdmsh,asdmet,ientt,1.0);
          }
        }else{
          quent = metqua<MFT,3,tdim_c,iquaf>(msh,asdmsh,asdmet,ientt,1.0);
        }
        // METRIS_ASSERT_MSG(tdim_c < msh.idim || quent <= 1 + 1.0e-15,
        //   "GT 1 quality with tdim = {} gdim = {} quent = {}",tdim_c,msh.idim,quent)
        //if(quent > 1 + 1.0e-15){
        //  printf("## DEBUG QUENT > 1 tdim = {} tdim_c = {} gdim = {} quent {} dif {}\n",
        //    tdim,tdim_c,msh.idim,quent,quent-1);
        //  exit(1);
        //}
      //}catch(...){
      //  *iinva = true;
      //} // Ignore exceptions in this context.
      if(lquae != NULL) (*lquae)[ientt] = quent;
      qtot += quent;
      (*qmin) = MIN(*qmin,quent);
      (*qmax) = MAX(*qmax,quent);
      CPRINTF3(" - getmetquamesh ientt {} dim {} qual = {}\n",ientt,tdim,quent);
    }
  }}CT_FOR1(tdim_c);
  }}CT_FOR1(ideg);

  *qavg = qtot / ncnt;

  //msh.met.setSpace(metspac0);
  return pow(qtot,1.0/msh.param->opt_pnorm);
}

// template double getmetquamesh<MetricFieldAnalytical>(
//          Mesh<MetricFieldAnalytical> &msh, int tdim, AsDeg asdmsh, AsDeg asdmet,
//          bool *iinva, double *qmin, double *qmax, double *qavg,
//          dblAr1 *lquae);
// template double getmetquamesh<MetricFieldFE>(
//          Mesh<MetricFieldFE> &msh, int tdim, AsDeg asdmsh, AsDeg asdmet,
//          bool *iinva, double *qmin, double *qmax, double *qavg,
//          dblAr1 *lquae);


#define EXPAND_TEMPLATE(r, SEQ) \
  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ), BOOST_PP_SEQ_ELEM(1, SEQ))
#define MFT_SEQ    (MetricFieldFE)(MetricFieldAnalytical)
#define QUAFUN_SEQ (QuaFun::Distortion)(QuaFun::Unit)(QuaFun::SizeShape)(QuaFun::StepDistance)

#define INSTANTIATE(MFT_VAL, QUAFUN_VAL) \
template double getmetquamesh< MFT_VAL , QUAFUN_VAL >( \
    Mesh<MFT_VAL>& msh, int tdim, AsDeg asdmsh, AsDeg asdmet, \
    bool* iinva, double* qmin, double* qmax, double* qavg, dblAr1* lquae);

BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE, (MFT_SEQ)(QUAFUN_SEQ))
#undef INSTANTIATE

#undef QUAFUN_SEQ
#undef MFT_SEQ
#undef EXPAND_TEMPLATE

} // End namespace
