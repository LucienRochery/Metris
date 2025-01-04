//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../quality/low_metqua.hxx"
#include "../quality/msh_metqua.hxx"
#include "../Mesh/Mesh.hxx"
#include "../CT_loop.hxx"


namespace Metris{


// Version with default "conformity error" settings that should be called except perhaps in specific cases
template <class MFT, int tdim, AsDeg asdmet>
double getmetquamesh(Mesh<MFT> &msh, bool *iinva, double *qmin, double *qmax, double *qavg, dblAr1 *lquae){
  return getmetquamesh<MFT,tdim,asdmet>(msh,-1,2,1.0,iinva,qmin,qmax,qavg,lquae);
}
template double getmetquamesh<MetricFieldAnalytical, 2, AsDeg::P1>(Mesh<MetricFieldAnalytical> &msh,
                                                               bool* iinva, double* qmin, double *qmax, double *qavg, dblAr1 *lquae);
template double getmetquamesh<MetricFieldFE        , 2, AsDeg::P1>(Mesh<MetricFieldFE        > &msh,
                                                               bool* iinva, double* qmin, double *qmax, double *qavg, dblAr1 *lquae);
template double getmetquamesh<MetricFieldAnalytical, 2, AsDeg::Pk>(Mesh<MetricFieldAnalytical> &msh,
                                                               bool* iinva, double* qmin, double *qmax, double *qavg, dblAr1 *lquae);
template double getmetquamesh<MetricFieldFE        , 2, AsDeg::Pk>(Mesh<MetricFieldFE        > &msh,
                                                               bool* iinva, double* qmin, double *qmax, double *qavg, dblAr1 *lquae);

template double getmetquamesh<MetricFieldAnalytical, 3, AsDeg::P1>(Mesh<MetricFieldAnalytical> &msh,
                                                               bool* iinva, double* qmin, double *qmax, double *qavg, dblAr1 *lquae);
template double getmetquamesh<MetricFieldFE        , 3, AsDeg::P1>(Mesh<MetricFieldFE        > &msh,
                                                               bool* iinva, double* qmin, double *qmax, double *qavg, dblAr1 *lquae);
template double getmetquamesh<MetricFieldAnalytical, 3, AsDeg::Pk>(Mesh<MetricFieldAnalytical> &msh,
                                                               bool* iinva, double* qmin, double *qmax, double *qavg, dblAr1 *lquae);
template double getmetquamesh<MetricFieldFE        , 3, AsDeg::Pk>(Mesh<MetricFieldFE        > &msh,
                                                               bool* iinva, double* qmin, double *qmax, double *qavg, dblAr1 *lquae);


template <class MFT, int tdim, AsDeg asdmet>
double getmetquamesh(Mesh<MFT> &msh, int power, int pnorm, double difto,
  bool *iinva, double *qmin, double *qmax, double *qavg, dblAr1 *lquae){

  static_assert(tdim == 2 || tdim == 3);

  MetSpace metspac0 = msh.met.getSpace();
  msh.met.setSpace(MetSpace::Log);

  int nentt = tdim == 2 ? msh.nface : msh.nelem;
  intAr2 &ent2poi = tdim == 2 ? msh.fac2poi : msh.tet2poi;

  double qtot = 0;
  *qmin = 1.0e30;
  *qmax = 0;

  int ncnt = 0;

  if(lquae != NULL){
    lquae->allocate(nentt);
    lquae->set_n(nentt);
  }

  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
    for(int ientt = 0; ientt < nentt; ientt++){
      if(isdeadent(ientt,ent2poi)) continue;
      ncnt ++;
      double quent = 0;
      try{
        if(msh.idim == 2){
          if constexpr(tdim == 2){
            quent = metqua<MFT,2,tdim>(msh,AsDeg::Pk,asdmet,
                                       ientt,power,pnorm,difto);
          }else{
            METRIS_THROW(WArgExcept());
          }
        }else{
          quent = metqua<MFT,3,tdim>(msh,AsDeg::Pk,asdmet,
                                     ientt,power,pnorm,difto);
        }
      }catch(...){
        *iinva = true; 
      } // Ignore exceptions in this context. 
      if(lquae != NULL) (*lquae)[ientt] = quent;
      qtot += quent;
      (*qmin) = MIN(*qmin,quent);
      (*qmax) = MAX(*qmax,quent);
    }
  }}CT_FOR1(ideg);

  *qavg = qtot / ncnt;

  msh.met.setSpace(metspac0);
  return pow(qtot,1.0/pnorm);
}

template double getmetquamesh<MetricFieldAnalytical, 2, AsDeg::P1>(Mesh<MetricFieldAnalytical>& msh, int power, int pnorm, double difto,
                                                               bool* iinva, double* qmin, double *qmax, double *qtot, dblAr1 *lquae);
template double getmetquamesh<MetricFieldFE        , 2, AsDeg::P1>(Mesh<MetricFieldFE        >& msh, int power, int pnorm, double difto,
                                                               bool* iinva, double* qmin, double *qmax, double *qtot, dblAr1 *lquae);
template double getmetquamesh<MetricFieldAnalytical, 2, AsDeg::Pk>(Mesh<MetricFieldAnalytical>& msh, int power, int pnorm, double difto,
                                                               bool* iinva, double* qmin, double *qmax, double *qtot, dblAr1 *lquae);
template double getmetquamesh<MetricFieldFE        , 2, AsDeg::Pk>(Mesh<MetricFieldFE        >& msh, int power, int pnorm, double difto,
                                                               bool* iinva, double* qmin, double *qmax, double *qtot, dblAr1 *lquae);

template double getmetquamesh<MetricFieldAnalytical, 3, AsDeg::P1>(Mesh<MetricFieldAnalytical>& msh, int power, int pnorm, double difto,
                                                               bool* iinva, double* qmin, double *qmax, double *qtot, dblAr1 *lquae);
template double getmetquamesh<MetricFieldFE        , 3, AsDeg::P1>(Mesh<MetricFieldFE        >& msh, int power, int pnorm, double difto,
                                                               bool* iinva, double* qmin, double *qmax, double *qtot, dblAr1 *lquae);
template double getmetquamesh<MetricFieldAnalytical, 3, AsDeg::Pk>(Mesh<MetricFieldAnalytical>& msh, int power, int pnorm, double difto,
                                                               bool* iinva, double* qmin, double *qmax, double *qtot, dblAr1 *lquae);
template double getmetquamesh<MetricFieldFE        , 3, AsDeg::Pk>(Mesh<MetricFieldFE        >& msh, int power, int pnorm, double difto,
                                                               bool* iinva, double* qmin, double *qmax, double *qtot, dblAr1 *lquae);
} // End namespace
