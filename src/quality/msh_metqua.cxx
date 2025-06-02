//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../quality/low_metqua.hxx"
#include "../quality/msh_metqua.hxx"
#include "../Mesh/Mesh.hxx"
#include "../utils/CT_loop.hxx"
#include "../utils/mprintf.hxx"


namespace Metris{


// Version with default "conformity error" settings that should be called except 
// perhaps in specific cases
template <class MFT>
double getmetquamesh(Mesh<MFT> &msh, int tdim, AsDeg asdmsh, AsDeg asdmet,
                     bool *iinva, double *qmin, double *qmax, double *qavg, 
                     dblAr1 *lquae){
  
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
          if constexpr(tdim_c == 2){
            quent = metqua<MFT,2,tdim_c>(msh,asdmsh,asdmet,ientt,1.0);
          }else{
            METRIS_THROW(WArgExcept());
          }
        }else{
          quent = metqua<MFT,3,tdim_c>(msh,asdmsh,asdmet,ientt,1.0);
        }
        if(quent > 1){
          printf("## DEBUG QUENT > 1 tdim = %d tdim_c = %d gdim = %d\n",
            tdim,tdim_c,msh.idim);
          exit(1);
        }
      //}catch(...){
      //  *iinva = true; 
      //} // Ignore exceptions in this context. 
      if(lquae != NULL) (*lquae)[ientt] = quent;
      qtot += quent;
      (*qmin) = MIN(*qmin,quent);
      (*qmax) = MAX(*qmax,quent);
      CPRINTF3(" - getmetquamesh ientt %d dim %d qual = %e\n",ientt,tdim,quent);
    }
  }}CT_FOR1(tdim_c);
  }}CT_FOR1(ideg);

  *qavg = qtot / ncnt;

  //msh.met.setSpace(metspac0);
  return pow(qtot,1.0/msh.param->opt_pnorm);
}
template double getmetquamesh<MetricFieldAnalytical>(
         Mesh<MetricFieldAnalytical> &msh, int tdim, AsDeg asdmsh, AsDeg asdmet,
         bool *iinva, double *qmin, double *qmax, double *qavg, 
         dblAr1 *lquae);
template double getmetquamesh<MetricFieldFE>(
         Mesh<MetricFieldFE> &msh, int tdim, AsDeg asdmsh, AsDeg asdmet,
         bool *iinva, double *qmin, double *qmax, double *qavg, 
         dblAr1 *lquae);

} // End namespace
