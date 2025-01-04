//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "low_ccoef.hxx"

#include "types.hxx"
#include "aux_utils.hxx"
#include "codegen_ccoef.hxx"
#include "Mesh/MeshBase.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "linalg/det.hxx"
#include "low_eval.hxx"
#include "low_geo.hxx"
#include "codegen_lag2bez.hxx"
//#include "aux_utils.hxx"
//#include "msh_structs.hxx"
//#include "Mesh/Mesh.hxx"


//#include <inc_hana.hxx>
//#include <boost/hana.hpp>

namespace Metris{
 


template<int gdim, int tdim, int ideg>
void getccoef(const MeshBase &msh, int ientt, double *nrmal, double *ccoef){
  static_assert(gdim == 2 || gdim == 3);
  const intAr2 &ent2poi = tdim == 2 ? msh.fac2poi : msh.tet2poi;

  // Constexpr prevents compilation of non-followed branches
  // ccoef_genbez do not exist for ideg == 1 
  if constexpr(ideg == 1 || ideg >= 4 || gdim != tdim){
    METRIS_ASSERT(nrmal != NULL);
    ccoef_eval<gdim,tdim,ideg>(msh.getBasis(),ent2poi,msh.coord,ientt,nrmal,ccoef);
  }else{
    static_assert(gdim == tdim);
    if(msh.getBasis() == FEBasis::Lagrange){
      ccoef_eval<gdim,gdim,ideg>(msh.getBasis(),ent2poi,msh.coord,ientt,nrmal,ccoef);
    }else{
      if constexpr(gdim == 2){
        ccoef_genbez2<ideg>(ent2poi,msh.coord,ientt,ccoef);
      }else{
        ccoef_genbez3<ideg>(ent2poi,msh.coord,ientt,ccoef);
      }
    }
  }

}

#define BOOST_PP_LOCAL_MACRO(n)\
template void getccoef<2,2,n>(const MeshBase &msh, int ientt, double *nrmal, double *ccoef);\
template void getccoef<3,2,n>(const MeshBase &msh, int ientt, double *nrmal, double *ccoef);\
template void getccoef<3,3,n>(const MeshBase &msh, int ientt, double *nrmal, double *ccoef);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

// Additionally returns whether mesh valid: prefer this over manually looping over
template<int gdim, int tdim, int ideg>
void getsclccoef(const MeshBase &msh, int ientt, double *nrmal, 
                 double *ccoef, bool *iinva){
  getccoef<gdim,tdim,ideg>(msh,ientt,nrmal,ccoef);
  constexpr int jdeg = tdim * (ideg - 1);
  constexpr int ncoef = tdim == 2 ? facnpps[jdeg] 
                                  : tetnpps[jdeg];
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  double meas = 
    getmeasentP1<gdim,tdim>(ent2poi[ientt],msh.coord,msh.param->vtol,NULL,iinva);
  constexpr int fact = ifact<tdim>();

  for(int icoef = 0; icoef < ncoef; icoef++){
    ccoef[icoef] /= (abs(meas) * fact);
    if(ccoef[icoef] >= msh.param->jtol) continue;
    *iinva = true;
  }
}
#define BOOST_PP_LOCAL_MACRO(n)\
template void getsclccoef<2,2,n>(const MeshBase &msh, int ientt, \
                              double *nrmal, double *ccoef, bool *iinva);\
template void getsclccoef<3,2,n>(const MeshBase &msh, int ientt, \
                              double *nrmal, double *ccoef, bool *iinva);\
template void getsclccoef<3,3,n>(const MeshBase &msh, int ientt, \
                              double *nrmal, double *ccoef, bool *iinva);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


template<int gdim, int tdim, int ideg>
void ccoef_eval(FEBasis ibasis, const intAr2& ent2poi, const dblAr2& coord, int ientt, double *nrmal, double* ccoef){
  static_assert(gdim == 2 || gdim == 3);

  // Get control coeffs by evaluating the Jacobian at the nodes
  double dum[gdim], jmat[tdim*gdim];
  constexpr int jdeg = tdim*(ideg-1);
  constexpr int nppj = tdim == 2 ? facnpps[jdeg] : tetnpps[jdeg];
  double rwrk[nppj];
  constexpr auto ordent = ORDELT(gdim);
  constexpr auto eval = tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg> ;

  double nrloc[3];
  //if constexpr(tdim == 2 && gdim == 3){
    // NO !
    //getnorfacP1(ent2poi[ientt],coord,nrmal);
    //double nrm = sqrt(getnrml2<3>(nrmal));
    //if(nrm < 1.0e-16) METRIS_THROW_MSG(GeomExcept(), "0 face P1 normal "<<nrm);
    //for(int ii = 0; ii < 3 ;ii++) nrmal[ii] /= nrm;
  //}

  for(int irnk = 0; irnk < nppj; irnk++){
    double bary[gdim+1];
    for(int ii = 0; ii < gdim + 1; ii++) 
      bary[ii] = ordent[jdeg][irnk][ii] / (double) jdeg;
    eval(coord, ent2poi[ientt],ibasis, DifVar::Bary, DifVar::None, bary, dum, jmat, NULL);

    if constexpr(tdim == 2 && gdim == 3){
      vecprod(&jmat[0],&jmat[3],nrloc);
      double det = sqrt(getnrml2<3>(nrloc));
      double sg  = getprdl2<3>(nrmal,nrloc);
      if(sg < 0) det = -det;
      rwrk[irnk] = det;

    }else{
      rwrk[irnk] = detmat<gdim>(jmat);
    }

  }

  constexpr auto lag2bez = tdim == 2 ? lag2bez2<jdeg,1> : lag2bez3<jdeg,1>;
  // Convert to Bézier field
  dblAr2 tmp(nppj,1,ccoef);
  lag2bez(&ccoef_eval_lfld::lfld[0],dblAr2(nppj,1,rwrk),tmp);
}


#define BOOST_PP_LOCAL_MACRO(n)\
template void ccoef_eval<2,2,n>(FEBasis ibasis, const intAr2& ent2poi, const dblAr2& coord, int ientt, double *nrmal, double* ccoef);\
template void ccoef_eval<3,2,n>(FEBasis ibasis, const intAr2& ent2poi, const dblAr2& coord, int ientt, double *nrmal, double* ccoef);\
template void ccoef_eval<3,3,n>(FEBasis ibasis, const intAr2& ent2poi, const dblAr2& coord, int ientt, double *nrmal, double* ccoef);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



} // End namespace
