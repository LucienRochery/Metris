//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "low_ccoef.hxx"
#include "codegen_ccoef.hxx"
#include "codegen_ccoef_d.hxx"

#include "types.hxx"
#include "utils/aux_misc.hxx"
#include "utils/fact_pow.hxx"
#include "Mesh/MeshBase.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "linalg/det.hxx"
#include "low_eval.hxx"
#include "low_geo.hxx"
#include "codegen_lag2bez.hxx"
//#include "utils/aux_misc.hxx"
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
  if constexpr(ideg == 1){
    bool iflat;
    ccoef[0] = ifact<tdim>()*getmeasentP1<gdim,tdim>(msh, ent2poi[ientt], nrmal, &iflat);
  }else if constexpr(ideg >= 4 || gdim != tdim){
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
  constexpr int ncoef = tdim == 2 ? getnnod2(jdeg) 
                                  : getnnod3(jdeg);
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  double meas = getmeasentP1<gdim,tdim>(msh,ent2poi[ientt],nrmal,iinva);
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



template<int idim, int ideg>
void getccoef_dcoord(const MeshBase &msh, int ientt, int icoor, double *ccoef, dblAr2& d_ccoef){
  METRIS_ENFORCE_MSG(msh.getBasis() == FEBasis::Bezier, 
    "control coefficient derivatives not implemented for Lagrange meshes");

  if(ccoef != NULL) getccoef<idim,idim,ideg>(msh,ientt,NULL,ccoef);

  if constexpr(idim == 2){
    d_ccoef_genbez2<ideg>(msh.fac2poi,msh.coord,ientt,icoor,d_ccoef);
  }else{
    d_ccoef_genbez3<ideg>(msh.tet2poi,msh.coord,ientt,icoor,d_ccoef);
  }
}
#define BOOST_PP_LOCAL_MACRO(n)\
template void getccoef_dcoord<2,n>(const MeshBase &msh, int ientt, int icoor, double *ccoef, dblAr2& d_ccoef);\
template void getccoef_dcoord<3,n>(const MeshBase &msh, int ientt, int icoor, double *ccoef, dblAr2& d_ccoef);
#define BOOST_PP_LOCAL_LIMITS     (2, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


template<int idim, int ideg>
void getccoef_dpoint(const MeshBase &msh, int ientt, int inode, double *ccoef, dblAr2& d_ccoef){
  METRIS_ENFORCE_MSG(msh.getBasis() == FEBasis::Bezier, 
    "control coefficient derivatives not implemented for Lagrange meshes");


  if constexpr(ideg == 1){
    int ifac = ifact<idim>();
    if(ccoef != NULL){
      bool iflat;
      const intAr2& ent2poi = msh.ent2poi(idim);
      ccoef[0] = ifac*getmeasentP1<idim,idim>(msh, ent2poi[ientt], NULL, &iflat);
    } 
    if constexpr(idim == 2){
      // def (x y)^orth = (y -x)
      vdiff_perp(msh.coord[msh.fac2poi(ientt,(inode+1)%3)], 
                 msh.coord[msh.fac2poi(ientt,(inode+2)%3)], d_ccoef[0]);
    }else{
      METRIS_THROW_MSG(TODOExcept(), "Control coeff derivatives per point not implemented in 3D")
    }
  }else{
    if(ccoef != NULL) getccoef<idim,idim,ideg>(msh,ientt,NULL,ccoef);
    if constexpr(idim == 2){
      d_pt_ccoef_genbez2<ideg>(msh.fac2poi,msh.coord,ientt,inode,d_ccoef);
    }else{
      METRIS_THROW_MSG(TODOExcept(), "Control coeff derivatives per point not implemented in 3D")
    }
  }
}
#define BOOST_PP_LOCAL_MACRO(n)\
template void getccoef_dpoint<2,n>(const MeshBase &msh, int ientt, int inode, double *ccoef, dblAr2& d_ccoef);\
template void getccoef_dpoint<3,n>(const MeshBase &msh, int ientt, int inode, double *ccoef, dblAr2& d_ccoef);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




template<int gdim, int tdim, int ideg>
void ccoef_eval(FEBasis ibasis, const intAr2& ent2poi, const dblAr2& coord, 
                int ientt, double *nrmal, double* ccoef){
  static_assert(gdim == 2 || gdim == 3);

  // Get control coeffs by evaluating the Jacobian at the nodes
  double dum[gdim], jmat[tdim*gdim];
  constexpr int jdeg = tdim*(ideg-1);
  constexpr int nppj = getnnode(tdim,jdeg);
  double rwrk[nppj];
  constexpr auto ordent = ORDELT(tdim);
  constexpr auto eval = tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg> ;

  double nrloc[3];

  for(int irnk = 0; irnk < nppj; irnk++){
    double bary[gdim+1];
    for(int ii = 0; ii < tdim + 1; ii++) 
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
