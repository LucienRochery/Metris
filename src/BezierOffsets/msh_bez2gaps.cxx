//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_bez2gaps.hxx"

#include "ho_constants.hxx"
#include "Mesh/Mesh.hxx"
#include "aux_topo.hxx"
#include "aux_exceptions.hxx"
#include "Mesh/MeshBase.hxx"

#include <boost/preprocessor/iteration/local.hpp>


namespace Metris{

template <int ndimn, int ideg>
void eltgapsbezconv(MeshBase &msh, int ielem, int isig, bool confirm){
  static_assert(ndimn == 2 || ndimn == 3);
  if constexpr(ideg == 1) return;
  if(!confirm) METRIS_THROW_MSG(WArgExcept(),"Probably don't want to call this");

  METRIS_ASSERT(ielem >= 0);

  constexpr int nedgl = (ndimn*(ndimn+1))/2;
  constexpr int nvert = ndimn + 1;
  constexpr int nfacl = ndimn == 2 ? 1 : 4;

  constexpr auto lnoed = []() -> auto {
    if constexpr(ndimn == 2) return lnoed2;
    else                     return lnoed3;
  }();
  constexpr auto ordelt = []() -> auto {
    if constexpr(ndimn == 2) return ordfac.s;
    else                     return ordtet.s;
  }();

  const intAr2& ent2poi = ndimn == 2 ? msh.fac2poi : msh.tet2poi;

  METRIS_ASSERT(!isdeadent(ielem,ent2poi));

  if constexpr(ideg > 1){
    for(int iedg = 0; iedg < nedgl; iedg++){
      int irnk0 = nvert + iedg * (ideg - 1);
      // Edge already done
      if(msh.poi2tag[0][ent2poi(ielem,irnk0)] >= msh.tag[0]) continue;
      int irnk1 = nvert + (iedg + 1) * (ideg - 1);

      // One or the other will be the first node of any tet's this shared edge
      msh.poi2tag[0][ent2poi[ielem][irnk0  ]] = msh.tag[0];
      msh.poi2tag[0][ent2poi[ielem][irnk1-1]] = msh.tag[0];


      int ip1 = ent2poi(ielem,lnoed[iedg][0]);
      int ip2 = ent2poi(ielem,lnoed[iedg][1]);

      for(int inode = irnk0; inode < irnk1; inode++){
        int ipoin = ent2poi(ielem,inode);

        int ii1 = ordelt[ideg][inode][lnoed[iedg][0]];
        int ii2 = ordelt[ideg][inode][lnoed[iedg][1]];

        for(int kk = 0; kk < ndimn; kk++){
          msh.coord[ipoin][kk] = msh.coord[ipoin][kk] + 
                      isig*( ii1*msh.coord[ip1][kk] + ii2*msh.coord[ip2][kk] ) / (double)ideg;
        }
      }
    }
  }
  if constexpr(ideg > 2){
    if(ndimn == 2) METRIS_THROW_MSG(TODOExcept(), "Implement face interior Deltas for dimension 2");

    for(int ifac = 0; ifac < nfacl; ifac++){
      int irnk0 = nvert + nedgl * (ideg - 1) 
                    + ifac * facnpps[ideg-2];
      if(msh.poi2tag(0,irnk0) >= msh.tag[0]) continue;
      int irnk1 = nvert + nedgl * (ideg - 1) 
                    + (ifac + 1) * facnpps[ideg-2];
  
      int ip1 = ent2poi(ielem,lnofa3[ifac][0]);
      int ip2 = ent2poi(ielem,lnofa3[ifac][1]);
      int ip3 = ent2poi(ielem,lnofa3[ifac][2]);
  
      for(int inode = irnk0; inode < irnk1; inode++){
        int ipoin = ent2poi(ielem,inode);
        msh.poi2tag(0,ipoin) = msh.tag[0];
  
        int ii1 = ordelt[ideg][inode][lnofa3[ifac][0]];
        int ii2 = ordelt[ideg][inode][lnofa3[ifac][1]];
        int ii3 = ordelt[ideg][inode][lnofa3[ifac][2]];
  
        for(int kk = 0; kk < ndimn; kk++){
          msh.coord[ipoin][kk] = msh.coord[ipoin][kk] +
            isig*( ii1*msh.coord[ip1][kk] + ii2*msh.coord[ip2][kk] + ii3*msh.coord[ip3][kk] ) / (double)ideg;
        }
      }
    }
  }
  if constexpr(ideg > 3){
    METRIS_THROW_MSG(TODOExcept(),"Quickly implement face Delta computation in curveMeshOffsets")
  }

}

#define BOOST_PP_LOCAL_MACRO(n)\
template void eltgapsbezconv<2, n>(MeshBase &msh, int ielem, int isig, bool confirm);\
template void eltgapsbezconv<3, n>(MeshBase &msh, int ielem, int isig, bool confirm);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


template <int ndimn, int ideg>
void eltgaps2bez(MeshBase &msh, int ielem){
  if constexpr(ideg == 1) return;
  msh.tag[0]++;
  eltgapsbezconv<ndimn,ideg>(msh,ielem,1,true);
}


template <int ndimn, int ideg>
void eltbez2gaps(MeshBase &msh, int ielem){
  if constexpr(ideg == 1) return;
  msh.tag[0]++;
  eltgapsbezconv<ndimn,ideg>(msh,ielem,-1,true);
}


template <int ideg>
void gaps2bez(MeshBase &msh){
  if constexpr(ideg == 1) return;
  msh.tag[0]++;

  if(msh.idim == 2){
    for(int iface = 0; iface < msh.nface; iface++){
      if(isdeadent(iface,msh.fac2poi)) continue;
      eltgapsbezconv<2,ideg>(msh,iface,1,true);
    }
  }else if(msh.idim == 3){
    for(int ielem = 0; ielem < msh.nelem; ielem++){
      if(isdeadent(ielem,msh.tet2poi)) continue;
      eltgapsbezconv<3,ideg>(msh,ielem,1,true);
    }
  }else{
    METRIS_THROW_MSG(WArgExcept(), "gaps2bez unsupported dim "<<msh.idim);
  }
}


template <int ideg>
void bez2gaps(MeshBase &msh){
  if constexpr(ideg == 1) return;
  msh.tag[0]++;
  if(msh.idim == 2){
    for(int iface = 0; iface < msh.nface; iface++){
      if(isdeadent(iface,msh.fac2poi)) continue;
      eltgapsbezconv<2,ideg>(msh,iface,-1,true);
    }
  }else if(msh.idim == 3){
    for(int ielem = 0; ielem < msh.nelem; ielem++){
      if(isdeadent(ielem,msh.tet2poi)) continue;
      eltgapsbezconv<3,ideg>(msh,ielem,-1,true);
    }
  }else{
    METRIS_THROW_MSG(WArgExcept(), "gaps2bez unsupported dim "<<msh.idim);
  }
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void bez2gaps<n>(MeshBase &msh);\
template void gaps2bez<n>(MeshBase &msh);\
template void eltbez2gaps<2, n>(MeshBase &msh, int ielem);\
template void eltbez2gaps<3, n>(MeshBase &msh, int ielem);\
template void eltgaps2bez<2, n>(MeshBase &msh, int ielem);\
template void eltgaps2bez<3, n>(MeshBase &msh, int ielem);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


} // End namespace
