//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../MetricField/msh_explogmet.hxx"
#include "../linalg/explogmet.hxx"
#include "../Mesh/MeshBase.hxx"
#include "../Mesh/Mesh.hxx"

namespace Metris{


template<int ndimn>
void setLogMetMesh0(const MeshBase &msh, dblAr2 &metfld){
  static_assert(ndimn == 2 || ndimn == 3);
  for(int ipoin = 0; ipoin < msh.npoin ;ipoin++){
    if(msh.poi2ent(ipoin,0) < 0) continue;
    try{
      getlogmet_inp<ndimn,double>(metfld[ipoin]);
    }catch(const MetrisExcept& e){
      printf("## EXCEPTION RAISED ON POINT %d COORD ",ipoin);
      dblAr1(msh.idim,msh.coord[ipoin]).print();
      for(int ii = 0; ii < MIN(10,msh.npoin); ii++){
        printf(" point %d met = ",ii);
        dblAr1((msh.idim*(msh.idim+1))/2,metfld[ii]).print();
      }
      throw(e);
    }
  }
}

template<int ndimn>
void setExpMetMesh0(const MeshBase &msh, dblAr2 &metfld){
  static_assert(ndimn == 2 || ndimn == 3);
  for(int ipoin = 0; ipoin < msh.npoin ;ipoin++){
    if(msh.poi2ent(ipoin,0) < 0) continue;
    try{
      getexpmet_inp<ndimn,double>(metfld[ipoin]);
    }catch(const MetrisExcept& e){
      printf("## EXCEPTION RAISED ON POINT %d COORD ",ipoin);
      dblAr1(msh.idim,msh.coord[ipoin]).print();
      for(int ii = 0; ii < MIN(10,msh.npoin); ii++){
        printf(" point %d met = ",ii);
        dblAr1((msh.idim*(msh.idim+1))/2,metfld[ii]).print();
      }
      throw(e);
    }
  }
}

template void setLogMetMesh0<2>(const MeshBase &msh, dblAr2 &metfld);
template void setLogMetMesh0<3>(const MeshBase &msh, dblAr2 &metfld);

template void setExpMetMesh0<2>(const MeshBase &msh, dblAr2 &metfld);
template void setExpMetMesh0<3>(const MeshBase &msh, dblAr2 &metfld);

} // End namespace
