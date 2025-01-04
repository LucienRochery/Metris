//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
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

  #ifndef NDEBUG
  double metcpy[6];
  #endif

  for(int ipoin = 0; ipoin < msh.npoin ;ipoin++){
    if(msh.poi2ent(ipoin,0) < 0) continue;
    #ifndef NDEBUG
    int nnmet = (ndimn*(ndimn+1))/2;
    for(int ii = 0; ii < nnmet; ii++){
      metcpy[ii] = metfld[ipoin][ii];
    }
    double metcpy[6];
    #endif
    try{
      getlogmet_inp<ndimn,double>(metfld[ipoin]);
    }catch(const MetrisExcept &e){
      std::cout<<"logmat failed for ipoin = "<<ipoin<<" poi2ent = "<<msh.poi2ent(ipoin,0)<<"\n";
      #ifndef NDEBUG
      printf("input was: ");
      dblAr1(nnmet,metcpy).print();
      #endif
      Mesh<MetricFieldAnalytical> *msh2 = (Mesh<MetricFieldAnalytical> *) (&msh);
      int ieleg = msh.poi2ent(ipoin,0);
      double metl[6];
      msh2->met.getMetPhys(AsDeg::Pk,DifVar::None,MetSpace::Exp,&ieleg,msh.coord[ipoin],metl,NULL);
      std::cout<<"recomputed metric: ";
      for(int ii = 0 ; ii < 6; ii++) printf("%20.16e ",metl[ii]); 
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
    }catch(const MetrisExcept &e){
      std::cout<<"expmat failed for ipoin = "<<ipoin<<" poi2ent = "<<msh.poi2ent(ipoin,0)<<"\n";
      throw(e);
    }
  }
}

template void setLogMetMesh0<2>(const MeshBase &msh, dblAr2 &metfld);
template void setLogMetMesh0<3>(const MeshBase &msh, dblAr2 &metfld);

template void setExpMetMesh0<2>(const MeshBase &msh, dblAr2 &metfld);
template void setExpMetMesh0<3>(const MeshBase &msh, dblAr2 &metfld);

} // End namespace
