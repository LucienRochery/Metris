//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_METQUA_D__
#define __METRIS_LOW_METQUA_D__


#include "quafun.hxx"
#include "../Mesh/MeshFwd.hxx"


namespace Metris{

enum class FEBasis;
enum class DifVar;
enum class AsDeg;

// Integrated quality 
template <class MetricFieldType, int gdim,
          QuaFun iquaf = QuaFun::Distortion, typename ftype = double>
ftype d_metqua(Mesh<MetricFieldType> &msh,
               AsDeg asdmsh, AsDeg asdmet, 
               int ielem, int power, 
               int ivar,
               FEBasis dofbas, DifVar idifmet, 
               ftype*__restrict__ dquael, ftype*__restrict__ hquael, 
               int pnorm = 1, double difto = 1);

template <class MetricFieldType, int gdim,
          QuaFun iquaf = QuaFun::Distortion, typename ftype = double>
ftype d_metqua0(Mesh<MetricFieldType> &msh,
                AsDeg asdmsh, AsDeg asdmet, 
                const int *ent2poi, int power, 
                int ivar,
                FEBasis dofbas, DifVar idifmet, 
                ftype*__restrict__ dquael, ftype*__restrict__ hquael, 
                int pnorm = 1, double difto = 1);


// Integrated quality 
template <class MetricFieldType, int gdim, int mshdeg, AsDeg asdmet, typename ftype = double>
ftype D_metqua(Mesh<MetricFieldType> &msh, int ielem, int power, 
               FEBasis dofbas, DifVar idifmet, 
               ftype*__restrict__ dquael, ftype*__restrict__ hquael, 
               int pnorm = 1, double difto = 1);

template <class MetricFieldType, int gdim, int mshdeg, AsDeg asdmet, typename ftype = double>
ftype D_metqua(Mesh<MetricFieldType> &msh, const int *ent2poi, int power, 
               FEBasis dofbas, DifVar idifmet, 
               ftype*__restrict__ dquael, ftype*__restrict__ hquael, 
               int pnorm = 1, double difto = 1);

// Differentiated wrt all nodes and coordinates 
template <class MFT, int gdim, int mshdeg, AsDeg asdmet, typename ftype>
ftype D_quafun_distortion(Mesh<MFT> &msh, 
                  const int* ent2poi,
                  const double*__restrict__ bary, 
                  int power, 
                  FEBasis dofbas, 
                  DifVar idifmet, 
                  ftype*__restrict__ dquael, 
                  ftype*__restrict__ hquael);







template<class MFT, int gdim, int tdim, class ftype=double>
std::function<ftype(Mesh<MFT>&,AsDeg,AsDeg,int,int,int,FEBasis,DifVar,
                    ftype*__restrict__,ftype*__restrict__,int,double)>
get_d_quafun(QuaFun iquaf){
  if(iquaf == QuaFun::Distortion){
    return d_metqua<MFT,gdim,QuaFun::Distortion,ftype>;
  }else if(iquaf == QuaFun::Unit){
    return d_metqua<MFT,gdim,QuaFun::Unit,ftype>;
  }else{
    METRIS_THROW_MSG(TODOExcept(),"cf quafun_")
  }
}



//template <class MetricFieldType, int ideg, int ivar, typename ftype = double>
//struct metqua3_d{
//	metqua3_d() = delete;
//	metqua3_d(const Mesh<MetricFieldType> &msh, int ielem, int power,
//						const double*  __restrict__ dmetvar, 
//						const double*  __restrict__ dpoivar, 
//						SANS::SurrealS<nvar,ftype> *qutet);
//};


//template <class MetricFieldType, int ideg>
//struct metqua3_shell_d{
//	metqua3_shell_d() = delete;
//	metqua3_shell_d(const Mesh<MetricFieldType> &msh, int ipoin, int iele0, int power, int mshell, 
//	                int* __restrict__ nshell, int* __restrict__ lshell,
//	                const double* __restrict__ dmetvar, 
//	                const double* __restrict__ dpoivar, 
//	                SANS::SurrealS<nvar,double>*  qushe);
//};

} // End namespace

#endif
