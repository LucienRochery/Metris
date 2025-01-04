//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../MetricField/MetricField.hxx"
#include "../Mesh/MeshBase.hxx"


namespace Metris{


template<int gdim, int tdim, int mshdeg, int tardeg>
void MetricFieldFE::getMetNodes(int ientt, double *metnod) const{
  METRIS_ENFORCE_MSG(ibasis != FEBasis::Undefined, "Metric was not initialized.");

  static_assert((tdim <= 3 && "## gen_nodes_degelev TOPO DIM > 3 NOT IMPLEMENTED \n"));
  METRIS_ASSERT(gdim == msh.idim);

  constexpr int nnmet = (gdim*(gdim+1))/2;
  METRIS_ENFORCE_MSG(nnmet == getnnmet(), "getMetNodes called with wrong gdim");

  if(this->ibasis != FEBasis::Lagrange) METRIS_THROW_MSG(TODOExcept(), 
    "Implement getMedNodes for Bézier (FE)")

  constexpr int npptar = (tdim == 1) ? edgnpps[tardeg]
                       :((tdim == 2) ? facnpps[tardeg]
                                     : tetnpps[tardeg]);

  constexpr auto ordent = ORDELT(tdim);

  const intAr2& ent2poi = (tdim == 1) ? msh.edg2poi
                        :((tdim == 2) ? msh.fac2poi
                                      : msh.tet2poi);

  double bary[tdim+1];
  //if(this->ibasis == FEBasis::Lagrange){
  for(int irnk = 0; irnk < npptar; irnk++){
    for(int ii = 0; ii < tdim+1 ;ii++) {
      bary[ii] = ordent[tardeg][irnk][ii] / (1.0*tardeg);
    }
    if constexpr(tdim == 1) eval1<nnmet,mshdeg>(rfld,ent2poi[ientt],this->ibasis,DifVar::None,DifVar::None,bary,&metnod[nnmet*irnk],NULL,NULL);
    else if     (tdim == 2) eval2<nnmet,mshdeg>(rfld,ent2poi[ientt],this->ibasis,DifVar::None,DifVar::None,bary,&metnod[nnmet*irnk],NULL,NULL);
    else                    eval3<nnmet,mshdeg>(rfld,ent2poi[ientt],this->ibasis,DifVar::None,DifVar::None,bary,&metnod[nnmet*irnk],NULL,NULL);

    #ifndef NDEBUG
      double nrm0 = 0;
      for(int ii = 0; ii < nnmet; ii++) nrm0+=abs(metnod[nnmet*irnk+ii]);
      METRIS_ASSERT(nrm0 < 1.0e10);
    #endif
  }
}



template<int gdim, int tdim, int mshdeg, int tardeg>
void MetricFieldAnalytical::getMetNodes(int ientt, double *metnod) const{
  METRIS_ENFORCE_MSG(ibasis != FEBasis::Undefined, "Metric was not initialized.");
  static_assert((tdim <= 3 && "## gen_nodes_degelev TOPO DIM > 3 NOT IMPLEMENTED \n"));

  constexpr int nnmet = (gdim*(gdim+1))/2;
  METRIS_ASSERT(gdim == msh.idim);

  if(this->ibasis != FEBasis::Lagrange) METRIS_THROW_MSG(TODOExcept(), 
    "Implement getMedNodes for Bézier (Analytical)")
    
  constexpr int npptar = (tdim == 1) ? edgnpps[tardeg]
                       :((tdim == 2) ? facnpps[tardeg]
                                     : tetnpps[tardeg]);


  constexpr auto ordent = ORDELT(tdim);
  //constexpr auto ordent = []() -> auto {
  //  if constexpr(tdim == 1) return ordedg.s;
  //  else if     (tdim == 2) return ordfac.s;
  //  else                    return ordtet.s;
  //}();

  const intAr2& ent2poi = (tdim == 1) ? msh.edg2poi
                        :((tdim == 2) ? msh.fac2poi
                                      : msh.tet2poi);

  constexpr auto eval = tdim == 1 ? eval1<gdim,mshdeg> 
                      : tdim == 2 ? eval2<gdim,mshdeg> 
                                  : eval3<gdim,mshdeg> ;

  double bary[tdim+1];
  double coop[gdim];
  for(int irnk = 0; irnk < npptar; irnk++){
    for(int ii = 0; ii < tdim+1 ;ii++){
      bary[ii] = ordent[tardeg][irnk][ii] / (1.0*tardeg);
    }
    eval(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                                 DifVar::None,bary,coop,NULL,NULL);
    anamet(NULL,coop,scale,0,&metnod[nnmet*irnk],NULL);
  }

  // Convert to log if log format expected
  if(this->ispace == MetSpace::Log){
    for(int irnk = 0; irnk < npptar; irnk++){
      getlogmet_inp<gdim,double>(&metnod[nnmet*irnk]);
    }
  }
}



#include <src/MetricField/MetricField_getMetNodes.ixx>


}// End namespace