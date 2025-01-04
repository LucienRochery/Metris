//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../MetricField/MetricField.hxx"
#include "../MetricField/msh_checkmet.hxx"

#include "../Mesh/MeshMetric.hxx"

#include "../linalg/eigen.hxx"


namespace Metris{


void checkMet(const MeshMetric<MetricFieldFE>& msh){
  if(msh.met.getSpace() == MetSpace::Log){
    std::cout<<"checkMet skip (Log) \n";
    return;
  }
      
  CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
    std::cout<<"-- checkMet call start: throw on negative eigval\n";
    double eigval[gdim];
    double eigvec[gdim*gdim];
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      if(msh.poi2ent(ipoin,0) < 0) continue;
      geteigsym<gdim,double>(msh.met[ipoin],eigval,eigvec);
      for(int ii = 0; ii < gdim; ii++) METRIS_ENFORCE(eigval[ii] > 1.0e-16);
    }
    std::cout<<"-- checkMet end all pass\n";
  }}CT_FOR1(gdim);
}


} // end namespace
