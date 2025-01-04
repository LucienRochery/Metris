//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../MetrisRunner/MetrisRunner.hxx"

#include "../Mesh/Mesh.hxx"

#include "../API/MetrisAPI.hxx"
#include "../metris_options.hxx"
#include "../aux_exceptions.hxx"
#include "../msh_checktopo.hxx"

namespace Metris{


MetrisRunner::MetrisRunner(int argc, char** argv, bool isilent) : 
opt(argc,argv),
param(opt){
  METRIS_ENFORCE_MSG(!(opt.count("met") && opt.count("anamet")),"Contradictory options: -back or -met and -anamet");
  constructorCommon(NULL,NULL);
}

MetrisRunner::MetrisRunner(MetrisAPI *data_front, MetrisAPI *data_back, 
                           MetrisParameters &param_): 
opt(),
param(param_){
  constructorCommon(data_front,data_back);

}

void MetrisRunner::constructorCommon(MetrisAPI *data_front, MetrisAPI *data_back){

  //hookedAPI = NULL;

  if(param.iverb >= 1){

    std::cout<<"\n\n"
    "Metris: high-order metric-based non-manifold tetrahedral remesher\n"
    "Copyright (C) 2023-2024, Massachusetts Institute of Technology\n"
    "Licensed under The GNU Lesser General Public License, version 2.1\n"
    "See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n";

    #ifndef NDEBUG
    std::cout<<"Debug build.\n";
    #endif

    std::cout<<"\n\n";
    
    // Here go Metris prints
    if(param.dbgfull){
      printf("\n\n##################################################\n");
      printf("### FULL DEBUG -> VERY EXPENSIVE ! -dbgfull option\n");
      printf("##################################################\n\n\n");
    }
  }
  if(param.iverb >= 1){
    const char *OCCrev;
    int eg_imajor, eg_iminor;
    EG_revision(&eg_imajor, &eg_iminor, &OCCrev);
    printf("Compiled with EGADS version %d.%d\n",eg_imajor,eg_iminor);
    printf("              OCC revision: %s\n\n",OCCrev);
  }
  
  if(data_front == NULL && data_back == NULL && !param.inpBack && !param.inpMesh) METRIS_THROW_MSG(WArgExcept(),
    "No meshes provided either through data or files");

  if(data_back == NULL && !param.inpBack){
    if(data_front){
      data_back = data_front; 
      data_front = NULL;
    }
  }

  if(data_front == data_back) data_front = NULL;


  this->metricFE = !param.anaMet;
  if(!metricFE){
    if(param.iverb >= 2){
      printf(" - Initialization as MetricFieldAnalytical\n");
    }
    msh_g = (MeshBase *) new Mesh<MetricFieldAnalytical>;
    iniMetris<MetricFieldAnalytical>(data_front,data_back);
  }else{
    if(param.iverb >= 2){
      printf(" - Initialization as metricFE\n");
    }
    msh_g = (MeshBase *) new Mesh<MetricFieldFE        >;
    iniMetris<MetricFieldFE        >(data_front,data_back);
  }

  if(param.dbgfull) check_topo(*msh_g, msh_g->nbpoi, msh_g->npoin, msh_g->nedge, msh_g->nface, msh_g->nelem,0);
}



template<class MetricFieldType>
void MetrisRunner::iniMetris(MetrisAPI *data_front, MetrisAPI *data_back){
  
  Mesh<MetricFieldType> &msh = *( (Mesh<MetricFieldType>*) msh_g);

  //int idim = 0;
  //int idimf, idimb;
  //if(data_front != NULL){
  //  idimf = data_front -> idim;
  //  idim = idimf;
  //}
  //if(data_back != NULL){
  //  idimb = data_back -> idim;
  //  idim = idimb;
  //}
  //if(data_front != NULL && data_back != NULL) METRIS_ASSERT(idimf == idimb);
  //METRIS_ASSERT(idim == 2 || idim == 3);

  //bak.set_gdim(idim);
  //msh.set_gdim(idim);

  bak.initialize(data_back, param);
  bak.setBasis(FEBasis::Bezier);
  bak.met.setSpace(MetSpace::Log);
  bak.met.setBasis(FEBasis::Lagrange);


  msh.initialize(data_front, bak, param);
  msh.setBasis(FEBasis::Bezier);
  msh.met.setSpace(MetSpace::Exp);
  msh.met.setBasis(FEBasis::Lagrange);

  if(param.dbgfull) check_topo(msh);
  

  //set_array_debugids<MetricFieldType>();
}

template void MetrisRunner::iniMetris<MetricFieldFE>(
  MetrisAPI *data_front, MetrisAPI *data_back);
template void MetrisRunner::iniMetris<MetricFieldAnalytical>(
  MetrisAPI *data_front, MetrisAPI *data_back);


template<class MetricFieldType>
void MetrisRunner::set_array_debugids(){
  //Add a debug id as needed
//  msh_g->edg2fac.set_dbgid(1);
}

template void MetrisRunner::set_array_debugids<MetricFieldFE>();
template void MetrisRunner::set_array_debugids<MetricFieldAnalytical>();

} // end namespace






