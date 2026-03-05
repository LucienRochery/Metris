//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../MetrisRunner/MetrisRunner.hxx"

#include "../Mesh/Mesh.hxx"

#include "../API/MetrisAPI.hxx"
#include "../metris_options.hxx"
#include "../aux_exceptions.hxx"
#include "../msh_checktopo.hxx"
#include "../utils/mprintf.hxx"
#include <filesystem>
namespace Metris{


MetrisRunner::MetrisRunner(int argc, char** argv) :
opt(argc,argv),
param_(opt),
param(&param_){
  METRIS_ENFORCE_MSG(!(opt.count("met") && opt.count("anamet")),"Contradictory options: -back or -met and -anamet");
  GETVDEPTH(param);
  if(DOPRINTS1()){
    fmt::print("Call: ");
    for(int ii = 1; ii < argc; ii++){
      fmt::print(" {} ",argv[ii]);
    }
    fmt::print("\n");
  }
  constructorCommon(NULL,NULL);
}

MetrisRunner::MetrisRunner(MetrisAPI *data_front, MetrisAPI *data_back,
                           MetrisParameters &param__):
opt(),
param_(param__),
param(&param_){
  MetrisParameters param_default;

  GETVDEPTH(param);
  CPRINTF1("Call: \n");
  if(DOPRINTS1()) ((MetrisParametersData&)param_default).printDifference(param__,"default",param->logFile);
  if(param->inpMesh){
    CPRINTF1("-- Input Mesh = {}\n", param->meshFileName);
  }else if(data_front){
    CPRINTF1("-- Input Mesh from data\n");
  }
  if(param->inpCAD){
    CPRINTF1("-- Input CAD  = {}\n", param->cadFileName);
  }
  if(param->inpBack){
    CPRINTF1("-- Input Back = {}\n", param->backFileName);
  }else if(data_back){
    CPRINTF1("-- Input Back from data\n");
  }
  if(param->inpMet){
    CPRINTF1("-- Input Met  = {}\n", param->metFileName);
  }
  constructorCommon(data_front,data_back);

}

void MetrisRunner::constructorCommon(MetrisAPI *data_front, MetrisAPI *data_back){

  param_.checkParameters();

  GETVDEPTH((&param_));

  if(DOPRINTS1()){

    MPRINTF("\n\n"
    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    "Metris: high-order metric-based non-manifold tetrahedral remesher\n"
    "Copyright (C) 2023-2025, Massachusetts Institute of Technology\n"
    "Licensed under The GNU Lesser General Public License, version 2.1\n"
    "See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n"
    );

    MPRINTF("Metris Git repository URL " METRIS_GIT_URL "\n");
    MPRINTF("commit SHA1 {} ",METRIS_GIT_COMMIT_HASH);

    #ifndef NDEBUG
    MPRINTF("\nDebug build.\n");
    #endif

    bool METRIS_USE_PETSC = false, use_absl=false;


    #ifdef METRIS_USE_PETSC
    METRIS_USE_PETSC = true;
    #endif

    #ifdef USE_ABSL
    use_absl = true;
    #endif

    if(METRIS_USE_PETSC || use_absl){
      MPRINTF("Compiled with libraries ");
      if(METRIS_USE_PETSC) fmt::print(LOGFILE__,"petsc ");
      if(use_absl) fmt::print(LOGFILE__,"absl");
      fmt::print(LOGFILE__,"\n");
    }

    if(param_.dbgfull){
      fmt::print(LOGFILE__,"\n\n");
      MPRINTF("##################################################\n");
      MPRINTF("### FULL DEBUG -> VERY EXPENSIVE ! -dbgfull option\n");
      MPRINTF("##################################################\n");
      fmt::print(LOGFILE__,"\n\n");
    }

    const char *OCCrev;
    int eg_imajor, eg_iminor;
    EG_revision(&eg_imajor, &eg_iminor, &OCCrev);
    fmt::print(LOGFILE__,"\n");
    MPRINTF("\nCompiled with EGADS version {}.{}\n",eg_imajor,eg_iminor);
    MPRINTF("              OCC revision: {}\n\n",OCCrev);
    MPRINTF("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
  }

  if(data_front == NULL && data_back == NULL && !param_.inpBack && !param_.inpMesh) METRIS_THROW_MSG(
    "No meshes provided either through data or files");

  if(data_back == NULL && !param_.inpBack){
    if(data_front){
      data_back = data_front;
      data_front = NULL;
    }
  }

  if(data_front == data_back) data_front = NULL;

  this->metricFE = !param_.anaMet;
  if(!metricFE){
    CPRINTF2(" - Initialization as MetricFieldAnalytical\n");
    msh_g = (MeshBase *) new Mesh<MetricFieldAnalytical>;
    iniMetris<MetricFieldAnalytical>(data_front,data_back);
  }else{
    CPRINTF2(" - Initialization as metricFE\n");
    msh_g = (MeshBase *) new Mesh<MetricFieldFE        >;
    iniMetris<MetricFieldFE        >(data_front,data_back);
  }

  if(param_.dbgfull)
    check_topo(*msh_g, msh_g->nbpoi, msh_g->npoin, msh_g->nedge, msh_g->nface, msh_g->nelem,0);

#ifdef OUTPUTTIMEANDUNITINFO
  printUnit = false;
  std::string outputTimeUnitFile = "timeAndUnitInfo.txt";
  bool fileExists = std::filesystem::exists(outputTimeUnitFile);

  foutputTimeUnit.open(outputTimeUnitFile, std::ios::app);
  METRIS_ASSERT_MSG(foutputTimeUnit.good(), "Error opening file: " + outputTimeUnitFile);

  if (!fileExists){

    foutputTimeUnit << std::setw(30) << "Unit %"
                    << std::setw(30) << "Time (s)"
                    << std::endl;
  }
#endif
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

  bak.initialize(data_back, param_);
  bak.setBasis(FEBasis::Bezier);
  bak.met.setSpace(MetSpace::Log);
  bak.met.setBasis(FEBasis::Lagrange);


  msh.initialize(data_front, bak, param_);
  msh.setBasis(FEBasis::Bezier);
  msh.met.setSpace(MetSpace::Exp);
  msh.met.setBasis(FEBasis::Lagrange);

  if(param_.dbgfull) check_topo(msh,0);


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






