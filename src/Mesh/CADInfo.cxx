//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Mesh/CADInfo.hxx"
#include "../MetrisRunner/MetrisParameters.hxx" 
#include "../Mesh/MeshBase.hxx"
#include "../Boundary/msh_inisurf.hxx"
#include "../aux_utils.hxx"
#include "../mprintf.hxx"
#include "../io_libmeshb.hxx"

namespace Metris{



//CADInfo::CADInfo(ego EGADS_model_) : EGADS_model(EGADS_model_) {
//  // We can leave EGADS_context = NULL
//  iniEGADSModel();
//}

void CADInfo::iniEGADSModel(){
  int ierro;
  ego geom;
  int oclass,mtype,nbody,*dum;
  ego *bodies;
  ierro = EG_getTopology(EGADS_model,&geom,&oclass,&mtype,NULL,&nbody,&bodies,&dum);
  if(ierro != 0){
    print_EGADS_error("EG_getTopology",ierro);
    METRIS_THROW(TopoExcept());
  }
  if(nbody == 0) METRIS_THROW_MSG(TopoExcept(),"CAD has nbody = "<<nbody);
  if(nbody  > 1) METRIS_THROW_MSG(TopoExcept(),"> 1 BODIES NOT SUPPORTED YET ");

  body = bodies[0];

  ego *buff; 
  ierro = EG_getBodyTopos(body,NULL,FACE,&ncadfa,&buff);
  if(ierro != 0){
    print_EGADS_error("EG_getBodyTopos (FACE)",ierro);
    METRIS_THROW(TopoExcept());
  }
  if(ncadfa == 0){
    printf("WARNING: Body with no faces !\n");
  }else{
    //printf("  body has %d faces \n",ncadfa);
  }
  // Non-owning set (EGADS frees buff?)
  cad2fac.set_buffer(ncadfa,buff); 
  cad2fac.set_n(ncadfa);


  ierro = EG_getBodyTopos(body,NULL,EDGE,&ncaded,&buff);
  if(ierro != 0){
    print_EGADS_error("EG_getBodyTopos (EDGE)",ierro);
    METRIS_THROW(TopoExcept());
  }
  if(ncaded == 0){
    printf("WARNING: Body with no edges !\n");
  }else{
    //printf("  body has %d edges \n",ncaded);
  }
  // Non-owning set (EGADS frees buff?)
  cad2edg.set_buffer(ncaded,buff); 
  cad2edg.set_n(ncaded);


  ierro = EG_getBodyTopos(body,NULL,LOOP,&ncadlp,&buff);
  if(ierro != 0){
    print_EGADS_error("EG_getBodyTopos (LOOP)",ierro);
    METRIS_THROW(TopoExcept());
  }
  if(ncaded == 0){
    printf("WARNING: Body with no loops !\n");
  }else{
    //printf("  body has %d edges \n",ncaded);
  }
  // Non-owning set (EGADS frees buff?)
  cad2lop.set_buffer(ncadlp,buff); 
  cad2lop.set_n(ncadlp);

  ierro = EG_getBodyTopos(body,NULL,NODE,&ncadno,&buff);
  if(ierro != 0){
    print_EGADS_error("EG_getBodyTopos (NODE)",ierro);
    METRIS_THROW(TopoExcept());
  }
  if(ncadno == 0){
    printf("WARNING: Body with no nodes !\n");
  }else{
    //printf("  body has %d nodes \n",ncadno);
  }
  cad2nod.set_buffer(ncadno,buff); 
  cad2nod.set_n(ncadno);
  ialloc = true;
}


void CADInfo::setModel(ego EGADS_context_, ego EGADS_model_){
  //if(EGADS_context_ != NULL){
    EGADS_context = EGADS_context_;
    EGADS_model   = EGADS_model_;
  //}else{
  //  // Create new context and hard-copy EGADS_model_ 
  //  int ierro = EG_open(&EGADS_context);
  //  if(ierro != EGADS_SUCCESS){
  //    print_EGADS_error("CADInfo::setModel: EG_open",ierro);
  //    METRIS_THROW(TopoExcept());
  //  }
  //  ierro = EG_contextCopy(EGADS_context, EGADS_model_, &EGADS_model);
  //  if(ierro != EGADS_SUCCESS){
  //    print_EGADS_error("CADInfo::setModel: EG_contextCopy",ierro);
  //    METRIS_THROW(TopoExcept());
  //  }
  //  printf("Debug hard-copied EGADS_model\n");
  //}
  iniEGADSModel();
}


void CADInfo::setModel(size_t nbyte, char* stream){
//  METRIS_THROW_MSG(TODOExcept(), "Fix CAD stream\n");
  printf("## DEBUG STREAM AS INT:\n");
  int *ptr = (int*)stream;
  for(int ii = 0; ii < MIN(10,nbyte / sizeof(int)); ii++){
    printf("%d : %d \n",ii,*ptr);
    ptr++;
  }
  int *nbodyptr = (int*) (&stream[0] + 3*sizeof(int) + 6*sizeof(double));
  printf("## DEBUG NBODY FROM STREAM %d \n",*nbodyptr);

  int ierro = EG_open(&EGADS_context);
  if(ierro != 0){
    print_EGADS_error("EG_open",ierro);
    METRIS_THROW(TopoExcept());
  }
  ierro = EG_importModel(EGADS_context, nbyte, stream, &EGADS_model);
  if(ierro != 0){
    print_EGADS_error("EG_importModel",ierro);
    METRIS_THROW(TopoExcept());
  }


  ego geom;
  int oclass,mtype,nbody,*dum;
  ego *bodies;
  nbody = 0;
  ierro = EG_getTopology(EGADS_model,&geom,&oclass,&mtype,NULL,&nbody,&bodies,&dum);
  printf("## DEBUG READING CHECK NBODY = %d \n",nbody);


  //size_t nbyte2;
  //char* stream2;
  //ierro = EG_exportModel(EGADS_model, &nbyte2, &stream2);
  //if(ierro != 0){
  //  print_EGADS_error("EG_exportModel",ierro);
  //  METRIS_THROW_MSG(TopoExcept(),"Failed to export model to stream2.");
  //}
  //printf(" - Stream2 of size %db \n",nbyte2);


  iniEGADSModel();
}


void CADInfo::iniCADLink(const MetrisParameters &param, MeshBase &msh, int nbpo0){
  GETVDEPTH(msh);

  if(EGADS_model == NULL){

    if(param.inpCAD){ 
    /* -------------- CAD File handling -------------- */
    // Throw out exceptions as these are not fatal. 
    METRIS_TRY0(
      CPRINTF1("-- Read CAD file %s and project.\n",param.cadFileName.c_str());
      int ierro = EG_open(&EGADS_context);
      if(ierro != 0){
        print_EGADS_error("EG_open",ierro);
        METRIS_THROW(TopoExcept());
      }

      CPRINTF2(" - Start reading CAD file.\n");
      int bitFlag = 0; 
      ierro = EG_loadModel(EGADS_context,bitFlag,param.cadFileName.c_str(),&EGADS_model);
      if(ierro != 0){
        print_EGADS_error("EG_loadModel",ierro);
        METRIS_THROW_MSG(WArgExcept(),"CAD Projection will not be available");
      }
      CPRINTF2(" - Done reading CAD file.\n");

      //printf("## Remove this \n");
      //size_t nbyte;
      //char* stream;
      //ierro = EG_exportModel(EGADS_model, &nbyte, &stream);
      //if(ierro != 0){
      //  print_EGADS_error("EG_exportModel",ierro);
      //  METRIS_THROW_MSG(TopoExcept(),"Failed to export model to stream.");
      //}
      //printf("Got nbyte %zu \n", nbyte);
      //ego EGADS_model2;
      //ierro = EG_importModel(EGADS_context, nbyte, stream, &EGADS_model2);
      //if(ierro != 0){
      //  print_EGADS_error("EG_importModel",ierro);
      //  METRIS_THROW(TopoExcept());
      //}

      //ego geom;
      ////int oclass,mtype,nbody,*dum;
      //int oclass;
      //int mtype;
      //int nbody;
      //int *dum;
      //ego *bodies;
      //ierro = EG_getTopology(EGADS_model,&geom,&oclass,&mtype,NULL,&nbody,&bodies,&dum);
      //if(ierro != 0){
      //  print_EGADS_error("EG_getTopology",ierro);
      //  METRIS_THROW(TopoExcept());
      //}
      //if(nbody == 0) METRIS_THROW_MSG(TopoExcept(),"EMPTY EGADS MODEL");
      //if(nbody  > 1) METRIS_THROW_MSG(TopoExcept(),"> 1 BODIES NOT SUPPORTED YET ");

      //wait();



    )
    }

  }

  // The EGADS_model is still NULL if initially NULL and !param.inpCAD
  if(EGADS_model == NULL){
    ncadfa = -1;
    ncaded = -1;
    ncadno = 0;

    for(int iface = 0; iface < msh.nface; iface++){
      if(isdeadent(iface,msh.fac2poi)) continue;
      int iref = msh.fac2ref[iface];
      if(iref < 0) METRIS_THROW_MSG(TopoExcept(),"Even without CAD: give faces refs!! iface = "<<iface<<" iref = "<<iref);
      if(iref > ncadfa) ncadfa = iref;
    }
    for(int iedge = 0; iedge < msh.nedge; iedge++){
      if(isdeadent(iedge,msh.edg2poi)) continue;
      int iref = msh.edg2ref[iedge];
      if(iref < 0) METRIS_THROW_MSG(TopoExcept(),"Even without CAD: give edges refs! iedge = !"<<iedge<<" iref = "<<iref);
      if(iref > ncaded) ncaded = iref;
    }
    for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
      int ityp = msh.bpo2ibi(ibpoi,1);
      if(ityp == 0) ncadno++;
    }

    // Refs are 0 - n-1
    ncadfa++;
    ncaded++; 

    if(param.iverb >= 1) printf("-- Counted refs node = %d edge = %d triangle = %d \n",ncadno,ncaded,ncadfa);

  }else{

    iniEGADSModel();
    prjMeshPoints(msh, nbpo0);
    if(DOPRINTS2()) writeMesh("inisurf.meshb",msh);

  }
}





CADInfo::~CADInfo(){
  free();
}

void CADInfo::free(){
  EGADS_context = NULL;
  EGADS_model = NULL; 
  ncadno = 0;
  ncaded = 0;
  ncadfa = 0;
  ncadlp = 0;

  //// Note these are non-owning, EGADS manages actual allocs. 
  //if(ialloc){
  //  EG_free(cad2nod);
  //  cad2nod.free();
  //  EG_free(cad2edg);
  //  cad2edg.free();
  //  EG_free(cad2fac);
  //  cad2fac.free();
  //}

  ialloc = false;
}

CADInfo& CADInfo::operator=(const CADInfo &inp){
  free();

  ialloc = false;

  EGADS_context = inp.EGADS_context;
  EGADS_model   = inp.EGADS_model;
  ncadno = inp.ncadno;
  ncaded = inp.ncaded;
  ncadfa = inp.ncadfa;
  ncadlp = inp.ncadlp;

  cad2nod = inp.cad2nod;
  cad2edg = inp.cad2edg;
  cad2fac = inp.cad2fac;
  cad2lop = inp.cad2lop;

  //cfa2tag.allocate(METRIS_MAXTAGS, ncadfa, true);
  //ced2tag.allocate(METRIS_MAXTAGS, ncaded, true);
  //cno2tag.allocate(METRIS_MAXTAGS, ncadno, true);

  //inp.cfa2tag.copyTo(cfa2tag, METRIS_MAXTAGS);
  //inp.ced2tag.copyTo(ced2tag, METRIS_MAXTAGS);
  //inp.cno2tag.copyTo(cno2tag, METRIS_MAXTAGS);

  return *this;
}


}//end namespace
