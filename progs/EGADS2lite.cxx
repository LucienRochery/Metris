//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include <src/MetrisRunner/MetrisParameters.hxx>
#include <src/metris_options.hxx>
#include <src/aux_exceptions.hxx>
#include <egads.h>
#include <src/aux_EGADSprinterr.hxx>
#include <fstream>
#include "../libs/libmeshb.hxx"
#include <string>
#include <algorithm>

using namespace Metris;

std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

int main(int argc, char** argv){
  MetrisOptions opt(argc, argv);
  MetrisParameters param(opt);


  const char *OCCrev;
  int eg_imajor, eg_iminor;
  EG_revision(&eg_imajor, &eg_iminor, &OCCrev);
  printf("Compiled with EGADS version %d.%d\n",eg_imajor,eg_iminor);
  printf("              OCC revision: %s\n\n",OCCrev);

  ego EGADS_context, EGADS_model;
  int ierro = EG_open(&EGADS_context);
  if(ierro != 0){
    print_EGADS_error("EG_open",ierro);
    METRIS_THROW_MSG(WArgExcept(),"EGADS context could not be initialized.");
  }

  size_t nbyte;
  char* stream;

  if(param.cadFileName.find(".legads") != std::string::npos){
    printf(" -- EGADSlite to EGADS \n");
    printf("  - Only read, load and get topo test \n");

    return 0;
  }


  std::cout<<" - Start reading .egads file "<<param.cadFileName<<"\n";
  int bitFlag = 0; 
  ierro = EG_loadModel(EGADS_context,bitFlag,param.cadFileName.c_str(),&EGADS_model);
  if(ierro != 0){
    print_EGADS_error("EG_loadModel",ierro);
    METRIS_THROW_MSG(WArgExcept(),"EGADS file could not be read.");
  }

  ierro = EG_exportModel(EGADS_model, &nbyte, &stream);
  if(ierro != 0){
    print_EGADS_error("EG_exportModel",ierro);
    METRIS_THROW_MSG(TopoExcept(),"Failed to export model to stream.");
  }

  printf("Stream of size %db \n",nbyte);

  std::string fname = param.cadFileName;
  //std::replace( fname.begin(), fname.end(), ".egads", ".legads"); 
  fname = ReplaceAll(fname, ".egads", ".legads");

  FILE *fp = fopen(fname.c_str(),"wb");
  METRIS_ENFORCE( fp != nullptr );
  printf("Opened filed %s \n",fname.c_str());
  fwrite(stream,sizeof(char),nbyte,fp);
  printf("Wrote byte flow to file\n");
  fclose(fp);



  fname = "toto.meshb";
  int64_t libIdx = GmfOpenMesh(fname.c_str(), GmfWrite, 3, 2);
  if(!libIdx)METRIS_THROW_MSG(WArgExcept(),
    "FILE COULDNT BE OPENED OR WRONG VERSION name = "<<fname);
  GmfWriteByteFlow(libIdx, stream, (int) nbyte);
  GmfCloseMesh(libIdx);



  EG_free(stream);

  
  return 0;
}