//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php
#ifndef __METRIS_CAD_INFO__
#define __METRIS_CAD_INFO__


#include "../types.hxx"
#include <egads.h>


namespace Metris{

struct MetrisParameters;
class MeshBase; 

class CADInfo{
public:
  CADInfo(){
    EGADS_context = NULL;
    EGADS_model   = NULL;
    body          = NULL;
    ialloc        = false;
  }
  //CADInfo(ego EGADS_model_);

  ~CADInfo();
  
  // Return whether useable 
  bool operator()() const {return EGADS_model != NULL;}

  void free();

  CADInfo &operator=(const CADInfo &src);

  //Get the EGADS body and populate ncadfa/ed/no & cad2fac/edg/nod 
  void iniEGADSModel();
  // EGADS_context_ can be null: in that case, a new one is created, and 
  // the EGADS_model_ is hard-copied. 
  // Otherwise, EGADS_model_ is simply referenced, as well as the context. 
  void setModel(ego EGADS_context_, ego EGADS_model_); // this calls iniEGADSModel
  // egadslite:
  void setModel(size_t nbyte, char* stream); // this calls iniEGADSModel
  void iniCADLink(const MetrisParameters &param, MeshBase &msh, int nbpo0); 

public:
  ego EGADS_model;
  ego body;
  int ncadfa, ncaded, ncadno, ncadlp;
  egoAr1 cad2fac, cad2edg, cad2nod, cad2lop; 
  bool ialloc; // to manage copies 

protected:
  ego EGADS_context;
};


}//end namespace

#endif
