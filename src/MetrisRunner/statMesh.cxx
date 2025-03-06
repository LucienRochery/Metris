//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "MetrisRunner.hxx"

#include "../Mesh/Mesh.hxx"

#include "../msh_lenedg.hxx"
#include "../aux_histogram.hxx"
#include "../quality/msh_metqua.hxx"
#include "../low_ccoef.hxx"
#include "../utils/mprintf.hxx"


namespace Metris{

void MetrisRunner::statMesh(MeshStat* stat){
  if(this->metricFE){
    statMesh0<MetricFieldFE>(stat);
  }else{
    statMesh0<MetricFieldAnalytical>(stat);
  }
}

template<class MFT>
void MetrisRunner::statMesh0(MeshStat* stat){


  Mesh<MFT> &msh = *( (Mesh<MFT>*) msh_g );
  GETVDEPTH(msh);

  if(!DOPRINTS1() && stat == NULL) return;
  
  msh.cleanup();

  intAr2 ilned;
  ilned.set_n(0);
  dblAr1 rlned;
  dblAr1 lenbds = {1.0/sqrt(2), sqrt(2)};
  getLengthEdges<MFT>(msh,ilned,rlned,LenTyp::Quad);

  if(stat != NULL){
    stat->minlen = 1.0e30;
    stat->maxlen = -1.0e30;
    stat->avglen = 0.0;
    int nunit = 0;
    for(int ii = 0; ii < ilned.get_n(); ii++){
      double len = rlned[ii];
      if(len >= 1.0/sqrt(2) && len <= sqrt(2)) nunit++;
      stat->avglen += len;
      stat->minlen = MIN(len, stat->minlen);
      stat->maxlen = MAX(len, stat->maxlen);
    }
    stat->avglen /= ilned.get_n();
    stat->pctunit = 100*((double)nunit) / ilned.get_n();
  }

  if(DOPRINTS1()){
    print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length (quadrature)");
  
    if(DOPRINTS2())
      getLengthEdges<MFT>(msh,ilned,rlned,LenTyp::GeoSiz);
      print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length (geometric)");{
    }
  }



  //getLengthEdges<MFT>(msh,ilned,rlned,LenTyp::GeoSiz);
  //print_histogram(rlned,IntrpTyp::Linear,lenbds,"l","Edge length (geom interp)");

  double qmin, qmax, qavg;
  bool iinva;
  dblAr1 lquae;
  dblAr1 dum = {1.0e-8, 0.1};
  if(msh.idim == 2){
    getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
  }else{
    getmetquamesh<MFT,3,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
  }
  if(stat != NULL){
    stat->minqua = qmin;
    stat->maxqua = qmax;
    stat->avgqua = qavg;
  }
  print_histogram(msh,lquae,IntrpTyp::Geometric,dum,"q","Element quality (As P1)");

  if(msh.curdeg > 1){
    if(msh.idim == 2){
      getmetquamesh<MFT,2,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    }else{
      getmetquamesh<MFT,3,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    }
    if(stat != NULL){//overwrite
      stat->minqua = qmin;
      stat->maxqua = qmax;
      stat->avgqua = qavg;
    }
    print_histogram(msh,lquae,IntrpTyp::Geometric,dum,"q","Element quality (As Pk)");

    int tdim = msh.get_tdim();
    if(tdim == msh.idim){
      int nentt = msh.nentt(tdim);
      intAr2 &ent2poi = msh.ent2poi(tdim);

      int jdeg = msh.idim * (msh.curdeg - 1);
      int ncoef = tdim == 2 ? getnnod2(jdeg) : getnnod3(jdeg);
      //dblAr1 ccoef(ncoef);
      //ccoef.set_n(ncoef);

      lquae.allocate(nentt*ncoef);
      lquae.set_n(nentt*ncoef);

      bool iinva;
      for(int ientt = 0; ientt < nentt; ientt++){
        if(isdeadent(ientt,ent2poi)) continue;
        CT_FOR0_INC(2,3,idim){if(idim == tdim){
          CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
            getsclccoef<idim,idim,ideg>(msh,ientt,NULL,&lquae[ientt*ncoef],&iinva);
          }}CT_FOR1(ideg);
        }}CT_FOR1(idim);
      }// for ientt
      
      dum[0] = msh.param->jtol;
      dum[1] = 1;
      print_histogram(msh,lquae,IntrpTyp::Geometric,dum,"J","Scaled Jacobian");

    }// if tdim 

  }
}

//#define BOOST_PP_LOCAL_MACRO(n)
template void MetrisRunner::statMesh0<MetricFieldAnalytical>(MeshStat *stat);
template void MetrisRunner::statMesh0<MetricFieldFE        >(MeshStat *stat);
//#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
//#include BOOST_PP_LOCAL_ITERATE()


}