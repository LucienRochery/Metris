//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../MetrisRunner/MetrisRunner.hxx"

#include "../Mesh/Mesh.hxx"

#include "../msh_lenedg.hxx"
#include "../aux_histogram.hxx"
#include "../quality/msh_metqua.hxx"
#include "../low_ccoef.hxx"


namespace Metris{

void MetrisRunner::statMesh(){
  if(this->metricFE){
    statMesh0<MetricFieldFE>();
  }else{
    statMesh0<MetricFieldAnalytical>();
  }
}

template<class MFT>
void MetrisRunner::statMesh0(){
  Mesh<MFT> &msh = *( (Mesh<MFT>*) msh_g );
  
  msh.cleanup();


  intAr2 ilned;
  ilned.set_n(0);
  dblAr1 rlned;
  dblAr1 lenbds = {1.0/sqrt(2), sqrt(2)};
  getLengthEdges<MFT>(msh,ilned,rlned,LenTyp::Quad);
  print_histogram(rlned,IntrpTyp::Linear,lenbds,"l","Edge length (quadrature)");

  if(param.iverb >= 2)
    getLengthEdges<MFT>(msh,ilned,rlned,LenTyp::GeoSiz);
    print_histogram(rlned,IntrpTyp::Linear,lenbds,"l","Edge length (geometric)");{
  }



  //getLengthEdges<MFT>(msh,ilned,rlned,LenTyp::GeoSiz);
  //print_histogram(rlned,IntrpTyp::Linear,lenbds,"l","Edge length (geom interp)");

  double qmin, qmax, qavg;
  bool iinva;
  dblAr1 lquae,dum = {0.1, 0.9};
  if(msh.idim == 2){
    getmetquamesh<MFT,2,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
  }else{
    getmetquamesh<MFT,3,AsDeg::P1>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
  }
  print_histogram(lquae,IntrpTyp::Geometric,dum,"q","Element quality (As P1)");

  if(msh.curdeg > 1){
    if(msh.idim == 2){
      getmetquamesh<MFT,2,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    }else{
      getmetquamesh<MFT,3,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
    }
    print_histogram(lquae,IntrpTyp::Geometric,dum,"q","Element quality (As Pk)");

    int tdim = msh.get_tdim();
    if(tdim == msh.idim){
      int nentt = msh.nentt(tdim);
      intAr2 &ent2poi = msh.ent2poi(tdim);

      int jdeg = msh.idim * (msh.curdeg - 1);
      int ncoef = tdim == 2 ? facnpps[jdeg] : tetnpps[jdeg];
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
      print_histogram(lquae,IntrpTyp::Geometric,dum,"J","Scaled Jacobian");

    }// if tdim 

  }
}

//#define BOOST_PP_LOCAL_MACRO(n)
template void MetrisRunner::statMesh0<MetricFieldAnalytical>();
template void MetrisRunner::statMesh0<MetricFieldFE        >();
//#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
//#include BOOST_PP_LOCAL_ITERATE()


}