//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "MetrisRunner.hxx"

#include "../Mesh/Mesh.hxx"

#include "../msh_lenedg.hxx"
#include "../aux_histogram.hxx"
#include "../quality/msh_metqua.hxx"
#include "../low_ccoef.hxx"
#include "../utils/mprintf.hxx"

#include "../SolutionField/SolutionField.hxx"
#include "../SolutionField/minInterpError.hxx"
#include "../SolutionField/interpError.hxx"

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
  GETVDEPTH(msh.param);

  if(!DOPRINTS1() && stat == NULL) return;

  const int tdim = msh.get_tdim();
  
  msh.cleanup();
  msh.met.setSpace(MetSpace::Exp);

  static const dblAr1 qbnds = {1.0e-8, 0.1};
  static const dblAr1 jbnds = {msh.param->jtol, 1};

  intAr2 ilned, ilned_bdry;
  dblAr1 rlned, rlned_bdry;
  getLengthEdges<MFT>(msh,tdim  ,-1,ilned     ,rlned     ,LenTyp::GeoSiz);
  getLengthEdges<MFT>(msh,tdim-1,-1,ilned_bdry,rlned_bdry,LenTyp::GeoSiz);

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


    stat->minlen_bdry = 1.0e30;
    stat->maxlen_bdry = -1.0e30;
    stat->avglen_bdry = 0.0;
    nunit = 0;
    for(int ii = 0; ii < ilned_bdry.get_n(); ii++){
      double len = rlned_bdry[ii];
      if(len >= 1.0/sqrt(2) && len <= sqrt(2)) nunit++;
      stat->avglen_bdry += len;
      stat->minlen_bdry = MIN(len, stat->minlen_bdry);
      stat->maxlen_bdry = MAX(len, stat->maxlen_bdry);
    }
    stat->avglen_bdry /= ilned_bdry.get_n();
    stat->pctunit_bdry = 100*((double)nunit) / ilned_bdry.get_n();
  }


  if(DOPRINTS1()){
    dblAr1 lenbds = {1.0/sqrt(2), sqrt(2)};
    print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length (geometric)");
    print_histogram(msh,rlned_bdry,IntrpTyp::Linear,lenbds,"l","Edge length (geo, bdry)");

    if(DOPRINTS3()){
      // This is very expensive to compute
      getLengthEdges<MFT>(msh,tdim,-1,ilned,rlned,LenTyp::Quad);
      print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length (quadrature)");
    }
  }


  if(param->anaSol && msh.idim == 2){
    SolutionFieldAnalytical sol(msh);
    sol.setAnalyticalSolution(param->ianasol);
    double errGlo = -1;
    //CT_FOR0_INC(2,3,idim){if(idim == msh.idim){
    CT_FOR0_INC(1,2,ideg){if(ideg == msh.curdeg){
    CT_FOR0_INC(1,2,pnorm){if(pnorm == param->intp_pnorm){
    CT_FOR0_INC(1,2,pdeg){if(pdeg == param->intp_pdeg){
      errGlo = interpErrGlo<2,ideg,pdeg,pnorm,true>(sol);
    }}CT_FOR1(pdeg);
    }}CT_FOR1(pnorm);
    }}CT_FOR1(ideg);
    //}}CT_FOR1(idim);
    CPRINTF1("-- Analytical solution interpolation error pnorm %d pdeg %d: %e\n",
             param->intp_pdeg, param->intp_pnorm, errGlo);
  } 


  double qmin, qmax, qavg;
  double qmin_bdry, qmax_bdry, qavg_bdry;
  bool iinva, iinva_bdry;
  dblAr1 lquae, lquae_bdry;
  getmetquamesh<MFT>(msh,tdim,AsDeg::P1,AsDeg::P1,
                     &iinva,&qmin,&qmax,&qavg,&lquae);

  double qua_surf_wt_normal = msh.param->qua_surf_wt_normal;
  if(tdim >= 3){
    msh.param->qua_surf_wt_normal = 0;
    getmetquamesh<MFT>(msh,tdim-1,AsDeg::P1,AsDeg::P1,
                       &iinva_bdry,&qmin_bdry,&qmax_bdry,&qavg_bdry,&lquae_bdry);
    msh.param->qua_surf_wt_normal = qua_surf_wt_normal;
  }
  if(stat != NULL){
    stat->minqua = qmin;
    stat->maxqua = qmax;
    stat->avgqua = qavg;

    stat->minqua_bdry = qmin_bdry;
    stat->maxqua_bdry = qmax_bdry;
    stat->avgqua_bdry = qavg_bdry;
  }
  if(DOPRINTS1()){
    CPRINTF1(" - Quality (as P1     ) min = %15.7e \n",qmin);
    CPRINTF1("                        max = %15.7e \n",qmax);
    CPRINTF1("                        avg = %15.7e \n",qavg);
    if(tdim >= 3){
      CPRINTF1(" - Quality (as P1 bdry) min = %15.7e \n",qmin_bdry);
      CPRINTF1("                        max = %15.7e \n",qmax_bdry);
      CPRINTF1("                        avg = %15.7e \n",qavg_bdry);
    }
  }
  if(DOPRINTS2()){
    print_histogram(msh,lquae,IntrpTyp::Geometric,qbnds,"q","Element quality (as P1)");
    if(tdim >= 3) print_histogram(msh,lquae_bdry,IntrpTyp::Geometric,qbnds,"q","Element quality (as P1 bdry)");
  }

  if(msh.curdeg > 1){
    getmetquamesh<MFT>(msh,tdim,AsDeg::Pk,AsDeg::Pk,
                       &iinva,&qmin,&qmax,&qavg,&lquae);
    msh.param->qua_surf_wt_normal = 0;
    if(tdim >= 3){
      getmetquamesh<MFT>(msh,tdim-1,AsDeg::Pk,AsDeg::Pk,
                         &iinva_bdry,&qmin_bdry,&qmax_bdry,&qavg_bdry,&lquae_bdry);
    }
    msh.param->qua_surf_wt_normal = qua_surf_wt_normal;
    if(stat != NULL){//overwrite
      stat->minqua = qmin;
      stat->maxqua = qmax;
      stat->avgqua = qavg;

      stat->minqua_bdry = qmin_bdry;
      stat->maxqua_bdry = qmax_bdry;
      stat->avgqua_bdry = qavg_bdry;
    }

    if(DOPRINTS2()){
      print_histogram(msh,lquae,IntrpTyp::Geometric,qbnds,"q","Element quality (as Pk)");
      print_histogram(msh,lquae_bdry,IntrpTyp::Geometric,qbnds,"q","Element quality (as Pk bdry)");
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
        
        print_histogram(msh,lquae,IntrpTyp::Geometric,jbnds,"J","Scaled Jacobian");
      }// if tdim 
    }// if DOPRINTS2
  }
}

//#define BOOST_PP_LOCAL_MACRO(n)
template void MetrisRunner::statMesh0<MetricFieldAnalytical>(MeshStat *stat);
template void MetrisRunner::statMesh0<MetricFieldFE        >(MeshStat *stat);
//#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
//#include BOOST_PP_LOCAL_ITERATE()


}