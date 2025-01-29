//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "low_gaps.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../Mesh/Mesh.hxx"
#include "../SANS/Surreal/SurrealS.h"
#include "../aux_exceptions.hxx"
#include "../ho_constants.hxx"
#include "../aux_utils.hxx"

#include "../linalg/eigen.hxx"
#include "../linalg/det.hxx"
#include "../linalg/utils.hxx"
#include "../linalg/matprods.hxx"

#include <boost/preprocessor/iteration/local.hpp>



namespace Metris{


/*
Only for degree 2 for now. 
- ientt element hosting edge
- tdim topological dimension of element
- dmet is a globally sized matrix:
  d2met(ipoin,j) = j-th component of the physical derivatives of the metric at ipoin
  The matrix can be size 0 if the metric is known analytically.

The offset does not depend on the element the edge is seen from: 

Delta_k = -1/2\sum_{st}  l_s l_t T_{skt}^M 
Enforce symmetry by doing 
Delta_k = -1/4 sum_st l_s l_t (T_{skt} + T_{tks})

CORRECTION: not l, but M^{+1/2}l !  obviously since originally N, and N = O(1) 

This way, no need to set constraints on metric derivatives tensor. 
*/

template<class MFT, int gdim, int ideg>
void getBezOffsetsEdge(Mesh<MFT> &msh, 
  int tdim, const int* ent2pol, int iedgl, double* offsets){
  METRIS_ASSERT(ideg == msh.curdeg || ideg == 1);
  METRIS_ASSERT(tdim == 2 || tdim == 3);

  double dedg0[gdim],dedg1[gdim];



  int nedgl = (tdim*(tdim+1))/2;
  const intAr2 lnoed(nedgl,2,tdim == 2 ? lnoed2[0] : lnoed3[0]);

  constexpr int nnmet = (gdim*(gdim+1))/2;
  double met[nnmet], dmet[gdim*nnmet];
  int ipoi1 = ent2pol[lnoed(iedgl,0)];
  int ipoi2 = ent2pol[lnoed(iedgl,1)];
  for(int ii = 0; ii < gdim; ii++){ 
    dedg0[ii] = msh.coord(ipoi2,ii) - msh.coord(ipoi1,ii);
  }

  #ifndef NDEBUG
  if(msh.param->iverb >= 4){
    printf("-- getBezOffsetsEdge debug points:\n");
    printf("ipoi1 = %d : ",ipoi1);
    dblAr1(gdim,msh.coord[ipoi1]).print();
    printf("ipoi2 = %d : ",ipoi2);
    dblAr1(gdim,msh.coord[ipoi2]).print();
  }
  #endif

  // If metric field analytical, we can get derivative just from edge
  // Otherwise, need the same dim element. 
  double bary[gdim+1] = {0};
  if constexpr (std::is_same<MFT,MetricFieldAnalytical>::value){
    int edg2pol[edgnpps[ideg]];
    edg2pol[0] = ipoi1; 
    edg2pol[1] = ipoi2; 
    int idx0 = tdim + 1 + iedgl*(ideg-1);
    for(int ii = 0; ii < ideg-1; ii++)  edg2pol[2+ii] = ent2pol[idx0+ii];

    //msh.met.template getMetBary<AsDeg::Pk,AsDeg::Pk>(DifVar::Phys,MetSpace::Exp,ent2pol,tdim,bary,met,dmet);

    double bary[2] = {0.5,0.5};
    msh.met.getMetBary(AsDeg::P1,DifVar::Phys,
                       MetSpace::Exp,
                       edg2pol,1,bary,met,dmet);
  }else{
    METRIS_ASSERT(tdim == gdim);
    bary[lnoed(iedgl,0)] = 0.5;
    bary[lnoed(iedgl,1)] = 0.5;
    msh.met.getMetBary(AsDeg::P1,DifVar::Phys,MetSpace::Exp,
                       ent2pol,tdim,bary,met,dmet);
  }
  #ifndef NDEBUG
    if(std::isnan(met[0]) || std::isnan(met[1]) || std::isnan(met[2])){
      printf("## DEBUG tdim %d iedgl %d ent2pol: ",tdim,iedgl);
      intAr1(tdim+1,ent2pol).print();
      printf("NaN met = ");
      dblAr1(nnmet,met).print();
      printf("entity mets dump:\n");
      for(int ii = 0; ii < tdim + 1; ii++){
        int ipdbg = ent2pol[ii];
        printf("%d : ",ipdbg);
        dblAr1(nnmet,msh.met[ipdbg]).print();
      }
      printf("bary = ");
      dblAr1(tdim+1,bary).print();
    }
  #endif

  //double len = getlenedg<gdim>(dedg0,met);
  //// if sqrt(u^T M u) = l, then u / l is unit in M: sqrt(u^T/l  M u/l) = 1
  //double scale = len;
  //// Hence we rescale the offsets by scale at the end. 
  //// We could also have rescaled the metric. 
  //// There is no perfect approach?

  // Compute M^{1/2} (no diff)
  double met_p12[nnmet]; 
  double eigval[gdim], eigvec[gdim*gdim];
  geteigsym<gdim, double>(met,eigval,eigvec);
  for(int ii = 0; ii < gdim; ii++){
    eigval[ii] = sqrt(eigval[ii]); 
  }
  eig2met<gdim, double>(eigval, eigvec, met_p12);
  symXvec<gdim>(met_p12, dedg0, dedg1);


  #ifndef NDEBUG
  if(msh.param->iverb >= 4){
    //printf("scale = %f \n",scale);
    printf("dedg0 : ");
    dblAr1(gdim,dedg0).print();
    printf("dedg1 : ");
    dblAr1(gdim,dedg1).print();
    printf("-- getBezOffsetsEdge debug print met:\n");
    dblAr1(nnmet,met).print();
    printf("-- getBezOffsetsEdge debug print d1met:\n");
    dblAr1(nnmet,&dmet[0]).print();
    printf("-- getBezOffsetsEdge debug print d2met:\n");
    dblAr1(nnmet,&dmet[nnmet]).print();
    printf("-- getBezOffsetsEdge debug print met^{1/2}:\n");
    dblAr1(nnmet,met_p12).print();
  }
  #endif

  //#ifndef NDEBUG
  //  bool havenan = false;
  //  for(int ii = 0; ii < nnmet; ii++){
  //    if(isnan(met[ii])) havenan = true;
  //    for(int jj = 0; jj < gdim; jj++){
  //      if(isnan(dmet[ii*gdim+jj])) havenan = true;
  //    }
  //  }

  //  if(havenan){
  //    printf("## NAN METRICS IN getBezOffsetsEdge\n");
  //    printf("gdim = %d ideg = %d tdim = %d iedgl = %d edg2pol:",gdim,ideg,tdim,iedgl);
  //    intAr1(edgnpps[ideg],edg2pol).print();
  //    printf("met = ");
  //    dblAr1(nnmet,met).print();
  //    printf("dmet = ");
  //    dblAr1(gdim*nnmet,dmet).print();
  //  }
  //#endif

  typedef SANS::SurrealS<gdim,double> doubleS; 
  doubleS metS[nnmet];
  getmet_dbl2SurS<gdim,gdim>(met,dmet,metS);

  if(msh.idbg[0] > 0){
    printf("Metric: ");
    for(int ii = 0; ii < nnmet; ii++) printf(" %f ",metS[ii].value());
    printf("\n");
    for(int jj = 0; jj < gdim; jj++){
      printf("d%d: ",jj);
      for(int ii = 0; ii < nnmet; ii++) printf(" %f ",metS[ii].deriv(jj));
      printf("\n");
    }
    printf("\n");
  }


  // Get M^{-1/2} into met12_m12
  doubleS met_m12[nnmet]; 
  doubleS eigvalS[gdim], eigvecS[gdim*gdim];
  geteigsym<gdim, doubleS>(metS,eigvalS,eigvecS);
  for(int ii = 0; ii < gdim; ii++){
    eigvalS[ii] = 1.0 / sqrt(eigvalS[ii]); 
  }
  eig2met<gdim, doubleS>(eigvalS, eigvecS, met_m12);

  #ifndef NDEBUG
  if(msh.param->iverb >= 4){
    printf("-- getBezOffsetsEdge debug print met^{-1/2}:\n");
    printf("M^{-1/2}: ");
    for(int ii = 0; ii < nnmet; ii++) printf(" %f ",met_m12[ii].value());
    printf("\n");
    for(int jj = 0; jj < gdim; jj++){
      printf("d%d: ",jj);
      for(int ii = 0; ii < nnmet; ii++) printf(" %f ",met_m12[ii].deriv(jj));
      printf("\n");
    }
    printf("\n");
  }
  #endif


  // Tensor T_ijk = sum_s g_si d_s g_kj 
  // Compute C_ijk = T_ijk + T_kji has all the right symmetries. 
  // with g = M^{-1/2}
  double tens[gdim*nnmet];
  for(int ii = 0; ii < gdim; ii++){
    for(int jj = 0; jj < gdim; jj++){
      for(int kk = jj; kk < gdim; kk++){
        tens[ii*nnmet + sym2idx(jj,kk)] = 0;
        for(int tt = 0; tt < gdim; tt++){
          //tens[ii*nnmet + sym2idx(jj,kk)] += 
          //  (met_m12[sym2idx(ii,tt)].value()*met_m12[sym2idx(jj,kk)].deriv(tt)
          //+  met_m12[sym2idx(kk,tt)].value()*met_m12[sym2idx(jj,ii)].deriv(tt))/2;
          tens[ii*nnmet + sym2idx(jj,kk)] += 
            met_m12[sym2idx(ii,tt)].value()*met_m12[sym2idx(jj,kk)].deriv(tt);
        }
      }
    }
  }
  #ifndef NDEBUG
  if(msh.param->iverb >= 4){
    printf("-- getBezOffsetsEdge debug print tensor T^M:\n");
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = 0; jj < gdim; jj++){
        for(int kk = jj; kk < gdim; kk++){
          printf("%d%d%d = %f\n",ii,jj,kk,tens[ii*nnmet + sym2idx(jj,kk)]);
        }
      }
    }
  }
  #endif
  


  for(int kk = 0; kk < gdim; kk++){
    offsets[kk] = 0;
    for(int tt = 0; tt < gdim; tt++){
      for(int ss = 0; ss < gdim; ss++){
        offsets[kk] += -dedg1[ss]*dedg1[tt]*tens[ss*nnmet + sym2idx(kk,tt)]/2;//*scale;
      }
    }
  }


  return;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template void getBezOffsetsEdge<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh, \
    int tdim, const int* ent2poi, int iedgl, double* offsets);\
template void getBezOffsetsEdge<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh, \
    int tdim, const int* ent2poi, int iedgl, double* offsets);\
template void getBezOffsetsEdge<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh, \
    int tdim, const int* ent2poi, int iedgl, double* offsets);\
template void getBezOffsetsEdge<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh, \
    int tdim, const int* ent2poi, int iedgl, double* offsets);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()
// (1, METRIS_MAX_DEG)

// srmet = M^{1/2}
// stJ0tR holds scale * J0^T * R^T
// In fact we don't really need the scale explicitely. 
template<int gdim>
void scalrotJ0(const MeshBase &msh, int ielem  ,
                const double* __restrict__ srmet  ,
                double* __restrict__ scale  ,
                double* __restrict__ stJ0tR  ){

  static_assert(gdim == 2 || gdim == 3);


  const intAr2& ent2poi = gdim == 2 ? msh.fac2poi : msh.tet2poi;

  // Used in both
  double jmat[gdim][gdim];
  for(int j = 0; j < gdim; j++){
    for(int i = 0; i < gdim; i++){
      jmat[j][i] = msh.coord[ent2poi[ielem][j+1]][i]
                 - msh.coord[ent2poi[ielem][  0]][i];
    }
    //jmat[0][i] = msh.coord[msh.tet2poi(ielem,1)][i]
    //           - msh.coord[msh.tet2poi(ielem,0)][i];
    //jmat[1][i] = msh.coord[msh.tet2poi(ielem,2)][i]
    //           - msh.coord[msh.tet2poi(ielem,0)][i];
    //jmat[2][i] = msh.coord[msh.tet2poi(ielem,3)][i]
    //           - msh.coord[msh.tet2poi(ielem,0)][i];
  }
  double detJ1 = detmat<gdim>(jmat[0]);
  double detJ0 = gdim == 2 ? sqrt(3)/2 : 1 / sqrt(2);
  printf("Not sqrt(3)??\n");
  wait();
  

  /*Get orientation*/

  // J^1^T M^{1/2} = scaling*J_0^T R^T
  matXsym<gdim>(jmat[0],srmet,stJ0tR);


  stJ0tR[2*0+0] = jmat[0][0]*srmet[0] + jmat[0][1]*srmet[1];
  printf("Recomputed %15.7e\n",stJ0tR[2*0+0]);



  /*
  Scale. If discrete metric, interpolate directly from front vertices as log-met, then expmet
  Otherwise, get 
  */
  //double bary[gdim+1] = {1.0/(gdim+1),1.0/(gdim+1),1.0/(gdim+1)};
  double detM12 = detsym<gdim>(srmet);
  if(detM12 < 1.0e-16) METRIS_THROW_MSG(GeomExcept(),"SINGULAR METRIC " << detM12
  <<" log-metric "<<srmet[0]<<" "<<srmet[1]<<" "<<srmet[2]<<" "
                  <<srmet[3]<<" "<<srmet[4]<<" "<<srmet[5]<<" ")

  *scale = pow(detJ1 * detM12 / detJ0,1.0/gdim);
}


template void scalrotJ0<2>(const MeshBase &msh, int ielem  ,
                           const double* __restrict__ srmet  ,
                           double* __restrict__ scale  ,
                           double* __restrict__ rotJ0  );
template void scalrotJ0<3>(const MeshBase &msh, int ielem  ,
                           const double* __restrict__ srmet  ,
                           double* __restrict__ scale  ,
                           double* __restrict__ rotJ0  );



} // End namespace
