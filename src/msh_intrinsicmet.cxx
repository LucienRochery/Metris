//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "msh_intrinsicmet.hxx"

////#include "msh_metric.hxx"
#include "low_geo.hxx"
#include "linalg/explogmet.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "io_libmeshb.hxx"

#include "../libs/lplib3.h"

namespace Metris{


//#include "test_subs.hxx"

// Metric field computed as Lagrange field
// TODO: ilagmet and ilogmet parameters: backmesh should be Bézier log
// We can simply to a Lag2Bez routine call on the metric field too
template<class MetricFieldType, int ideg>
void getMetMesh(const MetrisParameters &param, MeshMetric<MetricFieldType> &msh){

  double hmin = param.hmin;
  double hmax = param.hmax;
  double lbdmin = 1.0 / sqrt(hmax);
  double lbdmax = 1.0 / sqrt(hmin);

  int nproc = param.nproc;

	FEBasis ibas0 = msh.met.getBasis();

	msh.met.forceBasisFlag(FEBasis::Lagrange);
	msh.met.forceSpaceFlag(MetSpace::Log);


  int nthread = GetNumberOfCores();
  if(nthread <= 0){
    if(param.iverb >= 1)printf("## WARNING: LPlib function GetNumberOfCores() returned negative threads. Set to default %d.\n",METRIS_MAXTAGS);
    nthread = METRIS_MAXTAGS;
  }else{
    if(param.iverb >= 3) printf("-- LPlib found ncore = %d \n",nthread);
    if(nthread > METRIS_MAXTAGS){
      if(param.iverb >= 1)printf("## WARNING: must verify nthread <= METRIS_MAXTAGS = %d. Increase in metris_constants.hxx.\n",METRIS_MAXTAGS);
      nthread = METRIS_MAXTAGS;
    }
  }
  if(nproc > 0) nthread = MIN(nthread, nproc);
  if(param.iverb >= 1) printf("Running intrinsic metric with nproc = %d \n",nthread);


  int tdim = msh.get_tdim();
  METRIS_ENFORCE(tdim == 2 || tdim == 3);

  int nentt = msh.nentt(tdim); 
  const intAr2 &ent2poi = msh.ent2poi(tdim); 


  int64_t LibIdx = InitParallel(nthread);
  int LP_elt = NewType(LibIdx, nentt);
  int LP_poi = NewType(LibIdx, msh.npoin);
  float LP_stat[2];

  // We don't need to add all points as deps, only vertices. 
  // The others create subset dependencies. 
  BeginDependency(LibIdx, LP_elt, LP_poi);
  for(int ientt = 0; ientt < nentt; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;
    for(int ii = 0; ii < tdim + 1; ii++) AddDependency(LibIdx, ientt+1, ent2poi(ientt,ii)+1);
  }
  EndDependency(LibIdx, LP_stat);
  //printf("LP stat %f %f \n",LP_stat[0],LP_stat[1]);


	// Placeholder
  msh.rwork.allocate(msh.npoin);
  msh.rwork.set_n(msh.npoin);
	for(int ipoin = 0; ipoin < msh.npoin; ipoin++) msh.rwork[ipoin] = 1.0;
	

	msh.tag[0]++;
	int poitag = msh.tag[0];


	void (*metcomp2_LPlib)(int,int,int,MeshMetric<MetricFieldType>*,double,double)
	=[] (int ipoi0, int ipoi1, int ithrd, MeshMetric<MetricFieldType> *msh, double lbdmin, double lbdmax){
		int nnmet = (msh->idim*(msh->idim+1))/2;
		for(int ipoin = ipoi0 - 1; ipoin < ipoi1; ipoin++){
      if(msh->poi2ent(ipoin,0) < 0) continue;
			for(int jj = 0; jj < nnmet; jj++) msh->met(ipoin,jj) /= msh->rwork[ipoin];
		}
		// Control sizes here if provided (hmin hmax)
	};

	float acc; 

  CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
    CT_FOR0_INC(2,c_gdim,tdim_){if(tdim_ == tdim){
      acc = LaunchParallelMultiArg(LibIdx, LP_elt, LP_poi, (void*)getMetMesh0_lplib<MetricFieldType,gdim,tdim_,ideg>, 
                                   2, &msh, poitag);
   }}CT_FOR1(tdim_);
  }}CT_FOR1(gdim);

  if(param.iverb >= 3) printf("Intrinsic metric accel 1 = %f \n",acc);
 	acc = LaunchParallelMultiArg(LibIdx, LP_poi, 0, (void*)metcomp2_LPlib, 
  	                            3, &msh, lbdmin, lbdmax);
  if(param.iverb >= 3) printf("Intrinsic metric accel 2 = %f \n",acc);

	if(ibas0 == FEBasis::Bezier) msh.met.setBasis(FEBasis::Bezier);
	
	msh.met.setSpace(MetSpace::Exp);


	StopParallel(LibIdx);
}
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void getMetMesh< MetricFieldAnalytical , n>(const MetrisParameters &param, MeshMetric<MetricFieldAnalytical> &msh);\
template void getMetMesh< MetricFieldFE         , n>(const MetrisParameters &param, MeshMetric<MetricFieldFE        > &msh);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




template<class MetricFieldType, int gdim, int tdim, int ideg>
void getMetMesh0_lplib(int ient0, int ient1,int ithread, MeshMetric<MetricFieldType> *msh_, int poitag){
  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim == 2 || tdim == 3);
  static_assert(tdim <= gdim);


  MeshMetric<MetricFieldType> &msh = *msh_;
  constexpr int nnmet = (gdim*(gdim + 1))/2;
  double bary[tdim+1],metl[nnmet],meas0;
  intAr2 &ent2poi = tdim == 2 ? msh.fac2poi : msh.tet2poi;
  constexpr auto ordent = ORDELT(tdim);

  constexpr int npps = tdim == 2 ? facnpps[ideg] : tetnpps[ideg];

  bool iflat;

  for(int ientt = ient0 - 1; ientt < ient1; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;

    double norfac[3];
    if constexpr(tdim == 2 && gdim == 3) getnorfacP1(ent2poi[ientt], msh.coord, norfac);


    meas0 = getmeasentP1<gdim,tdim>(msh,ent2poi[ientt],norfac,&iflat);

    //METRIS_ENFORCE_MSG(!iflat,"Element "<<ientt<<" tdimn = "<<tdim<<" is flat"); // Not a program failure -> not an assert
    if(iflat){
      printf("Element %d gdim = %d tdim = %d is flat \n",ientt,gdim,tdim);
      printf("vertices: ");
      intAr1(tdim + 1, ent2poi[ientt]).print();
      writeMesh("debug_flat",msh);

      meas0 = getmeasentP1<gdim,tdim>(msh, ent2poi[ientt], norfac, &iflat);
      METRIS_THROW(GeomExcept());
    }

    for(int iver = 0; iver < npps; iver++){ 
      int ipoin = ent2poi(ientt,iver);
      METRIS_ASSERT(ipoin >= 0 && ipoin < msh.npoin);

      for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = ordent[ideg][iver][ii] / (double) ideg;

      if(msh.poi2tag(0,ipoin) < poitag){
        msh.poi2tag(0,ipoin) = poitag;
        msh.rwork[ipoin]    = meas0;
        METRIS_ENFORCE((!getintmetxi<gdim,tdim,ideg>(msh.coord,ent2poi[ientt],
                                               msh.getBasis(),bary,msh.met[ipoin])));
        getlogmet_inp<gdim>(msh.met[ipoin]);
        for(int jj = 0; jj < nnmet; jj++) msh.met(ipoin,jj) *= meas0;
      }else{
        msh.rwork[ipoin] += meas0;
        METRIS_ENFORCE((!getintmetxi<gdim,tdim,ideg>(msh.coord,ent2poi[ientt],
                                                 msh.getBasis(),bary,metl) ));
        getlogmet_inp<gdim>(metl);
        for(int jj = 0; jj < nnmet ; jj++) msh.met(ipoin,jj) += metl[jj]*meas0;
      }
    }
  }
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void getMetMesh0_lplib< MetricFieldAnalytical , 2, 2, n>(int ient0, int ient1,int ithread, MeshMetric<MetricFieldAnalytical> *msh_, int poitag);\
template void getMetMesh0_lplib< MetricFieldAnalytical , 3, 2, n>(int ient0, int ient1,int ithread, MeshMetric<MetricFieldAnalytical> *msh_, int poitag);\
template void getMetMesh0_lplib< MetricFieldAnalytical , 3, 3, n>(int ient0, int ient1,int ithread, MeshMetric<MetricFieldAnalytical> *msh_, int poitag);\
template void getMetMesh0_lplib< MetricFieldFE         , 2, 2, n>(int ient0, int ient1,int ithread, MeshMetric<MetricFieldFE        > *msh_, int poitag);\
template void getMetMesh0_lplib< MetricFieldFE         , 3, 2, n>(int ient0, int ient1,int ithread, MeshMetric<MetricFieldFE        > *msh_, int poitag);\
template void getMetMesh0_lplib< MetricFieldFE         , 3, 3, n>(int ient0, int ient1,int ithread, MeshMetric<MetricFieldFE        > *msh_, int poitag);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




} // End namespace




