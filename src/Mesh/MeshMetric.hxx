//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MESH_METRIC__
#define __METRIS_MESH_METRIC__

#include "MeshBase.hxx"
#include "../MetricField/MetricField.hxx"

#include "../CT_loop.hxx"


namespace Metris{


template <class MFT> class MeshMetric;

template<int idim, int ideg,class MetricFieldType>
double getDomainVolume0(MeshMetric<MetricFieldType> &msh);


template<class MetricFieldType>
class MeshMetric : public MeshBase{
 

public: 
  MetricClass metricClass() const override { return met.metricClass(); }


  // // true if back mesh supplied from file, false otherwise
  // // this is used outside of IO in deg_elevate to know whether
  // // to interpolate from new or to localize in back
	// bool hasbak;
  // // Analytical metric 
  // int ianamet;
  // void (*anamet)(void* ctx, double *crd, int idif1, double *met, double *dmet);
  // //template<int ndimn>
  // //void (*anametS)(void *ctx, double *crd, SANS::SurrealS<ndimn,double> metS);
  // dblAr2 met;
  MetricFieldType met;

  MeshClass meshClass() const override { return MeshClass::MeshMetric; }

  MeshMetric() = delete;

	MeshMetric(int nipwk_=1, int niewk_=1, int nifwk_=1, int nitwk_=1, int nrpwk_=1) :
	 MeshBase(nipwk_, niewk_, nifwk_, nitwk_, nrpwk_), met(*this){}


	~MeshMetric(){}

  void set_npoin(int npoin, bool skipallocf = false) override; 
  void set_nentt(int tdimn, int nentt, bool skipallocf = false) override; 

	//void copyConstants(const MeshBase &msh, int MAX_DEG = METRIS_MAX_DEG){
	//	MeshBase::copyConstants(msh,MAX_DEG); 
	//}
	//void readConstants(int64_t libIdx, int usrMinDeg){
	//	MeshBase::readConstants(libIdx, usrMinDeg); 
	//}

	//void allocate(){allocate(getSysMem());}
	//void allocate(unsigned long long int mem);
  /*
	void setMetBasis(FEBasis ibasis){
		METRIS_ASSERT_MSG(ibasis != FEBasis::Undefined, "Undefined basis passed to setMetBasis");
		if(ibasis == met.getBasis()) return;

		int nwork = met.getnnmet()*npoin;
		double *rwork = getrwork(nwork);
		if(rwork == NULL) METRIS_THROW_MSG(DMemExcept(), "Unable to recover rwork array");

		met.setBasis(ibasis,nwork,rwork);
	}
  */

	
	double getDomainVolume(){
    double ret = -1;
    CT_FOR0_INC(2,3,gdim){if(gdim == this->idim){
      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){
        if(ideg == this->curdeg) ret = getDomainVolume0<gdim,ideg>(*this);
      }CT_FOR1(ideg);
    }}CT_FOR1(gdim);
    return ret;
	}


};


} // End namespace

#endif
