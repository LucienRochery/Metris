//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../MetricField/MetricField.hxx"
#include "../Mesh/MeshBase.hxx"
#include "../Mesh/MeshMetric.hxx"

#include "../msh_anamet.hxx"
#include "../aux_utils.hxx"
#include "../low_eval.hxx"
#include "../linalg/utils.hxx"
#include "../linalg/explogmet.hxx"

#include <cmath>

namespace Metris{




MetricFieldAnalytical::MetricFieldAnalytical(MeshMetric<MetricFieldAnalytical> &msh) 
:MetricFieldFE(msh){
	ianamet = -1;
	anamet = NULL;
	scale = 1.0;
  #ifdef DEBUG_ARRAYS_FULL
  rfld.set_narray(&(msh.npoin));
  #endif
}


void MetricFieldAnalytical::setAnalyticalMetric(int ianamet_){
	int idim = msh.idim;
	METRIS_ASSERT(idim == 2 || idim == 3);
	
 	if(ianamet_ <= 0 || ianamet_ > MAX_ANAMET_DEFINED(idim ) )
 	  METRIS_THROW_MSG(WArgExcept(),"Invalid index: 1 - "<<MAX_ANAMET_DEFINED(idim )<<" accepted");
	
	this->ianamet = ianamet_;
  this->anamet  = (idim == 2 ? __ANAMET2D[this->ianamet-1] : __ANAMET3D[this->ianamet-1]);
}

void MetricFieldAnalytical::setAnalyticalMetric(anamet_proto anamet_ptr){
  METRIS_ASSERT(anamet_ptr != NULL);
  
  this->ianamet = -1;
  this->anamet  = anamet_ptr;
}
 
MetricFieldAnalytical &MetricFieldAnalytical::operator=(const MetricFieldAnalytical& inp){
	MetricFieldFE::operator=(inp);
	this->ianamet = inp.ianamet;
	this->anamet  = inp.anamet;
	this->scale   = inp.scale;
	return *this;
}


void MetricFieldAnalytical::normalize(double coeff){
  // Easy !
  scale = coeff;
}



// --- 
void MetricFieldAnalytical::getMetPhys(AsDeg metdeg,
                                       DifVar idiff, MetSpace tarspac, 
                                       int *ieleg, 
                                       const double*__restrict__ coop, 
                                       double*__restrict__ metl, 
                                       double*__restrict__ dmet, 
                                       int ithread) {
    CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
      getMetPhys0<gdim>(idiff,tarspac,coop,metl,dmet);
    }}CT_FOR1(gdim);
}




// --- 
void MetricFieldAnalytical::getMetBary(AsDeg asdmet,
                                       DifVar idiff, MetSpace tarspac, 
                                       const int*__restrict__ ent2pol, 
                                       int tdimn, 
                                       const double*__restrict__  bary,  
                                       double*__restrict__ metl, 
                                       double*__restrict__ dmet) {

  CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
    if(asdmet == AsDeg::P1){

      getMetBary0<gdim,1>(idiff,tarspac,ent2pol,
                          tdimn,bary,metl,dmet);

    }else{

      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
        getMetBary0<gdim,ideg>(idiff,tarspac,ent2pol,
                                    tdimn,bary,metl,dmet);
      }}CT_FOR1(ideg);

    }
  }}CT_FOR1(gdim);
}


// --- 
void MetricFieldAnalytical::
  getMetFullinfo(AsDeg asdmet,DifVar idiff, MetSpace tarspac, 
                 const int*__restrict__ ent2pol, int tdimn, 
                 const double*__restrict__  bary,  
                 const double*__restrict__  coop,  
                 double*__restrict__ metl,         
                 double*__restrict__ dmet) {
  getMetPhys(asdmet,idiff,tarspac,NULL,coop,metl,dmet);
}



template<int gdim, int ideg>
void MetricFieldAnalytical::
  getMetBary0(DifVar idiff, MetSpace tarspac, 
              const int*__restrict ent2pol, int tdimn, 
              const double*__restrict__ bary, 
                    double*__restrict__ metl, 
                    double*__restrict__ dmet) {
	METRIS_ASSERT(gdim == msh.idim);
	METRIS_ENFORCE_MSG(idiff != DifVar::Bary, "Implement Barycentric derivatives cf reverse of MetricFieldAnalytical.")
  
	double coop[gdim];

  if(tdimn == 1){
    eval1<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),DifVar::None,DifVar::None,bary,coop,NULL,NULL);
  }else if(tdimn == 2){
    eval2<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),DifVar::None,DifVar::None,bary,coop,NULL,NULL);
  }else{
    METRIS_ASSERT(gdim == 3); // Or we could just make tdimn a template parameter and avoid linker errors here by using constexpr if 
    eval3<3,ideg>(msh.coord,ent2pol,msh.getBasis(),DifVar::None,DifVar::None,bary,coop,NULL,NULL);
  }

	getMetPhys0<gdim>(idiff,tarspac,coop,metl,dmet);

}


template<int gdim>
void MetricFieldAnalytical::
  getMetPhys0(DifVar idiff, MetSpace tarspac, 
              const double*__restrict__  coop,  
              double*__restrict__ metl, 
              double*__restrict__ dmet){
	METRIS_ASSERT(gdim == msh.idim);
	int idifa = 0;
	if(idiff != DifVar::None) idifa = 1;

	anamet(NULL,coop,scale,idifa,metl,dmet);

	constexpr int nnmet = (gdim*(gdim+1))/2;
	
	if(tarspac != MetSpace::Exp){
    // This is undesireable but we rather prepare for it. 
    // Let's log that it happened.
    this->nspace_miss++;

    if(idiff != DifVar::None){
      SANS::SurrealS<gdim,double> metS[nnmet];
      getmet_dbl2SurS<gdim,gdim>(metl,dmet,metS);
      getspacmet_inp<gdim,SANS::SurrealS<gdim,double>>(metS, tarspac);
      getmet_SurS2dbl<gdim,gdim>(metS,metl,dmet);
    }else{
      #ifndef NDEBUG
      for(int ii = 0; ii < nnmet; ii++){
        if(std::isnan(metl[ii])){
          printf("NaN analytical metric at ");
          dblAr1(gdim,coop).print();
          METRIS_THROW(GeomExcept());
        }
      }
      #endif
      getspacmet_inp<gdim,double>(metl, tarspac);
    }
  }
}


#define BOOST_PP_LOCAL_MACRO(ideg)\
template void MetricFieldAnalytical::getMetBary0<2,ideg>\
     (DifVar idiff, MetSpace tarspac, const int*__restrict ent2poi, int tdimn,\
                            const double*__restrict__ bary, \
                                  double*__restrict__ metl, double*__restrict__ dmet) ;\
template void MetricFieldAnalytical::getMetBary0<3,ideg>\
     (DifVar idiff, MetSpace tarspac, const int*__restrict ent2poi, int tdimn,\
                            const double*__restrict__ bary, \
                                  double*__restrict__ metl, double*__restrict__ dmet) ;
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


template void MetricFieldAnalytical::getMetPhys0<2>
                            (DifVar idiff, MetSpace tarspac, 
                            const double*__restrict__ coop, 
                                  double*__restrict__ metl, double*__restrict__ dmet) ;
template void MetricFieldAnalytical::getMetPhys0<3>
                            (DifVar idiff, MetSpace tarspac, 
                            const double*__restrict__ coop, 
                                  double*__restrict__ metl, double*__restrict__ dmet) ;


} // End namespace
