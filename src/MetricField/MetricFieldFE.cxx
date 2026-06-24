//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../MetricField/MetricField.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../MetricField/msh_explogmet.hxx"
#include "../Mesh/MeshBase.hxx"
#include "../msh_lag2bez.hxx"
#include "../low_eval.hxx"
#include "../io_libmeshb.hxx"
#include "../msh_structs.hxx"
#include "../linalg/invmat.hxx"
#include "../linalg/explogmet.hxx"
#include "../linalg/utils.hxx"
#include "../linalg/symidx.hxx"

#include "../utils/CT_loop.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"
#include "../utils/aux_misc.hxx"

#include "../Localization/msh_localization.hxx"


namespace Metris{

MetricFieldFE::MetricFieldFE(MeshBase &msh_) : msh(msh_){
  ispace = MetSpace::Exp;
  ibasis = FEBasis::Undefined;
  nspace_miss = 0;
  rfld.set_n(msh.npoin);
}


void MetricFieldFE::setSpace(MetSpace ispacn, bool iforce){
  #ifndef NDEBUG
  if(msh.meshClass() == MeshClass::MeshBack && !iforce && ispacn != this->ispace){
    METRIS_THROW_MSG("setSpace called on back mesh! Should not be done.")
  }
  #endif
  METRIS_ASSERT(ispacn != MetSpace::Undefined);
  if(ispacn == MetSpace::Log) setLog();
  else                        setExp();
}

int MetricFieldFE::getnnmet()const{
  return (msh.idim*(msh.idim+1))/2;
}

void MetricFieldFE::allocate(){
  int mpoin = msh.mpoin;
  METRIS_ASSERT(mpoin > 0);
  int idim = msh.idim;
  METRIS_ASSERT(idim == 2 || idim == 3);
  int nnmet = getnnmet();
  this->rfld.allocate(mpoin,nnmet);
}

MetricFieldFE &MetricFieldFE::operator=(const MetricFieldFE& inp){
  METRIS_ASSERT(inp.msh.npoin == msh.npoin);
  METRIS_ASSERT(inp.msh.idim == msh.idim);
  METRIS_ASSERT(msh.idim == 2 || msh.idim == 3);
  this->ispace = inp.ispace;
  this->ibasis = inp.ibasis;
  //this->npoin  = inp.npoin;
  //this->nnmet  = inp.nnmet;
  allocate();

  inp.rfld.copyTo(rfld,msh.npoin);
  return *this;
}


void MetricFieldFE::setLog(){

  if(ibasis == FEBasis::Undefined){
    GETVDEPTH(msh.param);
    CPRINTF1("# WARNING: Metric field basis not defined, possible surface or line mesh: skip setLog()\n");
    return;
  }

  METRIS_ASSERT(msh.idim == 2 || msh.idim == 3);

  if(this->ispace == MetSpace::Log) return;

	if(msh.idim == 2){
		setLogMetMesh0<2>(msh,this->rfld);
	}else{
		setLogMetMesh0<3>(msh,this->rfld);
	}

  this->ispace = MetSpace::Log;
}

void MetricFieldFE::setExp(){

  if(ibasis == FEBasis::Undefined){
    GETVDEPTH(msh.param);
    CPRINTF1("# WARNING: Metric field basis not defined, possible surface or line mesh: skip setExp()\n");
    return;
  }

  METRIS_ASSERT(msh.idim == 2 || msh.idim == 3);

  if(this->ispace == MetSpace::Exp) return;

	if(msh.idim == 2){
		setExpMetMesh0<2>(msh,this->rfld);
	}else{
		setExpMetMesh0<3>(msh,this->rfld);
	}

  this->ispace = MetSpace::Exp;
}



void MetricFieldFE::setLagrange(){
  if(this->ibasis == FEBasis::Lagrange) return;

  CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){
    if(ideg == this->msh.curdeg){
      if(this->msh.idim == 2){
        setFieldLagrange<ideg,3>(this->msh,this->rfld);
      }else{
        setFieldLagrange<ideg,6>(this->msh,this->rfld);
      }
    }
  }CT_FOR1(ideg);

  this->ibasis = FEBasis::Lagrange;
  return;
}

void MetricFieldFE::setBezier(){
  if(this->ibasis == FEBasis::Bezier) return;

  CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){
    if(ideg == this->msh.curdeg){
      if(this->msh.idim == 2){
        setFieldBezier<ideg,3>(this->msh,this->rfld);
      }else{
        setFieldBezier<ideg,6>(this->msh,this->rfld);
      }
    }
  }CT_FOR1(ideg);

  this->ibasis = FEBasis::Bezier;
  return;
}




void MetricFieldFE::readMetricFile(std::string inpname){

  int metdim;
  int64_t libIdx = MetrisOpenMeshFile<GmfRead>(inpname.c_str(), &metdim);

  METRIS_ENFORCE_MSG(metdim == msh.idim, "Metric and mesh dimension must agree");

  GETVDEPTH(msh.param);

  int ilag = 1 - GmfStatKwd(libIdx, GmfBezierBasis);
  if(ilag == 1){
    CPRINTF1("Metric in Lagrange format.\n");
    this->ibasis = FEBasis::Lagrange;
  }else{
    CPRINTF1("Metric in Bézier format.\n");
    this->ibasis = FEBasis::Bezier;
  }

  this->ispace = MetSpace::Exp;
  int ncom = GmfStatKwd(libIdx, GmfComments);
  if(ncom) CPRINTF1("Metric file contains {} comments\n", ncom);

  char buff[257];
  GmfGotoKwd(libIdx, GmfComments);
  for(int icom = 0; icom < ncom; icom++){
    GmfGetLin(libIdx, GmfComments, buff);
    std::string strbuf(buff);
    if(strbuf.find("log") != std::string::npos){
      CPRINTF1("log-Metric supplied.\n");
      this->ispace = MetSpace::Log;
    }else{
      PRINTF("## UNKNOWN COMMENT LINE: {}\n", strbuf);
      METRIS_THROW_MSG(
        "Known bug: Comments works in binary file but not plaintext after transmesh.\n")
    }
  }


  int nsolf,ltyp[GmfMaxTyp],szfls;
  int npoif = GmfStatKwd(libIdx,GmfSolAtVertices, &nsolf, &szfls, ltyp);
  METRIS_ENFORCE_MSG(npoif == msh.npoin,"Metric file npoin = {} does not agree with mesh = {}", npoif, msh.npoin);

  METRIS_ENFORCE(nsolf == 1);
  METRIS_ENFORCE(ltyp[0] == GmfSymMat);

  int nnmet = getnnmet();
  METRIS_ASSERT_MSG(nnmet == 3 || nnmet == 6," nnmet = {}", nnmet);

  // Possible reallocation
  this->rfld.allocate(msh.mpoin,nnmet);
  this->rfld.set_n(msh.npoin);

  GmfGotoKwd(libIdx, GmfSolAtVertices);

  GmfGetBlock(libIdx, GmfSolAtVertices, 1, msh.npoin, 0, NULL, NULL,
    GmfDoubleVec, nnmet, &this->rfld[0][0], &this->rfld[msh.npoin-1][0]);

  GmfCloseMesh(libIdx);
}

void MetricFieldFE::writeMetricFile(std::string outname, bool iprefix){

  GETVDEPTH(msh.param);

  // For now, always force writing as exp metric.
  const MetSpace outspac = MetSpace::Exp;

  if(ibasis == FEBasis::Undefined){
    CPRINTF1("# WARNING: Metric field basis not defined, possible surface or line mesh: skip writeMetricFile = {}\n", outname);
    return;
  }

  dblAr2 buffer;
  if(outspac == this->ispace){
    buffer = this->rfld;
  }else{
    buffer.allocate(this->rfld.get_n(),this->rfld.get_stride());
    buffer.set_n(this->rfld.get_n());
    for(int ii = 0; ii < this->rfld.get_n(); ii++){
      for(int jj = 0; jj < this->rfld.get_stride(); jj++){
        buffer(ii,jj) = this->rfld(ii,jj);
      }
    }
    if(outspac == MetSpace::Log){
      if(msh.idim == 2){
        setLogMetMesh0<2>(msh, buffer);
      }else{
        setLogMetMesh0<3>(msh, buffer);
      }
    }else{
      if(msh.idim == 2){
        setExpMetMesh0<2>(msh, buffer);
      }else{
        setExpMetMesh0<3>(msh, buffer);
      }
    }
  }


  std::string metName = correctExtension_solb(outname);
  int idim = msh.idim;
  METRIS_ASSERT(idim == 2 || idim == 3);

  metName = iprefix ? msh.param->outmPrefix + metName : metName;

  int nnmet = getnnmet();

  int szfld = GmfSymMat;

  int64_t libIdx;

  libIdx = MetrisOpenMeshFile<GmfWrite>(metName, idim);
  if(this->ibasis == FEBasis::Bezier) GmfSetKwd( libIdx, GmfBezierBasis, 1);
  if(outspac == MetSpace::Log ){
    GmfSetKwd(libIdx,GmfComments,1);
    const char buf[8] = "log";
    GmfSetLin(libIdx,GmfComments,buf);
  }

  CPRINTF1("-- Write file file {} npoin {}\n",metName,msh.npoin);


  GmfSetKwd(libIdx, GmfSolAtVertices, msh.npoin, 1, &szfld);
  GmfSetBlock(libIdx, GmfSolAtVertices, 1, msh.npoin, 0, NULL, NULL,
              GmfDoubleVec, nnmet, buffer[0], buffer[msh.npoin-1]);


  GmfCloseMesh( libIdx );
}



// coeff is size multiplier
void MetricFieldFE::normalize(double coeff){
  METRIS_ASSERT(coeff > 0.0);
  double size = 1/coeff/coeff;
  int nnmet = getnnmet();
  if(ispace == MetSpace::Exp){
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      if(msh.isdeadpoint(ipoin)) continue;
      for(int ii = 0; ii < nnmet; ii++){
        rfld(ipoin,ii) *= size;
      }
    }
  }else{
    double lsize = log(size);
    // Identity commutes with everyone so exp(log(M) + c I)
    // = M * exp(c)I
    // hence replace c with log(c) and add identity to the metrics.
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      if(msh.isdeadpoint(ipoin)) continue;
      for(int ii = 0; ii < msh.idim; ii++){
        rfld[ipoin][sym2idx(ii,ii)] += lsize;
      }
    }
  }
}


void MetricFieldFE::correctMetric(){
  //const MetrisParameters *param = msh.param;

  double eigval[3], eigvec[9];
  for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
    INCVDEPTH(msh.param);
    if(msh.isdeadpoint(ipoin)) continue;
    if(msh.idim == 2){
      geteigsym<2,double>(rfld[ipoin],eigval,eigvec);
    }else{
      geteigsym<3,double>(rfld[ipoin],eigval,eigvec);
    }
    for(int ii = 0; ii < msh.idim; ii++){
      if(eigval[ii] < 1.0e-32)
        PRINTF("## CORRECTED SMALL OR NEGATIVE EIGVALS: {}\n",dblAr1(msh.idim,eigval));
      eigval[ii] = abs(eigval[ii]);
    }
    if(msh.idim == 2){
      eig2met<2,double>(eigval,eigvec,rfld[ipoin]);
    }else{
      eig2met<3,double>(eigval,eigvec,rfld[ipoin]);
    }
  }

}



// ----
void MetricFieldFE::getMetBary(AsDeg asdmet,
                               DifVar idiff,  MetSpace tarspac,
                               const int* __restrict__ ent2pol, int tdimn,
                               const double*__restrict__  bary,
                               double*__restrict__ metl,
                               double*__restrict__ dmet) {
  METRIS_ENFORCE_MSG(ibasis != FEBasis::Undefined, "Metric was not initialized.");

  // barycentric is defined but not physical. Probably a mistake either way
  if(msh.idim != tdimn && idiff != DifVar::None)
    METRIS_THROW_MSG( "Differing tdim = {} and gdim = {} are you sure about idiff?", tdimn, msh.idim)



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

  if(std::isnan(metl[0])){
    GETVDEPTH(msh.param);
    PRINTF("## DEBUG NAN METRIC IN GETMETBARY\n");
    PRINTF("bary = {}\n",dblAr1(msh.idim+1,bary));
    for(int inode = 0; inode < getnnode(msh.idim,msh.curdeg); inode++){
      PRINTF("elt node {} met = {}\n",ent2pol[inode],
             dblAr1((msh.idim*(msh.idim+1))/2, this->rfld[ent2pol[inode]]));
    }
  }

}




// ----
void MetricFieldFE::getMetFullinfo(AsDeg asdmet,
                                   DifVar idiff, MetSpace tarspac,
                                   const int*__restrict__ ent2pol,
                                   int tdimn,
                                   const double*__restrict__  bary,
                                   [[maybe_unused]] const double*__restrict__  coop,
                                   double*__restrict__ metl,
                                   double*__restrict__ dmet) {
  METRIS_ENFORCE_MSG(ibasis != FEBasis::Undefined,
                     "Metric was not initialized.");
  getMetBary(asdmet,idiff,tarspac,ent2pol,tdimn,bary,metl,dmet);
}






template<int gdim, int ideg>
void MetricFieldFE::getMetBary0( DifVar idiff,  MetSpace tarspac,
                                const int*__restrict__ ent2pol,
                                int tdimn,
                                const double*__restrict__ bary,
                                      double*__restrict__ metl,
                                      double*__restrict__ dmet) {
  constexpr int nnmet = (gdim*(gdim+1))/2;

  // METRIS_ASSERT_MSG(this->ispace == MetSpace::Log,
  //                   "Do not call getMetBary on an Exp met field");

  double *dmet0 = dmet;
  RoutineWorkMemory<double> tmp(this->rwrkmem);
  DifVar idife = DifVar::None;
  if(idiff == DifVar::Phys){
    dmet0 = tmp.allocate(nnmet*gdim);
    idife = DifVar::Bary;
  }

  // Get M(xi) and dM / dxi at bary
  if(tdimn == 1){
    eval1<nnmet,ideg>(rfld,ent2pol,this->ibasis,idife,DifVar::None,bary,metl,dmet0,NULL);
  }else if(tdimn == 2){
    eval2<nnmet,ideg>(rfld,ent2pol,this->ibasis,idife,DifVar::None,bary,metl,dmet0,NULL);
  }else{
    eval3<nnmet,ideg>(rfld,ent2pol,this->ibasis,idife,DifVar::None,bary,metl,dmet0,NULL);
  }

  if(tdimn > gdim) METRIS_THROW_MSG( "getMetBary0 topo dim > gdim");


  if(idiff == DifVar::None){
    // If not returning now, we'll do this with SurrealS expmet
    if(tarspac != this->ispace) getspacmet_inp<gdim,double>(metl,tarspac);
    return;
  }


  if(idiff == DifVar::Phys){
    double coop[gdim];
    double jmat3[gdim][gdim];
    // Get dF_K / dxi at bary

    if(tdimn == 1){
      eval1<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),DifVar::Bary,
                         DifVar::None,bary,coop,jmat3[0],NULL);
    }else if(tdimn == 2){
      if constexpr(gdim >= 2)
        eval2<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),DifVar::Bary,
                           DifVar::None,bary,coop,jmat3[0],NULL);
    }else{
      if constexpr(gdim >= 3)
        eval3<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),DifVar::Bary,
                           DifVar::None,bary,coop,jmat3[0],NULL);
    }

    // Get dM / dX at X (coop)
    // 1. Invert dF_K/dxi
    METRIS_ENFORCE(!invmat<gdim>(jmat3[0]));
    // 2. Get dM/dxi (dF_K/dxi)^{-1}
    int kk = 0;
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = 0; jj < nnmet; jj++){
        dmet[kk] = jmat3[ii][0]*dmet0[nnmet*0+jj];
        for(int ll = 1; ll < gdim; ll++){
          dmet[kk] += jmat3[ii][ll]*dmet0[nnmet*ll+jj];
        }
        kk++;
      }
    }
  }

  if(tarspac != this->ispace){
    // This is undesireable but we rather prepare for it.
    // Let's log that it happened.
    this->nspace_miss++;

    SANS::SurrealS<gdim,double> metS[nnmet];
    getmet_dbl2SurS<gdim,gdim>(metl,dmet,metS);
    getspacmet_inp<gdim,SANS::SurrealS<gdim,double>>(metS, tarspac);
    getmet_SurS2dbl<gdim,gdim>(metS,metl,dmet);
  }


}
#define BOOST_PP_LOCAL_MACRO(ideg)\
template void MetricFieldFE::getMetBary0<2, ideg>\
                            ( DifVar idiff,  MetSpace tarspac, \
                                const int*__restrict__ ent2pol,\
                                int tdimn,\
                                const double*__restrict__ bary,\
                                      double*__restrict__ metl,\
                                      double*__restrict__ dmet);\
template void MetricFieldFE::getMetBary0<3, ideg>\
                            ( DifVar idiff,  MetSpace tarspac, \
                                const int*__restrict__ ent2pol,\
                                int tdimn,\
                                const double*__restrict__ bary,\
                                      double*__restrict__ metl,\
                                      double*__restrict__ dmet);
#define BOOST_PP_LOCAL_LIMITS (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



} // End namespace
