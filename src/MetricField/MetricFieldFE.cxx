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

#include <array>
#include <utility>


namespace Metris{

MetricFieldFE::MetricFieldFE(MeshBase &msh_) : msh(msh_){
  ispace = MetSpace::Exp;
  ibasis = FEBasis::Undefined;
  nspace_miss = 0;
  rfld.set_n(msh.npoin);
}

MetricFieldFE::MetricFieldFE(const MetricFieldFE &inp) : msh(inp.msh), nspace_miss(0){
  inp.rfld.copyTo(rfld);
  inp.elementToCoefficient.copyTo(elementToCoefficient);
  inp.localCoefficients.copyTo(localCoefficients);
  layout = inp.layout;
  metricDegree = inp.metricDegree;
  ispace = inp.ispace;
  ibasis = inp.ibasis;
}

void MetricFieldFE::setSpace(MetSpace ispacn, bool iforce){
  #ifndef NDEBUG
  if(msh.meshClass() == MeshClass::MeshBack && !iforce && ispacn != this->ispace){
    METRIS_THROW_MSG("setSpace called on back mesh! Should not be done.")
  }
  #endif
  if(layout == MetricFieldLayout::Discontinuous){
    METRIS_ENFORCE(ispacn == MetSpace::Log || ispacn == MetSpace::Exp);
    if(ispacn == ispace) return;
    // Commit only after every coefficient has been converted successfully.
    dblAr2 converted;
    localCoefficients.copyTo(converted);
    for(int icoeff = 0; icoeff < converted.get_n(); icoeff++){
      if(ispacn == MetSpace::Log) getlogmet_inp<2,double>(converted[icoeff]);
      else getexpmet_frechet<2>(converted[icoeff],nullptr,0,converted[icoeff],nullptr);
      for(int k = 0; k < 3; k++) METRIS_ENFORCE(std::isfinite(converted(icoeff,k)));
    }
    localCoefficients = std::move(converted);
    ispace = ispacn;
    return;
  }
  METRIS_ASSERT(ispacn != MetSpace::Undefined);
  if(ispacn == MetSpace::Log) setLog();
  else                        setExp();
}

int MetricFieldFE::getnnmet()const{
  return (msh.idim*(msh.idim+1))/2;
}

void MetricFieldFE::allocate(){
  requireContinuous("Point-indexed metric allocation");
  int mpoin = msh.mpoin;
  METRIS_ASSERT(mpoin > 0);
  int idim = msh.idim;
  METRIS_ASSERT(idim == 2 || idim == 3);
  int nnmet = getnnmet();
  this->rfld.allocate(mpoin,nnmet);
}

MetricFieldFE &MetricFieldFE::operator=(const MetricFieldFE& inp){
  if(this == &inp) return *this;
  if(inp.layout == MetricFieldLayout::Discontinuous){
    METRIS_ENFORCE_MSG(metricClass() == MetricClass::MetricFieldFE,
                      "An analytical metric cannot use discontinuous FE storage");
    METRIS_ENFORCE_MSG(msh.idim == inp.msh.idim && msh.nface == inp.msh.nface
                      && msh.curdeg == inp.msh.curdeg && msh.npoin == inp.msh.npoin,
                      "DG metric assignment requires matching mesh layouts");
    for(int ientt = 0; ientt < msh.nface; ientt++){
      const bool dead = msh.fac2poi(ientt,0) < 0;
      METRIS_ENFORCE_MSG(dead == (inp.msh.fac2poi(ientt,0) < 0),
                        "DG metric assignment requires matching live entities");
      if(dead) continue;
      for(int node = 0; node < getnnode(2,msh.curdeg); node++){
        METRIS_ENFORCE_MSG(msh.fac2poi(ientt,node) == inp.msh.fac2poi(ientt,node),
                          "DG metric assignment requires matching entity ordering");
      }
    }
    inp.elementToCoefficient.copyTo(elementToCoefficient);
    inp.localCoefficients.copyTo(localCoefficients);
    layout = inp.layout;
    metricDegree = inp.metricDegree;
    ispace = inp.ispace;
    ibasis = inp.ibasis;
    return *this;
  }
  METRIS_ASSERT(inp.msh.npoin == msh.npoin);
  METRIS_ASSERT(inp.msh.idim == msh.idim);
  METRIS_ASSERT(msh.idim == 2 || msh.idim == 3);
  this->ispace = inp.ispace;
  this->ibasis = inp.ibasis;
  layout = MetricFieldLayout::Continuous;
  metricDegree = 0;
  elementToCoefficient.free();
  localCoefficients.free();
  //this->npoin  = inp.npoin;
  //this->nnmet  = inp.nnmet;
  allocate();

  inp.rfld.copyTo(rfld,msh.npoin);
  return *this;
}

int MetricFieldFE::getDegree() const{
  return layout == MetricFieldLayout::Continuous ? msh.curdeg : metricDegree;
}

int MetricFieldFE::getCoefficientCount() const{
  return layout == MetricFieldLayout::Continuous ? msh.npoin : localCoefficients.get_n();
}

void MetricFieldFE::initializeDiscontinuous(int degree){
  requireContinuous("Discontinuous metric initialization");
  METRIS_ENFORCE_MSG(metricClass() == MetricClass::MetricFieldFE,
                    "Discontinuous storage requires an FE metric");
  METRIS_ENFORCE_MSG(msh.idim == 2 && msh.get_tdim() == 2 && msh.curdeg == 2,
                    "Discontinuous metric storage requires planar Q2 triangles");
  METRIS_ENFORCE_MSG(degree == 1, "Discontinuous metric degree must be one");
  int count = 0;
  elementToCoefficient.allocate(msh.nface,3);
  elementToCoefficient.set_n(msh.nface);
  elementToCoefficient.fill(-1);
  for(int ientt = 0; ientt < msh.nface; ientt++){
    if(msh.fac2poi(ientt,0) < 0) continue;
    for(int node = 0; node < 3; node++) elementToCoefficient(ientt,node) = count++;
  }
  METRIS_ENFORCE_MSG(count > 0, "Discontinuous metric requires live triangles");
  localCoefficients.allocate(count,getnnmet());
  localCoefficients.set_n(count);
  localCoefficients.fill(0.);
  metricDegree = degree;
  ibasis = FEBasis::Lagrange;
  ispace = MetSpace::Log;
  layout = MetricFieldLayout::Discontinuous;
}

int MetricFieldFE::getElementCoefficientIndex(int ientt, int localNode) const{
  const int dim = msh.get_tdim();
  METRIS_ENFORCE(dim == 2 || dim == 3);
  METRIS_ENFORCE(ientt >= 0 && ientt < msh.nentt(dim));
  METRIS_ENFORCE(msh.ent2poi(dim)(ientt,0) >= 0);
  METRIS_ENFORCE(localNode >= 0 && localNode < getnnode(dim,getDegree()));
  if(layout == MetricFieldLayout::Continuous) return msh.ent2poi(dim)(ientt,localNode);
  METRIS_ENFORCE_MSG(msh.idim == 2 && dim == 2 && msh.curdeg == 2
                    && msh.nface == elementToCoefficient.get_n(),
                    "Mesh layout does not match DG metric storage");
  const int coefficient = elementToCoefficient(ientt,localNode);
  METRIS_ENFORCE(coefficient >= 0 && coefficient < localCoefficients.get_n());
  return coefficient;
}

const double *MetricFieldFE::elementCoefficient(int ientt, int localNode) const{
  const int index = getElementCoefficientIndex(ientt,localNode);
  return layout == MetricFieldLayout::Continuous ? rfld[index] : localCoefficients[index];
}

double *MetricFieldFE::elementCoefficient(int ientt, int localNode){
  const int index = getElementCoefficientIndex(ientt,localNode);
  return layout == MetricFieldLayout::Continuous ? rfld[index] : localCoefficients[index];
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
  requireContinuous("Vertex metric file input");

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
  requireContinuous("Vertex metric file output");

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
  if(layout == MetricFieldLayout::Discontinuous){
    METRIS_ENFORCE(std::isfinite(coeff) && coeff > 0.);
    const double shift = -2.*std::log(coeff);
    const double factor = std::exp(shift);
    METRIS_ENFORCE(ispace == MetSpace::Log || ispace == MetSpace::Exp);
    if(ispace == MetSpace::Exp) METRIS_ENFORCE(std::isfinite(factor) && factor > 0.);
    dblAr2 normalized;
    localCoefficients.copyTo(normalized);
    for(int icoeff = 0; icoeff < normalized.get_n(); icoeff++){
      if(ispace == MetSpace::Log){
        for(int k = 0; k < 2; k++) normalized(icoeff,sym2idx(k,k)) += shift;
      }else{
        for(int k = 0; k < 3; k++) normalized(icoeff,k) *= factor;
      }
      for(int k = 0; k < 3; k++) METRIS_ENFORCE(std::isfinite(normalized(icoeff,k)));
    }
    localCoefficients = std::move(normalized);
    return;
  }
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
  requireContinuous("Metric correction");
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



template<int gdim, int metdeg, int geomdeg>
void MetricFieldFE::getLogMetOnEntity0(int ientt, DifVar idiff,
                                      const double *bary, double *logm,
                                      double *dlogm) const{
  constexpr int nnmet = gdim*(gdim+1)/2;
  constexpr int nnode = getnnode(gdim,metdeg);
  std::array<double,nnode*nnmet> coefficients;
  std::array<int,nnode> nodes;
  for(int node = 0; node < nnode; node++){
    nodes[node] = node;
    const double *value = elementCoefficient(ientt,node);
    for(int k = 0; k < nnmet; k++) coefficients[node*nnmet+k] = value[k];
    if(ispace == MetSpace::Exp)
      getlogmet_inp<gdim,double>(&coefficients[node*nnmet]);
  }
  const dblAr2 field(nnode,nnmet,coefficients.data());
  double referenceDerivatives[gdim*nnmet];
  const DifVar referenceDiff = idiff == DifVar::None ? DifVar::None : DifVar::Bary;
  if constexpr(gdim == 2)
    eval2<nnmet,metdeg>(field,nodes.data(),ibasis,referenceDiff,DifVar::None,
                       bary,logm,referenceDerivatives,nullptr);
  else
    eval3<nnmet,metdeg>(field,nodes.data(),ibasis,referenceDiff,DifVar::None,
                       bary,logm,referenceDerivatives,nullptr);
  if(idiff == DifVar::None) return;
  if(idiff == DifVar::Bary){
    for(int k = 0; k < gdim*nnmet; k++) dlogm[k] = referenceDerivatives[k];
    return;
  }

  // eval stores J^T: row r contains derivatives with respect to xi_r.
  // Its inverse therefore maps reference gradients to physical gradients.
  double point[gdim], inverseTranspose[gdim*gdim];
  const int *geometryNodes = msh.ent2poi(gdim)[ientt];
  if constexpr(gdim == 2)
    eval2<gdim,geomdeg>(msh.coord,geometryNodes,msh.getBasis(),DifVar::Bary,
                      DifVar::None,bary,point,inverseTranspose,nullptr);
  else
    eval3<gdim,geomdeg>(msh.coord,geometryNodes,msh.getBasis(),DifVar::Bary,
                      DifVar::None,bary,point,inverseTranspose,nullptr);
  METRIS_ENFORCE_MSG(!invmat<gdim>(inverseTranspose),
                    "Singular geometry Jacobian in entity metric evaluation");
  for(int r = 0; r < gdim; r++) for(int k = 0; k < nnmet; k++){
    dlogm[r*nnmet+k] = 0.;
    for(int s = 0; s < gdim; s++)
      dlogm[r*nnmet+k] += inverseTranspose[r*gdim+s]*referenceDerivatives[s*nnmet+k];
  }
}

void MetricFieldFE::getLogMetOnEntity(int ientt, DifVar idiff,
                                     const double *bary, double *logm,
                                     double *dlogm) const{
  METRIS_ENFORCE_MSG(metricClass() == MetricClass::MetricFieldFE,
                    "Entity-local log-metric evaluation requires an FE metric");
  METRIS_ENFORCE(msh.idim == msh.get_tdim());
  METRIS_ENFORCE(msh.idim == 2 || msh.idim == 3);
  METRIS_ENFORCE(msh.curdeg >= 1 && msh.curdeg <= METRIS_MAX_DEG);
  METRIS_ENFORCE(ibasis == FEBasis::Lagrange || ibasis == FEBasis::Bezier);
  METRIS_ENFORCE(ispace == MetSpace::Log || ispace == MetSpace::Exp);
  METRIS_ENFORCE(idiff == DifVar::None || idiff == DifVar::Bary || idiff == DifVar::Phys);
  METRIS_ENFORCE(bary != nullptr && logm != nullptr);
  METRIS_ENFORCE(idiff == DifVar::None || dlogm != nullptr);
  getElementCoefficientIndex(ientt,0);
  if(layout == MetricFieldLayout::Discontinuous){
    METRIS_ENFORCE(metricDegree == 1);
    getLogMetOnEntity0<2,1,2>(ientt,idiff,bary,logm,dlogm);
    return;
  }
  CT_FOR0_INC(1,METRIS_MAX_DEG,degree){if(msh.curdeg == degree){
    if(msh.idim == 2)
      getLogMetOnEntity0<2,degree,degree>(ientt,idiff,bary,logm,dlogm);
    else
      getLogMetOnEntity0<3,degree,degree>(ientt,idiff,bary,logm,dlogm);
  }}CT_FOR1(degree);
}

void MetricFieldFE::getMetOnEntity(int ientt, DifVar idiff, const double *bary,
                                 double *met, double *dmet) const{
  METRIS_ENFORCE(met != nullptr);
  METRIS_ENFORCE(idiff == DifVar::None || dmet != nullptr);
  double logm[6], dlogm[18];
  getLogMetOnEntity(ientt,idiff,bary,logm,dlogm);
  const int nderiv = idiff == DifVar::None ? 0 : msh.idim;
  if(msh.idim == 2) getexpmet_frechet<2>(logm,dlogm,nderiv,met,dmet);
  else getexpmet_frechet<3>(logm,dlogm,nderiv,met,dmet);
}

// ----
void MetricFieldFE::getMetBary(AsDeg asdmet,
                               DifVar idiff,  MetSpace tarspac,
                               const int* __restrict__ ent2pol, int tdimn,
                               const double*__restrict__  bary,
                               double*__restrict__ metl,
                               double*__restrict__ dmet) {
  requireContinuous("Geometry-node-indexed metric evaluation");
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
  constexpr int max_nnode = getnnode(gdim,ideg);

  METRIS_ENFORCE(tarspac == MetSpace::Log || tarspac == MetSpace::Exp);
  METRIS_ENFORCE(this->ispace == MetSpace::Log
                 || this->ispace == MetSpace::Exp);

  double *dmet0 = dmet;
  RoutineWorkMemory<double> tmp(this->rwrkmem);
  DifVar idife = DifVar::None;
  if(idiff == DifVar::Phys){
    dmet0 = tmp.allocate(nnmet*gdim);
    idife = DifVar::Bary;
  }

  if(tdimn > gdim) METRIS_THROW_MSG( "getMetBary0 topo dim > gdim");

  const auto evaluate_log_metric = [&](const auto &field,
                                       const int *field_nodes){
    if(tdimn == 1){
      eval1<nnmet,ideg>(field,field_nodes,this->ibasis,idife,
                        DifVar::None,bary,metl,dmet0,NULL);
    }else if(tdimn == 2){
      eval2<nnmet,ideg>(field,field_nodes,this->ibasis,idife,
                        DifVar::None,bary,metl,dmet0,NULL);
    }else{
      eval3<nnmet,ideg>(field,field_nodes,this->ibasis,idife,
                        DifVar::None,bary,metl,dmet0,NULL);
    }
  };

  // Metric interpolation is always log-Euclidean. If the field happens to be
  // stored in exponential (SPD) space, convert each local FE coefficient to
  // log space before applying the basis functions. Converting only the final
  // interpolated matrix would instead compute log(sum_i phi_i M_i), which is
  // not log-Euclidean interpolation.
  if(this->ispace == MetSpace::Log){
    evaluate_log_metric(rfld,ent2pol);
  }else{
    const int nnode = getnnode(tdimn,ideg);
    std::array<double,max_nnode*nnmet> log_coefficients;
    std::array<int,max_nnode> local_nodes;
    for(int inode = 0; inode < nnode; inode++){
      local_nodes[inode] = inode;
      for(int imet = 0; imet < nnmet; imet++){
        log_coefficients[inode*nnmet + imet]
            = rfld(ent2pol[inode],imet);
      }
      getlogmet_inp<gdim,double>(
          &log_coefficients[inode*nnmet]);
    }
    const dblAr2 log_field(nnode,nnmet,log_coefficients.data());
    evaluate_log_metric(log_field,local_nodes.data());
  }


  if(idiff == DifVar::None){
    // If not returning now, we'll do this with SurrealS expmet
    if(tarspac == MetSpace::Exp) getexpmet_inp<gdim,double>(metl);
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

  if(tarspac == MetSpace::Exp){
    SANS::SurrealS<gdim,double> metS[nnmet];
    getmet_dbl2SurS<gdim,gdim>(metl,dmet,metS);
    getexpmet_inp<gdim,SANS::SurrealS<gdim,double>>(metS);
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
