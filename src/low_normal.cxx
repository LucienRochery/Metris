//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_normal.hxx"
#include "low_eval.hxx"
#include "low_geo.hxx"
#include "linalg/det.hxx"
#include "Mesh/MeshBase.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "utils/mprintf.hxx"
#include "utils/CT_loop.hxx"
#include "io_libmeshb.hxx"
#include "metris_constants.hxx"

namespace Metris{



void getnorfacP1(const int *fac2pol, const dblAr2 &coord, double *nrmal){
  METRIS_ASSERT(coord.get_stride() == 3);
  vecprod_vdif(coord[fac2pol[1]],coord[fac2pol[0]],
               coord[fac2pol[2]],coord[fac2pol[0]],nrmal);
}

void getnorfac(const MeshBase &msh, int iface, 
               const double *bary, AsDeg asdmsh, double *nrmal){

  METRIS_ASSERT(msh.idim == 3);

  if(asdmsh == AsDeg::P1 || msh.curdeg == 1){
    getnorfacP1(msh.fac2poi[iface], msh.coord, nrmal);
    return;
  }

  double coop[3], jmat[2][3];
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
    eval2<3, ideg>(msh.coord, msh.fac2poi[iface], msh.getBasis(), 
                   DifVar::Bary, DifVar::None, 
                   bary, coop, jmat[0], NULL);
    vecprod(jmat[0],jmat[1],nrmal);
  }}CT_FOR1(ideg);
  return;
}

static int warning_print_CADnor = 0;
// To compute the CAD normal, the safest is to average the vertex normals.
// This is because taking the average of the (u,v)'s can send us just about
// anywhere.
int getnorfacCAD(const MeshBase &msh, int iface, double *nrmal){
  GETVDEPTH(msh.param);
  bool oneOK = false;
  for(int ii = 0; ii < 3; ii++) nrmal[ii] = 0;
  for(int iver = 0; iver < 3; iver++){
    int ipoin = msh.fac2poi(iface,iver);
    int ibpoi = msh.poi2ebp(ipoin,2,iface,-1);
    CPRINTF2(" - getnorfacCAD iface %d iver %d ipoin %d ibpoi %d\n", iface, iver, ipoin, ibpoi);
    METRIS_ASSERT_MSG(ibpoi >= 0," boundary point "<<ipoin<<" has no ibpoi");

    double dum[3];
    if(getnorpoiCAD2(msh,ibpoi,dum)){
      if(warning_print_CADnor++ < 10){
        CPRINTF2(" # ibpoi %d ipoin %d skipped, possible singularity\n",ibpoi,ipoin);
        if(warning_print_CADnor >= 10) CPRINTF2(" # suppressing print\n");
      }
      continue;
    }

    if(DOPRINTS3()){
      CPRINTF3(" - ipoin %d ibpoi %d (u,v) = %e %e +nor ",ipoin,ibpoi,
               msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1));
      dblAr1(3,dum).print();
    }
    oneOK = true;
    for(int ii = 0; ii < 3; ii++) nrmal[ii] += dum[ii];
  }

  METRIS_ASSERT(oneOK);

  if(oneOK) return 0;
  return 1;
}


// Return outgoing normal of edge (2D only)
int getnorpoiCAD1(const MeshBase &msh, int ipoin, std::map<ego,int> &edgorient, 
                  double *norpoi){
  METRIS_ASSERT(msh.idim == 2);

  int ibpoi = msh.poi2bpo[ipoin];
  METRIS_ASSERT(ibpoi >= 0);

  int itype = msh.bpo2ibi(ibpoi,1);
  if(itype == 1){
    int iedge = msh.bpo2ibi(ibpoi,2);
    int iref  = msh.edg2ref[iedge];
    ego obj = msh.CAD.cad2edg[iref];

    double result[18];
    int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
    if(ierro != 0) return ierro;

    double *du = &result[3];

    int isens = edgorient[obj];

    norpoi[0] =  isens*du[1];
    norpoi[1] = -isens*du[0];
  }else{
    // Else is corner 
    norpoi[0] = 0;
    norpoi[1] = 0;
    do{
      ibpoi = msh.bpo2ibi(ibpoi,3);
      if(ibpoi < 0) break;

      int iedge = msh.bpo2ibi(ibpoi,2);
      int iref  = msh.edg2ref[iedge];
      ego obj   = msh.CAD.cad2edg[iref];

      double result[18];
      int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
      if(ierro != 0) return ierro;

      double *du = &result[3];

      int isens = edgorient[obj];

      norpoi[0] +=  isens*du[1];
      norpoi[1] += -isens*du[0];
    }while(ibpoi >= 0);
  }

  double nrm = msh.idim == 2 ? getnrml2<2>(norpoi) : getnrml2<3>(norpoi);
  METRIS_ENFORCE(nrm >= 1.0e-32);
  nrm = 1.0/sqrt(nrm);
  for(int ii = 0; ii < msh.idim; ii++) norpoi[ii] *= nrm;

  return 0;
}


// Return outgoing normal of face (3D only)
int getnorpoiCAD2(const MeshBase &msh, int ibpoi, double *norpoi){

  METRIS_ASSERT(ibpoi >= 0);
  METRIS_ASSERT(msh.bpo2ibi(ibpoi,1) == 2);
  METRIS_ASSERT(msh.CAD());

  int iface = msh.bpo2ibi(ibpoi,2);
  METRIS_ASSERT(iface >= 0);
  int iref  = msh.fac2ref[iface];
  METRIS_ASSERT(iref >= 0);
  ego obj   = msh.CAD.cad2fac[iref];
  METRIS_ASSERT(obj != NULL);

  int mtype = obj->mtype;
  METRIS_ASSERT(mtype == 1 || mtype == -1);

  double result[18];
  int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
  if(ierro != 0) return ierro;
  double *du = &result[3];
  double *dv = &result[6];
  
  vecprod(du,dv,norpoi);

  if(normalize_vec<3>(norpoi)) return 1;

  for(int ii = 0; ii < 3; ii++) norpoi[ii] *= mtype;

  return 0;
}



template <int ideg>
void getnorballref(MeshBase &msh, const intAr1 &lball, int iref, double* norpoi){
  // Discrete 
  for(int ii = 0; ii < 3; ii++) norpoi[ii] = 0;
  double norfac[3];


  for(int iface : lball){
    INCVDEPTH(msh.param);

    METRIS_ASSERT(!isdeadent(iface,msh.fac2poi));

    int iref2 = msh.fac2ref[iface];
    if(iref2 != iref && iref >= 0) continue;

    if constexpr(ideg == 1){
      getnorfacP1(msh.fac2poi[iface],msh.coord,norfac);
      CPRINTF1(" - iface %d vertices %d %d %d normal %e %e %e \n",iface,
        msh.fac2poi(iface,0),msh.fac2poi(iface,1),msh.fac2poi(iface,2),
        norfac[0], norfac[1], norfac[2]);
    }else{
      METRIS_THROW_MSG(TODOExcept(),"Implement normal computation HO");
    }
    // Note the normal is already area weighted.
    for(int ii = 0; ii < 3; ii++) norpoi[ii] += norfac[ii];
  }
  if(normalize_vec<3>(norpoi) != 0){
    writeMesh("debug_zero_normal.meshb",msh);
    printf("Normal norm = %e ball size %d\n",sqrt(getnrml2<3>(norpoi)),lball.get_n());
    for(int iface : lball){
      int iref2 = msh.fac2ref[iface];
      if(iref2 != iref && iref >= 0) continue;
      getnorfacP1(msh.fac2poi[iface],msh.coord,norfac);
      printf(" - iface %d vertices %d %d %d normal %e %e %e norm %e \n",iface,
        msh.fac2poi(iface,0),msh.fac2poi(iface,1),msh.fac2poi(iface,2),
        norfac[0], norfac[1], norfac[2], sqrt(getnrml2<3>(norfac)));
    }
    METRIS_THROW_MSG(GeomExcept(),"normal vanishes");
  }
  METRIS_ENFORCE(normalize_vec<3>(norpoi) == 0);
}


// Compute normal of point ipoin using CAD
// iref can be provided as a constraint. If < 0, use all faces, otherwise only
// matching iref.
template <int ideg>
void getnorpoiref(const MeshBase &msh, int ipoin, int iref, double* norpoi){
  METRIS_ASSERT(msh.idim == 3);

  METRIS_ASSERT(msh.CAD());

  // Actually it's free when called from some cavity callers
  //if(msh.CAD()) METRIS_ASSERT(nball == 0); // We don't need this, bpos give us all

  for(int ii = 0; ii < 3; ii++) norpoi[ii] = 0;
  double norfac[3];


  // Face point -> ref or not we can get the unique normal 
  if(msh.CAD()){

    double result[18];
    double *du,*dv;
    double nrm;


    // Whether tdimp 2 or less, we can do this loop, it'll have 1 iter if tdimp == 2 !
    // This is mainly because of periodic surface
    for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0 ; ibpoi = msh.bpo2ibi(ibpoi,3)){
      int bdim = msh.bpo2ibi(ibpoi,1);
      if(bdim != 2) continue;

      int ientt = msh.bpo2ibi(ibpoi,2);

      int iref2 = msh.fac2ref[ientt];
      METRIS_ASSERT(iref2 >= 0);
      if(iref2 != iref && iref >= 0) continue;

      ego obj = msh.CAD.cad2fac[iref2];

      int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
      METRIS_ASSERT(ierro == 0);

      du = &result[3];
      dv = &result[6];
  
      vecprod(du,dv,norfac);
  
      nrm = sqrt(getnrml2<3>(norfac));
      METRIS_ASSERT_MSG(nrm > 1.0e-16, "nrm = "<<nrm);
  
      for(int ii = 0; ii < 3; ii++) norpoi[ii] += norfac[ii] / nrm;
    }
    
    nrm = sqrt(getnrml2<3>(norpoi));
    METRIS_ASSERT_MSG(nrm > 1.0e-16, "(2) nrm = "<<nrm);

    for(int ii = 0; ii < 3; ii++) norpoi[ii] /= nrm;

    return;
  }



}

#define BOOST_PP_LOCAL_MACRO(n)\
template void getnorballref<n>(MeshBase &msh, const intAr1& lball, int iref, double* norpoi);\
template void getnorpoiref<n>(const MeshBase &msh, int ipoin, int iref, double* norpoi);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



}// namespace
