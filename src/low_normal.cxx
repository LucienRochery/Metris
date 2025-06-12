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

// To compute the CAD normal, the safest is to average the vertex normals.
// This is because taking the average of the (u,v)'s can send us just about
// anywhere.
int getnorfacCAD(const MeshBase &msh, int iface, double *nrmal){
  GETVDEPTH(msh.param);
  bool oneOK = false;
  for(int ii = 0; ii < 3; ii++) nrmal[ii] = 0;
  for(int iver = 0; iver < 3; iver++){
    INCVDEPTH(msh.param);
    int ipoin = msh.fac2poi(iface,iver);
    int ibpoi = msh.poi2ebp(ipoin,2,iface,-1);
    CPRINTF2(" - getnorfacCAD iface %d iver %d ipoin %d ibpoi %d\n", iface, iver, ipoin, ibpoi);
    METRIS_ASSERT_MSG(ibpoi >= 0," boundary point "<<ipoin<<" has no ibpoi");

    double dum[3];
    static int warning_print_CADnor = 0;
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


// Compute normal of point ipoin using CAD or discretely
// iref can be provided as a constraint. If < 0, use all faces, otherwise only
// matching iref.
void getnorpoiref(const MeshBase &msh, int ipoin, int iref, double* norpoi){
  METRIS_ASSERT(msh.idim == 3);

  METRIS_ASSERT(msh.CAD());

  // Actually it's free when called from some cavity callers
  //if(msh.CAD()) METRIS_ASSERT(nball == 0); // We don't need this, bpos give us all

  for(int ii = 0; ii < 3; ii++) norpoi[ii] = 0;
  double norfac[3];


  // Face point -> ref or not we can get the unique normal 

  double result[18];
  double *du,*dv;
  double nrm;


  // Whether tdimp 2 or less, we can do this loop, it'll have 1 iter if tdimp == 2 !
  // This is mainly because of periodic surface
  for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0 ; ibpoi = msh.bpo2ibi(ibpoi,3)){
    int bdim = msh.bpo2ibi(ibpoi,1);
    if(bdim != 2) continue;

    int iface = msh.bpo2ibi(ibpoi,2);

    int iref2 = msh.fac2ref[iface];
    METRIS_ASSERT(iref2 >= 0);
    if(iref2 != iref && iref >= 0) continue;

    if(msh.CAD()){
      ego obj = msh.CAD.cad2fac[iref2];

      int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
      METRIS_ASSERT(ierro == 0);

      du = &result[3];
      dv = &result[6];

      vecprod(du,dv,norfac);
    }else{
      double bary[3] = {0, 0, 0};
      // bary is unused if ideg == 1
      if(msh.curdeg > 1){
        int iver = msh.template getverfac<1>(iface, ipoin);
        METRIS_ASSERT(iver >= 0);
        bary[iver] = 1;
      }
      getnorfac(msh, iface, bary, AsDeg::Pk, norfac);
    }
    
    if(normalize_vec<3>(norfac)) METRIS_THROW_MSG(GeomExcept(), "norfac nrm = "<<nrm);

    for(int ii = 0; ii < 3; ii++) norpoi[ii] += norfac[ii];
  }
  
  if(normalize_vec<3>(norpoi)) METRIS_THROW_MSG(GeomExcept(), "norpoi nrm = "<<nrm);

}

#define BOOST_PP_LOCAL_MACRO(n)\
template void getnorballref<n>(MeshBase &msh, const intAr1& lball, int iref, double* norpoi);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



int gettanpoiref(const MeshBase &msh, int ipoin, int iref, double* tanpoi){
  METRIS_ASSERT(iref >= 0);

  const int gdim = msh.idim;

  int tdim  = msh.poi2ent(ipoin,1);
  if(tdim > 1) return 1;

  if(!msh.CAD()){

    int iedg1 = msh.poi2ent(ipoin,0);
    int iver1 = msh.template getveredg<1>(iedg1, ipoin);
    int iref1 = msh.edg2ref[iedg1];

    int iedg2 = msh.edg2edg(iedg1, 1-iver1);
    int iver2 = msh.template getveredg<1>(iedg2, ipoin);
    int iref2 = msh.edg2ref[iedg2];

    if(iref >= 0 && iref1 != iref && iref2 != iref) return 2;

    // If only one to be considered, make sure it's the first
    bool skipref2 = false;
    if(iref >= 0 && iref1 != iref){
      skipref2 = true;

      int itmp = iedg1;
      iedg1 = iedg2;
      iedg2 = itmp;

      itmp = iref1;
      iref1 = iref2;
      iref2 = itmp;

      itmp = iver1;
      iver1 = iver2;
      iver2 = itmp;
    }

    double tane1[3], tane2[3];
    if(msh.curdeg == 1){
      for(int ii = 0; ii < gdim; ii ++) tane1[ii] = msh.coord(ipoin,ii)
                                                  - msh.coord(msh.edg2poi(iedg1,1-iver1),ii);
      if(!skipref2){
        for(int ii = 0; ii < gdim; ii ++) tane2[ii] = msh.coord(msh.edg2poi(iedg2,1-iver2),ii)
                                                    - msh.coord(ipoin,ii);
      }
    }else{
      double bary[2] = {0, 0}, dum[3];
      bary[iver1] = 1;
      CT_FOR0_INC(2,3,gdim_c){if(gdim_c == gdim){
      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
        eval1<gdim_c,ideg>(msh.coord, msh.edg2poi[iedg1], msh.getBasis(), 
                           DifVar::Bary, DifVar::None, 
                           bary, dum, tane1, NULL);
        if(!skipref2){
          bary[iver1] = 0;
          bary[iver2] = 1;
          eval1<gdim_c,ideg>(msh.coord, msh.edg2poi[iedg2], msh.getBasis(), 
                             DifVar::Bary, DifVar::None, 
                             bary, dum, tane2, NULL);
        }
      }}CT_FOR1(ideg);
      }}CT_FOR1(gdim_c);
    }

    if(!skipref2){
      // Average normalized tans weighted by inverse of the edge lengths. 
      // A proxy for this is the tangent norm. 
      // So we're doing the following average (of the unnormalized):
      // tane1 / norm(tane1)^2 + tane2 / norm(tane2)^2
      double nrm1 = gdim == 2 ? getnrml2<2>(tane1) : getnrml2<3>(tane1);
      double nrm2 = gdim == 2 ? getnrml2<2>(tane2) : getnrml2<3>(tane2);

      // vecNrmTol is for the squared norm as we have here
      METRIS_ASSERT(nrm1 > Constants::vecNrmTol);
      METRIS_ASSERT(nrm2 > Constants::vecNrmTol);

      for(int ii = 0; ii < gdim; ii++) 
        tanpoi[ii] = tane1[ii] / nrm1 + tane2[ii] / nrm2;
    }else{
      for(int ii = 0; ii < gdim; ii++) tanpoi[ii] = tane1[ii];
    }

  }else{

    if(iref >= 0){
      // Check loop case and throw error in debug mode. 
      #ifndef NDEBUG
      if(msh.getpoitdim(ipoin) == 0){
        int nref = 0;
        for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
          tdim = msh.bpo2ibi(ibpoi, 1);
          if(tdim != 1) continue;
          int ientt = msh.bpo2ibi(ibpoi,2);
          int iref2 = msh.edg2ref[ientt];
          if(iref2 == iref) nref++;
        }
        if(nref > 2) METRIS_THROW_MSG(TopoExcept(), 
          " ## "<<nref<<" edges of same ref at point "<<ipoin)
        if(nref == 2)METRIS_THROW_MSG(TODOExcept(), 
          " ## Implement case "<<nref<<" edges of same ref at point "<<ipoin)
      }
      #endif

      ego obj = msh.CAD.cad2edg[iref];
      int ibpoi = msh.poi2ebp(ipoin,1,-1,iref);
      double result[18];
      int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
      if(ierro != EGADS_SUCCESS) return 10+ierro;

      for(int ii = 0; ii < gdim; ii++) tanpoi[ii] = result[3+ii];
    }else{
      // Simply average all normalized tangents at the point. 
      // If point is dim 1, we can get iref from ent2poi and do a single eval. 
      // There would only be one entry anyways
      double result[18];
      for(int ii = 0; ii < gdim; ii++) tanpoi[ii] = 0;
      for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
        tdim = msh.bpo2ibi(ibpoi, 1);
        if(tdim == 0) continue;
        if(tdim  > 1) break;
        int iedge = msh.bpo2ibi(ibpoi, 2);
        int iref2 = msh.edg2ref[iedge];
        ego obj = msh.CAD.cad2edg[iref2];
        int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
        if(ierro != EGADS_SUCCESS) return 10+ierro;

        ierro = gdim == 2 ? normalize_vec<2>(&result[3])
                          : normalize_vec<3>(&result[3]);
        if(ierro != 0) return 5;

        for(int ii = 0; ii < gdim; ii++) tanpoi[ii] += result[3+ii];
      }
    }

  }

  int ierro = gdim == 2 ? normalize_vec<2>(tanpoi)
                        : normalize_vec<3>(tanpoi);
  if(ierro != 0) return 3;


  return 0;
}


}// namespace
