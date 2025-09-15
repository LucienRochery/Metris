//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "normal.hxx"
#include "misc.hxx"

#include "../low_eval.hxx"
#include "../linalg/det.hxx"
#include "../Mesh/MeshBase.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"
#include "../utils/CT_loop.hxx"
#include "../io_libmeshb.hxx"
#include "../metris_constants.hxx"

namespace Metris{



void getnorfacP1(const int *fac2pol, const dblAr2 &coord, double *nrmal){
  METRIS_ASSERT(coord.get_stride() == 3);
  vecprod_vdif(coord[fac2pol[1]],coord[fac2pol[0]],
               coord[fac2pol[2]],coord[fac2pol[0]],nrmal);
}

void getnorfac(const MeshBase&__restrict__ msh, int iface, 
               const double*__restrict__ bary, AsDeg asdmsh, double*__restrict__ nrmal){
  return getnorfac(msh, msh.fac2poi[iface], bary, asdmsh, nrmal);
}


void getnorfac(const MeshBase&__restrict__ msh, const int*__restrict__ fac2pol, 
               const double*__restrict__ bary, AsDeg asdmsh, double*__restrict__ nrmal){

  METRIS_ASSERT(msh.idim == 3);

  if(asdmsh == AsDeg::P1 || msh.curdeg == 1){
    getnorfacP1(fac2pol, msh.coord, nrmal);
    return;
  }

  double coop[3], jmat[2][3];
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
    eval2<3, ideg>(msh.coord, fac2pol, msh.getBasis(), 
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
  int nod2bpo[3];
  nod2bpo[0] = msh.poi2ebp(msh.fac2poi(iface,0), 2, iface, msh.fac2ref[iface]);
  nod2bpo[1] = msh.poi2ebp(msh.fac2poi(iface,1), 2, iface, msh.fac2ref[iface]);
  nod2bpo[2] = msh.poi2ebp(msh.fac2poi(iface,2), 2, iface, msh.fac2ref[iface]);

  return getnorfacCAD(msh, msh.fac2poi[iface], nod2bpo, nrmal);
}

int getnorfacCAD(const MeshBase&__restrict__ msh, const int*__restrict__ fac2pol,
                 const int*__restrict__ nod2bpo, double*__restrict__ nrmal){
  GETVDEPTH(msh.param);
  bool oneOK = false;
  for(int ii = 0; ii < 3; ii++) nrmal[ii] = 0;
  for(int iver = 0; iver < 3; iver++){
    INCVDEPTH(msh.param);
    int ibpoi = nod2bpo[iver];
    if(ibpoi < 0) continue;

    int ipoin = fac2pol[iver];
    CPRINTF3(" - getnorfacCAD iver {} ipoin {} ibpoi {}\n", iver, ipoin, ibpoi);

    double dum[3];
    static int warning_print_CADnor = 0;
    if(getnorpoiCAD2(msh,ibpoi,dum)){
      if(warning_print_CADnor++ < 10){
        CPRINTF3(" # ibpoi {} ipoin {} skipped, possible singularity\n",ibpoi,ipoin);
        if(warning_print_CADnor >= 10) CPRINTF2(" # suppressing print\n");
      }
      continue;
    }

    CPRINTF3(" - ipoin {} ibpoi {} (u,v) = {} {} +nor {} {} {}\n",ipoin,ibpoi,
             msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1),
             dum[0],dum[1],dum[2]);
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
  GETVDEPTH(msh.param);
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
      CPRINTF1(" - iface {} vertices {} {} {} normal {} {} {} \n",iface,
        msh.fac2poi(iface,0),msh.fac2poi(iface,1),msh.fac2poi(iface,2),
        norfac[0], norfac[1], norfac[2]);
    }else{
      METRIS_THROW_MSG("TODO: Implement normal computation HO");
    }
    // Note the normal is already area weighted.
    for(int ii = 0; ii < 3; ii++) norpoi[ii] += norfac[ii];
  }
  if(normalize_vec<3>(norpoi) != 0){
    writeMesh("debug_zero_normal.meshb",msh);
    MPRINTF("Normal norm = {} ball size {}\n",sqrt(getnrml2<3>(norpoi)),lball.get_n());
    for(int iface : lball){
      int iref2 = msh.fac2ref[iface];
      if(iref2 != iref && iref >= 0) continue;
      getnorfacP1(msh.fac2poi[iface],msh.coord,norfac);
      MPRINTF(" - iface {} vertices {} {} {} normal {} {} {} norm {} \n",iface,
        msh.fac2poi(iface,0),msh.fac2poi(iface,1),msh.fac2poi(iface,2),
        norfac[0], norfac[1], norfac[2], sqrt(getnrml2<3>(norfac)));
    }
    METRIS_THROW_MSG("normal vanishes");
  }
  METRIS_ENFORCE(normalize_vec<3>(norpoi) == 0);
}


// Compute normal of point ipoin using CAD or discretely
// iref can be provided as a constraint. If < 0, use all faces, otherwise only
// matching iref.
void getnorpoiref(const MeshBase &msh, int ipoin, int iref, double* norpoi){
  METRIS_ASSERT(msh.idim == 3);
  GETVDEPTH(msh.param);

  // Actually it's free when called from some cavity callers
  //if(msh.CAD()) METRIS_ASSERT(nball == 0); // We don't need this, bpos give us all

  for(int ii = 0; ii < 3; ii++) norpoi[ii] = 0;
  double norfac[3];


  // Face point -> ref or not we can get the unique normal 

  double result[18];
  double *du,*dv;


  // Whether tdimp 2 or less, we can do this loop, it'll have 1 iter if tdimp == 2 !
  // This is mainly because of periodic surface
  int pdim = -1;
  ego obj;
  for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0 ; ibpoi = msh.bpo2ibi(ibpoi,3)){
    int bdim = msh.bpo2ibi(ibpoi,1);
    if(pdim < 0) pdim = bdim;
    if(bdim != 2) continue;

    int iface = msh.bpo2ibi(ibpoi,2);

    int iref2 = msh.fac2ref[iface];
    METRIS_ASSERT(iref2 >= 0);
    if(iref2 != iref && iref >= 0) continue;

    if(msh.CAD()){
      obj = msh.CAD.cad2fac[iref2];

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
    
    if(normalize_vec<3>(norfac)){
      // Attempt recovery by using a nearby point on the triangle.
      bool irecovered = false;
      if(msh.CAD() && pdim < 2){
        double lbpo[3];
        int iver = msh.template getverfac<1>(iface, ipoin);
        METRIS_ASSERT(iver >= 0);
        for(int ii = 0; ii < 3; ii++){
          int ipoi2 = msh.fac2poi(iface, ii);
          lbpo[ii] = msh.poi2ebp(ipoi2, 2, iface, -1);
          METRIS_ASSERT(lbpo[ii] >= 0);
        }
        double uv[2];
        for(int ii = 0; ii < 2; ii++){
          uv[ii] = 0.9*msh.bpo2rbi(lbpo[iver], ii) 
                 + 0.1*msh.bpo2rbi(lbpo[(iver+1)%3], ii)
                 + 0.1*msh.bpo2rbi(lbpo[(iver+2)%3], ii);
        }
        int ierro = EG_evaluate(obj, uv, result);
        METRIS_ASSERT(ierro == 0);

        du = &result[3];
        dv = &result[6];
        vecprod(du,dv,norfac);

        if(!normalize_vec<3>(norfac)){
          irecovered = true;
          CPRINTF3(" - recovered normal for ipoin {}, iref {} iface {} "
                    "iver {} uv = {} {}\n", ipoin, iref, iface, iver,
                    uv[0], uv[1]);
        }
      }

      if(!irecovered){
        MPRINTF("norfac vanishes in getnorpoiref for ipoin {}, iref {}\n", 
              ipoin, iref);
        MPRINTF("iface = {} nodes {}\n",iface,
                intAr1(getnnode(msh.curdeg,2),msh.fac2poi[iface]));

        MPRINTF("Using msh.CAD() = {}\n",msh.CAD());
        MPRINTF("result = {}\n",dblAr1(18,result));

        METRIS_THROW_MSG("## norfac vanishes");
      }
    }

    for(int ii = 0; ii < 3; ii++) norpoi[ii] += norfac[ii];
  }
  
  if(normalize_vec<3>(norpoi)) METRIS_THROW_MSG("## norpoi vanishes");

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
    if(iref >= 0 && iedg2 < 0 && iref1 != iref) return 3;

    int iver2 = -1, iref2 = iref;
    bool skipref2 = false;
    if(iedg2 >= 0){
      iver2 = msh.template getveredg<1>(iedg2, ipoin);
      iref2 = msh.edg2ref[iedg2];
    }else{
      skipref2 = true;
    }

    if(iref >= 0 && iref1 != iref && iref2 != iref) return 2;

    // If only one to be considered, make sure it's the first
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
        METRIS_ASSERT_MSG(nref <= 2, 
          " ## {} edges of same ref at point {}", nref, ipoin)
        if(nref == 2) METRIS_THROW_MSG( 
          "TODO: Implement case {} edges of same ref at point {}", nref, ipoin)
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

// Normalized dotprod difference to 1 of CAD/elt normals accumulated over the nodes
// Prefer calling this over the callee.
template<int ideg>
double getnordev(const MeshBase&__restrict__ msh, int iface, const double*__restrict__ norfac){
  constexpr int nnode = getnnode(2, ideg);
  int nod2bpo[nnode];
  int iref = msh.fac2ref[iface];
  for(int ii = 0; ii < nnode; ii++) nod2bpo[ii] = msh.poi2ebp(msh.fac2poi(iface,ii), 2, iface, iref);
  return getnordev<ideg>(msh, msh.fac2poi[iface], nod2bpo, norfac);
}
#define BOOST_PP_LOCAL_MACRO(n)\
template double getnordev<n>(const MeshBase&__restrict__ msh, int iface, const double*__restrict__ norfac_);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



// Normalized dotprod difference to 1 of CAD/elt normals accumulated over the nodes
// nod2bpo holds the ibpois for the element nodes. 
template<int ideg>
double getnordev(const MeshBase&__restrict__ msh, const int*__restrict__ fac2pol, 
                 const int*__restrict__ nod2bpo, const double*__restrict__ norfac_){
  METRIS_ASSERT(msh.idim == 3);
  GETVDEPTH(msh.param);

  CPRINTF1("-- START getnordev fac2pol = {} {} {} nod2bpo = {} {} {}\n",fac2pol[0],fac2pol[1],fac2pol[2],
           nod2bpo[0],nod2bpo[1],nod2bpo[2]);
  if(DOPRINTS2()){
    if(nod2bpo[0] >= 0) CPRINTF2(" - nod2bpo[0] = ibpoi {} : {}\n",nod2bpo[0],intAr1(nibi,msh.bpo2ibi[nod2bpo[0]]));
    if(nod2bpo[1] >= 0) CPRINTF2(" - nod2bpo[1] = ibpoi {} : {}\n",nod2bpo[1],intAr1(nibi,msh.bpo2ibi[nod2bpo[1]]));
    if(nod2bpo[2] >= 0) CPRINTF2(" - nod2bpo[2] = ibpoi {} : {}\n",nod2bpo[2],intAr1(nibi,msh.bpo2ibi[nod2bpo[2]]));
  }

  constexpr int gdim = 3;
  constexpr int tdim = 2;
  constexpr int nnode = getnnode(tdim, ideg);

  METRIS_ASSERT(nod2bpo != NULL);
  // At least one of the vertices (not nodes) should have valid ibpoi.
  METRIS_ASSERT(nod2bpo[0] >= 0 || nod2bpo[1] >= 0 || nod2bpo[2] >= 0);
  METRIS_ASSERT(nod2bpo[0] < msh.nbpoi || nod2bpo[1]  < msh.nbpoi || nod2bpo[2]  < msh.nbpoi);
  const int iref = nod2bpo[0] >= 0 ? msh.fac2ref[msh.bpo2ibi(nod2bpo[0],2)] :
                   nod2bpo[1] >= 0 ? msh.fac2ref[msh.bpo2ibi(nod2bpo[1],2)] :
                                     msh.fac2ref[msh.bpo2ibi(nod2bpo[2],2)] ;
  METRIS_ASSERT(iref >= 0 && iref < msh.CAD.ncadfa);
  #ifndef NDEBUG
  for(int inode = 0; inode < nnode; inode++){
    if(nod2bpo[inode] < 0) continue;
    METRIS_ASSERT_MSG(iref == msh.fac2ref[msh.bpo2ibi(nod2bpo[inode],2)] , 
      "iref {}, inode {}, fac2poi {}, ibpoi {} : {},fac2ref = {}",
      iref, inode, fac2pol[inode], nod2bpo[inode], intAr1(nibi,msh.bpo2ibi[nod2bpo[inode]]),
      msh.fac2ref[msh.bpo2ibi(nod2bpo[inode],2)]);
  }
  #endif

  if(!msh.CAD()){
    static int nwarnprt = 0;
    if(nwarnprt++ < 10) MPRINTF("\n\n ## IMPLEMENT DISCRETE NORMALS IN GETNORDEV\n\n");
    return 0.0;
  }

  double result[18];
  double norCAD[gdim], norelt[gdim];
  double *du = &result[3];
  double *dv = &result[6];
  const ego obj  = msh.CAD.cad2fac[iref];

  // Even if we use CAD normals at all vertices, we can compute this one just once.
  if(ideg == 1 && norfac_ == NULL){
    getnorfacP1(fac2pol, msh.coord, norelt);
    if(normalize_vec<gdim>(norelt)){
      writeMesh("debug_ibpoi",msh);
      MPRINTF("norfac vanished face nodes {}\n", intAr1(nnode, fac2pol));
      for(int ii = 0; ii < gdim; ii++) MPRINTF("{}: {:23.15e}\n",ii,norelt[ii]);
      METRIS_THROW_MSG("Normal (elt) vanishes");
    }
  }else if(ideg == 1 && norfac_ != NULL){
    for(int ii = 0; ii < gdim; ii++) norelt[ii] = norfac_[ii];
    if(normalize_vec<gdim>(norelt)){
      MPRINTF("## PROVIDED NORFAC VANISHES\n");
      METRIS_THROW_MSG("Normal (elt) vanishes");
    }
  }


  double nordev = 0;
  int nsum = 0;
  for(int inode = 0; inode < nnode; inode++){
    INCVDEPTH(msh.param);
    int ipoin = fac2pol[inode];
    int ibpoi = nod2bpo[inode];
    if(ibpoi < 0){
      CPRINTF1(" - inode {} ipoin {} ibpoi {} skipped\n",inode,ipoin,ibpoi);
      continue;
    }

    int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
    METRIS_ENFORCE_MSG(ierro == 0, "metqua0 EG_evaluate error {}", ierro);
    vecprod(du,dv,norCAD);
    if(normalize_vec<gdim>(norCAD)){
      // legitimate if e.g. cone tip.
      if(msh.getpoitdim(ipoin) != 0){
        #if 0
        MPRINTF("ipoin = {} point dim = {}",ipoin,msh.getpoitdim(ipoin));
        MPRINTF("using ibpoi = {} : ",ibpoi);
        intAr1(nibi, msh.bpo2ibi[ibpoi]).print();
        int  idbgdim = msh.bpo2ibi(ibpoi,1);
        int  idbgent = msh.bpo2ibi(ibpoi,2);
        MPRINTF(" entity ref = {} \n", msh.ent2ref(idbgdim)[idbgent]);
        MPRINTF("get du = ");
        dblAr1(gdim, du).print();
        MPRINTF("get dv = ");
        dblAr1(gdim, dv).print();
        MPRINTF("vecprod = ");
        dblAr1(gdim, norCAD).print();

        MPRINTF("(u,v) = {:24.15e} {:24.15e}\n", msh.bpo2rbi(ibpoi,0), msh.bpo2rbi(ibpoi,1));
        MPRINTF("eval coop {:24.15e} {:24.15e} {:24.15e} \n",
          result[0],result[1],result[2]);

        double nrm = getnrml2<gdim>(norCAD);
        MPRINTF("nrm = {:24.15e}\n",nrm);
        #endif
        METRIS_THROW_MSG("Normal (CAD) vanishes at ipoin {}", ipoin);
      }
      nordev += 0;
      continue;
    }

    if constexpr (ideg > 1){
      constexpr auto ordelt = ORDELT(tdim);
      double bary[tdim+1];
      for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = ordelt[ideg][inode][ii]/((double) (ideg));
      getnorfac(msh, fac2pol, bary, AsDeg::Pk, norelt);
      if(normalize_vec<gdim>(norelt)){
        writeMesh("debug_ibpoi",msh);
        MPRINTF("norfac vanished face node {} point {} nodes {}\n",inode,ipoin,
                intAr1(nnode, fac2pol));
        for(int ii = 0; ii < gdim; ii++) MPRINTF("{}: {:23.15e}\n",ii,norelt[ii]);
        METRIS_THROW_MSG("Normal (elt) vanishes");
      }
    }

    CPRINTF2(" - norCAD = {} {} {} , norfac = {} {} {}\n",norCAD[0],norCAD[1],norCAD[2],
             norelt[0],norelt[1],norelt[2]);
    //#ifndef NDEBUG
    //if(DOPRINTS3()){
    //  for(int ii = 0; ii < gdim; ii++){
    //    norsfac[inode][ii] = norelt[ii];
    //    norsCAD[inode][ii] = norCAD[ii];
    //  }
    //}
    //#endif
    double dtprd = getprdl2<gdim>(norelt, norCAD);
    double tmp = 1 - abs(dtprd);
    METRIS_ASSERT_MSG(tmp >= -1.0e-15, 
      "Normalized normals dotprod > ? tmp = {:e} (<0?)\n"
      "norelt = {}, norCAD = {}",tmp,dblAr1(3,norelt), dblAr1(3,norCAD));
    // Value can be slightly negative due to roundoff
    tmp = MAX(tmp, 0.0);
    nordev += abs(tmp);
    nsum++;
    CPRINTF1(" - face {} {} {} inode = {} ibpoi = {} local dev = {}\n",fac2pol[0],fac2pol[1],fac2pol[2],inode,ibpoi,tmp);
    CPRINTF2(" - ibpoi info = {} {} {} , ref = {}, t/(u,v) = {} {} \n", 
             msh.bpo2ibi(ibpoi,0), msh.bpo2ibi(ibpoi,1), 
             msh.bpo2ibi(ibpoi,2), 
             msh.bpo2ibi(ibpoi,1) > 0 ? msh.ent2ref(msh.bpo2ibi(ibpoi,1))[msh.bpo2ibi(ibpoi,2)] : -1,
             msh.bpo2rbi(ibpoi,0), msh.bpo2rbi(ibpoi,1));
  }
  CPRINTF1("-- END getnordev got nordeg = {}\n",nordev/nsum);
  //if(DOPRINTS3()){
  //  int facsave[nnode]; 
  //  int ifac0 = 0;
  //  for(int ii = 0; ii < inode; ii++) facsave[ii] = msh.fac2poi(ifac0,ii);
  //  for(int ii = 0; ii < inode; ii++) msh.fac2poi(ifac0,ii) = fac2pol[ii];
  //  writeMesh("debug_nordev.meshb",msh,msh.nedge,ifac0,msh.nelem,msh.nedge,ifac0+1,msh.nelem);
  //  for(int ii = 0; ii < inode; ii++) msh.fac2poi(ifac0,ii) = facsave[ii];
//
  //  dblWrkAr1 norpoifac = msh.get_rwork(3*msh.npoin);
  //  dblWrkAr1 norpoiCAD = msh.get_rwork(3*msh.npoin);
  //  for(int inode = 0; inode < nnode; inode++){
  //    for(int ii = 0; ii < 3; ii++) norpoifac[3*fac2pol[inode]+ii] = norsfac[inode][ii];
  //    for(int ii = 0; ii < 3; ii++) norpoiCAD[3*fac2pol[inode]+ii] = norsCAD[inode][ii];
  //  }
  //  writeField("debug_norfac.solb", msh, SolTyp::CG, norpoifac.get_array(), 3);
  //  writeField("debug_norCAD.solb", msh, SolTyp::CG, norpoiCAD.get_array(), 3);
  //}
  return nordev/nsum;
}
#define BOOST_PP_LOCAL_MACRO(n)\
template double getnordev<n>(const MeshBase&__restrict__ msh, const int*__restrict__ fac2pol, const int*__restrict__ nod2bpo, const double*__restrict__ norfac_);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



}// namespace
