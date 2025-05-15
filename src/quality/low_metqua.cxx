//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_metqua.hxx"

#include "../ho_constants.hxx"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "../aux_topo.hxx"
#include "../low_geo.hxx"
#include "../low_normal.hxx"
#include "../io_libmeshb.hxx"
#include "../Mesh/Mesh.hxx"
#include "../linalg/det.hxx"

#include "../utils/aux_pp_inc.hxx"


namespace Metris{



template <class MFT, int gdim, int tdim, QuaFun iquaf, typename ftype>
ftype metqua(Mesh<MFT> &msh, AsDeg asdmsh, AsDeg asdmet,
             int ientt, ftype difto){
	static_assert(tdim==2 || tdim==3);

  const intAr2 &ent2poi = msh.ent2poi(tdim);

  METRIS_ASSERT(!isdeadent(ientt, ent2poi));

  static_assert(gdim==2 || gdim==3);

  double bary[tdim+1];

  const int pnorm = msh.param->opt_pnorm;

  ftype qutet; 
  double nordev = 0;
  bool do_nordev = tdim == 2 && gdim == 3 
    && msh.CAD()
    && abs(msh.param->qua_surf_wt_normal) > 1.0e-9*abs(msh.param->qua_surf_wt_quality);
  if(tdim == 1 && gdim >= 2) 
    METRIS_THROW_MSG(TODOExcept(), "TODO: Edge quality with normal dev")


  // Performance impact should be zero
  constexpr auto quafun_xi = get_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
  constexpr auto ordelt = ORDELT(tdim);

  //// Compute normal at the nodes. This is then used to interpolate a normal
  //// within the element. Fewer EG_evaluate calls needed and more robust as 
  //// (u,v) interpolation followed by evaluation is not necessarily very stable.
  //double norfld[getnnode(tdim,METRIS_MAX_DEG)][gdim];
  //if(do_nordev){
  //  double result[18];
  //  double *du = &result[3];
  //  double *dv = &result[6];
  //  const int nnode = getnnode(tdim, asdmsh == AsDeg::P1 ? 1 : ideg);
  //  const int iref = msh.fac2ref[ientt];
  //  const ego obj  = msh.cad2fac[iref];

  //  for(int inode = 0; inode < nnode; inode++){
  //    int ipoin = ent2poi(ientt, inode);
  //    int ibpoi = msh.poi2ebp(ent2poi[ientt], tdim, ientt, iref);
  //    METRIS_ASSERT(ibpoi >= 0 && ibpoi < msh.nbpoi);
  //    int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
  //    METRIS_ENFORCE_MSG(ierro == 0, "metqua0 EG_evaluate error " << ierro);
  //    vecprod(du,dv,norfld[inode]);
  //    if(normalize_vec<gdim>(norfld[inode])) METRIS_THROW_MSG(GeomExcept(), "Normal vanishes");
  //  }
  //}
  const int ideg = msh.curdeg;

  // Accumulate normal error at the nodes (depending on asdmsh)
  if(do_nordev){
    double result[18];
    double norCAD[gdim], norelt[gdim];
    double *du = &result[3];
    double *dv = &result[6];
    const int ideg_eff = asdmsh == AsDeg::P1 ? 1 : ideg;

    const int nnode = getnnode(tdim, ideg_eff);
    const int iref = msh.fac2ref[ientt];
    const ego obj  = msh.CAD.cad2fac[iref];

    // Even if we use CAD normals at all vertices, we can compute this one just once.
    if(ideg == 1){
      getnorfacP1(ent2poi[ientt], msh.coord, norelt);
      if(normalize_vec<gdim>(norelt)){
        writeMesh("debug_ibpoi",msh);
        printf("norelt vanished ientt = %d nodes ",ientt);
        intAr1(nnode, ent2poi[ientt]).print();
        METRIS_THROW_MSG(GeomExcept(), "Normal (elt) vanishes");
      }
    }

    for(int inode = 0; inode < nnode; inode++){
      int ipoin = ent2poi(ientt, inode);
      int ibpoi = msh.poi2ebp(ipoin, tdim, ientt, iref);
      if(ibpoi < 0){
        for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
          printf("ibpoi %d : ",ibpoi);
          intAr1(nibi, msh.bpo2ibi[ibpoi]).print();
        }
        writeMesh("debug_ibpoi",msh);
      }
      METRIS_ASSERT_MSG(ibpoi >= 0 && ibpoi < msh.nbpoi, 
        "iface = "<<ientt<<" iref  "<<iref<<" inode = "<<inode
        <<" ipoin = "<<ipoin<<" ibpoi = "<< ibpoi);
      int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
      METRIS_ENFORCE_MSG(ierro == 0, "metqua0 EG_evaluate error " << ierro);
      vecprod(du,dv,norCAD);
      if(normalize_vec<gdim>(norCAD)){
        // legitimate if e.g. cone tip.
        if(msh.getpoitdim(ipoin) != 0){
          printf("ipoin = %d point dim = %d",ipoin,msh.getpoitdim(ipoin));
          printf("using ibpoi = %d : ",ibpoi);
          intAr1(nibi, msh.bpo2ibi[ibpoi]).print();
          printf("get du = ");
          dblAr1(gdim, du).print();
          printf("get dv = ");
          dblAr1(gdim, dv).print();
          printf("vecprod = ");
          dblAr1(gdim, norCAD).print();

          METRIS_THROW_MSG(GeomExcept(), "Normal (CAD) vanishes at ipoin "<<ipoin);
        }
        nordev += 0;
        continue;
      }

      if(ideg > 1){
        for(int ii = 0; ii < tdim + 1; ii++)
          bary[ii] = ordelt[ideg_eff][inode][ii]/((double) (ideg_eff));
        getnorfac(msh, ientt, bary, asdmsh, norelt);
        if(normalize_vec<gdim>(norelt)){
          writeMesh("debug_ibpoi",msh);
          printf("norelt vanished ientt = %d node %d point %d nodes ",ientt,inode,ipoin);
          intAr1(nnode, ent2poi[ientt]).print();
          METRIS_THROW_MSG(GeomExcept(), "Normal (elt) vanishes");
        }
      }

      double dtprd = getprdl2<gdim>(norelt, norCAD);
      double tmp = 1-dtprd;
      METRIS_ASSERT(tmp >= 0);
      nordev += tmp*tmp;
    }
    nordev /= nnode;
    nordev = sqrt(nordev);
  }

  if(asdmet == AsDeg::Pk && ideg > 1){
    qutet = 0.0;
    nordev= 0.0;

    const int idegj = SMOO_DEGJ(ideg);
    const int nnodj = tdim == 2 ? getnnod2(idegj) : getnnod3(idegj);



    for(int iquad = 0; iquad < nnodj; iquad++){
      for(int ii = 0; ii < tdim + 1; ii++){
        bary[ii] = ordelt[idegj][iquad][ii]/((double) (idegj));
      }
      ftype qua0 = quafun_xi(msh,asdmet,asdmsh,ent2poi[ientt],bary,NULL);
      qutet += pow(abs(qua0 - difto),pnorm) / nnodj;

    }

  }else{
    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim  + 1);
    #ifndef NDEBUG
      try{
    #endif
    qutet = quafun_xi(msh,asdmet,asdmsh,ent2poi[ientt],bary,NULL);
    #ifndef NDEBUG
      }catch(const MetrisExcept &e){
        printf("## metqua ent2pol \n");
        intAr1(getnnode(tdim,ideg), ent2poi[ientt]).print();
        throw(e);
      }
    #endif
    qutet = pow(abs(qutet - difto),pnorm);
  }

  if(do_nordev){
    METRIS_ASSERT(msh.param->qua_surf_wt_quality >= 0);
    METRIS_ASSERT(msh.param->qua_surf_wt_normal  >= 0);
    qutet = msh.param->qua_surf_wt_quality*qutet 
          + msh.param->qua_surf_wt_normal*pow(nordev, pnorm); // for homogeneity
  }

  return qutet;
}

#define EXPAND_TEMPLATE(r,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ),\
                              BOOST_PP_SEQ_ELEM(2, SEQ))
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical) 
#define QUAFUN_SEQ (QuaFun::Distortion)(QuaFun::Unit)


#define INSTANTIATE(MFT_VAL,QUAFUN,FTYPE)\
template FTYPE metqua< MFT_VAL , 2, 2, QUAFUN, FTYPE>\
            (Mesh<MFT_VAL> &msh, AsDeg, AsDeg, \
             int ielem, FTYPE difto);\
template FTYPE metqua< MFT_VAL , 3, 2, QUAFUN, FTYPE>\
            (Mesh<MFT_VAL> &msh, AsDeg, AsDeg, \
             int ielem, FTYPE difto);\
template FTYPE metqua< MFT_VAL , 3, 3, QUAFUN, FTYPE>\
            (Mesh<MFT_VAL> &msh, AsDeg, AsDeg, \
             int ielem, FTYPE difto);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,(MFT_SEQ)(QUAFUN_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE



#undef EXPAND_TEMPLATE
#undef MFT_SEQ // note these two could go into headers 

} // End namespace
