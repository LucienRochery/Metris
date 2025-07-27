//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_metqua.hxx"

#include "../ho_constants.hxx"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "../aux_topo.hxx"
#include "../low_geo/misc.hxx"
#include "../low_geo/normal.hxx"
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

  ftype qutet = 0; 
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
  const int ideg_eff = asdmsh == AsDeg::P1 ? 1 : ideg;
  const int nnode = getnnode(tdim, ideg_eff);
  constexpr int nnmet = (gdim*(gdim+1))/2;

  // Accumulate normal error at the nodes (depending on asdmsh)
  if(do_nordev){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      nordev = getnordev<ideg>(msh, ientt);
    }}CT_FOR1(ideg);
  }

  if(ideg_eff > 1){
    qutet = 0.0;
    nordev= 0.0;

    //const int idegj = SMOO_DEGJ(ideg);
    //const int nnodj = tdim == 2 ? getnnod2(idegj) : getnnod3(idegj);


    METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);

    //double metl[nnmet];

    for(int iquad = 0; iquad < nnode; iquad++){
      for(int ii = 0; ii < tdim + 1; ii++){
        bary[ii] = ordelt[ideg][iquad][ii]/((double) (ideg));
      }
      const int ipoin = ent2poi(ientt, iquad);
      ftype qua0 = quafun_xi(msh,asdmet,asdmsh,ent2poi[ientt],bary,msh.met[ipoin]);
      // About 6x slower if MetSpace::Log: leave an assert here. 
      //if(msh.met.getSpace() == MetSpace::Exp){
      //  qua0 = quafun_xi(msh,asdmet,asdmsh,ent2poi[ientt],bary,msh.met[ipoin]);
      //}else{
      //  for(int ii = 0; ii < nnmet; ii++) metl[ii] = msh.met(ipoin,ii);
      //  getexpmet_inp<gdim>(metl);
      //  qua0 = quafun_xi(msh,asdmet,asdmsh,ent2poi[ientt],bary,metl);
      //}

      // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
      // time to run pow() here even if pnorm = 2 or 1.
      qua0 = abs(qua0 - difto);
      if(pnorm == 2){
        qua0 *= qua0;
      }else if(pnorm > 2){
        qua0 = pow(qua0, pnorm);
      }
      qutet += qua0 / nnode;
      
      //qutet += pow(abs(qua0 - difto),pnorm) / nnode;

    }

  }else{

    #if 0
      // Value doesn't matter. (P1)
      for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim  + 1);
      METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
      for(int iquad = 0; iquad < tdim + 1; iquad++){
        int ipoin = ent2poi(ientt, iquad);
        ftype qua0 = quafun_xi(msh,asdmet,asdmsh,ent2poi[ientt],bary,msh.met[ipoin]);
        // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
        // time to run pow() here even if pnorm = 2 or 1.
        qua0 = abs(qua0 - difto);
        if(pnorm == 2){
          qua0 *= qua0;
        }else if(pnorm > 2){
          qua0 = pow(qua0, pnorm);
        }
        qutet += qua0 / (tdim + 1);
      }
    
    #elif 0
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
      // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
      // time to run pow() here even if pnorm = 2 or 1.
      qutet = abs(qutet - difto);
      if(pnorm == 2){
        qutet *= qutet;
      }else if(pnorm > 2){
        qutet = pow(qutet, pnorm);
      }

    #else
      METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
      double met[nnmet];
      for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim  + 1);
      for(int jj = 0; jj < nnmet; jj++) met[jj] = 0;
      for(int ii = 0; ii < tdim + 1; ii++){
        int ipoin = ent2poi(ientt, ii);
        for(int jj = 0; jj < nnmet; jj++){
          met[jj] += msh.met(ipoin,jj) / (tdim + 1);
        }
      }
      qutet = quafun_xi(msh,AsDeg::P1,asdmet,ent2poi[ientt],bary,met);
      // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
      // time to run pow() here even if pnorm = 2 or 1.
      qutet = abs(qutet - difto);
      if(pnorm == 2){
        qutet *= qutet;
      }else if(pnorm > 2){
        qutet = pow(qutet, pnorm);
      }
    #endif

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
