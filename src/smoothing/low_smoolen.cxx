//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "low_smoolen.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"
#include "../low_geo/lenedg.hxx"
#include "../low_geo/measure.hxx"
#include "../low_topo.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../Adaptation/Insertion/EdgeSeed.hxx"
#include "../io_libmeshb.hxx"

namespace Metris{

// TODO: if used, this routine's preprocessing can be moved out.
// Return 0 if nothing done, > 0 error, -1 moved.
template<class MFT>
int movePointCavLen(Mesh<MFT>& msh, const MshCavity &cav,
                   int miter, int ithrd1){

  GETVDEPTH(msh.param);

  const int tdim = msh.getpoitdim(cav.ipins);
  if(tdim == 0) return 0;

  const intAr1 &lcent = cav.lcent(tdim);
  intWrkAr1 lpoin = msh.get_iwork(lcent.get_n());
  poi2poi(msh, cav.ipins, tdim, lcent, lpoin.get_array(), ithrd1);

  // Now we have our list of cavity boundary points.
  // The smoothing is a simple weighted average
  return smoopoilen(msh, cav.ipins, lpoin.get_array(), NULL, NULL, NULL, miter);
}

template
int movePointCavLen<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
  const MshCavity &cav, int miter, int ithrd1);
template
int movePointCavLen<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
  const MshCavity &cav, int miter, int ithrd1);



// TODO: if used, this routine's preprocessing can be moved out.
// Return 0 if nothing done, > 0 error, -1 moved.
template<class MFT>
int movePointBallLen(Mesh<MFT>& msh, int ipmov,
                    int miter, double *qlen0, double *qlen1, int ithrd1){

  GETVDEPTH(msh.param);
  const int tdim = msh.getpoitdim(ipmov);
  if(tdim == 0) return 0;

  intWrkAr1 lbedg = msh.get_iwork(10);
  intWrkAr1 lbfac = msh.get_iwork(30);
  intWrkAr1 lbtet = msh.get_iwork(msh.get_tdim() >= 3 ? 100 : 0);

  int iopen;
  int ierro = ball(msh, ipmov, lbedg.get_array(), lbfac.get_array(), lbtet.get_array(),
                   &iopen, false, ithrd1);
  if(ierro != 0) return ierro;

  intWrkAr1 lpoin = msh.get_iwork(100);
  intAr1& lbent = tdim == 1 ? lbedg.get_array()
                : tdim == 2 ? lbfac.get_array()
                            : lbtet.get_array();
  lpoin.set_n(0);
  poi2poi(msh, ipmov, tdim, lbent, lpoin.get_array(), ithrd1);

  // The smoothing is a simple weighted average
  return smoopoilen(msh, ipmov, lpoin.get_array(), 
                    &(lbedg.get_array()), &(lbfac.get_array()), &(lbtet.get_array()),
                    miter, qlen0, qlen1);
}

template
int movePointBallLen<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
   int ipmov, int miter, double *qlen0, double *qlen1, int ithrd1);
template
int movePointBallLen<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
   int ipmov, int miter, double *qlen0, double *qlen1, int ithrd1);



// Return 0 if nothing done, > 0 error, -1 moved.
// Optionally provide list of elements lbedg, lbfac, lbtet to control validity of. 
// This is not necessary when moving within a cavity, for example,
// as we follow with cavity correction. 
template<class MFT>
int smoopoilen(Mesh<MFT>& msh, int ipmov, 
               const intAr1 &lpoin, 
               const intAr1* lbedg, const intAr1* lbfac, const intAr1* lbtet,
               int miter,
               double *qlen0, double *qlen1){

  const double damp = 0.1;
  const double nordev = -1; // improve in the future by getting initial ball nordev?

  GETVDEPTH(msh.param);

  const int tdimp = msh.getpoitdim(ipmov);
  if(tdimp == 0) return 0;

  const int iseed = msh.poi2ent(ipmov,0);
  const int iref = msh.ent2ref(tdimp)[iseed];

  CPRINTF1("-- START smoopoilen ipmov {} lpoin.n = {}\n",ipmov,lpoin.get_n());


  ego obj = NULL;
  int ibmov = -1;
  bool use_uv = false;
  if(msh.isboundary_tdim(tdimp)){
    ibmov = msh.poi2ebp(ipmov,tdimp,iseed,iref);
    METRIS_ASSERT(ibmov >= 0);
    if(msh.CAD()){
      obj = tdimp == 1 ? msh.CAD.cad2edg[iref]
                       : msh.CAD.cad2fac[iref] ;
      use_uv = true;
    }
  }
  int ndof = use_uv ? tdimp : msh.idim;
  int ipdof = use_uv ? ibmov : ipmov;
  dblAr2 &mshdof = use_uv ? msh.bpo2rbi : msh.coord;

  

  const int npoil = lpoin.get_n();
  dblWrkAr1 rpoin = msh.get_rwork(npoil);
  int edg2pol[2] = {ipmov, -1};
  double sz[2], result[18];
  double coord_opt[3], uv_opt[2], met_opt[6];
  double coord0[3], uv0[2], met0[6];
  const int nnmet = (msh.idim*(msh.idim + 1)) / 2;
  for(int ii = 0; ii < msh.idim; ii++) coord0[ii] = coord_opt[ii] = msh.coord(ipmov,ii);
  if(use_uv) 
    for(int ii = 0; ii < 2; ii++) uv0[ii] = uv_opt[ii] = msh.bpo2rbi(ibmov,ii);
  for(int ii = 0; ii < nnmet   ; ii++) met0[ii] = met_opt[ii] = msh.met(ipmov,ii);
  
  // Initialize optimal cost at initial cost, check no regression
  double wt_ini = 0, avglen = 0, minlen = 1.0e30, maxlen = -1;
  double avglen_opt, minlen_opt, maxlen_opt;
  double lenqua = -1, lenqua_opt = -1;
  for(int ii = 0; ii < npoil; ii++){
    INCVDEPTH(msh.param);
    int ipoin = lpoin[ii];
    edg2pol[1] = ipoin;
    double len = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,edg2pol,sz) 
                               : getlenedg_geosz<MFT,3,1>(msh,edg2pol,sz);
    avglen += len;
    minlen = MIN(minlen, len);
    maxlen = MAX(maxlen, len);
    // Lengths above 1 are weighted more strongly, such that ipins is drawn towards those points
    // Lengths below 1 are weighted weakly, to repel ipins
    double fact = len < 1 ? pow(len, 4) : pow(len, 4);
    wt_ini += fact;
    double quaed = len < 1.0 ? 1.0 - len 
                             : 1.0 - 1.0 / len;
    lenqua = MAX(quaed, lenqua);
    CPRINTF1(" - initial ipoin {} len {} weight {}\n",ipoin,len,fact);
  }
  avglen /= npoil;
  avglen_opt = avglen;
  minlen_opt = minlen;
  maxlen_opt = maxlen;
  CPRINTF1(" - initial avglen = {} min {} max {} cost = {} qua {}\n",avglen,minlen,maxlen,wt_ini,lenqua);
  lenqua_opt = lenqua;
  if(qlen0 != NULL) *qlen0 = lenqua;

  bool one_update = false;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);
    double wttot = 0;
    avglen = 0;
    minlen = 1.0e30;
    maxlen = -1;
    lenqua = -1;
    for(int ii = 0; ii < npoil; ii++){
      INCVDEPTH(msh.param);
      int ipoin = lpoin[ii];
      edg2pol[1] = ipoin;
      double len = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,edg2pol,sz) 
                                 : getlenedg_geosz<MFT,3,1>(msh,edg2pol,sz);
      avglen += len;
      minlen = MIN(minlen, len);
      maxlen = MAX(maxlen, len);
      // Lengths above 1 are weighted more strongly, such that ipins is drawn towards those points
      // Lengths below 1 are weighted weakly, to repel ipins
      rpoin[ii] = len < 1 ? pow(len, 8) : pow(len, 4);
      wttot += rpoin[ii];
      double quaed = len < 1.0 ? 1.0 - len 
                               : 1.0 - 1.0 / len;
      lenqua = MAX(quaed, lenqua);
      CPRINTF1(" - iter {} ipoin {} len {} weight {}\n",niter,ipoin,len,rpoin[ii]);
    }
    avglen /= npoil;

    bool ivalid = true;
    CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
      if(!ivalid) CT_CONTINUE(gdim);
      CT_FOR0_INC(2,gdim,tdim){if(tdim <= msh.get_tdim()){
        if(!ivalid) CT_CONTINUE(tdim);
        const intAr1* lbent = tdim == 1 ? lbedg
                            : tdim == 2 ? lbfac
                                        : lbtet;
        if(lbent == NULL){
          CPRINTF1(" # no lbent for tdim {}\n",tdim);
          CT_CONTINUE(tdim);
        }
        for(int ientt : *lbent){
          bool ieval;
          ieval = isvalideltP1<gdim,tdim>(msh, ientt, NULL, NULL, nordev);
          CPRINTF2(" - tdim {} ball elt {} valid = {}\n", tdim, ientt, ieval);
          if(ieval) continue;
          CPRINTF1(" - dim {} element {} became invalid -> reject iteration\n",tdim,ientt);
          ivalid = false;
          break;
        }
      }}CT_FOR1(tdim);
    }}CT_FOR1(gdim);

    CPRINTF1(" - valid {} update new lenqua {:.2e} opt = {:.2e}\n",ivalid,lenqua,lenqua_opt);

    if(lenqua_opt > lenqua && ivalid){
      CPRINTF1(" - iter {} new optimum cost {} -> {}\n",niter, lenqua_opt, lenqua);
      lenqua_opt = lenqua;
      minlen_opt = minlen;
      maxlen_opt = maxlen;
      avglen_opt = avglen;
      for(int ii = 0; ii < msh.idim; ii++) coord_opt[ii] = msh.coord(ipmov,ii);
      if(use_uv) for(int ii = 0; ii < 2       ; ii++) uv_opt[ii]    = msh.bpo2rbi(ibmov,ii);
      for(int ii = 0; ii < nnmet   ; ii++) met_opt[ii]   = msh.met(ipmov,ii);
      one_update = true;
    }

    if(niter == miter - 1) break;

    double pdofn[3];
    for(int ii = 0; ii < ndof; ii++) pdofn[ii] = 0;
    for(int ii = 0; ii < npoil; ii++){
      int ipoin = lpoin[ii];
      int idof = ipoin;
      if(use_uv){
        idof = msh.poi2ebp(ipoin,tdimp,iseed,iref);
        #ifndef NDEBUG
        if(idof < 0){
          PRINTF("# idof = {} with ipoin {} tdimp {} iseed {} iref {}\n",
                 idof, ipoin, tdimp, iseed, iref);
          writeMesh("debug_seed",msh);
        }
        #endif
      }
      for(int jj = 0; jj < ndof; jj++)
        pdofn[jj] += mshdof(idof,jj) * rpoin[ii] / wttot;
    }
    CPRINTF1(" - iter {} avglen = {} min {} max {} cost {} lenqua {} new dof {}\n",niter,
             avglen,minlen,maxlen,wttot,lenqua,dblAr1(ndof,pdofn));
    for(int ii = 0; ii < ndof; ii++) mshdof(ipdof,ii) = (1-damp)*mshdof(ipdof,ii)
                                                         + damp *pdofn[ii];
    if(use_uv){
      // In this case, mshdof was bpo2rbi and we need to call EG_evaluate
      int ierro = EG_evaluate(obj, mshdof[ipdof], result);
      if(ierro != 0){
        CPRINTF1(" # EG_evaluate error {}\n",ierro);
        return 3;
      }
      for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipmov,ii) = result[ii];
    }

    int ierro = msh.interpMetBack(ipmov);
    if(ierro != 0){
      CPRINTF1(" # interpMetBack failed ierro = {} \n",ierro);
      break;
    }
  }

  if(!one_update){
    CPRINTF1("-- END movePointCavLen: no update done\n");
    if(use_uv){
      for(int ii = 0; ii < 2; ii++) msh.bpo2rbi(ibmov,ii) = uv0[ii];
    }
    for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipmov,ii) = coord0[ii];
    for(int ii = 0; ii < nnmet; ii++)    msh.met(ipmov,ii)   = met0[ii];
    return 0;
  }

  if(use_uv){
    for(int ii = 0; ii < 2; ii++) msh.bpo2rbi(ibmov,ii) = uv_opt[ii];
  }
  for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipmov,ii) = coord_opt[ii];
  for(int ii = 0; ii < nnmet; ii++)    msh.met(ipmov,ii)   = met_opt[ii];

  if(qlen1 != NULL) *qlen1 = lenqua_opt;

  if(qlen1 != NULL) METRIS_ASSERT(*qlen1 < *qlen0);

  CPRINTF1("-- END movePointCavLen: ipoin {} lenqua {:.2e} avglen {:.2f} minlen {:.2f} maxlen {:.2f}",
          ipmov, lenqua_opt, avglen_opt, minlen_opt, maxlen_opt);
  if(qlen0 != NULL && qlen1 != NULL && DOPRINTS1()) PRINTF(" qlen0 = {:.2e} qlen1 = {:.2e}\n", *qlen0, *qlen1);
  else if(DOPRINTS1()) PRINTF("\n");

  return -1;
}

template
int smoopoilen<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
              int ipmov, const intAr1 &lpoin, 
              const intAr1* lbedg, const intAr1* lbfac, const intAr1* lbtet, 
              int miter, 
              double *qlen0, double *qlen1);
template
int smoopoilen<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
              int ipmov, const intAr1 &lpoin,
              const intAr1* lbedg, const intAr1* lbfac, const intAr1* lbtet, 
              int miter, 
              double *qlen0, double *qlen1);
} // namespace Metris