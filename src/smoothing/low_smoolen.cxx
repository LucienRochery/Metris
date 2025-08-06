//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "low_smoolen.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../utils/mprintf.hxx"
#include "../low_geo/lenedg.hxx"

namespace Metris{

template<class MFT>
int smoopoilen(Mesh<MFT>& msh, int ipmov, const intAr1 &lpoin, int miter, int tdimp, int iseed){
  GETVDEPTH(msh.param);

  if(tdimp < 0) tdimp = msh.getpoitdim(ipmov);
  if(iseed < 0) iseed = msh.poi2ent(ipmov,0);

  if(tdimp != msh.get_tdim()){
    CPRINTF1(" # smoopoilen only implemented for interior points\n");
    return 3;
  }

  int iref = msh.ent2ref(tdimp)[iseed];
  
  const double damp = 0.1;

  
  const int npoil = lpoin.get_n();
  dblWrkAr1 rpoin = msh.get_rwork(npoil);
  int edg2pol[2] = {ipmov, -1};
  double sz[2];
  double coord_opt[3], met_opt[6];
  const int nnmet = (msh.idim*(msh.idim + 1)) / 2;
  for(int ii = 0; ii < msh.idim; ii++) coord_opt[ii] = msh.coord(ipmov,ii);
  for(int ii = 0; ii < nnmet; ii++) met_opt[ii] = msh.met(ipmov,ii);
  
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
      double fact = len < 1 ? pow(len, 8) : pow(len, 4);
      rpoin[ii] = fact;
      wttot += rpoin[ii];
      double quaed = len < 1.0 ? 1.0 - len 
                               : 1.0 - 1.0 / len;
      lenqua = MAX(quaed, lenqua);
      CPRINTF1(" - iter {} ipoin {} len {} weight {}\n",niter,ipoin,len,rpoin[ii]);
    }
    avglen /= npoil;

    if(lenqua_opt > lenqua){
      CPRINTF1(" - iter {} new optimum cost {} -> {}\n",niter, lenqua_opt, lenqua);
      lenqua_opt = lenqua;
      minlen_opt = minlen;
      maxlen_opt = maxlen;
      avglen_opt = avglen;
      for(int ii = 0; ii < msh.idim; ii++) coord_opt[ii] = msh.coord(ipmov,ii);
      for(int ii = 0; ii < nnmet; ii++) met_opt[ii] = msh.met(ipmov,ii);
    }

    double coorn[3];
    for(int ii = 0; ii < msh.idim; ii++) coorn[ii] = 0;
    for(int ii = 0; ii < npoil; ii++){
      int ipoin = lpoin[ii];
      for(int jj = 0; jj < msh.idim; jj++){
        coorn[jj] += msh.coord(ipoin,jj) * rpoin[ii] / wttot;
      }
    }
    CPRINTF1(" - iter {} avglen = {} min {} max {} cost {} lenqua {}\n",niter,
             avglen,minlen,maxlen,wttot,lenqua);
    for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipmov,ii) = (1-damp)*msh.coord(ipmov,ii)
                                                                 + damp * coorn[ii];

    int ierro = msh.interpMetBack(ipmov,tdimp,iseed,iref,NULL);
    if(ierro != 0){
      CPRINTF1(" # interpMetBack failed ierro = {} \n",ierro);
      return 2;
    }
  }


  for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipmov,ii) = coord_opt[ii];
  for(int ii = 0; ii < nnmet; ii++)    msh.met(ipmov,ii) = met_opt[ii];

  CPRINTF1("-- END movePointCavLen: lenqua {} avglen {} minlen {} maxlen {}\n",
          lenqua_opt, avglen_opt, minlen_opt, maxlen_opt);

  return 0;
}

template
int smoopoilen<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, int ipmov, const intAr1 &lpoin, int miter, int tdimp, int iseed);
template
int smoopoilen<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, int ipmov, const intAr1 &lpoin, int miter, int tdimp, int iseed);
} // namespace Metris