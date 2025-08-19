//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

/*
Routine for "direct" smoothing as P1. From each (facet, metric) pair, generate remaining vertex to be unit. Then average over ball.
Simplest possible approach.
*/

#include "low_smoolen.hxx"


#include "../Mesh/Mesh.hxx" 
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../smoothing/msh_smooball.hxx"
#include "../smoothing/low_smooballdiff.hxx"

#include "../low_geo/lenedg.hxx"

#include "../aux_topo.hxx"
#include "../utils/aux_timer.hxx"
#include "../low_topo.hxx"
#include "../utils/mprintf.hxx"
#include "../quality/low_metqua.hxx"
#include "../io_libmeshb.hxx"
#include "../msh_checktopo.hxx"

#include "libs/lplib3/lplib3.h"



namespace Metris{


template<class MFT>
double smoothMeshLength(Mesh<MFT> &msh, int tdim, int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);

  if(tdim == 1){
    CPRINTF1(" # smoothMeshLength tdim = 1 not implemented as it requires nordev to ensure validity\n");
    return 0;
  }

  //printf("## DEBUG SET MAX prints\n");
  //msh.param->iverb = 3;
  //msh.param->ivdepth = 15;
  //wait();
  

  const int nedgl = (tdim*(tdim+1))/2;
  const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                             tdim == 2 ? lnoed2[0] : lnoed3[0]);

  const int nentt = msh.nentt(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim); 

  CPRINTF1("-- START smoothMeshLength tdim = {}\n",tdim);



  const int miter = msh.param->opt_smoo_niter;
  //const double maxwt = 20.0;
  //const double qrthr = 2.0;
  const double tolavg = msh.param->opt_smoo_tol;
  const double tolmax = msh.param->opt_smoo_tol;

  dblAr1 work;
  if(msh.param->iflag2 != 0){
    MPRINTF("\n\n##WARNING EXPERIMENTAL SMOOTHING FUNCTION 2\n");
  }

  // 1 -> no maximum quality increase allowed 
  //const double maxinc_worst = 1.00;

  const int nnmet = (msh.idim*(msh.idim+1))/2;

  METRIS_ENFORCE(msh.param->opt_power < 0); // Otherwise rework the mins / maxs


  const int mball = 100;
  intAr1 lball(mball);

  msh.tag[ithrd1]++;

  const int medge = tdim == 2 ? msh.nface : msh.nelem;
  HshTab_I2I ledge;
  ledge.reserve(medge);

  double t0s = get_wall_time();
  int nedge_tot = 0;
  for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
    INCVDEPTH(msh.param);
    if(isdeadent(ientt,ent2poi)) continue;
    for(int ied = 0; ied < nedgl; ied++){
      INCVDEPTH(msh.param);

      // Check edge already seen
      int ip1 = ent2poi(ientt, lnoed(ied,0));
      int ip2 = ent2poi(ientt, lnoed(ied,1));
      auto key = stup2(ip1,ip2);
      if(ledge.find(key) != ledge.end()) continue;

      nedge_tot++;
      
      ledge[key] = ientt;
    }// for ied
  }// for ientt
  double t1s = get_wall_time();
  CPRINTF1(" - init time {:.2e}s nlong = {}\n",t1s-t0s,(int)ledge.size());


  intWrkAr1 lpoin = msh.get_iwork(100);
  double noper = 0;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);

    double qmin = 1.0e30, qmax = -1.0e30, qavg = 0.0;
    int imax = -1;
    int navg = 0;
    for(auto iedge : ledge){
      int ip1 = std::get<0>(iedge.first);
      int ip2 = std::get<1>(iedge.first);

      double sz[2];
      int edg2pol[2] = {ip1, ip2};
      double len = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,edg2pol,sz)
                                 : getlenedg_geosz<MFT,3,1>(msh,edg2pol,sz);

      double quaed = len < 1.0 ? 1.0 - len 
                               : 1.0 - 1.0 / len;

      qmin = MIN(qmin, quaed);
      qmax = MAX(qmax, quaed);
      qavg += quaed;
    }

    qavg /= ledge.size();
    double t0 = get_wall_time();
    CPRINTF1(" - smoo iter {:3} init {:10.6e} < q < {:10.6e} (at {}), avg = {:10.6e}\n",
             niter,qmin,qmax,imax,qavg,msh.param->opt_pnorm);
    //if(iverb >= 2 && qmax >= 1e10){
    //  printf("## HIGH QMAX mshdeg = {} \n",msh.curdeg);
    //  std::string fname = "qmax"+std::to_string(imax);
    //  writeMesh(fname,msh);
    //  //wait();
    //} 

    int nsucc = 0;
    int nmov  = 0;

    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      if(msh.poi2ent(ipoin, 0) < 0) continue;
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
      INCVDEPTH(msh.param);

      if(msh.getpoitdim(ipoin) != tdim){
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];
        continue;
      }

      double qnrm0, qnrm1;
      int ierro = movePointBallLen(msh, ipoin, 10, &qnrm0, &qnrm1, ithrd2);

        //PRINTF("## DEBUG PAUSE HERE ierro = {} \n",ierro);
        //writeMesh("debug_movepointball",msh);
        //wait();

      if(ierro < 0){
        nsucc++;
        CPRINTF1(" - success smoothing ipoin {} q max {:10.6e} -> {:10.6e}\n",ipoin,qnrm0,qnrm1);


        bool imov = false;
        // qnrm1 should be < qnrm0 for there to be progress 
        if(qnrm0 - qnrm1 > tolavg) imov = true;
        //// idem qmax 
        //if(qmax0 - qmax1 > tolmax) imov = true;
        if(imov){
          nmov ++;
          // We have no choice but to recompute as movePointBallLen
          // uses lowest-dim, but we need to tag all.
          ierro = poi2poi(msh, ipoin, -1, lpoin.get_array(), NULL, ithrd2);
          METRIS_ENFORCE(ierro == 0);
          for(int ipoi2 : lpoin.get_array()){
            msh.poi2tag(ithrd1,ipoi2) = msh.tag[ithrd1] - 1; // reactivate
          }
        }else{
          msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1]; // deactivate
        }
      }else{
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1]; // deactivate
      }
    }

    double t1 = get_wall_time();
    CPRINTF1(" - Iteration end time = {:.2e}s nsuccess = {} nmov = {} \n",
                          t1-t0,nsucc,nmov);
    noper += nmov;
    if(nmov == 0) break;
  } // end for niter

  return noper / (double) nentt;
}


template double smoothMeshLength<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh,
                                 int tdim, int ithrd1, int ithrd2);
template double smoothMeshLength<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh,
                                 int tdim, int ithrd1, int ithrd2);


} // end namespace
