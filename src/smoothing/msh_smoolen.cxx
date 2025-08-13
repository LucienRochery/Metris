//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

/*
Routine for "direct" smoothing as P1. From each (facet, metric) pair, generate remaining vertex to be unit. Then average over ball.
Simplest possible approach.
*/

#include "../Mesh/Mesh.hxx" 
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../smoothing/msh_smooball.hxx"
#include "../smoothing/low_smooballdiff.hxx"

#include "../aux_topo.hxx"
#include "../utils/aux_timer.hxx"
#include "../low_topo.hxx"
#include "../utils/mprintf.hxx"
#include "../quality/low_metqua.hxx"
#include "../io_libmeshb.hxx"

#include "libs/lplib3/lplib3.h"



namespace Metris{


// idim: gdim = tdim
template<class MFT>
double smoothMeshLength(Mesh<MFT> &msh,int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);

  constexpr int tdim = idim;
  //constexpr int gdim = idim;
  if(ideg > idim + 1){
    PRINTF("## SMOOTHING DISABLED FOR DEGREE {} and dim {} \n",ideg,idim);
    return -1.0;
  }


  int nentt = msh.nentt(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim); 
  const intAr2 &ent2ent = msh.ent2ent(tdim); 



  // Eventually move all constants to MetrisParameters 
  // L2 conformity error from 0 to 1 
  const double difto = 1.0;
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

  constexpr int nnmet = (idim*(idim+1))/2;

  METRIS_ENFORCE(msh.param->opt_power < 0); // Otherwise rework the mins / maxs
  // Otherwise not only edge nodes 
  METRIS_ENFORCE(ideg <= tdim + 1); 


  const int mball = 100;
  intAr1 lball(mball);

  msh.tag[ithrd1]++;

  auto quafun = get_quafun<MFT,idim,idim>(iquaf);

  double noper = 0;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);

    double qmin = 1.0e30,qmax = -1.0e30, qnrm = 0.0;
    int imax = -1;
    int navg = 0;
    for(int ientt = 0; ientt < nentt; ientt++){

      if(isdeadent(ientt,ent2poi)) continue;

      navg++;

      double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ientt,difto);

      qnrm += quael;
      qmin = MIN(qmin,quael);
      if(qmax < quael){
        imax = ientt;
        qmax = quael;
      }
    }

    qnrm /= navg;
    double t0 = get_wall_time();
    CPRINTF1(" - smoo iter {:3} init {:10.6e} < q < {:10.6e} (at {}), avg = {:10.6e} " 
                   "(p = {})\n",niter,qmin,qmax,imax,qnrm,msh.param->opt_pnorm);
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

      int ib = msh.poi2bpo[ipoin];
      if(ib >= 0) continue;

      int ientt = getpoient(msh, ipoin, tdim);  
      int iver = tdim == 2 ? msh.template getverfac<ideg>(ientt, ipoin)
                           : msh.template getvertet<ideg>(ientt, ipoin);

      METRIS_ASSERT(iver >= 0);


      CPRINTF1(" - smoo pt {} seed elt {} \n", ipoin,ientt);
      int ierro = 0;

      if(iver < tdim+1){
        int iopen;
        // Vertex case 
        if constexpr (idim == 2){
          intAr1 dum; 
          bool imani = false;
          ierro = ball2(msh,ipoin,ientt,lball,dum,&iopen,&imani,ithrd2);
          METRIS_ASSERT(imani == true);
        }else{
          ierro = ball3(msh,ipoin,ientt,lball,&iopen,ithrd2);
        }
        METRIS_ASSERT(iopen == 0);
      }else{
        // HO node
        int nppe = getnnod1(ideg) - 2;
        int ied = (iver - (tdim + 1)) / nppe;
        if constexpr (tdim == 2){
          lball.set_n(0);
          lball.stack(ientt);
          METRIS_ASSERT(ied < 4);

          int ifnei = ent2ent(ientt,ied);
          // nm not impl yet and bdry should never happen
          METRIS_ASSERT(ifnei >= 0);

          lball.stack(ifnei);
        }else{
          int ip1 = msh.tet2poi(ientt, lnoed3[ied][0]);
          int ip2 = msh.tet2poi(ientt, lnoed3[ied][1]);
          static intAr1 dum;
          int iopen;
          shell3(msh, ip1, ip2, ientt, lball, dum, &iopen);
          METRIS_ASSERT(iopen < 0);
        }

      }
      METRIS_ASSERT(ierro == 0); 

      double coor0[idim];
      double met0[nnmet];
      for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
      for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipoin,ii);
      double qnrm0, qmax0, qnrm1, qmax1;
      ierro = movePointCavLen(msh, cav, tdim, ientt, miter, ithrd1);
      try{
        //ierro = smooballdirect<MFT,idim,ideg>(msh,ipoin,lball,qball,
        //                       &qnrm0,&qmax0,&qnrm1,&qmax1,
        //                       qpower,qpnorm,difto,maxwt,inorm,iverb,ithrd2);
        if(msh.param->iflag2 == 0){
          ierro = smooballdiff<MFT,idim,ideg>(msh,ipoin,lball,
                                 &qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
        }else{
          ierro = smooballdiff_luksan<MFT,idim,ideg>(msh,ipoin,lball,
                                     &qnrm0,&qmax0,&qnrm1,&qmax1,work,iquaf);
        }
        if(qmax1 > qmax){
          CPRINTF1(" - reject move, worst above last worst {:15.7e} > {:15.7e}\n", 
                   qmax1, qmax);
          for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
          for(int ii = 0; ii < nnmet;ii++) msh.met(ipoin,ii)   =  met0[ii];
          ierro = 1;
        }
      }catch(const MetrisExcept &e){
        PRINTF("## FAILED  smooballdirect\n");
        writeMesh("smooth_error.meshb",msh);
        throw(e);
      }
      if(ierro == 0){
        nsucc++;
        CPRINTF1(" - success smoothing {} q avg {:10.6e} -> {:10.6e} " 
                 "max {:10.6e} -> {:10.6e}\n",ipoin,qnrm0,qnrm1,qmax0,qmax1);

        bool imov = false;
        // qnrm1 should be < qnrm0 for there to be progress 
        if(qnrm0 - qnrm1 > tolavg) imov = true;
        // idem qmax 
        if(qmax0 - qmax1 > tolmax) imov = true;
        if(imov){
          nmov ++;
          for(int iele2 : lball){
            for(int ii = 0; ii < idim+1; ii++){
              int ipoi2 = ent2poi(iele2,ii);
              msh.poi2tag(ithrd1,ipoi2) = msh.tag[ithrd1] - 1; // reactivate
            }
          }
        }else{
          msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1]; // deactivate
        }
      }else{
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1]; // deactivate
      }

    } // for ipoin // for ientt

    double t1 = get_wall_time();
    CPRINTF1(" - Iteration end time = {:.2e}s nsuccess = {} nmov = {} \n",
                          t1-t0,nsucc,nmov);
    noper += nmov;
    if(nmov == 0) break;
  } // end for niter

  return noper / (double) nentt;
}


template double smoothMeshLength<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh,
                                        int ithrd1, int ithrd2);
template double smoothMeshLength<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh,
                                        int ithrd1, int ithrd2);


} // end namespace
