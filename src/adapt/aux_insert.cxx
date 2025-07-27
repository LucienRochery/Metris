//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "aux_insert.hxx"
#include "low_insert.hxx" // for error codes

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../utils/mprintf.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../msh_structs.hxx"
#include "../low_topo.hxx"
#include "../low_geo/normal.hxx"
#include "../low_geo/measure.hxx"
#include "../io_libmeshb.hxx"
#include "../linalg/det.hxx"
#include "../low_geo/lenedg.hxx"

#include "../msh_checktopo.hxx"


namespace Metris{


template<class MFT>
int aux_findCloseConstrained(Mesh<MFT>& msh, MshCavity &cav, 
                             int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);

  int nced0 = cav.lcedg.get_n();
  int ncfa0 = cav.lcfac.get_n();
  int ncte0 = cav.lctet.get_n();

  const int tdimm = msh.get_tdim();
  const intAr1 &lcent = cav.lcent(tdimm);
  const int ncen0 = tdimm == 1 ? nced0 : 
                    tdimm == 2 ? ncfa0 : ncte0;
  int edg2pol[2] = {cav.ipins, -1};
  double sz[2];

  const int pdimins = msh.getpoitdim(cav.ipins);

  int iice0 = 0, iice1 = ncen0;
  // tag points whose dist to ipins and balls have been computed:
  msh.tag[ithrd1]++; 
  for(int niter = 0; niter < 2; niter++){
    // Loop over current neighbourhood and check if points close.
    for(int iicen = iice0; iicen < iice1; iicen++){
      int icent = lcent[iicen];
      for(int iver = 0; iver < tdimm + 1; iver++){
        int ipoin = msh.ent2poi(tdimm)(icent,iver);
        if(!msh.poicstr[ipoin]) continue;
        if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
        // Only consider points of dimension <= of ipins
        if(msh.getpoitdim(ipoin) > pdimins) continue;

        // Points whose length has been computed are tagged itag.
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];
        edg2pol[1] = ipoin;
        double len = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,edg2pol,sz)
                                   : getlenedg_geosz<MFT,3,1>(msh,edg2pol,sz);
        if(len >= 1.0/sqrt(2)) continue;
        CPRINTF1(" # At least one close constrained point, stop here\n");
        msh.tag[ithrd1]++;
        return 1;
      }
    }
    if(niter == 1) break;

    // Next extend by point balls and restart
    for(int iicen = 0; iicen < ncen0; iicen++){
      int icent = lcent[iicen];
      for(int iver = 0; iver < tdimm + 1; iver++){
        int ipoin = msh.ent2poi(tdimm)(icent,iver);
        if(msh.poi2tag(ithrd1,ipoin) > msh.tag[ithrd1]) continue;
        // Points whose ball has been computed are tagged itag + 1.
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1] + 1;
        // Append ball
        int iopen;
        ball(msh, ipoin, cav.lcedg, cav.lcfac, cav.lctet, &iopen, true, ithrd2);
      }
    }

    iice0 = iice1;
    iice1 = lcent.get_n();
  }

  msh.tag[ithrd1]++;
  cav.lcedg.set_n(nced0);
  cav.lcfac.set_n(ncfa0);
  cav.lctet.set_n(ncte0);
  return 0;
}
template int aux_findCloseConstrained<MetricFieldAnalytical>(
  Mesh<MetricFieldAnalytical>& msh, MshCavity &cav, int ithrd1, int ithrd2);
template int aux_findCloseConstrained<MetricFieldFE       >(
  Mesh<MetricFieldFE       >& msh, MshCavity &cav, int ithrd1, int ithrd2);


template<class MFT>
int aux_bisecPointLen(Mesh<MFT> &msh, 
                      int tdim, int ientt, int iedl,
                      int ibins, int tdimp, int iseed, int iref, 
                      bool icollapse, 
                      const MshCavity &cav){
  GETVDEPTH(msh.param);
  const auto lnoed = tdim == 1 ? lnoed1 : 
                     tdim == 2 ? lnoed2 : lnoed3;
  const intAr2 &ent2poi = msh.ent2poi(tdim);

  ego obj = NULL;
  if(msh.CAD()){
    obj = tdimp == 1 && msh.isboundary_edges() ? msh.CAD.cad2edg[iref] :
          tdimp == 2 && msh.isboundary_faces() ? msh.CAD.cad2fac[iref] : NULL;
  }

  int ip1 = ent2poi(ientt,lnoed[iedl][0]);
  int ip2 = ent2poi(ientt,lnoed[iedl][1]);

  double bar1_min = 1.0e-6, bar1_max = 1 - 1.0e-6;
  double algnd[3];

  constexpr int mnode = getnnode(1,METRIS_MAX_DEG);
  const int nnode = getnnode(1,msh.curdeg);
  int edg2pol[mnode], edg2po2[2] = {cav.ipins, -1};
  int idx[4] = {0};
  int idx1[2];
  for(int inoed = 0; inoed <= msh.curdeg; inoed++){
    idx[lnoed[iedl][0]] = msh.curdeg - inoed;
    idx[lnoed[iedl][1]] = inoed;
    idx1[0] = msh.curdeg - inoed;
    idx1[1] = inoed;
    int inode_sup = mul2nod(tdim,idx);
    int inode_sub = mul2nod(1,idx1);
    edg2pol[inode_sub] = ent2poi(ientt,inode_sup);
  }
  CPRINTF2(" - edg2pol = ");
  if(DOPRINTS2()) intAr1(nnode,edg2pol).print();

  int ierro;
//restart_bisection:
  // Get bar1 s.t. new edges are not short. There can be other short edges, but 
  // not from splitting the parent edge. 
  bool fnd_len = false;
  double bar1_opt = -1, err_opt = 1.0e30, bar1;
  for(int ntry_len = 0; ntry_len < 10; ntry_len++){
    INCVDEPTH(msh.param);
    bar1 = (bar1_min + bar1_max) / 2;
    double bar2[2] = {bar1, 1 - bar1};

    // Evaluate ipins on CAD or element, also get algnd for interpMetBack 
    if(ibins >= 0 && msh.CAD()){
      int ib[2];
      // Correct ibs : attach to ref or edge/face as needed
      for(int ii = 0; ii < 2; ii++){
        ib[ii] = msh.poi2ebp(edg2pol[ii],tdimp,iseed,iref);
        METRIS_ASSERT(ib[ii] >= 0);
      }

      for(int ii = 0; ii < 2; ii++) msh.bpo2rbi(ibins,ii) = 
          bar1*msh.bpo2rbi[ib[0]][ii] + (1.0 - bar1)*msh.bpo2rbi[ib[1]][ii];

      CPRINTF1(" - boundary point new t/(u,v) = %f %f\n",
               msh.bpo2rbi(ibins,0),msh.bpo2rbi(ibins,1));

      double result[18];
      METRIS_ASSERT(obj != NULL);
      ierro = EG_evaluate(obj, msh.bpo2rbi[ibins], result);
      if(ierro != 0) return ierro;

      for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins,ii) = result[ii];

      if(tdimp == 1){
        for(int ii = 0; ii < msh.idim; ii++) algnd[ii] = result[3+ii];
      }else{
        vecprod(&result[3], &result[6], algnd);
      }
    }else if(ibins >= 0 && !msh.CAD()){ 
      METRIS_ASSERT(tdimp <= 2);
      // No reevaluation, but initialize algnd to edge tangent 
      CPRINTF1(" - discrete algnd initialization tdimp %d \n",tdimp);
      if(tdimp == 1){
        // To compute at higher degree, copy more vertices into ip
        MSH_DIM_DEG0(msh){
          eval1<gdim,ideg>(msh.coord, edg2pol,
                           msh.getBasis(), DifVar::Bary, DifVar::None,
                           bar2, msh.coord[cav.ipins], algnd, NULL);
        }MSH_DIM_DEG1();
      }else{
        MSH_DIM_DEG0(msh){
          eval1<gdim,ideg>(msh.coord, edg2pol,
                           msh.getBasis(), DifVar::None, DifVar::None,
                           bar2, msh.coord[cav.ipins], NULL, NULL);
        }MSH_DIM_DEG1();
        if(msh.idim == 3) getnorfacP1(ent2poi[ientt],msh.coord,algnd);
      }
    }else{
      MSH_DIM_DEG0(msh){
        eval1<gdim,ideg>(msh.coord, edg2pol,
                         msh.getBasis(), DifVar::None, DifVar::None,
                         bar2, msh.coord[cav.ipins], NULL, NULL);
      }MSH_DIM_DEG1();
    }

    ierro = msh.interpMetBack(cav.ipins, tdimp, iseed, iref, algnd);
    if(ierro != 0) return INS2D_ERR_INTERPMETBACK;

    //if(DOPRINTS3()){
    //  int ipnew = msh.newpoitopo(0);
    //  msh.newbpotopo(ipnew, 0, ipnew);
    //  const int nnmet = (msh.idim*(msh.idim+1))/2;
    //  for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipnew, ii) = msh.coord(cav.ipins, ii);
    //  for(int ii = 0; ii < nnmet; ii++) msh.met(ipnew, ii) = msh.met(cav.ipins, ii);
    //}


    double sz[2];
    edg2po2[1] = ip1;
    double len1 = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,edg2po2,sz)
                                : getlenedg_geosz<MFT,3,1>(msh,edg2po2,sz);
    edg2po2[1] = ip2;
    double len2 = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,edg2po2,sz)
                                : getlenedg_geosz<MFT,3,1>(msh,edg2po2,sz);

    CPRINTF1(" - %d bar1 = %e lens = %e %e valid %d %d (err = %e %e) dist %e %e sumlen %e err to sqrt(2) = %e\n",
              ntry_len,bar1,len1,len2,
              len1 > 1/sqrt(2), len2 > 1/sqrt(2),
              abs(len1 - 1/sqrt(2)), abs(len2 - 1/sqrt(2)), 
              abs(len1-1), abs(len2-1),
              len1+len2,
              abs(sqrt(2) - len1 - len2));

    if(len1 > len2){
      // make len1 shorter by increasing bar1 (pulling ipins towards ip1)
      bar1_min = bar1;
    }else{
      bar1_max = bar1;
    }

    bool ivalid = !icollapse ? len1 >= 1/sqrt(2) && len2 >= 1/sqrt(2) : true;

    if(ivalid){
      fnd_len = true;
      // Once a viable is found, it is possible length distance to 1 will 
      // make a couple of iterates not viable, so it's important to keep a viable
      // bar1. 
      double err = abs(abs(1-len1) - abs(1-len2));
      if(err < err_opt){
        err_opt = err;
        bar1_opt = bar1;
      }
      CPRINTF1(" - config error %e \n",err);
      if(err < 1.0e-2) break;
    }
  }// for ntry_len

  //if(DOPRINTS3()){
  //  writeMesh("debug_bisection",msh);
  //  msh.met.writeMetricFile("debug_bisection");
  //  for(int ii = npoi0; ii < msh.npoin; ii++){
  //    msh.killpoint(ii);
  //  }
  //  msh.set_npoin(npoi0);
  //  msh.set_nbpoi(nbpo0);
  //}

  if(!fnd_len) return 1;

  CPRINTF1("-- END aux_bisecPointLen w/ bar1 = %f\n",bar1_opt);
  return 0;
}

template int aux_bisecPointLen<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, 
                                                      int tdim, int ientt, int iedl,
                                                      int ibins, int tdimp, int iseed, int iref,
                                                      bool icollapse,
                                                      const MshCavity &cav);
template int aux_bisecPointLen<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, 
                                                      int tdim, int ientt, int iedl,
                                                      int ibins, int tdimp, int iseed, int iref,
                                                      bool icollapse,
                                                      const MshCavity &cav);



template<class MFT>
int aux_movePointCav(Mesh<MFT>& msh, MshCavity &cav, 
                     int tdimp, int iseed, int iref, double *algnd){
  GETVDEPTH(msh.param);
  int ierro = 0;

  // Interior and non-surface case, most straightforward
  if(tdimp == msh.idim){
    const intAr2 &ent2poi = msh.ent2poi(msh.get_tdim());
    int tdim = tdimp;

    double bary[4];
    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim + 1);

    double newp[3] = {0,0,0};
    double eval[3];
    //double metl[6];
    //const int nnmet = (msh.idim*(msh.idim+1))/2;
    double meast = 0;
    const intAr1& lcent = cav.lcent(tdim);
    for(int ientt : lcent){
      double wt;
      bool iflat;
      CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
        // gdim == tdim here
        if constexpr(gdim == 2){
          eval2<gdim,1>(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                        DifVar::None,bary,eval,NULL,NULL);
          wt = getmeasentP1<gdim,2>(msh, ent2poi[ientt], algnd, &iflat);
        }else{
          eval3<gdim,1>(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                        DifVar::None,bary,eval,NULL,NULL);
          wt = getmeasentP1<gdim,3>(msh, ent2poi[ientt], algnd, &iflat);
        }
      }}CT_FOR1(gdim);
      if(iflat) continue;

      //msh.met.getMetBary(AsDeg::P1, DifVar::None, msh.met.getSpace(), 
      //                   ent2poi[ientt], msh.get_tdim(), bary, 
      //                   metl, NULL);

      // For simply barycentre, use meas0
      // To skew towards the largest elements use meas0*meas0
      for(int ii = 0; ii < msh.idim;ii++) newp[ii] += wt*eval[ii];
      meast += wt;
    }
    METRIS_ASSERT(meast > 0);
    for(int ii = 0; ii < msh.idim;ii++){
      newp[ii] /= meast;
      msh.coord(cav.ipins,ii) = newp[ii];
    }

  }else{
    // Boundary case. 
    return 0;
  }// if tdimp == msh.get_tdim()

  // reinterp metric. This is always interior case, no need for ref of bdry dir
  ierro = msh.interpMetBack(cav.ipins,tdimp,iseed,iref,algnd);
  if(ierro != 0){
    CPRINTF1(" - interpMetBack failed ierro = %d \n",ierro);
    ierro = INS2D_ERR_INTERPMETBACK;
  }
  
  return ierro;
}

template
int aux_movePointCav<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, MshCavity &cav, 
                     int tdimp, int iseed, int iref, double *algnd);
template
int aux_movePointCav<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, MshCavity &cav, 
                     int tdimp, int iseed, int iref, double *algnd);

} // end namespace
