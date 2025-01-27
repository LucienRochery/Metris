//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_insert2D.hxx"
#include "low_insert2D.hxx"

#include "../low_lenedg.hxx"
#include "../aux_topo.hxx"
#include "../aux_timer.hxx"
#include "../mprintf.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../adapt/msh_swap2D.hxx"
#include "../BezierOffsets/low_gaps.hxx"
#include "../low_geo.hxx"
#include "../Mesh/Mesh.hxx"

#include <cmath>


namespace Metris{

// lpins work array sized dynamically by this routine ; it's an argument solely because this will be called several times, save on alloc
// also: as iterations go, fewer and fewer edges are long, no use allocating more than once to maximum needed size (first iter)
template<class MFT, int gdim, int ideg>
double insertLongEdges(Mesh<MFT> &msh, int *ninser, int ithrd1, int ithrd2, int ithrd3){
  GETVDEPTH(msh);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(ithrd1 != ithrd3);
  METRIS_ASSERT(ithrd2 != ithrd3);

  bool imovmet = true; 
  bool imovavg = true;
  bool type2 = false;

  // Swap norm -1: length-based. 
  //swapOptions swapOpt(100, -1, 0.0);
  swapOptions swapOpt(100, 0, 0.005);

  msh.met.setSpace(MetSpace::Log);


  if(msh.get_tdim() != 2) METRIS_THROW_MSG(TODOExcept(), 
    "Implement insertLongEdges on tdim != 2, got tdim = "<<msh.get_tdim());

  if(msh.bak != NULL) msh.bak->met.setSpace(MetSpace::Log);

  CPRINTF2("-- insertLongEdges start \n");
  #ifndef NDEBUG
  CPRINTF2(" - Note: improve by generating several points per edge. Generated but not used cf loop nn/2 \n");
  CPRINTF2(" - Note: improve by filtering point propositions \n");
  #endif
  int edg2pol[edgnpps[ideg]];


//  int npins = 0;
//  int mpins = lpins.size();

  double stat = 0;

  constexpr int tdim = gdim;

  constexpr auto getedgent = tdim == 2 ? getedgfac : getedgtet;

  const int nedgl = (tdim*(tdim+1))/2;

  int lnoed1[1][2] = {{0, 1}};
  const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                             tdim == 2 ? lnoed2[0] : lnoed3[0]);
  intAr2  &ent2poi = msh.ent2poi(tdim);
  intAr2r &ent2tag = msh.ent2tag(tdim);
  intAr2  &ent2ent = msh.ent2ent(tdim);
  //intAr2 &ent2ent = tdim == 1 ? msh.edg2edg : 
  //                  tdim == 2 ? msh.fac2fac : msh.tet2tet;
  // Ref of lower tdim entities
  //intAr1 &ent2ref = tdim == 2 ? msh.edg2ref : msh.fac2ref;


  //intAr2 &ent2tag = tdim == 1 ? msh.edg2tag : 
  //                  tdim == 2 ? msh.fac2tag : msh.tet2tag;


  const int mcaverr = CAV_ERR_NERROR;
  intAr1 lcaverr(mcaverr);
  lcaverr.set_n(mcaverr);

  const int minserr = 100;
  intAr1 linserr(minserr);
  linserr.set_n(minserr);

  const int miter1 = 30, miter2 = 30;

  double bar1[2] = {0.5,0.5};
  double coop[gdim], t0;
  int ierro;
  int ninser1 = 0;
  int ninser2 = 0;
  *ninser = 0;
  int niter1 = 0;
  // Untaged elements are to be considered 
  #ifndef GLOFRO
  msh.tag[ithrd1]++;
  #endif
  do{
    INCVDEPTH(msh);
    ninser1 = 0;

    msh.cleanup();

    //// New elements are untagged for ithrd1 -> no need to go over them 
    //t0 = get_wall_time();
    //double stat_swap = swap2D<MFT,gdim,ideg>(msh, swapOpt, ithrd2, ithrd3);
    //t2 = get_wall_time();
    //if(iverb >= 2){
    //  printf("   - Post inser swap time %f stat = %f \n",t2-t0,stat_swap);
    //}

    int niter2 = 0;
    do{
      ninser2 = 0;


      lcaverr.fill(0);
      linserr.fill(0);
      int nedgt = 0;
      int nerro = 0;
      int nent0 = msh.nentt(tdim);
      t0 = get_wall_time();
      for(int ientt = 0; ientt < nent0; ientt++){
        INCVDEPTH(msh);
        
        if(isdeadent(ientt,ent2poi)) continue;
        if(ent2tag(ithrd1,ientt) >= msh.tag[ithrd1]) continue;

        bool dobrk = false;
        // Whatever happens, this element will be tagged as inert. Indeed, 
        // it is either going to give rise to an insertion (thus disappear)
        // or it won't (thus become inert)
        // The only exception is if an insertion is rejected due to short edge, 
        // or other cavity extension routines, as neighbours influence this decision. 
        ent2tag(ithrd1,ientt) = msh.tag[ithrd1];
        for(int ied = 0; ied < nedgl; ied++){
          nedgt++;

          // Pretty accurate!
          double sz[2];
          double len = getlenedg_geosz<MFT,gdim,ideg>(msh,ientt,tdim,ied,sz);
          CPRINTF1(" - try ientt = %d ied = %d len = %f \n",ientt,ied,len);
          if(len < sqrt(2)) continue;

          edg2pol[0] = ent2poi(ientt,lnoed[ied][0]);
          edg2pol[1] = ent2poi(ientt,lnoed[ied][1]);
          int idx0 = tdim + 1 + ied*(ideg-1);
          for(int ii = 0; ii < ideg-1; ii++) edg2pol[2+ii] = ent2poi[ientt][idx0+ii];

          eval1<gdim,ideg>(msh.coord, edg2pol, msh.getBasis(), 
                           DifVar::None, DifVar::None, 
                           bar1, coop, NULL, NULL);

          if(DOPRINTS1()){
            CPRINTF1(" - enact ins ientt = %d ied = %d len = %f edg %d %d  coord = ",
                   ientt,ied,len,edg2pol[0],edg2pol[1]);
            dblAr1(msh.idim,coop).print();
            CPRINTF1(" - bary = ");
            dblAr1(3, bar1).print();
          }

          #if 0
          if(imovmet){
            // Do better than this, compute the Bézier offset for the metric and 
            // place the point here -> follow curvature (even if P1)
            double offset[gdim];
            getBezOffsetsEdge<MFT,gdim,ideg>(msh,tdim,ent2poi[ientt],ied,offset);

            double nrm1 = geterrl2<gdim>(msh.coord[edg2pol[0]],msh.coord[edg2pol[0]]);
            double nrm2 = getnrml2<gdim>(offset);

            // Maximum ratio -> dampen if more 
            const double ratlen = 0.1;

            double fac = 1;
            if(nrm2 > ratlen*ratlen*nrm1){
              fac = ratlen*sqrt(nrm1/nrm2);
            }

            for(int ii = 0; ii < gdim; ii++) coop[ii] += offset[ii]*fac;
            for(int ii = 0; ii < gdim; ii++) coop[ii] = 0.25 * msh.coord[edg2pol[0]][ii]
              + 0.5 * coop[ii] + 0.25 * msh.coord[edg2pol[1]][ii];
            METRIS_ASSERT(ideg == 1);
          }else if(imovavg){
            int ineig = ent2ent(ientt,ied);
            //double bary[tdim+1] = {1.0/(tdim+1)}; 
            //double eval[gdim];
            //for(int ii = 0; ii < gdim; ii++) coop[ii] = 0;
            //double wttot = 0;
            //for(int ient2 : {ientt, ineig}){
            //  if(ient2 < 0) continue;
            //  eval2<gdim,1>(msh.coord,msh.fac2poi[ient2],msh.getBasis(),
            //                DifVar::None,DifVar::None,
            //                bary,eval,NULL,NULL);
            //  double meas0 = getmeasentP1<gdim>(ent2poi[ient2], msh.coord);
            //  for(int jj = 0; jj < gdim; jj++) coop[jj] += meas0*eval[jj];
            //  wttot += meas0;
            //}
            //for(int jj = 0; jj < gdim; jj++) coop[jj] /= wttot;
            if(ineig >= 0){
              int ipoi1 = ent2poi(ientt,ied);
              int ie2 = getedgent(msh,ineig,edg2pol[0],edg2pol[1]);
              int ipoi2 = ent2poi(ineig,ie2);
              for(int jj = 0; jj < gdim; jj++){
                coop[jj] = 0.5*msh.coord(ipoi1,jj) + 0.5*msh.coord(ipoi2,jj);
              }
            }
          }
          #endif

          // Try insert point coop
          int nent00 = msh.nentt(tdim); 
          int itry = 0;
          do{
            ierro = insedgesurf(msh,ientt,ied,coop,bar1[0],
                                lcaverr,ithrd2,ithrd3);
            if(ierro <= 0) break;
            itry++;
            if(itry >= 1 + imovmet + imovavg) break;
            if(ierro == 0) CPRINTF1(" - After trying ierro = 0 \n");
            if(ierro > 0 && (itry == 0 && imovmet || imovavg)){
              CPRINTF1(" -> insedgesurf fail: try again w/ imovmet %d imovavg %d\n",imovmet,imovavg);
              if(DOPRINTS1()){
                CPRINTF1(" - initial ipins = ");
                dblAr1(gdim,coop).print();
              }
              if(imovmet && itry == 1){
                CPRINTF1(" -> do imovmet\n");
                // Do better than this, compute the Bézier offset for the metric and 
                // place the point here -> follow curvature (even if P1)
                double offset[gdim];
                getBezOffsetsEdge<MFT,gdim,ideg>(msh,tdim,ent2poi[ientt],ied,offset);

                double nrm1 = geterrl2<gdim>(msh.coord[edg2pol[0]],msh.coord[edg2pol[0]]);
                double nrm2 = getnrml2<gdim>(offset);

                // Maximum ratio -> dampen if more 
                const double ratlen = 0.1;

                double fac = 1;
                if(nrm2 > ratlen*ratlen*nrm1){
                  fac = ratlen*sqrt(nrm1/nrm2);
                }

                for(int ii = 0; ii < gdim; ii++) coop[ii] += offset[ii]*fac;
                for(int ii = 0; ii < gdim; ii++) coop[ii] = 0.25 * msh.coord[edg2pol[0]][ii]
                  + 0.5 * coop[ii] + 0.25 * msh.coord[edg2pol[1]][ii];
              }else if(imovavg){
                // This can put the point outside ! Depends on the boundary
                int ineig = ent2ent(ientt,ied);
                //double bary[tdim+1] = {1.0/(tdim+1)}; 
                //double eval[gdim];
                //for(int ii = 0; ii < gdim; ii++) coop[ii] = 0;
                //double wttot = 0;
                //for(int ient2 : {ientt, ineig}){
                //  if(ient2 < 0) continue;
                //  eval2<gdim,1>(msh.coord,msh.fac2poi[ient2],msh.getBasis(),
                //                DifVar::None,DifVar::None,
                //                bary,eval,NULL,NULL);
                //  double meas0 = getmeasentP1<gdim>(ent2poi[ient2], msh.coord);
                //  for(int jj = 0; jj < gdim; jj++) coop[jj] += meas0*eval[jj];
                //  wttot += meas0;
                //}
                //for(int jj = 0; jj < gdim; jj++) coop[jj] /= wttot;
                if(ineig >= 0){
                  int ipoi1 = ent2poi(ientt,ied);
                  int ie2 = getedgent(msh,ineig,edg2pol[0],edg2pol[1]);
                  int ipoi2 = ent2poi(ineig,ie2);
                  CPRINTF1(" -> do imovavg use ineig %d ip1/2 %d %d\n",
                    ineig,ipoi1,ipoi2);
                  for(int jj = 0; jj < gdim; jj++){
                    coop[jj] = 0.5*msh.coord(ipoi1,jj) + 0.5*msh.coord(ipoi2,jj);
                  }
                }
              }
              if(DOPRINTS1()){
                CPRINTF1(" - final ipins = ");
                dblAr1(gdim,coop).print();
              }
            }
          }while(ierro != 0);

          if(ierro <= 0){
            ninser2++;
            dobrk = true;
            int nent11 = msh.nentt(tdim);
            for(int ientn = nent00; ientn < nent11; ientn++){
              for(int ii = 0; ii < tdim + 1 ; ii++){
                int ineig = ent2ent(ientn,ii);
                if(ineig < 0) continue;
                METRIS_ASSERT(!isdeadent(ineig,ent2poi));
                ent2tag(ithrd1,ineig) = msh.tag[ithrd1] - 1;
              }
            }
          }else{
            CPRINTF1(" - insertion failed ierro = %d \n",ierro);
            linserr[ierro - 1] ++;
            nerro++;
            // If cavity selection error, untag as inert, as neighbours may influence yet
            //if(ierro == INS2D_ERR_INCCAV2D 
            //|| ierro == INS2D_ERR_INCCAV2D2 
            //|| ierro == INS2D_ERR_SHORTEDG){
            //  ent2tag(ithrd1,ientt) = msh.tag[ithrd1] - 1;
            //}
          }

          if(dobrk) break;
        }

      }

      double t1 = get_wall_time();
      int ncallps = 1000*(int)((ninser2 / (t1-t0)) / 1000);
      CPRINTF1(" - Loop 1 end t = %f ninser %d = %d /s; nerro %d coll \n",
                t1-t0,ninser2,ncallps,nerro);
      if(DOPRINTS2() && nerro > 0){
        CPRINTF2(" - cavity ierro list:\n");
        for(int ii = 0; ii < mcaverr; ii++){
          if(lcaverr[ii] == 0) continue;
          CPRINTF2("   ierro = %d : %d \n",ii+1,lcaverr[ii]);
        }
        CPRINTF2(" - inspoi ierro list:\n");
        for(int ii = 0; ii < minserr; ii++){
          if(linserr[ii] == 0) continue;
          CPRINTF2("   ierro = %d : %d \n",ii+1,linserr[ii]);
        }
      }

      //if(msh.param->iverb == 4){
      //  printf("DEBUG KILL AFTER FIRST IVERB 4\n");
      //  exit(1);
      //}

      //if(ninser2 == 50){
      //  printf("## DEBUG SET IVERB = 4\n");
      //  msh.param->iverb = 4;
      //}

      // LOOP 2: over elements ; try element barycentres 
      if(type2){

        ninser2 = 0;
        lcaverr.fill(0);
        linserr.fill(0);
        nerro = 0;
        nent0 = msh.nentt(tdim);

        int ipdum = msh.newpoitopo(2,-1); 
        int edg2pol[2];
        nerro = 0;
        for(int ientt = 0; ientt < nent0; ientt++){
          if(isdeadent(ientt,ent2poi)) continue;
          if(ent2tag(ithrd1,ientt) >= msh.tag[ithrd1]) continue;

          double bar2[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
          eval1<gdim,ideg>(msh.coord, msh.fac2poi[ientt], msh.getBasis(), 
                           DifVar::None, DifVar::None, 
                           bar2, msh.coord[ipdum], NULL, NULL);
      
          //ierro = msh.interpMetBack(ipdum,tdim,ientt,
          //                          -1, NULL,
          //                          &msh.poi2bak(msh.fac2poi(ientt,0),tdim-1),
          //                          msh.met[ipdum]);
          ierro = msh.interpMetBack(ipdum,tdim,ientt,-1, NULL);
          
          if(ierro > 0){
            nerro++;
            continue;
          }
          edg2pol[0] = ipdum;
          bool oneshort = false;
          for(int ii = 0; ii < 3; ii++){
            double sz[2];
            edg2pol[1] = msh.fac2poi(ientt,ii);
            double len = getlenedg_geosz<MFT,gdim,1>(msh,edg2pol,sz);
            if(len < 1.0/sqrt(2)) oneshort = true;
            if(oneshort) break;
          }
          if(oneshort) continue;

          int nent00 = msh.nentt(tdim); 
          ierro = insfacsurf(msh,ientt,coop,lcaverr,ithrd2,ithrd3);
          if(ierro <= 0){
            ninser2++;
            int nent11 = msh.nentt(tdim);
            for(int ientn = nent00; ientn < nent11; ientn++){
              for(int ii = 0; ii < tdim + 1 ; ii++){
                int ineig = ent2ent(ientn,ii);
                if(ineig < 0) continue;
                METRIS_ASSERT(!isdeadent(ineig,ent2poi));
                ent2tag(ithrd1,ineig) = msh.tag[ithrd1] - 1;
              }
            }
          }else{
            linserr[ierro - 1] ++;
            nerro++;
          }
        }


        double t2 = get_wall_time();
        ncallps = 1000*(int)((ninser2 / (t2-t1)) / 1000);
          CPRINTF1("   - Loop 2 end t = %f ninser %d = %d /s; nerro %d coll \n",
                   t2-t1,ninser2,ncallps,nerro);
        if(DOPRINTS2() && nerro > 0){
          CPRINTF2("   - cavity ierro list:\n");
          for(int ii = 0; ii < mcaverr; ii++){
            if(lcaverr[ii] == 0) continue;
            CPRINTF2("     ierro = %d : %d \n",ii+1,lcaverr[ii]);
          }
          CPRINTF2("   - inspoi ierro list:\n");
          for(int ii = 0; ii < minserr; ii++){
            if(linserr[ii] == 0) continue;
            CPRINTF2("     ierro = %d : %d \n",ii+1,linserr[ii]);
          }
        }
      }

      if(nedgt == 0) stat = 0;
      else           stat = MAX(stat, (double)ninser2/(double)nedgt);

      ninser1 += ninser2;
    }while(ninser2 > 0 && niter2++ < miter2);

    *ninser += ninser1;
  }while(ninser1 > 0 && niter1++ < miter1);


  //int ipnew = msh.npoin;
  //msh.set_npoin(msh.npoin+1);
  //int ngood = 0;
  //int ngood2 = 0;
  //int nerro = 0;
  //double vmin = 1.0e30, vmax = -1.0e30;
  //for(int iface = 0; iface < msh.nface; iface++){

  //  if(isdeadent(iface,msh.fac2poi)) continue;

  //  double bar2[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
  //  eval1<gdim,ideg>(msh.coord, msh.fac2poi[iface], msh.getBasis(), 
  //                   DifVar::None, DifVar::None, 
  //                   bar1, msh.coord[ipnew], NULL, NULL);

  //  ierro = msh.interpMetBack(msh.coord[ipnew],iface,
  //                                  &msh.poi2bak[msh.fac2poi(iface,0)],
  //                                  msh.met[ipnew]);
  //  if(ierro > 0){
  //    nerro++;
  //    continue;
  //  }
  //  double nrmal[3];
  //  bool iflat;
  //  double vtol = 1.0e-12;
  //  double meas0 = getmeasentP1<2,2>(msh.fac2poi[iface], msh.coord, 
  //                                   vtol, &iflat, iverb-1);

  //  double det = detsym<2,double>(msh.met[ipnew]);

  //  double vol = meas0*sqrt(det);
  //  vmin = MIN(vol,vmin);
  //  vmax = MAX(vol,vmax);

  //  if(vol > 3*sqrt(3)/2){
  //    ngood2++;
  //  }

  //  bool oneshort = false;
  //  int edg2pol[2];
  //  edg2pol[0] = ipnew;
  //  for(int ii = 0; ii < 3; ii++){
  //    double sz[2];
  //    edg2pol[1] = msh.fac2poi(iface,ii);
  //    double len = getlenedg_geosz<MFT,gdim,1>(msh,edg2pol,sz);
  //    if(len < 1.0/sqrt(2)) oneshort = true;
  //    if(oneshort) break;
  //  }
  //  if(!oneshort) ngood++;
  //}
  //if(ngood > 0) printf("## DEBUG: NERRO = %d NGOOD1 = %d NGOOD2 = %d vmin = %f vmax = %f\n",
  //  nerro,ngood,ngood2,vmin,vmax);
  //msh.set_npoin(msh.npoin-1);




  return stat;
}

#if 0

// lpins work array sized dynamically by this routine ; it's an argument solely because this will be called several times, save on alloc
// also: as iterations go, fewer and fewer edges are long, no use allocating more than once to maximum needed size (first iter)
template<class MFT, int gdim, int ideg>
double insertLongEdges(Mesh<MFT> &msh,   int iverb, int ithrd1 ){

  if(msh.get_tdim() != 2) METRIS_THROW_MSG(TODOExcept(), 
    "Implement insertLongEdges on tdim != 2, got tdim = "<<msh.get_tdim());

  msh.met.setSpace(MetSpace::Log);
  if(msh.bak != NULL) msh.bak->met.setSpace(MetSpace::Log);

  if(iverb >= 1) printf("  -- insertLongEdges start \n");
  if(iverb >= 1) printf("   - Note: improve by generating several points per edge. Generated but not used cf loop nn/2 \n");
  if(iverb >= 1) printf("   - Note: improve by filtering point propositions \n");
  int edg2pol[edgnpps[ideg]];

//  int npins = 0;
//  int mpins = lpins.size();

  double stat = 0;

  int tdim = msh.get_tdim();

  const int nedgl = (tdim*(tdim+1))/2;

  int lnoed1[1][2] = {{0, 1}};
  const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                             tdim == 2 ? lnoed2[0] : lnoed3[0]);

  intAr2 &ent2poi = tdim == 1 ? msh.edg2poi : 
                    tdim == 2 ? msh.fac2poi : msh.tet2poi;

  //intAr2 &ent2ent = tdim == 1 ? msh.edg2edg : 
  //                  tdim == 2 ? msh.fac2fac : msh.tet2tet;
  // Ref of lower tdim entities
  //intAr1 &ent2ref = tdim == 2 ? msh.edg2ref : msh.fac2ref;


  //intAr2 &ent2tag = tdim == 1 ? msh.edg2tag : 
  //                  tdim == 2 ? msh.fac2tag : msh.tet2tag;


  int nentt = tdim == 1 ? msh.nedge : 
              tdim == 2 ? msh.nface : msh.nelem;

  const int mcaverr = CAV_ERR_NERROR;
  intAr1 lcaverr(mcaverr);
  lcaverr.set_n(mcaverr);

  const int minserr = 100;
  intAr1 linserr(minserr);
  linserr.set_n(minserr);


  int ierro;
  int ninser = 0;
  do{
    lcaverr.fill(0);
    linserr.fill(0);
    int nedgt = 0;
    ninser = 0;
    int nerro = 0;
    int ncall = 0;
    int nent0 = nentt;
    //msh.tag[ithrd1]++;
    double t0 = get_wall_time();
    for(int ientt = 0; ientt < nent0; ientt++){
      if(isdeadent(ientt,ent2poi)) continue;
      //ent2tag(ithrd1,ientt) = msh.tag[ithrd1];
      bool dobrk = false;
      for(int ied = 0; ied < nedgl; ied++){
        nedgt++;

        //int ient2 = msh.fac2fac(ientt,ied);
        //if(ient2 >= 0){
        //  if(ent2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
        //}

        //double len1 = getlenedg_quad<MFT,gdim,ideg>(msh,ientt,2,ied,nquad);

        // Pretty accurate!
        double sz[2];
        double len = getlenedg_geosz<MFT,gdim,ideg>(msh,ientt,tdim,ied,sz);

        if(len < sqrt(2)) continue;

        int nn = round(len) - 1;

        // If len < 1.5 but > sqrt(2)
        if(nn == 0) nn = 1;

        METRIS_ASSERT(len / (nn + 1) > 1/sqrt(2));




        edg2pol[0] = ent2poi(ientt,lnoed[ied][0]);
        edg2pol[1] = ent2poi(ientt,lnoed[ied][1]);
        int idx0 = tdim + 1 + ied*(ideg-1);
        for(int ii = 0; ii < ideg-1; ii++) edg2pol[2+ii] = ent2poi[ientt][idx0+ii];
        double bar1[2];
        double coop[gdim];

        //bool ibdry = false;
        //int ibent = -1;
        //if(tdim == 2){
        //  ibent = getedgglo(msh,edg2pol[0],edg2pol[1]);
        //  if(ibent >= 0) ibdry = true;
        //}else{
        //  METRIS_THROW_MSG(TODOExcept(),"3D surface. Note check whether face or edge");
        //}
        //METRIS_ASSERT(!(ibent < 0 && ent2ent(ientt,ied) < 0));

        // geometric length model is as follows
        // a = sz[0] / sz[1]
        // len = sz[1] int_0^1 a^t dt 
        // hence length between t0 and t1 is 
        // either (t1 - t0) if a = 1
        // or: len(t0->t1) = sz[1]*(a^t1 - a^t0)/ln(a)

        double a = sz[1] / sz[0];
        double at1;//, t1;
        double loga;

        const double tola = 1.0e-14;

        if(!(abs(a-1.0) < tola)){
          at1 = (a+nn) / (nn + 1.0); // a^t1 directly
          loga = log(a);
          //t1 = log((a+nn) / (nn + 1.0)) / log(a);
        }
        //for(int ii = 0; ii < nn; ii++){
        for(int ii = (nn-1)/2; ii < (nn-1)/2 + 1; ii++){
          double t;
          if(abs(a-1.0) < tola){
            t = (ii + 1.0) / (nn + 1.0);
          }else{
            //t = (ii + 1.0) * at1 - ii;
            t = log((ii + 1.0) * at1 - ii) / loga;
            //printf("Debug mode 2 \n");
          }
          bar1[0] = t;
          bar1[1] = 1.0 - t;

          //printf("Debug ii = %d t = %f \n",ii,t);
          //if(!(abs(a-1.0) < tola) && ii == 0) wait();
          ////wait();

          METRIS_ASSERT(t > 1.0/1000.0);


          //if(!msh.CAD() || !ibdry){
            eval1<gdim,ideg>(msh.coord, edg2pol, msh.getBasis(), 
                             DifVar::None, DifVar::None, 
                             bar1, coop, NULL, NULL);
          //}else if(tdim == 2){
          //  double param;
          //  int ib[2];
          //  int iref = ent2ref[ibent];
          //  METRIS_ASSERT(iref >= 0);

          //  for(int ii = 0; ii < 2; ii++){
          //    ib[ii] = msh.poi2bpo[edg2pol[ii]];

          //    if(msh.bpo2ibi[ib[ii]][1] < (tdim - 1)){ // Get ent specific
          //      ib[ii] = getent2bpo(msh,ib[ii],ibent,tdim-1);
          //    }else{ // Single for ref
          //      ib[ii] = getref2bpo(msh,ib[ii],iref,tdim-1);
          //    }
          //    METRIS_ASSERT(ib[ii] >= 0);
          //  }

          //  param = t * msh.bpo2rbi[ib[0]][0] + (1.0 - t) * msh.bpo2rbi[ib[1]][0];

          //  double result[18];
          //  ego obj = msh.CAD.cad2edg[iref];
          //  ierro = EG_evaluate(obj,&param,result);

          //  METRIS_ASSERT(ierro == 0);

          //  if(ierro != 0){
          //    nerro++;
          //    continue;
          //  }

          //  for(int ii = 0; ii < gdim; ii++) coop[ii] = result[ii];

          //}else{
          //  METRIS_THROW_MSG(TODOExcept(),"3D surface");
          //}




          // Try insert point coop
          ierro = insedgesurf(msh,ientt,ied,coop,bar1[0],iverb,lcaverr,ithrd1);
          ncall++;
          if(ierro <= 0) ninser++;

          if(ierro <= 0) dobrk = true;
          else{
            linserr[ierro - 1] ++;
            nerro++;
            #ifndef NDEBUG
            if(msh.param->dbgfull){
              printf("Wait error = %d \n",ierro);
              wait();
            }
            #endif
          }
        } 

        if(dobrk) break;
      }

    }

    double t1 = get_wall_time();
    int ncallps = 1000*(int)((ninser / (t1-t0)) / 1000);
    if(iverb >= 1){
      printf("   - Loop end ninser %d = %d /s; nerro %d coll \n",
        ninser,ncallps,nerro);
      if(iverb >= 2){
        if(nerro > 0){
          printf("   - cavity ierro list:\n");
          for(int ii = 0; ii < mcaverr; ii++){
            if(lcaverr[ii] == 0) continue;
            printf("     ierro = %d : %d \n",ii+1,lcaverr[ii]);
          }
          printf("   - inspoi ierro list:\n");
          for(int ii = 0; ii < minserr; ii++){
            if(linserr[ii] == 0) continue;
            printf("     ierro = %d : %d \n",ii+1,linserr[ii]);
          }
        }
      }
    }
    stat = MAX(stat, (double)ninser/(double)nedgt);

  }while(ninser > 0);

  return stat;
}
#endif

// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double insertLongEdges<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                                int* ninser, int ithrd1, int ithrd2, int ithrd3);\
template double insertLongEdges<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                                int* ninser, int ithrd1, int ithrd2, int ithrd3);\
template double insertLongEdges<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                                int* ninser, int ithrd1, int ithrd2, int ithrd3);\
template double insertLongEdges<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                                int* ninser, int ithrd1, int ithrd2, int ithrd3);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()





} // end namespace