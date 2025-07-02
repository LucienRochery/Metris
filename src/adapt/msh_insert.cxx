//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_insert.hxx"
#include "low_insert.hxx"

#include "../low_lenedg.hxx"
#include "../aux_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../utils/aux_timer.hxx"
#include "../utils/mprintf.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../adapt/msh_swap.hxx"
#include "../BezierOffsets/low_gaps.hxx"
#include "../low_geo.hxx"
#include "../Mesh/Mesh.hxx"
#include "../msh_checktopo.hxx"
#include "../aux_histogram.hxx"
#include "../msh_lenedg.hxx"

#include <cmath>


namespace Metris{

// This version tags elements and also edges via a hash table. 
// Initially, all elements are untagged (active) and the hash table is empty. 
// When an insertion fails on an edge, the edge is added to the hashtable. 
// When insertions fail on all edges of an element, the element becomes tagged (inactive)
// When an insertion is carried out, neighbouring edges are not untagged because
// the majority of rejections are due to error INS2D_ERR_SHORTEDG which is only
// made worse by an insertion. 
// Edges eliminated by other errors can be removed from the hashtable in a second
// time. 
// Note their errors might depend on greater than 1 neighbourhood; hence there 
// is no perfect solution other than to allow a full run at some point (in the end). 

// lpins work array sized dynamically by this routine ; it's an argument solely because this will be called several times, save on alloc
// also: as iterations go, fewer and fewer edges are long, no use allocating more than once to maximum needed size (first iter)
template<class MFT, int gdim, int ideg>
double insertLongEdges(Mesh<MFT> &msh, int *ninser, int ithrd1, int ithrd2, int ithrd3){
  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(ithrd1 != ithrd3);
  METRIS_ASSERT(ithrd2 != ithrd3);

  bool ishort_conservative = false;

  bool imovmet = true; 

  // Swap norm -1: length-based. 
  //swapOptions swapOpt(100, -1, 0.0);
  swapOptions swapOpt(100, 0, 0.005);

  //msh.met.setSpace(MetSpace::Log);
  msh.met.setSpace(MetSpace::Exp);

  CPRINTF2("-- insertLongEdges start \n");
  #ifndef NDEBUG
  CPRINTF2(" - Note: improve by generating several points per edge. Generated but not used cf loop nn/2 \n");
  CPRINTF2(" - Note: improve by filtering point propositions \n");
  #endif
  int edg2pol[getnnod1(ideg)];


  double stat = 0;

  const int tdim = msh.get_tdim();

  const int nedgl = (tdim*(tdim+1))/2;

  const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                             tdim == 2 ? lnoed2[0] : lnoed3[0]);
  intAr2  &ent2poi = msh.ent2poi(tdim);
  intAr2r &ent2tag = msh.ent2tag(tdim);
  intAr2  &ent2ent = msh.ent2ent(tdim);

  const int mcaverr = CAV_ERR_NERROR;
  intAr1 lcaverr(mcaverr);

  const int minserr = 100;
  intAr1 linserr(minserr);

  const int miter2 = 3;

  // Store errors per edge. Errors unused for now, only presence in hash table.
  HshTab_I2I ledro; 

  double bar1[2] = {0.5,0.5};
  double coop[gdim], t0;
  int ierro;
  *ninser = 0;

  msh.tag[ithrd1]++;
  
  MshCavity cav(100,100,1);
  CavWrkArrs work;

  int niter2 = 0;
  int ninser2 = 1; // otherwise doesnt enter loop
  for(niter2 = 0; niter2 < miter2 && ninser2 > 0; niter2++){
    INCVDEPTH(msh.param);
    ninser2 = 0;

    lcaverr.fill(0);
    linserr.fill(0);
    int nedgt = 0;
    int nlong = 0;
    int nskip = 0;
    int nerro = 0;

    int max_nnewp = -1;
    double avg_nnewp = 0;

    int nent0 = msh.nentt(tdim);
    t0 = get_wall_time();
    for(int ientt = 0; ientt < nent0; ientt++){
      INCVDEPTH(msh.param);
      
      if(isdeadent(ientt,ent2poi)) continue;
      if(ent2tag(ithrd1,ientt) >= msh.tag[ithrd1]) continue;

      bool dobrk = false;
      // Whatever happens, this element will be tagged as inert. Indeed, 
      // it is either going to give rise to an insertion (thus disappear)
      // or it won't (thus become inert)
      // The only exception is if an insertion is rejected due to short edge, 
      // or other cavity extension routines, as neighbours influence this decision. 
      ent2tag(ithrd1,ientt) = msh.tag[ithrd1];
      bool safe_continue = true;
      for(int ied = 0; ied < nedgl; ied++){
        // Forbid continues because of ledro management. Error must be set. 
        METRIS_ASSERT(safe_continue);
        safe_continue = false;

        GETVDEPTH(msh.param);
        nedgt++;

        int ip1 = ent2poi(ientt, lnoed(ied,0));
        int ip2 = ent2poi(ientt, lnoed(ied,1));
        auto key = stup2(ip1,ip2);
        auto tt  = ledro.find(key);
        if(tt != ledro.end()){
          CPRINTF2(" - ientt %d ied %d ip %d %d in hash tab -> skip\n",
                   ientt, ied, ip1, ip2);
          nskip++;
          safe_continue = true;
          continue;
        }


        double sz[2];
        double len = getlenedg_geosz<MFT,gdim,ideg>(msh,ientt,tdim,ied,sz);
        CPRINTF2(" - try ientt = %d ied = %d len = %f \n",ientt,ied,len);
        if(len < sqrt(2)){
          ledro[key] = -1;
          safe_continue = true;
          continue;
        }
        nlong++;

        int nnewp;
        if(ishort_conservative){
          // Make no edges shorter than 1
          // In this case, take the floor of length. 
          nnewp = floor(len);
        }else{
          // Allow full range [1/sqrt(2),sqrt(2)]
          bool ifnd = false;
          // Distance to 1 decreases until it increases. 
          for(nnewp = 1; nnewp < 1000; nnewp++){
            // Length is 1 / (nnewp + 1)
            double err1 = abs(len/(nnewp+1) - 1);
            double err2 = abs(len/(nnewp+2) - 1);
            if(err2 > err1){
              ifnd = true;
              break;
            }
          }
          METRIS_ASSERT_MSG(ifnd, "Infinite length? "<<len);
          METRIS_ASSERT_MSG(len/(nnewp+1) > 1/sqrt(2),
            "with len = "<<len<<" got nnewp = "<<nnewp<<" each is "
            <<len/nnewp);
        }
        avg_nnewp += nnewp;
        max_nnewp = MAX(max_nnewp, nnewp);

        edg2pol[0] = ent2poi(ientt,lnoed[ied][0]);
        edg2pol[1] = ent2poi(ientt,lnoed[ied][1]);


        CPRINTF2(" - No CAD link for this edge -> eval at element\n");
        int idx0 = tdim + 1 + ied*(ideg-1);
        for(int ii = 0; ii < ideg-1; ii++) edg2pol[2+ii] = ent2poi[ientt][idx0+ii];


        // Note: no CAD use as insertEdge handles CAD eval. This routine only 
        // generates points approximately. 
        eval1<gdim,ideg>(msh.coord, edg2pol, msh.getBasis(), 
                         DifVar::None, DifVar::None, 
                         bar1, coop, NULL, NULL);

        if(DOPRINTS1()){
          CPRINTF2(" - enact ins ientt = %d ied = %d len = %f edg %d %d coord = ",
                 ientt,ied,len,edg2pol[0],edg2pol[1]);
          dblAr1(msh.idim,coop).print();
        }

        // Try insert point coop
        int nent00 = msh.nentt(tdim); 
        int itry = 0;
        do{
          if(DOPRINTS2()){
            writeMesh("preins",msh);
          }

          ierro = insertEdge(msh,tdim,ientt,ied,coop,bar1[0],
                             cav,work,lcaverr,ithrd2,ithrd3);

          if(ierro <= 0) break;
          itry++;
          if(itry >= 1 + imovmet) break;
          if(ierro == 0) CPRINTF2(" - After trying ierro = 0 \n");
          if(ierro > 0 && (itry == 0 && imovmet)){
            CPRINTF2(" -> insertEdge fail: try again w/ imovmet %d\n",imovmet);
            if(DOPRINTS1()){
              CPRINTF2(" - initial ipins = ");
              dblAr1(gdim,coop).print();
            }
            if(imovmet && itry == 1){
              CPRINTF2(" -> do imovmet\n");
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
              for(int ii = 0; ii < gdim; ii++) coop[ii] = 0.25 * msh.coord(edg2pol[0],ii)
                + 0.5 * coop[ii] + 0.25 * msh.coord(edg2pol[1],ii);
            }
            if(DOPRINTS1()){
              CPRINTF2(" - final ipins = ");
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
          CPRINTF2(" - insertion failed ierro = %d \n",ierro);
          linserr[ierro - 1] ++;
          nerro++;
          ledro[key] = ierro;
        }


        if(msh.param->dbgfull) check_topo(msh,ithrd2);

        if(dobrk) break;

        safe_continue = true;
      }// for ied

    }// for ientt

    avg_nnewp /= nlong;
    double t1 = get_wall_time();
    int ncallps = 1000*(int)((ninser2 / (t1-t0)) / 1000);
    CPRINTF2(" - Loop 1 end t = %f nlong %d nskip %d ninser %d = %d /s; nerro %d; nnewp avg = %f, max = %d\n",
              t1-t0,nlong,nskip,ninser2,ncallps,nerro,avg_nnewp, max_nnewp);
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


    msh.cleanup();

    if(nedgt == 0) stat = MAX(stat, 0);
    else           stat = MAX(stat, (double)ninser2/(double)nedgt);

    *ninser += ninser2;
  }// for niter2


  return stat;
}

#if 0
//Same as next: did not insert several points at once, but computed how many might
//be made. Was still looping over elements, not edges. 

// lpins work array sized dynamically by this routine ; it's an argument solely because this will be called several times, save on alloc
// also: as iterations go, fewer and fewer edges are long, no use allocating more than once to maximum needed size (first iter)
template<class MFT, int gdim, int ideg>
double insertLongEdges(Mesh<MFT> &msh, int *ninser, int ithrd1, int ithrd2, int ithrd3){
  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(ithrd1 != ithrd3);
  METRIS_ASSERT(ithrd2 != ithrd3);

  bool ishort_conservative = false;

  bool imovmet = true; 

  // Swap norm -1: length-based. 
  //swapOptions swapOpt(100, -1, 0.0);
  swapOptions swapOpt(100, 0, 0.005);

  //msh.met.setSpace(MetSpace::Log);
  msh.met.setSpace(MetSpace::Exp);

  CPRINTF2("-- insertLongEdges start \n");
  #ifndef NDEBUG
  CPRINTF2(" - Note: improve by generating several points per edge. Generated but not used cf loop nn/2 \n");
  CPRINTF2(" - Note: improve by filtering point propositions \n");
  #endif
  int edg2pol[getnnod1(ideg)];


  double stat = 0;

  const int tdim = msh.get_tdim();

  const int nedgl = (tdim*(tdim+1))/2;

  int lnoed1[1][2] = {{0, 1}};
  const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                             tdim == 2 ? lnoed2[0] : lnoed3[0]);
  intAr2  &ent2poi = msh.ent2poi(tdim);
  intAr2r &ent2tag = msh.ent2tag(tdim);
  intAr2  &ent2ent = msh.ent2ent(tdim);

  const int mcaverr = CAV_ERR_NERROR;
  intAr1 lcaverr(mcaverr);

  const int minserr = 100;
  intAr1 linserr(minserr);

  const int miter2 = 3;


  double bar1[2] = {0.5,0.5};
  double coop[gdim], t0;
  int ierro;
  *ninser = 0;

  msh.tag[ithrd1]++;
  
  // At most one edge in an insertion that doesn't collapse a point.
  MshCavity cav(100,100,1);
  CavWrkArrs work;

  int niter2 = 0;
  int ninser2 = 1; // otherwise doesnt enter loop
  for(niter2 = 0; niter2 < miter2 && ninser2 > 0; niter2++){
    INCVDEPTH(msh.param);
    ninser2 = 0;

    lcaverr.fill(0);
    linserr.fill(0);
    int nedgt = 0;
    int nlong = 0;
    int nerro = 0;

    int max_nnewp = -1;
    double avg_nnewp = 0;

    int nent0 = msh.nentt(tdim);
    t0 = get_wall_time();
    for(int ientt = 0; ientt < nent0; ientt++){
      INCVDEPTH(msh.param);
      
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
        GETVDEPTH(msh.param);
        nedgt++;

        double sz[2];
        double len = getlenedg_geosz<MFT,gdim,ideg>(msh,ientt,tdim,ied,sz);
        CPRINTF2(" - try ientt = %d ied = %d len = %f \n",ientt,ied,len);
        if(len < sqrt(2)) continue;
        nlong++;

        int nnewp;
        if(ishort_conservative){
          // Make no edges shorter than 1
          // In this case, take the floor of length. 
          nnewp = floor(len);
        }else{
          // Allow full range [1/sqrt(2),sqrt(2)]
          bool ifnd = false;
          // Distance to 1 decreases until it increases. 
          for(nnewp = 1; nnewp < 1000; nnewp++){
            // Length is 1 / (nnewp + 1)
            double err1 = abs(len/(nnewp+1) - 1);
            double err2 = abs(len/(nnewp+2) - 1);
            if(err2 > err1){
              ifnd = true;
              break;
            }
          }
          METRIS_ASSERT_MSG(ifnd, "Infinite length? "<<len);
          METRIS_ASSERT_MSG(len/(nnewp+1) > 1/sqrt(2),
            "with len = "<<len<<" got nnewp = "<<nnewp<<" each is "
            <<len/nnewp);
        }
        avg_nnewp += nnewp;
        max_nnewp = MAX(max_nnewp, nnewp);

        edg2pol[0] = ent2poi(ientt,lnoed[ied][0]);
        edg2pol[1] = ent2poi(ientt,lnoed[ied][1]);


        CPRINTF2(" - No CAD link for this edge -> eval at element\n");
        int idx0 = tdim + 1 + ied*(ideg-1);
        for(int ii = 0; ii < ideg-1; ii++) edg2pol[2+ii] = ent2poi[ientt][idx0+ii];


        // Note: no CAD use as insertEdge handles CAD eval. This routine only 
        // generates points approximately. 
        eval1<gdim,ideg>(msh.coord, edg2pol, msh.getBasis(), 
                         DifVar::None, DifVar::None, 
                         bar1, coop, NULL, NULL);

        if(DOPRINTS1()){
          CPRINTF2(" - enact ins ientt = %d ied = %d len = %f edg %d %d coord = ",
                 ientt,ied,len,edg2pol[0],edg2pol[1]);
          dblAr1(msh.idim,coop).print();
        }

        // Try insert point coop
        int nent00 = msh.nentt(tdim); 
        int itry = 0;
        do{
          if(DOPRINTS2()){
            writeMesh("preins",msh);
          }

          ierro = insertEdge(msh,tdim,ientt,ied,coop,bar1[0],
                             cav,work,lcaverr,ithrd2,ithrd3);

          if(ierro <= 0) break;
          itry++;
          if(itry >= 1 + imovmet) break;
          if(ierro == 0) CPRINTF2(" - After trying ierro = 0 \n");
          if(ierro > 0 && (itry == 0 && imovmet)){
            CPRINTF2(" -> insertEdge fail: try again w/ imovmet %d\n",imovmet);
            if(DOPRINTS1()){
              CPRINTF2(" - initial ipins = ");
              dblAr1(gdim,coop).print();
            }
            if(imovmet && itry == 1){
              CPRINTF2(" -> do imovmet\n");
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
              for(int ii = 0; ii < gdim; ii++) coop[ii] = 0.25 * msh.coord(edg2pol[0],ii)
                + 0.5 * coop[ii] + 0.25 * msh.coord(edg2pol[1],ii);
            }
            if(DOPRINTS1()){
              CPRINTF2(" - final ipins = ");
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
          CPRINTF2(" - insertion failed ierro = %d \n",ierro);
          linserr[ierro - 1] ++;
          nerro++;
        }


        if(msh.param->dbgfull) check_topo(msh,ithrd2);

        if(dobrk) break;
      }

    }

    avg_nnewp /= nlong;
    double t1 = get_wall_time();
    int ncallps = 1000*(int)((ninser2 / (t1-t0)) / 1000);
    CPRINTF2(" - Loop 1 end t = %f nlong %d ninser %d = %d /s; nerro %d; nnewp avg = %f, max = %d\n",
              t1-t0,nlong,ninser2,ncallps,nerro,avg_nnewp, max_nnewp);
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


    msh.cleanup();

    if(nedgt == 0) stat = MAX(stat, 0);
    else           stat = MAX(stat, (double)ninser2/(double)nedgt);

    *ninser += ninser2;
  }// for niter2


  return stat;
}
#endif



#if 0 
// Old routine only inserted one point at a time


// lpins work array sized dynamically by this routine ; it's an argument solely because this will be called several times, save on alloc
// also: as iterations go, fewer and fewer edges are long, no use allocating more than once to maximum needed size (first iter)
template<class MFT, int gdim, int ideg>
double insertLongEdges(Mesh<MFT> &msh, int *ninser, int ithrd1, int ithrd2, int ithrd3){
  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(ithrd1 != ithrd3);
  METRIS_ASSERT(ithrd2 != ithrd3);


  bool imovmet = true; 
  static int nwarnprt1 = 0;

  // Swap norm -1: length-based. 
  //swapOptions swapOpt(100, -1, 0.0);
  swapOptions swapOpt(100, 0, 0.005);

  //msh.met.setSpace(MetSpace::Log);
  msh.met.setSpace(MetSpace::Exp);

  CPRINTF2("-- insertLongEdges start \n");
  #ifndef NDEBUG
  CPRINTF2(" - Note: improve by generating several points per edge. Generated but not used cf loop nn/2 \n");
  CPRINTF2(" - Note: improve by filtering point propositions \n");
  #endif
  int edg2pol[getnnod1(ideg)];


  double stat = 0;

  const int tdim = msh.get_tdim();

  const int nedgl = (tdim*(tdim+1))/2;

  int lnoed1[1][2] = {{0, 1}};
  const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                             tdim == 2 ? lnoed2[0] : lnoed3[0]);
  intAr2  &ent2poi = msh.ent2poi(tdim);
  intAr2r &ent2tag = msh.ent2tag(tdim);
  intAr2  &ent2ent = msh.ent2ent(tdim);

  const int mcaverr = CAV_ERR_NERROR;
  intAr1 lcaverr(mcaverr);

  const int minserr = 100;
  intAr1 linserr(minserr);

  const int miter2 = 3;


  double bar1[2] = {0.5,0.5};
  double coop[gdim], t0;
  int ierro;
  *ninser = 0;

  msh.tag[ithrd1]++;
  
  // At most one edge in an insertion that doesn't collapse a point.
  MshCavity cav(100,100,1);
  CavWrkArrs work;

  int niter2 = 0;
  int ninser2 = 1; // otherwise doesnt enter loop
  for(niter2 = 0; niter2 < miter2 && ninser2 > 0; niter2++){
    INCVDEPTH(msh.param);
    ninser2 = 0;

    lcaverr.fill(0);
    linserr.fill(0);
    int nedgt = 0;
    int nerro = 0;
    int nent0 = msh.nentt(tdim);
    t0 = get_wall_time();
    for(int ientt = 0; ientt < nent0; ientt++){
      INCVDEPTH(msh.param);
      
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
        GETVDEPTH(msh.param);
        nedgt++;

        double sz[2];
        double len = getlenedg_geosz<MFT,gdim,ideg>(msh,ientt,tdim,ied,sz);
        CPRINTF2(" - try ientt = %d ied = %d len = %f \n",ientt,ied,len);
        if(len < sqrt(2)) continue;

        edg2pol[0] = ent2poi(ientt,lnoed[ied][0]);
        edg2pol[1] = ent2poi(ientt,lnoed[ied][1]);

        bool iuse_CAD = false;

        if(!iuse_CAD){
          CPRINTF2(" - No CAD link for this edge -> eval at element\n");
          int idx0 = tdim + 1 + ied*(ideg-1);
          for(int ii = 0; ii < ideg-1; ii++) edg2pol[2+ii] = ent2poi[ientt][idx0+ii];

          eval1<gdim,ideg>(msh.coord, edg2pol, msh.getBasis(), 
                           DifVar::None, DifVar::None, 
                           bar1, coop, NULL, NULL);
        }else{
          METRIS_THROW_MSG(TODOExcept(),"Implement useCAD case in msh_insert. "
            "Note that this is done in low_insert.cxx")
        }// if(!iuse_CAD)

        if(DOPRINTS1()){
          CPRINTF2(" - enact ins ientt = %d ied = %d len = %f edg %d %d  coord = ",
                 ientt,ied,len,edg2pol[0],edg2pol[1]);
          dblAr1(msh.idim,coop).print();
        }

        // Try insert point coop
        int nent00 = msh.nentt(tdim); 
        int itry = 0;
        do{
          if(DOPRINTS2()){
            writeMesh("preins",msh);
          }

          ierro = insertEdge(msh,tdim,ientt,ied,coop,bar1[0],
                             cav,work,lcaverr,ithrd2,ithrd3);

          if(ierro <= 0) break;
          itry++;
          if(itry >= 1 + imovmet) break;
          if(ierro == 0) CPRINTF2(" - After trying ierro = 0 \n");
          if(ierro > 0 && (itry == 0 && imovmet)){
            CPRINTF2(" -> insertEdge fail: try again w/ imovmet %d\n",imovmet);
            if(DOPRINTS1()){
              CPRINTF2(" - initial ipins = ");
              dblAr1(gdim,coop).print();
            }
            if(imovmet && itry == 1){
              CPRINTF2(" -> do imovmet\n");
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
              for(int ii = 0; ii < gdim; ii++) coop[ii] = 0.25 * msh.coord(edg2pol[0],ii)
                + 0.5 * coop[ii] + 0.25 * msh.coord(edg2pol[1],ii);
            }
            if(DOPRINTS1()){
              CPRINTF2(" - final ipins = ");
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
          CPRINTF2(" - insertion failed ierro = %d \n",ierro);
          linserr[ierro - 1] ++;
          nerro++;
        }


        if(msh.param->dbgfull) check_topo(msh,ithrd2);

        if(dobrk) break;
      }

    }

    double t1 = get_wall_time();
    int ncallps = 1000*(int)((ninser2 / (t1-t0)) / 1000);
    CPRINTF2(" - Loop 1 end t = %f ninser %d = %d /s; nerro %d\n",
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


    msh.cleanup();

    if(nedgt == 0) stat = MAX(stat, 0);
    else           stat = MAX(stat, (double)ninser2/(double)nedgt);

    *ninser += ninser2;
  }// for niter2


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