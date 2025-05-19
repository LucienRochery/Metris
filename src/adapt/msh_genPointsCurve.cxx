//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_lineadapt.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../utils/aux_misc.hxx"
#include "../low_lenedg.hxx"
#include "../io_libmeshb.hxx"
#include "../low_geo.hxx"
#include "../utils/mprintf.hxx"


namespace Metris{


// crvlen is the length of the CAD curve computed in the metric field. 
// icor0 is seed corner on edge 
// lnewt: t coord of new points
// ledge: seed edge in front mesh that t should lie on
template<class MFT>
void genPointsCurve(Mesh<MFT>& msh, int iref, int icor0, double crvlen, 
                    const double range[2], dblAr1 &lnewt, intAr1 &ledge){
  GETVDEPTH(msh.param);
  const double tessltar = 0.1;
  // Get CAD parameter range
  ego CADed = msh.CAD.cad2edg[iref];
  int ierro;

  int nnewp = MAX(round(crvlen) - 1, 2);
  double tarlen = crvlen / nnewp;

  CPRINTF1(" - edge ref %d cor0 %d range = (%f,%f), length %f, creating %d points for tarlen %f\n",
                        iref,icor0,range[0],range[1],crvlen,nnewp,tarlen);



  // Unlike the prior version, let's not go by physical length, but by parametric
  // length. Hence the mesh edges will be shorter than needed but, if ctrl pts
  // are added, they will get closer to 1 in length. 
  // We know the actual curve length, this decides how many points we put. 
  // While less perfectionist than previous version, this should also be a lot 
  // more robust. 


  /* -- Background tessellation */ 
  // To make things even more simple, we'll do a tessellation of the edge prior
  // Tessellation is binary such that all edges should be lower in length than 
  // tessltar. Start from the ts of this edge, we already have metrics there. 
  // tentative sizes:
  dblAr1 tspoi(3*nnewp);     // t coord
  dblAr1 tsmet(3*nnewp);     // metric size
  intAr2 tsedg(3*nnewp, 2);  // tess edge connectivity
  intAr2 ted2ed(3*nnewp, 2); // edge neighbours 
  intAr1 tsbke(3*nnewp);     // original mesh edge this came from 
  tspoi.set_n(1);
  tsmet.set_n(1);
  tsedg.set_n(0);
  ted2ed.set_n(0);
  tsbke.set_n(0);

  int ibcr0 = msh.poi2ebp(icor0,1,-1,iref);
  METRIS_ASSERT(ibcr0 >= 0);
  int iedg0 = msh.bpo2ibi(ibcr0,2);
  METRIS_ASSERT(iedg0 >= 0);

  if(abs(msh.bpo2rbi(ibcr0,0) - range[0]) < abs(msh.bpo2rbi(ibcr0,0) - range[1])){
    tspoi[0] = range[0];
  }else{
    tspoi[0] = range[1];
  }

  CPRINTF1(" - init tess point 0 t = %f range[0] = %f range[1] = %f\n",
           tspoi[0],range[0],range[1]);
  CPRINTF1(" - ibcr0 = %d from icor0 = %d has bpo2rbi = %f %f\n",
            ibcr0,icor0,msh.bpo2rbi(ibcr0,0),msh.bpo2rbi(ibcr0,1));


  { // don't pollute namespace with frequent variable names
    MetSpace ispac0 = msh.met.getSpace();
    msh.met.setSpace(MetSpace::Exp);

    int iedge = iedg0;
    int ipprv = icor0;
    double result[18];
    ierro = EG_evaluate(CADed, msh.bpo2rbi[ibcr0], result);
    METRIS_ASSERT(ierro == 0);
    // Computing the length of the tangent will later let us compute an edge
    // length by geo interp as in the usual case. 
    tsmet[0] = msh.idim == 2 ? getlenedg<2>(&result[3],msh.met[icor0])
                             : getlenedg<3>(&result[3],msh.met[icor0]);
    tsmet[0] = tsmet[0]*tsmet[0];

    while(true){
      INCVDEPTH(msh.param);
      int iver = msh.template getverent<1>(iedge, 1, ipprv);
      METRIS_ASSERT(iver >= 0);

      int ipnxt = msh.edg2poi(iedge, 1 - iver);
      int ibnxt = msh.poi2ebp(ipnxt, 1, iedge, iref);
      double tval = msh.bpo2rbi(ibnxt,0);

      int ntpoi = tspoi.get_n();
      tspoi.stack(tval);


      #ifndef NDEBUG
        // Verify the back edge ts contain the new t
        double tedg[2];
        for(int ii = 0; ii < 2; ii++){
          int ipoin = msh.edg2poi(iedge,ii);
          int ibpoi = msh.poi2ebp(ipoin,1,iedge,iref);
          METRIS_ASSERT(ibpoi >= 0);
          tedg[ii] = msh.bpo2rbi(ibpoi,0);
        }        
        bool inbd = tval >= tedg[0] && tval <= tedg[1]
                 || tval <= tedg[0] && tval >= tedg[1];
        if(!inbd) MPRINTF(" (1) New t %f not in edge bounds %f %f \n",tval, 
                          tedg[0],tedg[1]);
        METRIS_ASSERT_MSG(inbd,"Ini tess point t not in edge t bounds.");
      #endif

      int ntedg = tsedg.get_n();
      tsedg.set_n(ntedg+1);
      tsedg(ntedg,0) = ntpoi-1;
      tsedg(ntedg,1) = ntpoi;
      tsbke.stack(iedge);

      ted2ed.set_n(ntedg+1);
      ted2ed(ntedg,0) = -1;
      ted2ed(ntedg,1) = -1;
      if(ntedg - 1 >= 0){
        ted2ed(ntedg-1,0) = ntedg;
        ted2ed(ntedg  ,1) = ntedg-1;
      }


      CPRINTF1(" - new tess pt %d t = %f new edge %d = %d %d orig mesh edge %d ",ntpoi,tval,
               ntedg,ntpoi-1,ntpoi,iedge);
      #ifndef NDEBUG
      if(DOPRINTS1()) printf(" ts = %f %f \n",tedg[0],tedg[1]);
      #else
      if(DOPRINTS1()) printf("\n");
      #endif

      ierro = EG_evaluate(CADed, &tval, result);
      METRIS_ASSERT(ierro == 0);
      // len of curve g = \int_t \sqrt{g(t)^T M(g(t)) g'(t)}
      // \simeq \sum_i dt_i \int_{t=t_i}^{t_i+dt_i} \sqrt{g(t_i)^T M(g(t_i)) g'(t_i)}
      //        = \sum_i \sqrt{dt_i^2 M_i}
      // with M_i = g(t_i)^T M(g(t_i)) g'(t_i) = len_M(g'(t_i))^2 <- hence squaring. 
      double met1 = msh.idim == 2 ? getlenedg<2>(&result[3],msh.met[ipnxt])
                                  : getlenedg<3>(&result[3],msh.met[ipnxt]);
      met1 = met1*met1;
      tsmet.stack(met1);

      // Break at corner
      ibnxt = msh.poi2bpo[ipnxt];
      if(msh.bpo2ibi(ibnxt,1) == 0) break;

      ipprv = ipnxt;
      iedge = msh.edg2edg(iedge,iver);
    }// while true
    CPRINTF1(" - Phase 1 tess npoin = %d\n",tspoi.get_n());
    msh.met.setSpace(ispac0);

    // Now go over tessellation edges, compute lengths and split if needed
    int itsed = 0;
    // Store tessellated edges lengths sum to adjust tarlen
    double lentes = 0;
    while(itsed < tsedg.get_n()){
      INCVDEPTH(msh.param);
      int itpo1 = tsedg(itsed,0);
      int itpo2 = tsedg(itsed,1);
      double sz1 = tsmet[itpo1];
      double sz2 = tsmet[itpo2];

      METRIS_ASSERT(itpo1 >= 0);
      METRIS_ASSERT(itpo2 >= 0);
      METRIS_ASSERT(!isdeadent(tsbke[itsed],msh.edg2poi));


      double met1 = sqrt(sz1*sz2);

      double len = sqrt(met1) * abs(tspoi[itpo1] - tspoi[itpo2]);

      CPRINTF1(" - edge %d / %d verts %d %d sz %f %f len %f mesh edge %d \n",
               itsed,tsedg.get_n(),itpo1,itpo2,sz1,sz2,len,tsbke[itsed]);

      if(len <= tessltar){
        itsed++;
        // Once an edge is done for, add its length to the total. 
        lentes += len;
        continue;
      }

      // Split edge. Could be more clever about this, using sz1 and sz2
      double tnewp = 0.5*tspoi[itpo1] + 0.5*tspoi[itpo2];

      int ipnew = tspoi.get_n();
      tspoi.stack(tnewp);
      tsmet.stack(met1);


      // Current edge gets itpo1 -> ipnew
      tsedg(itsed,1) = ipnew;
      int ienew = tsedg.get_n();
      // a new edge gets ipnew -> itpo2
      tsedg.set_n(ienew+1);
      ted2ed.set_n(ienew+1);
      tsbke.stack(tsbke[itsed]);

      tsedg(ienew,0) = ipnew;
      tsedg(ienew,1) = itpo2;

      #ifndef NDEBUG
        // Verify the back edge ts contain the new t
        double tedg[2];
        int iseed = tsbke[itsed];
        for(int ii = 0; ii < 2; ii++){
          int ipoin = msh.edg2poi(iseed,ii);
          int ibpoi = msh.poi2ebp(ipoin,1,iseed,iref);
          METRIS_ASSERT(ibpoi >= 0);
          tedg[ii] = msh.bpo2rbi(ibpoi,0);
        }        
        bool inbd = tnewp >= tedg[0] && tnewp <= tedg[1]
                 || tnewp <= tedg[0] && tnewp >= tedg[1];
        if(!inbd){
          MPRINTF("New t %f not in edge bounds %f %f \n",tnewp, 
                          tedg[0],tedg[1]);
          MPRINTF("Edge is %d vertices %d %d \n",iseed,msh.edg2poi(iseed,0)
                  ,msh.edg2poi(iseed,1));
          MPRINTF("itpo1 = %d itpo2 = %d tspoi = %f %f\n",itpo1,itpo2,
                   tspoi[itpo1],tspoi[itpo2]);
          writeMesh("debug_t",msh);
        }
        METRIS_ASSERT_MSG(inbd,"New tess point t not in edge t bounds.");
      #endif

      // update neighbours
      // the "outside" neighbours
      int inei1 = ted2ed(itsed,0); // adjacent itpo2: opposite ipnew in new edge

      // old edge opposite unchanged point has new edge
      ted2ed(itsed,0) = ienew;
      // new edge opposite unchanged point has old edge 
      ted2ed(ienew,1) = itsed;

      //int inei2 = ted2ed(itsed,1); // adjacent itpo1: opposite ipnew in old edge 
      ted2ed(ienew,0) = inei1;

      // The left neighbour to itsed does not change, nor does itsed change
      // hence no update w inei2
      // ted2ed(itsed,1) = inei2;

      // But we do need to update for the right one 
      if(inei1 >= 0){
        if(ted2ed(inei1,0) == itsed){
          ted2ed(inei1,0) = ienew;
        }else if(ted2ed(inei1,1) == itsed){
          ted2ed(inei1,1) = ienew;
        }else{
          printf("## inei1 does not point to old edge!\n");
          printf("itsed %d inei1 %d \n",itsed,inei1);
          for(int itse2 = 0; itse2 < tsedg.get_n(); itse2++){
            printf(" %d : %d %d neighbours %d %d \n",itse2,tsedg(itse2,0),tsedg(itse2,1),
              ted2ed(itse2,0),ted2ed(itse2,1));
          }
          wait();
          METRIS_THROW(TopoExcept());
        }
      }


    }// while itsed

    double adjusted_tarlen = lentes / (MAX(round(lentes) - 1,1));
    CPRINTF1(" - Phase 2 tess npoin = %d new length %f tarlen %f -> %f \n",
             tspoi.get_n(),lentes,tarlen,adjusted_tarlen);
    tarlen = adjusted_tarlen;
  }


  // Next we use the background tessellation to produce points at distance 
  // approx 1 from each other.
  // Begin by recomputing the length as given by the tessellation, for consistency
  // and avoiding any bad surprises at the end. 
  // Replace tsmet by accumulated length. 

  // There is no non linearity in computing the edge lengths here:
  // no CAD updates, no metric updates. Hence we converge in one tarlen correction. 
  for(int niter = 0; niter < 2; niter++){

    double lentot = 0;
    int ipprv = 0; // 0 is always the leftmost point as we always leave the old 
    // edge left and only change its 1-th entry
    int itsed = 0;
    double lastlen = 0;
    lnewt.set_n(0);
    ledge.set_n(0);
    while(true){
      INCVDEPTH(msh.param);
      
      if(itsed < 0){
        lastlen = lentot;
        CPRINTF1("-- reached end: lastlen = %e recall range %f %f \n",lastlen,
          range[0],range[1]);
        break;
      }

      // Tessellation vertices and sizes
      int itpo1 = tsedg(itsed,0);
      int itpo2 = tsedg(itsed,1);
      double sz1 = tsmet[itpo1];
      double sz2 = tsmet[itpo2];
      int iseed  = tsbke[itsed];

      double met1 = sqrt(sz1 * sz2);
      double len = sqrt(met1) * abs(tspoi[itpo1] - tspoi[itpo2]);

      double lento0 = lentot;
      lentot += len;

      CPRINTF1(" - edge %d / %d points %d %d sz %f %f len %f tot %f tar %f\n",itsed,tsedg.get_n(),
                itpo1,itpo2,sz1,sz2,len,lentot,tarlen);

      // In any case, whether new pt or not, we'll need to walk to the next edge. 
      if(itpo1 == ipprv){
        ipprv = itpo2;
        itsed = ted2ed(itsed,0);
      }else if(itpo2 == ipprv){
        ipprv = itpo1;
        itsed = ted2ed(itsed,1);
      }else{
        printf("## FAILED TO find ipprv = %d ",ipprv);
        METRIS_THROW(TopoExcept())
      }


      // If we're yet to go above tarlen, continue to next edge. 
      if(lentot < tarlen*1.99) continue;


      // Otherwise, we just passed a threshold, create new point. 
      // theta st (1-theta)*lent0 + theta*lentot = tarlen
      double theta = (tarlen - lento0) / len;
      theta = MAX(0,MIN(1,theta)); // for close to machine zero errors
      //METRIS_ASSERT_MSG(theta >= 0 && theta <= 1,"theta out of bounds: "
      //  <<theta<<" diff above 1 "<<theta - 1);

      // Set point at (1-theta) itpo1 + theta itpo2
      double tnewp = (1-theta)*tspoi[itpo1] + theta*tspoi[itpo2];
      lnewt.stack(tnewp);
      ledge.stack(iseed);

      #ifndef NDEBUG
        double tedg[2];
        for(int ii = 0; ii < 2; ii++){
          int ipoin = msh.edg2poi(iseed,ii);
          int ibpoi = msh.poi2ebp(ipoin,1,iseed,iref);
          METRIS_ASSERT(ibpoi >= 0);
          tedg[ii] = msh.bpo2rbi(ibpoi,0);
        }
        CPRINTF1(" -> new point %d w theta = %f t = %f seed edg = %d, tedg %f %f \n",
                 lnewt.get_n()-1,theta,tnewp,iseed,tedg[0],tedg[1]);
        METRIS_ASSERT_MSG(tnewp >= tedg[0] && tnewp <= tedg[1]
                       || tnewp <= tedg[0] && tnewp >= tedg[1],
                       "New point t not in edge t bounds. Using tspoi: "
                       <<tspoi[itpo1]<<" "<<tspoi[itpo2]<<
                       "tess points "<<itpo1<<" "<<itpo2);
      #else
        CPRINTF1(" -> new point %d w theta = %f t = %f seed edg = %d\n",
                 lnewt.get_n()-1,theta,tnewp,iseed);
      #endif
      

      // reset lentot to the leftover length from this tessellation edge. 
      // that's lentot - tarlen
      lentot = lentot - tarlen;
    }//while true

    // We created nnewp1 + 1 edges 
    // the last one is length lastlen, when it should be tarlen 
    // Hence we need to adjust tarlen by adding (lastlen - tarlen) / (nnewp1 + 1)
    int nnewp1 = lnewt.get_n();
    double adjusted_tarlen = tarlen + (lastlen - tarlen) / (nnewp1 + 1);
    
    CPRINTF1(" - generated %d / %d points along curve last len %f "
           " tarlen %f -> %f \n",nnewp1,nnewp,lastlen,tarlen,adjusted_tarlen);
    if(DOPRINTS2() && msh.param->dbgfull){
      double result[18];
      ego obj = msh.CAD.cad2edg[iref]; 
      int npoi0 = msh.npoin;
      for(int inewt = 0; inewt < nnewp1; inewt++){
        double tcur = lnewt[inewt];
        METRIS_ENFORCE(EG_evaluate(obj, &tcur, result) == EGADS_SUCCESS);
        int ipnew = msh.newpoitopo(0, -1);
        msh.newbpotopo(ipnew,0,ipnew);
        for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipnew,ii) = result[ii];
      }
      CPRINTF2(" - First gen point %d t = %f \n",npoi0,lnewt[0]);
      writeMesh("genPoints_ref"+std::to_string(iref),msh);
      for(int ipoin = npoi0; ipoin < msh.npoin;ipoin++){
        msh.killpoint(ipoin);
      }
    }
    tarlen = adjusted_tarlen;
  }



}



template void genPointsCurve<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, 
  int iref, int icor0, double crvlen, const double range[2], dblAr1 &lnewt, intAr1& ledge);
template void genPointsCurve<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, 
  int iref, int icor0, double crvlen, const double range[2], dblAr1 &lnewt, intAr1& ledge);


} //namespace Metris
