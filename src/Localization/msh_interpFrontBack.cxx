//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "msh_localization.hxx"
#include "msh_interpFrontBack.hxx"

#include "../msh_lag2bez.hxx"
#include "../utils/aux_misc.hxx"
#include "../linalg/det.hxx"
#include "../low_topo.hxx"
#include "../low_normal.hxx"
#include "../low_geo.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/aux_timer.hxx"
#include "../utils/mprintf.hxx"
#include "../io_libmeshb.hxx"
#include "../Boundary/low_projsurf.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

namespace Metris{




// Very rudimentary routine to be enhanced with a kd-tree in the future 
// Localize points, interpolate metric and store seed 
// Needs two slots in tag arrays 
// if ipoi0 > 0, it should be the number of P1 points ! These points already have seeds.
template<class MetricFieldType, int bdeg>
void interpFrontBack(Mesh<MetricFieldType> &msh, MeshBack &bak, int ipoi0){
  INCVDEPTH(msh.param);
  if(bak.getBasis() == FEBasis::Lagrange && bak.curdeg > 1) 
      METRIS_THROW_MSG(WArgExcept(), "Back should be in Bézier format!");

  //METRIS_ENFORCE_MSG(msh.idim == msh.get_tdim(), "Mesh is surface or line in plane.");

  METRIS_ENFORCE(METRIS_MAXTAGS >= 2);

  msh.setBasis(FEBasis::Lagrange);

  int gdim  = msh.idim;
  int nnmet = (gdim*(gdim+1))/2;
  int ierro;

  int ithread = 0;


  MetSpace ispac0 = msh.met.getSpace();
  FEBasis ibas0 = msh.met.getBasis();

  msh.met.forceSpaceFlag(MetSpace::Log);
  msh.met.forceBasisFlag(FEBasis::Lagrange);

  if(ispac0 == MetSpace::Exp){
    for(int ipoin = 0; ipoin < ipoi0; ipoin++){
      if(msh.idim == 2){
        getlogmet_inp<2,double>(msh.met[ipoin]);
      }else{
        getlogmet_inp<3,double>(msh.met[ipoin]);
      }
    }
  }
  if(ibas0 == FEBasis::Bezier && ipoi0 > 0){
    if(msh.idim == 2){
      setFieldLagrange<1,3>(msh,msh.met.rfld);
    }else{
      setFieldLagrange<1,6>(msh,msh.met.rfld);
    }
  }


  // Match corners quadratically 
  intAr1 lcorb(10), lcorf(10);
  lcorb.set_n(0);
  lcorf.set_n(0);
  for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
    if(msh.bpo2ibi(ibpoi,0) < 0 ) continue;
    if(msh.bpo2ibi(ibpoi,1) != 0) continue;
    int ipoin = msh.bpo2ibi(ibpoi,0);
    lcorf.stack(ipoin);
  }
  for(int ibpoi = 0; ibpoi < bak.nbpoi; ibpoi++){
    INCVDEPTH(msh.param);
    if(bak.bpo2ibi(ibpoi,0) < 0 ) continue;
    if(bak.bpo2ibi(ibpoi,1) != 0) continue;
    int ipoin = bak.bpo2ibi(ibpoi,0);
    lcorb.stack(ipoin);
    CPRINTF1(" - back corner point %d \n",ipoin);
  }
  METRIS_ASSERT_MSG(lcorb.get_n() == lcorf.get_n(),
    "Found back cor n = "<<lcorb.get_n()<<" front n = "<<lcorf.get_n());

  for(int ipoif : lcorf){
    INCVDEPTH(msh.param);
    double dmin = 1.0e30;
    int imin = -1;
    int ii = -1;
    for(int ipoib : lcorb){
      ii++;

      // Already matched
      if(ipoib < 0) continue;

      double dist = msh.idim == 2 ? geterrl2<2>(msh.coord[ipoif],bak.coord[ipoib])
                                  : geterrl2<3>(msh.coord[ipoif],bak.coord[ipoib]);
      if(dist < dmin){
        dmin = dist;
        imin = ii;
      }
    }
    int ipoib = lcorb[imin];
    CPRINTF1(" - Corner match %d (back) with %d (front) distance %e \n",
             ipoib,ipoif,sqrt(dmin));
    msh.poi2bak[ipoif] = ipoib;
    for(int ii = 0; ii < nnmet; ii++)
      msh.met(ipoif,ii) = bak.met(ipoib,ii);
    // Make unmatchable
    lcorb[imin] = -1;
  }
  CPRINTF1("-- Matched %d back <-> front corners\n",lcorf.get_n());








  // Will increase in size if needed 
  intAr1 lerro(100);
  lerro.set_n(0); 

  double t0 = get_wall_time();

  intAr1 lpfro;
  if(ipoi0 > 0){
    lpfro.allocate(ipoi0);
    lpfro.set_n(0);
    for(int ipoin = 0; ipoin < ipoi0; ipoin++){
      if(msh.poi2ent(ipoin,0) < 0) continue;
      lpfro.stack(ipoin);
    }
  }else{
    lpfro.allocate(lcorf.get_n());
    lcorf.copyTo(lpfro);
  }

  //if(ipoi0 != 0){
  //  // Case where we assume only HO points to initialize. We have valid seeds.
  //  for(int ipoin = ipoi0; ipoin < msh.npoin; ipoin++){

  //    INCVDEPTH(msh.param);

  //    if(msh.poi2ent(ipoin,0) < 0) continue;

  //    int pdim = msh.getpoitdim(ipoin);
  //    // Corner not to be localized but matched prior
  //    if(pdim == 0) continue;

  //    ierro = msh.interpMetBack(ipoin);
  //    if(ierro != 0) METRIS_THROW_MSG(TODOExcept(),"Implement bad stack to retry later");
  //  }
  //}else{
    // Case where no points (except corners) are seeded. 
    CPRINTF1("-- Frontal back -> front links\n");
    // Create a point front. Each point gathers ball but only of tdim >= its own, 
    // initializes neighbours' poi2bak using its own poi2bak and potentially edg2fac, fac2tet
    // Add thusly initialized neighbours to stack, pop self. 
    //intAr1 lpfro(lcorf.get_n());
    //lcorf.copyTo(lpfro);
    METRIS_ASSERT(lcorf.get_n() > 0);
    intAr1 lentt[3] = {10, 100, msh.idim == 3 ? 100 : 0};
    msh.tag[ithread]++;
    int ptag = msh.tag[ithread];
    if(DOPRINTS2()){
      CPRINTF2(" - Init front with %d points",lpfro.get_n());
      if(lpfro.get_n() <= 10){
        printf(": ");
        lpfro.print();
      }else{
        printf("\n");
      }
    }
    for(int ipoin : lpfro){
      msh.poi2tag(ithread,ipoin) = ptag;
    }
    // Initial front is successes.
    int nsucc = lpfro.get_n();
    while(lpfro.get_n() > 0){
      INCVDEPTH(msh.param);
      int ipseed = lpfro.pop();
      METRIS_ASSERT_MSG(ipseed >= 0 && ipseed < msh.npoin," Got ipseed = "<<ipseed);
      METRIS_ASSERT(msh.poi2ent(ipseed,0) >= 0);

      int psdim = msh.getpoitdim(ipseed);

      int iopen;
      {
      INCVDEPTH(msh.param)
      int ientt = msh.poi2ent(ipseed, 0);
      int tdime = msh.poi2ent(ipseed, 1);
      int iver = msh.getverent(ientt,tdime,ipseed);
      int iedl = -1;
      ierro = 0;
      lentt[0].set_n(0);
      lentt[1].set_n(0);
      lentt[2].set_n(0);
      if(iver < tdime+1){
        ierro = ball(msh, ipseed, lentt[0], lentt[1], lentt[2], &iopen, false, ithread);
      }else{
        bool doshell = false;
        if(msh.curdeg <= 2){
          METRIS_ASSERT(msh.curdeg == 2);
          int nppe = getnnod1(msh.curdeg) - 2;
          iedl = (iver - (tdime + 1)) / nppe;
          doshell = true;
        }else{

          int nppe = getnnod1(msh.curdeg) - 2;

          if(iver <= tdime + 1 + nppe * (tdime*(tdime+1))/2){
            iedl = (iver - (tdime + 1)) / nppe;
            doshell = true;
          }

          // Interior control point, unless volume, could be face...
          METRIS_THROW_MSG(TODOExcept(), "Implement P3+ case")
          lentt[tdime-1].stack(ientt);
        }
        if(doshell){
          METRIS_ASSERT(iedl >= 0);
          const int* lnoed = tdime == 1 ? lnoed1[0] :
                             tdime == 2 ? lnoed2[iedl] : lnoed3[iedl];
          int ipsed1 = msh.ent2poi(tdime)(ientt, lnoed[0]);
          int ipsed2 = msh.ent2poi(tdime)(ientt, lnoed[1]);
          shell(msh, ipsed1, ipsed2, tdime, ientt, lentt[0], lentt[1], lentt[2], &iopen);
        }
      }

      }
      METRIS_ASSERT(ierro == 0);

      CPRINTF1(" - front point %d ball nedge %d nface %d ntetr %d\n",ipseed,
               lentt[0].get_n(),lentt[1].get_n(),lentt[2].get_n());

      double algnd_[3];
      for(int tdim = 1; tdim <= msh.get_tdim(); tdim++){
        for(int ientt : lentt[tdim-1]){
          INCVDEPTH(msh.param);
          CPRINTF1(" - ball dim %d entity %d check for new points\n",tdim,ientt);
          METRIS_ASSERT_MSG(ientt >= 0 && ientt < msh.nentt(tdim),
            "front point "<<ipseed
            <<"\nseed edg "<<lentt[0]
            <<"\nseed fac "<<lentt[1]
            <<"\nseed tet "<<lentt[2] << "\n"
            <<"Error at tdim "<<tdim<<" ientt = "<<ientt);
          int nnode = getnnode(tdim,msh.curdeg);

          for(int inode = 0; inode < nnode; inode++){
            int ipoin = msh.ent2poi(tdim)(ientt,inode);
            if(ipoin < ipoi0) continue; // If we decide to use this with ipoi0 != 0
            if(msh.poi2tag(ithread,ipoin) >= ptag) continue;

            int pdim = msh.getpoitdim(ipoin);
            // Corner not to be localized but matched prior
            if(pdim == 0) continue;
            if(pdim < psdim){
              CPRINTF1(" - point dim %d < seed dim = %d -> skip\n",pdim, psdim);
              continue;
            }

            // Get the ref of the point. 
            int iref = msh.ent2ref(pdim)[msh.poi2ent(ipoin,0)];
            METRIS_ASSERT(iref >= 0);

            double *algnd = NULL;
            if(pdim < msh.idim){
              msh.get_algnd(ipoin, algnd_);
              algnd = algnd_;
            }


            ierro = msh.interpMetBack00(ipoin, pdim, iref, ipseed, algnd);
            if(ierro != 0){
              CPRINTF1(" # interpMetBack00 failed (bad seed), retry point later\n");
              continue;
            }

            nsucc++;

            CPRINTF1(" - added point %d dim %d to front\n",ipoin,pdim);
            lpfro.stack(ipoin);
            msh.poi2tag(ithread,ipoin) = ptag;
          }
        }
      }// for tdim

    }// while lpfro
    if(nsucc != msh.npoin){
      MPRINTF("## Failed %d points\n",msh.npoin - nsucc);
      for(int ipoin = ipoi0; ipoin < msh.npoin; ipoin++){
        if(msh.poi2tag(ithread,ipoin) >= ptag) continue;
        printf("Failed point %d \n",ipoin);
      }
      METRIS_THROW(GeomExcept());
    }
  //}

  double t1 = get_wall_time();
  CPRINTF1("-- Interp Back -> Front time %f pt/s %d nerror %d \n",t1-t0,
                                        (int)(msh.npoin/(t1-t0)),lerro.get_n());

  if(lerro.get_n() == 0) return;

  METRIS_THROW_MSG(TODOExcept(), "Error handling unchanged since poi2bak 2D");
  int tdim  = gdim; 

  intAr1 lball(100);
  int nnode = msh.nnode(tdim);
  intAr2 &ent2poi = msh.ent2poi(tdim);
  intAr2 &bak_ent2tag = bak.ent2tag(tdim);
  int nloop = 0;
  int nerro = 0;
  int nfix = 0;
  do{
    if(nloop++ > 10) METRIS_THROW(GeomExcept()); 

    nfix  = 0;
    nerro = 0;
    for(int ipoin : lerro){
      INCVDEPTH(msh.param);
      // Fixed previously 
      if(msh.poi2tag(0,ipoin) < msh.tag[0]) continue;
      // get ball and try using neighbours 
      int ientt = getpoient(msh,ipoin,tdim); 
      int iopen;
      bool imani;
      if(tdim == 2){
        intAr1 dum;
        ierro = ball2(msh,ipoin,ientt,lball,dum,&iopen,&imani,1);
      }else{
        ierro = ball3(msh,ipoin,ientt,lball,&iopen,1);
      }
      METRIS_ASSERT(ierro == 0);

      bak.tag[0]++;
      for(int iebal : lball){
        INCVDEPTH(msh.param);
        for(int ii = 0 ;ii < nnode; ii++){
          int ipoi2 = ent2poi(iebal,ii);
          // These are the points that have failed 
          if(msh.poi2tag(0,ipoi2) >= msh.tag[0]) continue;

          int ieleg = msh.poi2bak[ipoi2];
          // Skip any seeds that have been tried 
          if(bak_ent2tag(0,ieleg) >= bak.tag[0]) continue;

          double bary[4], coopr[3];
          for(int ii = 0; ii < tdim + 1; ii++)  bary[ii] = 1.0 / (tdim + 1);
          if(tdim == 2){
            // dummy tdim 
            METRIS_THROW_MSG(TODOExcept(), 
              "Error handling unchanged since poi2bak 2D");
            ierro = locMesh<2,2,bdeg>(bak,&ieleg,msh.coord[ipoin],
                                      msh.get_tdim(),NULL,-1,NULL,
                                      coopr,bary,1.0e-6,0,true);
          }else{
            // dummy tdim 
            METRIS_THROW_MSG(TODOExcept(), 
              "Error handling unchanged since poi2bak 2D");
            ierro = locMesh<3,2,bdeg>(bak,&ieleg,msh.coord[ipoin],
                                      msh.get_tdim(),NULL,-1,NULL,
                                      coopr,bary,1.0e-6,0,true);
          }

          if(ierro == 0){
            msh.poi2bak[ipoin] = ieleg;
            nfix++;
            msh.poi2tag(0,ipoin)--; // untag as invalid, could help a neighbour 
            int *ent2pol = tdim == 2 ? bak.fac2poi[ieleg] : bak.tet2poi[ieleg];
            bak.met.getMetBary(AsDeg::Pk,
                               DifVar::None,MetSpace::Log,
                               ent2pol,tdim,bary,msh.met[ipoin],NULL);
            goto nxtpoi;
          }
        }
      }

      nerro++;

      nxtpoi:
      continue;

    }

    double t2 = get_wall_time();
    printf("-- Interp Back -> Front phase 2 time %f nfix %d nerror %d \n",t2-t1,
      nfix, nerro);
  }while(nerro > 0 && nfix > 0);

  if(nerro == 0) return;

  printf("## ERROR EXIT DUMP ERROR POINTS \n");
  int ii = 0;
  for(int ipoin : lerro){
    if(msh.poi2tag(0,ipoin) < msh.tag[0]) continue;
    printf("%d : ipoin = %d \n",ii++,ipoin);
  }
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void interpFrontBack<MetricFieldAnalytical, n >(\
Mesh<MetricFieldAnalytical> &msh, MeshBack &bak, int ipoi0);\
template void interpFrontBack<MetricFieldFE        , n >(\
Mesh<MetricFieldFE        > &msh, MeshBack &bak, int ipoi0);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()






}