//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../cavity/msh_cavity.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../Mesh/Mesh.hxx"
#include "../low_geo.hxx"
#include "../low_topo.hxx"
#include "../low_ccoef.hxx"
#include "../aux_topo.hxx"
#include "../aux_utils.hxx"
#include "../CT_loop.hxx"
#include "../mprintf.hxx"


namespace Metris{

// lbad[i][0] element index
// lbad[i][1] element topo dim
template<class MFT, int ideg>
int correct_cavity_fast(Mesh<MFT> &msh,
                        MshCavity &cav, CavOprOpt &opts, 
                        int npoi0, int nedg0, int nfac0, int nele0,
                        intAr2 &lbad, 
                        CavWrkArrs &work, 
                        int ithread){
  METRIS_ASSERT(lbad.get_stride() == 2);
  int iret; 

  CT_FOR0_INC(2,3,idim){if(msh.idim == idim){
    iret = correct_cavity_fast0<MFT,idim,ideg>(msh,cav,opts,npoi0,nedg0,nfac0,nele0,lbad,work,ithread);
  }}CT_FOR1(idim);

  return iret;
}

template<class MFT, int gdim, int ideg>
int correct_cavity_fast0(Mesh<MFT> &msh, 
                         MshCavity &cav, CavOprOpt &opts, 
                         int npoi0, int nedg0, int nfac0, int nele0, 
                         intAr2 &lbad, 
                         CavWrkArrs &work,
                         int ithread){
  METRIS_ASSERT(lbad.get_stride() == 2);
  GETVDEPTH(msh);
  
  lbad.set_n(0);


  //dblAr2 lmeas(2,mmeas,rwrk.allocate(2*mmeas));
  //dblAr2 &lmeas = work.lmeas;
  //lmeas.set_n(0); 
  //int mmeas = MAX(msh.nelem - nele0, msh.nface - nfac0);
  //lmeas.allocate(mmeas);

  double ccoef[tetnpps[gdim*(ideg-1)]]; // Largest possible




  //// first the quick checks for P1 
  //CT_FOR0_INC(2,gdim,tdim){
  //  const intAr2 &ent2poi = tdim == 2 ? msh.fac2poi : msh.tet2poi;

  //  int nent0 = tdim == 2 ? nfac0 : nele0;
  //  int nentt = tdim == 2 ? msh.nface : msh.nelem;

  //  for(int ientt = nent0; ientt < nentt; ientt++){

  //    if constexpr(tdim == 2 && gdim == 3){
  //      int iref = msh.fac2ref[ientt];
  //      getnorballref<1>(msh,cav.lcfac,iref,nrmal);
  //    }

  //    bool iflat;
  //    double meas0 = getmeasentP1<gdim,tdim>(ent2poi[ientt], msh.coord, msh.param->vtol, nrmal, &iflat);
  //    if(iflat) goto isbad;
  //    if constexpr(ideg >= 2 && gdim == tdim){
  //      getccoef<gdim,ideg>(msh,ientt,ccoef);
  //      for(int ii = 0; ii < ncoef; ii++){
  //        if(ccoef[ii] < msh.param->jtol * meas0) goto isbad;
  //      }
  //    }
  //    continue;
  //    
  //    isbad:
  //    if(*nbad >= mbad) METRIS_THROW(DMemExcept());
  //    lbad[*nbad][0] = ientt;
  //    lbad[*nbad][1] = tdim;
  //    (*nbad)++;
  //
  //  }
  //}CT_FOR1(tdim);



  if constexpr(ideg > 1){


    CPRINTF1("-- correct_cavity_fast phase 1 : curve & project\n");
    // No HO curvature yet, just CAD projection
    if(!msh.CAD()){
      CPRINTF1(" - no CAD context, skip\n");
    }else{
      double result[18];

      // ... 
      METRIS_ASSERT(msh.getBasis() == FEBasis::Lagrange);

      msh.tag[ithread]++;
      for(int tdim = 1; tdim <= 2; tdim++){
        if(tdim == 1 && !msh.isboundary_edges()) break; // Becasue we start with 1: 2 cannot be bdry then
        if(tdim == 2 && !msh.isboundary_faces()) break; // Because there's nothing after, but basically a continue

        int nent0 = tdim == 1 ? nedg0 : nfac0;
        int nentt = tdim == 1 ? msh.nedge : msh.nface; 
        intAr2 &ent2poi = msh.ent2poi(tdim);
        intAr1 &ent2ref = msh.ent2ref(tdim); 

        egoAr1 &cad2ent = tdim == 1 ? msh.CAD.cad2edg : msh.CAD.cad2fac;

        int nnode = tdim == 1 ? edgnpps[ideg] : facnpps[ideg];


        for(int ientt = nent0; ientt < nentt; ientt++){

          int iref = ent2ref[ientt];
          METRIS_ASSERT(iref >= 0);

          ego obj = cad2ent[iref];

          double nrm0 = getepsent<gdim>(msh,tdim,ientt);

          //int ip1 = msh.edg2poi(ientt,0);
          //int ip2 = msh.edg2poi(ientt,1);
          //int iref = msh.edg2ref[ientt];
          //METRIS_ASSERT(iref >= 0);
          //int ib1 = getpoiref2edgbpo(msh,ip1,ientt);
          //int ib2 = getpoiref2edgbpo(msh,ip1,ientt);
          
          // Start with something very simple. The uvs are already interpolated on point creation. Just evaluate
          for(int ii = tdim+1; ii < nnode; ii++){
            int ipoin = ent2poi(ientt,ii);
            // Obviously only update new points !!
            if(ipoin < npoi0) continue;

            // This not only avoids duplicates within one tdim, but ensures only the lowest tdim is 
            // responsible for the evaluate. 
            if(msh.poi2tag(ithread,ipoin) >= msh.tag[ithread]) continue;
            msh.poi2tag(ithread,ipoin) = msh.tag[ithread];

            METRIS_ASSERT(msh.bpo2ibi(msh.poi2bpo[ipoin],1) == tdim); // Actually using mark this should be true. 
            int ibpoi = msh.poi2ebp(ipoin,tdim,ientt,-1);
            METRIS_ASSERT(ibpoi >= 0);

            int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
            if(ierro != 0) return 100 + ierro;

            //double err = geterrl2<gdim>(msh.coord[ipoin],result);

            //// Technically equivalent to an assert, but not really an assert (bound to be removed)
            //if(err > nrm0){
            //  // Try inv_evaluate. 
            //  ierro = EG_invEvaluate(obj, msh.coord[ipoin], msh.bpo2rbi[ibpoi], result);
            //  if(ierro != 0) return 100 + ierro;
            //  err = geterrl2<gdim>(msh.coord[ipoin],result);
            //  if(err > nrm0){
            //    #ifndef NDEBUG
            //      printf("## DEBUG high CAD gap = %f nrm0 = %f\n",err,nrm0);
            //      wait();
            //    #else
            //      return CAV_ERR_CADFAR;
            //    #endif
            //  }
            //  //METRIS_THROW_MSG(GeomExcept(), "Very large geometric gap? Manual check " << err);
            //}

            for(int ii = 0; ii < gdim; ii++) msh.coord(ipoin,ii) = result[ii];
          }

        } // end for ientt

      } // end for tdim


    } // end if egads context


    // Interpolate metric at new points 
    double algnd_[3];
    double *algnd;
    msh.tag[ithread]++;
    for(int tdim = 1; tdim <= 2; tdim++){
      const intAr2 &ent2poi = msh.ent2poi(tdim); 
      const intAr1 &ent2ref = msh.ent2ref(tdim);
      int nnode = tdim == 1 ? edgnpps[ideg] : facnpps[ideg];
      int nent0 = tdim == 1 ? nedg0 : nfac0;
      int nentt = msh.nentt(tdim);

      for(int ientt = nent0; ientt < nentt; ientt++){
        INCVDEPTH(msh);
        METRIS_ASSERT(!isdeadent(ientt,ent2poi));
        int iref = ent2ref[ientt];
        for(int ii = tdim+1; ii < nnode; ii++){
          int ipoin = ent2poi(ientt,ii);
          if(ipoin < npoi0) continue;
          if(msh.poi2tag(ithread,ipoin) >= msh.tag[ithread]) continue;
          msh.poi2tag(ithread,ipoin) = msh.tag[ithread];

          //if(tdim == 1) METRIS_THROW_MSG(TODOExcept(), "algnd in cav reinterp");

          CPRINTF1("- update HO pt %d interp seed %d dim %d \n",ipoin,ientt,tdim);
          algnd = NULL;
          if(tdim < msh.get_tdim() && msh.CAD()){
            int ibpoi = msh.poi2ebp(ipoin,tdim,ientt,-1);
            METRIS_ASSERT(ibpoi >= 0);

            ego obj = tdim == 1 ? msh.CAD.cad2edg[iref] : msh.CAD.cad2fac[iref];

            double result[18];
            int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
            if(ierro != 0){
              algnd = NULL;
              CPRINTF1("# EG_eval failed\n");
            }else{
              for(int ii = 0; ii < msh.idim; ii++) algnd_[ii] = result[3+ii];
              algnd = algnd_;
            }
          }
          msh.interpMetBack(ipoin, tdim, ientt, iref, algnd);
        }
      }// for ientt
    }// for tdim 
  }



  CPRINTF1("-- correct_cavity_fast phase %d : verify validity\n",1+ideg>1);

  //double quael;
  CT_FOR0_INC(2,gdim,tdim){

    double nrmal[3];

    int nent0 = tdim == 2 ? nfac0 : nele0;
    int nentt = tdim == 2 ? msh.nface : msh.nelem;

    const intAr2& ent2poi = msh.ent2poi(tdim);

    for(int ientt = nent0; ientt < nentt; ientt++){
      INCVDEPTH(msh);
      if constexpr(tdim == 2 && gdim == 3){
        getnorfacP1(msh.fac2poi[ientt], msh.coord, nrmal);
        //if(msh.CAD()){
        //  getnorfacCAD(msh, ientt, nrmal);
        //}else{
        //  getnorfacP1(msh.fac2poi[ientt], msh.coord, nrmal);
        //}
      }

      bool iflat;
      if constexpr(ideg == 1){
        double meas = getmeasentP1<gdim,tdim>(msh, ent2poi[ientt], nrmal, &iflat);
        if(DOPRINTS1()){
          CPRINTF1(" - %d tdim %d ientt %d meas %f iflat %d using normal ",
                   ientt-nent0,tdim,ientt,meas,iflat);
          dblAr1(gdim,nrmal).print();
        }
        if(DOPRINTS2()){
          if(iflat || meas < 0){
            getnorfacP1(msh.fac2poi[ientt], msh.coord, nrmal);
            normalize_vec<3>(nrmal);
            CPRINTF2(" - discrete normal = ");
            dblAr1(gdim,nrmal).print();
          }
        }
      }else{
        getsclccoef<gdim,tdim,ideg>(msh, ientt, nrmal, ccoef, &iflat);
      }
      if(iflat) goto isbad;

      continue;
      
      isbad:
      int nbad = lbad.get_n();
      lbad.inc_n();
      lbad[nbad][0] = ientt;
      lbad[nbad][1] = tdim;
  
    }
  }CT_FOR1(tdim);
  return 0;

}


#define BOOST_PP_LOCAL_MACRO(n)\
template int correct_cavity_fast<MetricFieldAnalytical,n>(Mesh<MetricFieldAnalytical> &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0, intAr2 &lbad,CavWrkArrs &work,int ithread);\
template int correct_cavity_fast<MetricFieldFE        ,n>(Mesh<MetricFieldFE        > &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0,intAr2 &lbad,CavWrkArrs &work,int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


#define BOOST_PP_LOCAL_MACRO(n)\
template int correct_cavity_fast0<MetricFieldAnalytical, 2, n>(Mesh<MetricFieldAnalytical> &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0, intAr2 &lbad, CavWrkArrs &work,int ithread);\
template int correct_cavity_fast0<MetricFieldFE        , 2, n>(Mesh<MetricFieldFE        > &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0, intAr2 &lbad, CavWrkArrs &work,int ithread);\
template int correct_cavity_fast0<MetricFieldAnalytical, 3, n>(Mesh<MetricFieldAnalytical> &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0, intAr2 &lbad, CavWrkArrs &work,int ithread);\
template int correct_cavity_fast0<MetricFieldFE        , 3, n>(Mesh<MetricFieldFE        > &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0, intAr2 &lbad, CavWrkArrs &work,int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


} // end namespace
