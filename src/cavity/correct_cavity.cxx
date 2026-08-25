//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../cavity/msh_cavity.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../Mesh/Mesh.hxx"
#include "../low_geo/measure.hxx"
#include "../low_geo/normal.hxx"
#include "../low_geo/validity.hxx"
#include "../low_topo.hxx"
#include "../low_geo/ccoef.hxx"
#include "../aux_topo.hxx"
#include "../linalg/det.hxx"
#include "../io_libmeshb.hxx"

#include "../utils/aux_misc.hxx"
#include "../utils/CT_loop.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"

namespace Metris{

// lbad[i][0] element index
// lbad[i][1] element topo dim
template<class MFT, int ideg>
int correct_cavity(Mesh<MFT> &msh,
                   MshCavity &cav, 
                   CavOprOpt &opts, 
                   int npoi0, int nedg0, int nfac0, int nele0,
                   intAr2 &lbad, 
                   CavWrkArrs &work, 
                   int ithread){
  METRIS_ASSERT(lbad.get_stride() == 2);
  int iret; 

  CT_FOR0_INC(2,3,idim){if(msh.idim == idim){
    iret = correct_cavity0<MFT,idim,ideg>(msh,cav,opts,npoi0,nedg0,nfac0,nele0,lbad,work,ithread);
  }}CT_FOR1(idim);

  return iret;
}

template<class MFT, int gdim, int ideg>
int correct_cavity0(Mesh<MFT> &msh, 
                    [[maybe_unused]] MshCavity &cav, 
                    [[maybe_unused]] CavOprOpt &opts, 
                    int npoi0, int nedg0, int nfac0, int nele0, 
                    intAr2 &lbad, 
                    [[maybe_unused]] CavWrkArrs &work,
                    int ithread){
  METRIS_ASSERT(lbad.get_stride() == 2);
  GETVDEPTH(msh.param);
  
  lbad.set_n(0);


  //dblAr2 lmeas(2,mmeas,rwrk.allocate(2*mmeas));
  //dblAr2 &lmeas = work.lmeas;
  //lmeas.set_n(0); 
  //int mmeas = MAX(msh.nelem - nele0, msh.nface - nfac0);
  //lmeas.allocate(mmeas);

  double ccoef[getnnod3(gdim*(ideg-1))]; // Largest possible




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
  //    METRIS_ASSERT(*nbad < mbad);
  //    lbad[*nbad][0] = ientt;
  //    lbad[*nbad][1] = tdim;
  //    (*nbad)++;
  //
  //  }
  //}CT_FOR1(tdim);



  int ptag0 = msh.tag[ithread];
  if constexpr(ideg > 1){

    if(DOPRINTS2()){
      MshCavity cav2(msh.nelem-nele0,msh.nedge-nedg0,msh.nface-nfac0);
      for(int ii = nele0; ii < msh.nelem; ii++) cav2.lctet.stack(ii);
      for(int ii = nfac0; ii < msh.nface; ii++) cav2.lcfac.stack(ii);
      for(int ii = nedg0; ii < msh.nedge; ii++) cav2.lcedg.stack(ii);
      cav2.ipins = cav.ipins;
      writeMeshCavity("cavity1",msh,cav2);
      cav2.print(msh);
    }

    CPRINTF1("-- correct_cavity phase 1 : curve & project\n");
    // No HO curvature yet, just CAD projection
    if(!msh.CAD()){
      CPRINTF1(" - no CAD context, skip\n");
    }else{
      double result[18];

      // ... 
      METRIS_ASSERT(msh.getBasis() == FEBasis::Lagrange);

      ptag0++;
      cav.maxtag = MAX(cav.maxtag, ptag0);
      for(int tdim = 1; tdim <= 2; tdim++){
        if(tdim == 1 && !msh.isboundary_edges()) break; // Becasue we start with 1: 2 cannot be bdry then
        if(tdim == 2 && !msh.isboundary_faces()) break; // Because there's nothing after, but basically a continue

        int nent0 = tdim == 1 ? nedg0 : nfac0;
        int nentt = tdim == 1 ? msh.nedge : msh.nface; 
        intAr2 &ent2poi = msh.ent2poi(tdim);
        intAr1 &ent2ref = msh.ent2ref(tdim); 

        egoAr1 &cad2ent = tdim == 1 ? msh.CAD.cad2edg : msh.CAD.cad2fac;

        int nnode = tdim == 1 ? getnnod1(ideg) : getnnod2(ideg);


        for(int ientt = nent0; ientt < nentt; ientt++){

          int iref = ent2ref[ientt];
          METRIS_ASSERT(iref >= 0);

          ego obj = cad2ent[iref];

          //double nrm0 = getepsent<gdim>(msh,tdim,ientt);

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
            if(msh.poi2tag(ithread,ipoin) >= ptag0) continue;
            msh.poi2tag(ithread,ipoin) = ptag0;

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
            //      printf("## DEBUG high CAD gap = {} nrm0 = {}\n",err,nrm0);
            //      wait();
            //    #else
            //      return CAV_ERR_CADFAR;
            //    #endif
            //  }
            //  //METRIS_THROW_MSG( "Very large geometric gap? Manual check " << err);
            //}

            for(int ii = 0; ii < gdim; ii++) msh.coord(ipoin,ii) = result[ii];
          }

        } // end for ientt

      } // end for tdim


    } // end if egads context


    // Interpolate metric at new points 
    ptag0++;
    cav.maxtag = MAX(cav.maxtag, ptag0);
    for(int tdim = 1; tdim <= 3; tdim++){
      int nent0 = tdim == 1 ? nedg0 
                : tdim == 2 ? nfac0 : nele0;
      int nentt = msh.nentt(tdim);
      CPRINTF1(" - Update HO points dim {} entities {} <= i < {}\n",tdim,nent0,nentt);
      for(int ientt = nent0; ientt < nentt; ientt++){
        INCVDEPTH(msh.param);
        METRIS_ASSERT(!isdeadent(ientt,msh.ent2poi(tdim)));
        //int iref = msh.ent2ref(tdim)[ientt];
        for(int ii = tdim+1; ii < getnnode(tdim,ideg); ii++){
          INCVDEPTH(msh.param);
          int ipoin = msh.ent2poi(tdim)(ientt,ii);
          CPRINTF1(" - ipoin = {} tag = {} <? {}\n",ipoin,msh.poi2tag(ithread,ipoin),ptag0);
          if(ipoin < npoi0) continue;
          if(msh.poi2tag(ithread,ipoin) >= ptag0) continue;
          msh.poi2tag(ithread,ipoin) = ptag0;

          CPRINTF1("- update HO pt {} interp seed {} dim {} \n",ipoin,ientt,tdim);
          if(msh.interpMetBack(ipoin) != 0) return CAV_ERR_INTERPMETBACK;
        }// for ii = tdim+1
      }// for ientt
    }// for tdim 
  }// if ideg > 1



  CPRINTF1("-- correct_cavity phase {} : verify validity\n",1+ideg>1);

  //double quael;
  CT_FOR0_INC(2,gdim,tdim){

    double nrmal[3];

    int nent0 = tdim == 2 ? nfac0 : nele0;
    int nentt = msh.nentt(tdim);
    
    //const intAr2& ent2poi = msh.ent2poi(tdim);

    for(int ientt = nent0; ientt < nentt; ientt++){
      INCVDEPTH(msh.param);
      if constexpr(tdim == 2 && gdim == 3){
        getnorfacP1(msh.fac2poi[ientt], msh.coord, nrmal);
      }

      bool iflat = false;
      ElementValidityResult validity;
      if constexpr(ideg == 1){
        double meas = getmeasentP1<gdim,tdim>(msh, ientt, nrmal, &iflat);
        if(DOPRINTS1()){
          if constexpr (tdim == 2){
            CPRINTF1(" - {} tdim {} ientt {} meas {} iflat {} using normal {}\n",
                     ientt-nent0,tdim,ientt,meas,iflat,dblAr1(gdim,nrmal));
          }else{
            CPRINTF1(" - {} tdim {} ientt {} meas {} iflat {}\n",
                     ientt-nent0,tdim,ientt,meas,iflat);
          }
        }
      }else if constexpr(tdim == gdim){
        validity = classify_element_validity<gdim,ideg>(msh,ientt);
        iflat = !validity.accepted_conservatively();
      }else{
        // Embedded surface elements require their own oriented polynomial
        // certificate. Retain the historical compatibility check until that
        // separate Phase 3 implementation is available.
        getsclccoef<gdim,tdim,ideg>(msh, ientt, nrmal, ccoef, &iflat);
      }
      if(iflat){
        CPRINTF1(" - tdim {} ientt {} invalid\n",tdim,ientt);
        if constexpr(ideg > 1 && tdim == gdim){
          CPRINTF1(" - validity status {} lower bound {} coefficient {} "
                   "witness {} sample {}\n",
                   element_validity_status_name(validity.status),
                   validity.normalized_lower_bound,
                   validity.lower_bound_coefficient_index,
                   validity.normalized_witness,
                   validity.witness_sample_index);
        }else if constexpr(ideg > 1){
          constexpr int idegj = tdim*(ideg-1);
          constexpr int nnodj = getnnode(tdim,idegj);
          CPRINTF2(" - ccoef = {} \n",dblAr1(nnodj,ccoef));
        }
      }
      if(iflat) goto isbad;

      continue;
      
      isbad:
      int nbad = lbad.get_n();
      lbad.inc_n();
      lbad(nbad,0) = ientt;
      lbad(nbad,1) = tdim;
  
    }
  }CT_FOR1(tdim);
  return 0;

}


#define BOOST_PP_LOCAL_MACRO(n)\
template int correct_cavity<MetricFieldAnalytical,n>(Mesh<MetricFieldAnalytical> &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0, intAr2 &lbad,CavWrkArrs &work,int ithread);\
template int correct_cavity<MetricFieldFE        ,n>(Mesh<MetricFieldFE        > &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0,intAr2 &lbad,CavWrkArrs &work,int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


#define BOOST_PP_LOCAL_MACRO(n)\
template int correct_cavity0<MetricFieldAnalytical, 2, n>(Mesh<MetricFieldAnalytical> &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0, intAr2 &lbad, CavWrkArrs &work,int ithread);\
template int correct_cavity0<MetricFieldFE        , 2, n>(Mesh<MetricFieldFE        > &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0, intAr2 &lbad, CavWrkArrs &work,int ithread);\
template int correct_cavity0<MetricFieldAnalytical, 3, n>(Mesh<MetricFieldAnalytical> &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0, intAr2 &lbad, CavWrkArrs &work,int ithread);\
template int correct_cavity0<MetricFieldFE        , 3, n>(Mesh<MetricFieldFE        > &msh,\
                         MshCavity &cav, CavOprOpt &opts, \
                            int npoi0, int nedg0, int nfac0, int nele0, intAr2 &lbad, CavWrkArrs &work,int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


} // end namespace
