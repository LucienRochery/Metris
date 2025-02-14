//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../MetrisRunner/MetrisRunner.hxx"

#include "../Mesh/Mesh.hxx"

#include "../API/MetrisAPI.hxx"
#include "../msh_degelev.hxx"
#include "../aux_topo.hxx"
#include "../Boundary/msh_inisurf.hxx"
#include "../LPopt/msh_maxccoef.hxx"
#include "../low_ccoef.hxx"
#include "../BezierOffsets/msh_curve_offsets.hxx"
#include "../aux_utils.hxx"
#include "../aux_timer.hxx"
#include "../mprintf.hxx"

#include "../Localization/msh_localization.hxx"

#include "../io_libmeshb.hxx"

#include "../LPopt/msh_maxccoef.hxx"


namespace Metris{


MetrisRunner::~MetrisRunner(){
  delete msh_g; // virtual destructor -> calls the appropriate derived one
}

int MetrisRunner::degElevate(){
  if(param_.usrTarDeg <= msh_g->curdeg) return 0;
  if(this->metricFE){
    degElevate0<MetricFieldFE>();
  }else{
    degElevate0<MetricFieldAnalytical>();
  }
  return 1;
}

template<class MFT>
void MetrisRunner::degElevate0(){
  GETVDEPTH((*this));

  int ithread = 0;
  bool useOptim = param_.curveType == 3;

  Mesh<MFT> &msh = *( (Mesh<MFT>*) msh_g );

  if(param_.inpBack) METRIS_THROW_MSG(TODOExcept(), 
    "Degree elevation with input back not implemented");
  
  double t1 = get_wall_time(); 
  

  int ideg0 = msh.curdeg; 
  int nbpo0 = msh.nbpoi;
  int npoi0 = msh.npoin;

  CT_FOR0_EXC(1,METRIS_MAX_DEG,ideg){
    CT_FOR0_INC(ideg+1,METRIS_MAX_DEG,tdeg){
      if(ideg == ideg0 && tdeg == param_.usrTarDeg){
        CPRINTF1("-- Degree elevation %d -> %d \n",ideg,tdeg);  
        deg_elevate<MFT,ideg,tdeg>(msh);
      }
    }CT_FOR1(tdeg);
  }CT_FOR1(ideg);
 

  CPRINTF1("-- Back metric interpolation back deg = %d\n",bak.curdeg);
     
  CT_FOR0_INC(1,METRIS_MAX_DEG,bdeg){if(bak.curdeg == bdeg){
    // It's Lagrange nodes that should be localized.
    msh.setBasis(FEBasis::Lagrange);
    interpFrontBack<MFT,bdeg>(msh,bak,npoi0);
  }}CT_FOR1(bdeg);

  if(DOPRINTS2()) writeMesh("interpBack",msh);
  if(DOPRINTS2()) msh.met.writeMetricFile("interpBack");




  int tdim = msh.get_tdim();
  int nentt = msh.nentt(tdim);
  const intAr2& ent2poi = msh.ent2poi(tdim);
  int jdeg = tdim * (msh.curdeg - 1);
  int ncoef = tdim == 2 ? facnpps[jdeg]
                        : tetnpps[jdeg];
  int nnode = tdim == 2 ? facnpps[msh.curdeg]
                        : tetnpps[msh.curdeg];
  double ccoef[tetnpps[3*(METRIS_MAX_DEG-1)]];


  #if 0
  if(param_.inpBack)  METRIS_THROW_MSG(TODOExcept(),
      "Implement back mesh update in case of ext file after degelev");
  if(bak.nelem > 0){
    for(int ipoin = npoi0+1; ipoin < msh.npoin; ipoin++){
      int ielem = getpoitet(msh,ipoin);
      if(ielem < 0 || ielem >= msh.nelem) METRIS_THROW_MSG(TopoExcept(),
        "Failed to find back element for (HO) ipoin = "<<ipoin);
      msh.poi2bak(ipoin,3-1) = ielem;
    }
  }
  if(bak.nface > 0){
    for(int ipoin = npoi0+1; ipoin < msh.npoin; ipoin++){
      int iface = getpoifac(msh,ipoin);
      //if(iface < 0 || iface >= msh.nface) METRIS_THROW_MSG(TopoExcept(),
      //  "Failed to find back face for (HO) ipoin = "<<ipoin);
      msh.poi2bak(ipoin,2-1) = iface;
    }
  }
  if(bak.nedge > 0){
    for(int ipoin = npoi0+1; ipoin < msh.npoin; ipoin++){
      int iedge = getpoiedg(msh,ipoin);
      //if(iedge < 0 || iedge >= msh.nedge){
      //  printf("## FAILED TO getpoiedg for ipoin %d got iedge = %d \n",ipoin,iedge);
      //  printf(" poi2ent = %d %d\n",msh.poi2ent(ipoin,0),msh.poi2ent(ipoin,1));

      //  int pdim = msh.poi2ent(ipoin,1);
      //  int ientt = msh.poi2ent(ipoin,0);
      //  intAr2& ent2poi = msh.ent2poi(pdim);
      //  auto entnpps = ENTNPPS(2);
      //  printf("element : ");
      //  intAr1(entnpps[msh.curdeg],ent2poi[ientt]).print();

      //   METRIS_THROW_MSG(TopoExcept(),
      //  "Failed to find back edge for (HO) ipoin = "<<ipoin);
      //}
      msh.poi2bak(ipoin,1-1) = iedge;
    }
  }
  #endif


  const OptDoF idofs = OptDoF::HO;
  dblAr2 coor0;
  if(msh.CAD() || param_.curveType > 0){
    if(idofs == OptDoF::HO){
      coor0.allocate(msh.npoin-npoi0,msh.idim);
      coor0.set_n(msh.npoin-npoi0);
      for(int ipoin = npoi0; ipoin < msh.npoin; ipoin++){
        for(int ii = 0; ii < msh.idim; ii++)
          coor0(ipoin-npoi0,ii) = msh.coord(ipoin,ii);
      }
    }else{
      coor0.allocate(msh.npoin,msh.idim);
      coor0.set_n(msh.npoin);
      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        for(int ii = 0; ii < msh.idim; ii++)
          coor0(ipoin,ii) = msh.coord(ipoin,ii);
      }
    }
  }

  if(msh.CAD()){
    // It's Lagrange nodes that should be projected.
    msh.setBasis(FEBasis::Lagrange);
    prjMeshPoints(msh,nbpo0,true,true);
    if(DOPRINTS2()) writeMesh("prjMesh0", msh);
  }



  // However both the optimizer (n order to use d_ccoef) and BezOffsets need
  // Bézier 
  msh.setBasis(FEBasis::Bezier);

  if(param_.curveType > 0){

    CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
      CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
        curveMeshOffsets<MFT,gdim,ideg>(msh,false);
      }}CT_FOR1(ideg);
    }}CT_FOR1(gdim);

    if(DOPRINTS2()){
      writeMesh("crv0", msh);
      dblAr1 lminc(nentt);
      lminc.set_n(nentt);
      for(int ientt = 0; ientt < nentt; ientt++){
        if(isdeadent(ientt,ent2poi)) continue;
        CT_FOR0_INC(2,3,idim){if(idim == tdim){
          CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
            bool iinva;
            getsclccoef<idim,idim,ideg>(msh,ientt,NULL,ccoef,&iinva);
            lminc[ientt] = 1.0e30;
            for(int ii = 0; ii < ncoef; ii++){
              lminc[ientt] = MIN(lminc[ientt], ccoef[ii]);
            }
          }}CT_FOR1(ideg);
        }}CT_FOR1(idim);
      }// for ientt

      writeField("crv0", msh, SolTyp::P0Elt, lminc, 1);
    }

  }

  // First: curve as much as possible while retaining validity (backtrack)
  double qbktr = 0.9;
  double rcurv = 1;
  int niter = 0;
  bool ifirst = true;


  while(true){

    if(tdim != msh.idim){
      printf("## UNTESTED SURFACE CORRECTION \n");
      printf("## GET NORMAL \n");
      exit(1);
      wait();
    }

    int nflat = 0;
    msh.tag[ithread]++;
    for(int ientt = 0; ientt < nentt; ientt++){
      if(isdeadent(ientt, ent2poi)) continue;
      bool iflat = false;
      CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
        CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
          getsclccoef<gdim,gdim,ideg>(msh,ientt,NULL,ccoef,&iflat);
        }}CT_FOR1(ideg);
      }}CT_FOR1(gdim);
      if(iflat){
        nflat++;
        for(int ii = 0; ii < nnode; ii++){
          int ipoin = ent2poi(ientt,ii);
          msh.poi2tag(ithread, ipoin) = msh.tag[ithread];
        }
      }
    }

    CPRINTF1(" - backtrack iter %d ninva = %d\n",niter,nflat);

    if(nflat == 0) break;

    ifirst = false;

    if(idofs == OptDoF::HO){
      for(int ipoin = npoi0; ipoin < msh.npoin; ipoin++){
        if(msh.poi2tag(ithread,ipoin) < msh.tag[ithread]) continue;
        for(int ii = 0; ii < msh.idim; ii++)
          msh.coord(ipoin,ii) = qbktr*msh.coord(ipoin,ii) 
                              + (1.0 - qbktr)*coor0(ipoin-npoi0,ii) ;
      }
    }else{
      for(int ipoin = npoi0; ipoin < msh.npoin; ipoin++){
        if(msh.poi2tag(ithread,ipoin) < msh.tag[ithread]) continue;
        for(int ii = 0; ii < msh.idim; ii++)
          msh.coord(ipoin,ii) = qbktr*msh.coord(ipoin,ii) 
                              + (1.0 - qbktr)*coor0(ipoin,ii) ;
      }
    }

    rcurv *= qbktr;

    if(niter++ > 500){
      printf("## 500 BACKTRACK ITERATIONS ? \n");
      wait();
    }
  }

  if(ifirst && !useOptim){
    CPRINTF1(" - initial curvature valid : return\n");
    return;
  }else{
    CPRINTF1(" - backtracked factor %15.7e \n",rcurv);
    if(param_.curveType == 4){
      printf(" - Exiting here\n");
      return;
    }
  }

  if(useOptim){
    int opt_niter0 = param_.opt_niter;
    if(opt_niter0 == 0) param_.opt_niter = 20;
    optimMesh();
    param_.opt_niter = opt_niter0;
    return;
  }

  // Proceed to correction

  double tt0 = get_wall_time();
  if(msh.curdeg == 2){

    if(DOPRINTS2()) writeMesh("prjMesh", msh);

    if(idofs == OptDoF::HO){
      for(int ipoin = npoi0; ipoin < msh.npoin; ipoin++){
        for(int ii = 0; ii < msh.idim; ii++){
          double tmp = msh.coord(ipoin,ii);
          msh.coord(ipoin,ii) = coor0(ipoin-npoi0,ii);
          coor0(ipoin-npoi0,ii) = tmp;
        }
      }
    }else{
      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        for(int ii = 0; ii < msh.idim; ii++){
          double tmp = msh.coord(ipoin,ii);
          msh.coord(ipoin,ii) = coor0(ipoin,ii);
          coor0(ipoin,ii) = tmp;
        }
      }
    }
  }
  double tt1 = get_wall_time();
  CPRINTF1(" - Done time = %f\n",tt1-tt0);

  if(msh.curdeg == 2){

    intAr1 idx_point(msh.npoin);
    idx_point.set_n(msh.npoin);
    int npopt = 0; 

    if(idofs == OptDoF::HO){

      for(int ipoin = 0; ipoin < npoi0; ipoin++) idx_point[ipoin] = -1;

      for(int ipoin = npoi0; ipoin < msh.npoin; ipoin++)
        idx_point[ipoin] = ipoin - npoi0;

      npopt = msh.npoin - npoi0;

    }else{

      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        // Skip corners
        int ibpoi = msh.poi2bpo[ipoin];
        if(ibpoi >= 0){
          int itype = msh.bpo2ibi(ibpoi,1);
          if(itype == 0) continue;
          if(ipoin < npoi0) continue; 
        }
        idx_point[ipoin] = npopt;

        // ipoin >= npopt always 
        for(int ii = 0; ii < msh.idim; ii++)
          coor0(npopt,ii) = coor0(ipoin,ii);
        npopt++;
      }

    }
    coor0.set_n(npopt);


    CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      METRIS_ASSERT(idofs == OptDoF::HO && npopt == msh.npoin - npoi0
                   || idofs == OptDoF::Full);

      CT_FOR0_INC(2,3,idim){if(idim == msh.idim){
        bezGapsLP<idim, idim, ideg>(msh, idx_point, coor0, 
                                    LPMethod::IPM, LPLib::alglib);
      }}CT_FOR1(idim);
    }}CT_FOR1(ideg);


    //printf("debug coord ip 6 %f %f \n", msh.coord(6,0), msh.coord(6,1));


    #if 0
    double jtol = msh.param->jtol;
    const int miter = 10;
    const double qstr = 0.9;
    for(int niter = 0; niter < miter; niter++){
      double min_ccoef = maximizeCcoef<2,2,2>(msh,OptDoF::HO, LPMethod::IPM, LPLib::alglib);
      // double min_ccoef = maximizeCcoef<2>(msh,OptDoF::HO, LPMethod::IPM);
      if(min_ccoef >= jtol) break;
      // faces if surface 
      for(int iedge = 0; iedge < msh.nedge; iedge++){
        if(isdeadent(iedge,msh.edg2poi)) continue;
        int ipoi1 = msh.edg2poi(iedge,0);
        int ipoi2 = msh.edg2poi(iedge,1);
        int ipoih = msh.edg2poi(iedge,2);
        for(int ii = 0; ii < msh.idim; ii++){
          msh.coord(ipoih,ii) = qstr * msh.coord(ipoih,ii)
          + (1.0 - qstr)*(msh.coord(ipoi1,ii) + msh.coord(ipoi2,ii))/2.0;
        }
      }
    }
    #endif
  }else{
    METRIS_THROW(TODOExcept());
  }

  if(msh.param->iverb >= 1){
    FEBasis ibas0 = msh.getBasis();
    msh.setBasis(FEBasis::Lagrange);
    writeMesh("degelev.mesh", msh);
    msh.setBasis(ibas0);
  }

  double t2 = get_wall_time(); 
  CPRINTF1("-- Degree elevation time = %f\n",t2-t1);
}

template void MetrisRunner::degElevate0<MetricFieldFE>();
template void MetrisRunner::degElevate0<MetricFieldAnalytical>();




void MetrisRunner::writeOutputs(){
  if(this->metricFE){
    writeOutputs0<MetricFieldFE>();
  }else{
    writeOutputs0<MetricFieldAnalytical>();
  }
}

template<class MFT>
void MetrisRunner::writeOutputs0(){
  GETVDEPTH((*this));

  Mesh<MFT> &msh = static_cast<Mesh<MFT>&>(*msh_g);


  if(param_.wrtMesh){
    std::string baseOutName;
    std::string effMeshFileName, effMetFileName;
    auto pos = param_.outmFileName.find(".mesh");
    if(pos == std::string::npos){
      baseOutName = param_.outmFileName;
      effMeshFileName = param_.outmFileName + ".meshb";
      effMetFileName = baseOutName + ".solb";
    }else{
      effMeshFileName = param_.outmFileName;
      baseOutName = param_.outmFileName.substr(0,pos);
      effMetFileName = baseOutName + ".solb";
    }
    writeMesh(effMeshFileName, msh, param->main_in_prefix);
    msh.met.writeMetricFile(effMetFileName, param->main_in_prefix);
  }// End sequential
}
template void MetrisRunner::writeOutputs0<MetricFieldFE>();
template void MetrisRunner::writeOutputs0<MetricFieldAnalytical>();












}//End namespace