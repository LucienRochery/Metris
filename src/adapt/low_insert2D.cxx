//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../adapt/low_insert2D.hxx"
#include "../adapt/low_increasecav.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../msh_structs.hxx"
#include "../low_topo.hxx"
#include "../mprintf.hxx"
#include "../low_geo.hxx"
#include "../io_libmeshb.hxx"


namespace Metris{

// Return 0 if done nothing, 1 if error, -1 if done swap
// bar1 is t along the edge with 1 if lnoed[iedl][0]
template<class MetricFieldType>
int insedgesurf(Mesh<MetricFieldType>& msh, int iface, int iedl, double *coop, 
               double bar1, intAr1 &lerro, int ithrd1, int ithrd2){
  GETVDEPTH(msh);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);

  int iret = 0;

  bool isellen = true;

  if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(), "Implement + tet nelem = "<<msh.nelem)

  const int mcfac = 100, mcedg = 1; 
  const int nwork = mcfac + mcedg;
  RoutineWorkMemory ipool(msh.iwrkmem);
  int *iwork = ipool.allocate(nwork);
  MshCavity cav(0,mcfac,mcedg,nwork,iwork);

  CavOprOpt opts;
  CavOprInfo info;
  CavWrkArrs work;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = true;
  opts.dryrun = false;
  opts.allow_remove_points = false;

  int mcavcorr = 5, ncavcorr;


  CPRINTF1("-- START insedgesurf iface = %d ied %d\n",iface,iedl);

  int ierro = 0, nprem;


  int ifac2 = msh.fac2fac(iface,iedl);
 
  int iedge = -1; 
  if(ifac2 <= -1){ // geometric edge 
    int ip1 = msh.fac2poi(iface,lnoed2[iedl][0]);
    int ip2 = msh.fac2poi(iface,lnoed2[iedl][1]);
    iedge = getedgglo(msh,ip1,ip2);

    METRIS_ASSERT(iedge >= 0);
    METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));
    cav.lcedg.stack(iedge); 
  }


  if(ifac2 < -1){ // Non manifold 
    shell2_nm(msh,iface,iedl,cav.lcfac);
  }else if(ifac2 == -1){ // only iface
    cav.lcfac.stack(iface); 
  }else{ // 2 faces 
    cav.lcfac.stack(iface); 
    cav.lcfac.stack(ifac2); 
  }

  int ibpoi = -1;


  CPRINTF1(" - Initial npoin = %d \n",msh.npoin);

  // Create the point, set info for localization 
  int iseed, tdimp, iref;
  ego obj;
  double algnd[2];
  if(iedge >= 0){
    cav.ipins = msh.newpoitopo(1,iedge);
    ibpoi = msh.newbpotopo(cav.ipins,1,iedge);
    iseed = iedge;
    tdimp = 1;
    iref = msh.edg2ref[iedge];
    if(msh.CAD())
      obj  = msh.CAD.cad2edg[iref];
  }else{
    cav.ipins = msh.newpoitopo(2,iface);
    if(msh.isboundary_faces()){
      ibpoi = msh.newbpotopo(cav.ipins,2,iface);
    }
    iseed = iface;
    tdimp = 2;
    iref  = msh.fac2ref[iface];
    if(msh.isboundary_faces() && msh.CAD())
      obj   = msh.CAD.cad2fac[iref];
  }
  METRIS_ASSERT(iref >= 0);
  METRIS_ASSERT(obj != NULL || tdimp == 2 && !msh.isboundary_faces());

  CPRINTF1(" - create ipins %d \n",cav.ipins);

  for(int ii = 0; ii < msh.idim; ii++) msh.coord[cav.ipins][ii] = coop[ii];


  // Evaluate ipins on CAD, also get algnd for interpMetBack 
  int ip[2] = {msh.fac2poi(iface,lnoed2[iedl][0]),
               msh.fac2poi(iface,lnoed2[iedl][1])};
  if(ibpoi >= 0 && msh.CAD()){
    int ib[2];
    // Correct ibs : attach to ref or edge/face as needed
    for(int ii = 0; ii < 2; ii++){
      ib[ii] = msh.poi2ebp(ip[ii],tdimp,iseed,iref);
      METRIS_ASSERT(ib[ii] >= 0);
    }

    for(int ii = 0; ii < 2; ii++) msh.bpo2rbi(ibpoi,ii) = 
        bar1 * msh.bpo2rbi[ib[0]][ii] + (1.0 - bar1) * msh.bpo2rbi[ib[1]][ii];


    if(tdimp == 1){
      double result[18];
      ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
      if(ierro != 0){
        ierro = INS2D_ERR_EGEVALUATE; 
        iret = ierro;
        goto cleanup;
      }
      if(DOPRINTS2()){
        CPRINTF2("EG_evaluate orig = ");
        dblAr1(msh.idim,msh.coord[cav.ipins]).print();
        CPRINTF2("new = ");
        dblAr1(msh.idim,result).print();
      }
      for(int ii = 0; ii < msh.idim; ii++) msh.coord[cav.ipins][ii] = result[ii];
      for(int ii = 0; ii < msh.idim; ii++) algnd[ii] = result[3+ii];
    }else{
      if(msh.idim == 3) METRIS_THROW_MSG(TODOExcept(), "CAD eval in insert face")
    }

    // If the point moved away from the initial cavity, increase it. 
    ierro = increase_cavity2D(msh,msh.coord[cav.ipins],
                              opts,cav,ithrd1);

    if(DOPRINTS2()){
      writeMeshCavity("insert_cavity0.meshb", msh, cav, ithrd2);
      CPRINTF2("increase_cavity2D after EGADS failed ipins = %d \n",cav.ipins);
    }
    if(ierro != 0){
      ierro = INS2D_ERR_INCCAV2D2;
      goto cleanup;
    } 
  }else if(!msh.CAD()){ 
    // No reevaluation, but initialize algnd to edge tangent 
    double dum[2];
    double bary[2] = {bar1, 1.0 - bar1};
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      eval1<2,ideg>(msh.coord, ip,
                    msh.getBasis(), DifVar::Bary, DifVar::None,
                    bary, dum, algnd, NULL);
    }}CT_FOR1(ideg);
  }



  #ifndef NDEBUG
  try{
  #endif

  ierro = msh.interpMetBack(cav.ipins, tdimp, iseed, iref, algnd);
  if(ierro != 0){
    ierro = INS2D_ERR_INTERPMETBACK;
    goto cleanup;
  }

  #ifndef NDEBUG
  }catch(const MetrisExcept& e){
    printf("## DEBUG interpMetBack error %d \n",ierro);
    writeMesh("debug_interpMetBack.meshb",msh);
    printf("ipins = %d \n",cav.ipins);
    dblAr1(msh.idim,msh.coord[cav.ipins]).print();
    writeMeshCavity("error_insert_cavity."+std::to_string(ncavcorr)+".meshb", 
                    msh,cav, ithrd2);
    throw(e);
  }
  #endif

  ncavcorr = 0;
  do{
   
    if(DOPRINTS2()) writeMeshCavity("insert_cavity0."+std::to_string(ncavcorr), 
                                  msh,cav, ithrd2);
    CPRINTF1(" - initial cavity size %d \n",cav.lcfac.get_n()); 
   
    increase_cavity_Delaunay(msh, cav, cav.ipins, ithrd1);
   
    CPRINTF1(" - +del cavity size %d \n",cav.lcfac.get_n()); 

    if(isellen){
      nprem = increase_cavity_lenedg(msh,cav,opts,cav.ipins,ithrd1,ithrd2); 
      if(nprem < 0){
        ierro = INS2D_ERR_SHORTEDG;
      }
      CPRINTF1(" - +len cavity size %d nprem = %d\n", cav.lcfac.get_n(),nprem); 
    }
    int ierr2 = increase_cavity2D(msh,msh.coord[cav.ipins],
                                  opts,cav,ithrd1);
    CPRINTF1(" - +val cavity size %d \n",cav.lcfac.get_n()); 
    if(ierr2 > 0 && ierro <= 0) ierro = INS2D_ERR_INCCAV2D;

    if(ierro <= 0) break; 

    if(tdimp < msh.get_tdim()){
      ierro = INS2D_ERR_BDRYNOCORR; 
      CPRINTF1(" - Cannot correct boundary point in insert2D\n");
      goto cleanup;
    }

    // Try relocate ipins to cavity barycenter but only if not boundary
    double bary[3] = {1.0/3,1.0/3,1.0/3}; 
    METRIS_ASSERT(msh.idim == 2);
    double newp[2] = {0,0};
    double eval[2];
    double meast = 0;
    for(int iface : cav.lcfac){
      eval2<2,1>(msh.coord,msh.fac2poi[iface],msh.getBasis(),DifVar::None,DifVar::None,
                 bary,eval,NULL,NULL);
      bool iflat;
      double wt = getmeasentP1<2,2>(msh, msh.fac2poi[iface], cav.nrmal, &iflat);
      // For simply barycentre, use meas0
      // To skew towards the largest elements use meas0 * meas0
      for(int ii = 0; ii < msh.idim;ii++) newp[ii] += wt * eval[ii];
      meast += wt;
    }
    for(int ii = 0; ii < msh.idim;ii++){
      newp[ii] /= meast;
      msh.coord[cav.ipins][ii] = newp[ii];
    }

    if(DOPRINTS2()) writeMeshCavity("insert_cavity1."+std::to_string(ncavcorr)+".meshb", 
                                  msh,cav, ithrd2);

    // reinterp metric. This is always interior case, no need for ref of bdry dir
    ierr2 = msh.interpMetBack(cav.ipins,tdimp,iseed,-1,NULL);
    if(ierr2 != 0){
      ierro = INS2D_ERR_INTERPMETBACK;
      break; 
    }
  }while(ierro > 0 && ncavcorr++ < mcavcorr);
  if(ierro > 0) goto cleanup;


  // In this case let's even just take the face normal
  double nrmal[3];
  if(msh.idim == 3 && (!msh.CAD() || !msh.param->dbgfull)){
    getnorfacP1(msh.fac2poi[iface],msh.coord,nrmal);
    cav.nrmal = nrmal;
  }




  CPRINTF1(" insert ipins =  %d \n",cav.ipins);

  #ifndef NDEBUG
    if(DOPRINTS2()){
      int *refold = (int *) malloc(msh.nface*sizeof(int));
      METRIS_ENFORCE(refold != NULL);
      for(int ii = 0; ii < msh.nface; ii++){
        refold[ii] = msh.fac2ref[ii];
        //msh.fac2ref[ii] = 1 + (msh.fac2tag(ithrd1,ii) == msh.tag[ithrd1]);
        msh.fac2ref[ii] = 1;
      }

      for(int ii = 0; ii < cav.lcfac.get_n(); ii++){
        msh.fac2ref[cav.lcfac[ii]] = 2;
      }
      // Add a corner at ipoin 
      int ipoin = msh.newpoitopo(-1,-1);
      int ibpoi = msh.newbpotopo(ipoin,0,ipoin);
      for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipoin,ii) = msh.coord[cav.ipins][ii] ;
      writeMesh("debug_insert0.meshb",msh);
      for(int ii = 0; ii < nibi; ii++) msh.bpo2ibi(ibpoi,ii)  = -1;
      msh.set_npoin(msh.npoin-1);
      msh.set_nbpoi(msh.nbpoi-1);


      for(int ii = 0; ii < msh.nface; ii++){
        msh.fac2ref[ii] = refold[ii];
      }
      free(refold);

    }
  #endif

  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MetricFieldType,ideg>(msh,cav,opts,work,info,ithrd1);
  }}CT_FOR1(ideg);

  if(ierro > 0){
    lerro[ierro-1]++;
  }


  if(info.done){
    CPRINTF1("-- END insedgesurf ipins = %d  \n",cav.ipins);

    #ifndef NDEBUG
      if(DOPRINTS2()) writeMesh("debug_insert1.meshb",msh);
    #endif
    return -1; // Return did op
  }

  cleanup:
  iret = ierro;

  if(ibpoi >= 0){
    for(int ii = 0; ii < nibi; ii++) msh.bpo2ibi(ibpoi,ii) = -1;
    msh.set_nbpoi(msh.nbpoi - 1);
  }
  msh.set_npoin(msh.npoin - 1);


  return iret;
}



template int insedgesurf<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                         int iface, int iedl, double *coop, double bar1, 
                         intAr1 &lerro, int ithrd1, int ithrd2);
template int insedgesurf<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                         int iface, int iedl, double *coop, double bar1, 
                         intAr1 &lerro, int ithrd1, int ithrd2);









// Return 0 if done nothing, 1 if error, -1 if done swap
// bar1 is t along the edge with 1 if lnoed[iedl][0]
template<class MetricFieldType>
int insfacsurf(Mesh<MetricFieldType>& msh, int iface, double *coop, 
               intAr1 &lerro, int ithrd1, int ithrd2){
  GETVDEPTH(msh);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);

  int iret = 0;

  bool isellen = true;

  if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(), "Implement + tet nelem = "<<msh.nelem)

  const int mcfac = 100, mcedg = 1; 
  const int nwork = mcfac + mcedg;
  RoutineWorkMemory ipool(msh.iwrkmem);
  int *iwork = ipool.allocate(nwork);
  MshCavity cav(0,mcfac,mcedg,nwork,iwork);

  CavOprOpt opts;
  CavOprInfo info;
  CavWrkArrs work;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = true;
  opts.allow_remove_points = true;
  opts.dryrun = false;
  opts.allow_remove_points = false;

  int mcavcorr = 5, ncavcorr;

  const int tdim = 2;


  CPRINTF1("-- START insfacsurf iface = %d \n",iface);

  int ierro = 0, nprem;

  cav.lcfac.stack(iface); 
 
  cav.ipins = msh.newpoitopo(2,iface);
  int ibpoi = -1;
  if(msh.isboundary_faces()){
    ibpoi = msh.newbpotopo(cav.ipins,2,iface);
  }

  if(msh.idim == 3) METRIS_THROW_MSG(TODOExcept(), 
                                           "Implement iref and algnd 3D insfac")

  for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins,ii) = coop[ii];

  ierro = msh.interpMetBack(cav.ipins,tdim,iface,-1,NULL);
  if(ierro != 0){
    ierro = INS2D_ERR_INTERPMETBACK;
    goto cleanup;
  }

  ncavcorr = 0;
  do{
    if(DOPRINTS2())
      writeMeshCavity("insert_cavity0."+std::to_string(ncavcorr)+".meshb", 
                                    msh,cav, ithrd2);

    CPRINTF1(" - initial cavity size %d \n",cav.lcfac.get_n()); 
    increase_cavity_Delaunay(msh, cav, cav.ipins, ithrd1);
    CPRINTF1(" - +del cavity size %d \n",cav.lcfac.get_n()); 

    if(isellen){
      nprem = increase_cavity_lenedg(msh,cav,opts,cav.ipins,ithrd1,ithrd2); 
      if(nprem < 0){
        ierro = INS2D_ERR_SHORTEDG;
      }
      CPRINTF1(" - +len cavity size %d nprem = %d\n",cav.lcfac.get_n(),nprem); 
    }
    int ierr2 = increase_cavity2D(msh,msh.coord[cav.ipins],
                                  opts,cav,ithrd1);
    CPRINTF1(" - +val cavity size %d \n",cav.lcfac.get_n()); 
    if(ierr2 > 0 && ierro <= 0) ierro = INS2D_ERR_INCCAV2D;
    if(ierro <= 0) break; 

    // Try relocate ipins to cavity barycenter
    double bary[3] = {1.0/3,1.0/3,1.0/3}; 
    METRIS_ASSERT(msh.idim == 2);
    bool iflat;
    double newp[2] = {0,0};
    double eval[2];
    double meast = 0;
    for(int iface : cav.lcfac){
      eval2<2,1>(msh.coord,msh.fac2poi[iface],msh.getBasis(),DifVar::None,DifVar::None,
                 bary,eval,NULL,NULL);
      double meas0 = getmeasentP1<2,2>(msh, msh.fac2poi[iface], cav.nrmal, &iflat);
      // For simply barycentre, use meas0
      // To skew towards the largest elements use meas0 * meas0
      double wt = meas0;
      for(int ii = 0; ii < msh.idim;ii++){
        newp[ii] += wt * eval[ii];
      }
      meast += wt;
    }
    for(int ii = 0; ii < msh.idim;ii++){
      newp[ii] /= meast;
      msh.coord[cav.ipins][ii] = newp[ii];
    }

    if(DOPRINTS2()) writeMeshCavity("insert_cavity1."+std::to_string(ncavcorr)+".meshb", 
                                  msh,cav, ithrd2);

    // reinterp metric 
    ierr2 = msh.interpMetBack(cav.ipins,tdim,iface,-1,NULL);
    if(ierr2 != 0){
      ierro = INS2D_ERR_INTERPMETBACK;
      break; 
    }
  }while(ierro > 0 && ncavcorr++ < mcavcorr);
  if(ierro > 0) goto cleanup;


  // In this case let's even just take the face normal
  double nrmal[3];
  if(msh.idim == 3 && (!msh.CAD() || !msh.param->dbgfull)){
    getnorfacP1(msh.fac2poi[iface],msh.coord,nrmal);
    cav.nrmal = nrmal;
  }


  // CAD link initialization including normal
  if(ibpoi >= 0 && msh.CAD()){

    int iref = msh.fac2ref[iface];
    METRIS_ASSERT(iref >= 0);


    if(msh.isboundary_faces()){
      ego obj = msh.CAD.cad2fac[iref];
      METRIS_ASSERT(obj != NULL);
      // Correct ibs : attach to ref or edge/face as needed
      int ib[3];
      for(int ii = 0; ii < 3; ii++){
        int ip = msh.fac2poi(iface,ii);
        ib[ii] = msh.poi2ebp(ip,2,iface,iref);
        METRIS_ASSERT(ib[ii] >= 0);
      }

      for(int ii = 0; ii < nrbi; ii++) msh.bpo2rbi(ibpoi,ii) = 
                                                     msh.bpo2rbi[ib[0]][ii]/3.0  
                                                   + msh.bpo2rbi[ib[1]][ii]/3.0 
                                                   + msh.bpo2rbi[ib[2]][ii]/3.0;

      double result[18];
      ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
      if(ierro != 0){
        ierro = INS2D_ERR_EGEVALUATE; 
        iret = ierro;
        goto cleanup;
      }
      if(DOPRINTS2()){
        CPRINTF2("EG_evaluate orig = ");
        dblAr1(msh.idim,msh.coord[cav.ipins]).print();
        CPRINTF2("new = ");
        dblAr1(msh.idim,result).print();
      }
      for(int ii = 0; ii < msh.idim; ii++) msh.coord[cav.ipins][ii] = result[ii];
    }

    // Use normal
    if(msh.idim == 3){
      METRIS_THROW_MSG(TODOExcept(), "Implement CAD normal")
    }


    // If the point moved away from the initial cavity, increase it. 

    ierro = increase_cavity2D(msh,msh.coord[cav.ipins],
                              opts,cav,ithrd1);


    if(DOPRINTS2()){
      writeMeshCavity("insert_cavity0.meshb", msh, cav, ithrd2);
      printf("increase_cavity2D after EGADS failed ipins = %d \n",cav.ipins);
    }
    if(ierro != 0){
      ierro = INS2D_ERR_INCCAV2D2;
      goto cleanup;
    } 


  }



  CPRINTF1(" insert ipins =  %d \n",cav.ipins);

  #ifndef NDEBUG
    if(DOPRINTS2()){
      int *refold = (int *) malloc(msh.nface*sizeof(int));
      METRIS_ENFORCE(refold != NULL);
      for(int ii = 0; ii < msh.nface; ii++){
        refold[ii] = msh.fac2ref[ii];
        //msh.fac2ref[ii] = 1 + (msh.fac2tag(ithrd1,ii) == msh.tag[ithrd1]);
        msh.fac2ref[ii] = 1;
      }

      for(int ii = 0; ii < cav.lcfac.get_n(); ii++){
        msh.fac2ref[cav.lcfac[ii]] = 2;
      }
      // Add a corner at ipoin 
      int ipoin = msh.newpoitopo(-1,-1);
      int ibpoi = msh.newbpotopo(ipoin,0,ipoin);
      for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipoin,ii) = msh.coord[cav.ipins][ii] ;
      writeMesh("debug_insert0.meshb",msh);
      for(int ii = 0; ii < nibi; ii++) msh.bpo2ibi(ibpoi,ii)  = -1;
      msh.set_npoin(msh.npoin-1);
      msh.set_nbpoi(msh.nbpoi-1);


      for(int ii = 0; ii < msh.nface; ii++){
        msh.fac2ref[ii] = refold[ii];
      }
      free(refold);

    }
  #endif

  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MetricFieldType,ideg>(msh,cav,opts,work,info,ithrd1);
  }}CT_FOR1(ideg);

  if(ierro > 0){
    lerro[ierro-1]++;
  }


  if(info.done){
    CPRINTF1("-- END insfacsurf ipins = %d  \n",cav.ipins);
    #ifndef NDEBUG
      if(DOPRINTS2()) writeMesh("debug_insert1.meshb",msh);
    #endif
    return -1; // Return did op
  }

  cleanup:
  iret = ierro;

  if(ibpoi >= 0){
    for(int ii = 0; ii < nibi; ii++) msh.bpo2ibi(ibpoi,ii) = -1;
    msh.set_nbpoi(msh.nbpoi - 1);
  }
  msh.set_npoin(msh.npoin - 1);


  return iret;
}

template int insfacsurf<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                         int iface, double *coop, 
                         intAr1 &lerro, int ithrd1, int ithrd2);
template int insfacsurf<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                         int iface, double *coop, 
                         intAr1 &lerro, int ithrd1, int ithrd2);


} // end namespace
