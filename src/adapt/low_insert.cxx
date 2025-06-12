//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_insert.hxx"
#include "low_increasecav.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../utils/mprintf.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../msh_structs.hxx"
#include "../low_topo.hxx"
#include "../low_normal.hxx"
#include "../low_geo.hxx"
#include "../io_libmeshb.hxx"
#include "../linalg/det.hxx"
#include "../low_lenedg.hxx"


namespace Metris{

// Return 0 if done nothing, 1 if error, -1 if done swap
// bar1 is t along the edge with 1 if lnoed[iedl][0]
template<class MFT>
int insertEdge(Mesh<MFT>& msh, 
               int tdim, int ientt, int iedl, 
               double *coop, double bar1, 
               MshCavity &cav, CavWrkArrs &work, 
               intAr1 &lerro, int ithrd1, int ithrd2){

  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(!isdeadent(ientt,msh.ent2poi(tdim)));

  const int nedgl = (tdim*(tdim+1))/2;
  METRIS_ASSERT(iedl >= 0 && iedl < nedgl);

  int nentt = msh.nentt(tdim);

  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const auto lnoed = tdim == 2 ? lnoed2 : lnoed3;
  METRIS_ASSERT(ientt >= 0 && ientt < nentt && !isdeadent(ientt, ent2poi));

  //if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(), "Implement + tet nelem = "<<msh.nelem)

  CavOprOpt opts;
  CavOprInfo info;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = false;
  opts.dryrun = false;
  opts.allow_remove_points = false; // good for an infinite loop
  opts.allow_remove_points_superdim = true; // For boundary
  opts.qmax_nec = -1;
  opts.qmax_suf = -1;
  opts.qmax_iff = -1;

  int mcavcorr = 5, ncavcorr;

  cav.lcedg.allocate(10);
  cav.lcfac.allocate(10);
  cav.lctet.allocate(10);
  cav.lcedg.set_n(0);
  cav.lcfac.set_n(0);
  cav.lctet.set_n(0);

  CPRINTF1("-- START insertEdge tdim = %d ientt = %d ied %d\n",tdim,ientt,iedl);

  int ierro = 0, nprem;
  int ip1 = ent2poi(ientt,lnoed[iedl][0]);
  int ip2 = ent2poi(ientt,lnoed[iedl][1]);

  int iopen;
  shell(msh,ip1,ip2,tdim,ientt,cav.lcedg,cav.lcfac,cav.lctet,&iopen);
  CPRINTF1(" - cavity seed nedge %d nface %d ntetr %d\n",
           cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());


  METRIS_ASSERT(cav.lcedg.get_n() > 0 || cav.lcfac.get_n() > 0 || cav.lctet.get_n() > 0);
  #ifndef NDEBUG
  if(msh.get_tdim() == 2) METRIS_ASSERT(cav.lcfac.get_n() > 0);
  if(msh.get_tdim() == 3) METRIS_ASSERT(cav.lctet.get_n() > 0);
  if(msh.param->dbgfull){
    for(int ientt : cav.lcent(msh.get_tdim())){
      METRIS_ASSERT_MSG(!isdeadent(ientt,msh.ent2poi(msh.get_tdim())),
        "Shell element dead");
    }
  }
  #endif  


  int tdimp = -1;
       if(cav.lcedg.get_n() > 0) tdimp = 1;
  else if(cav.lcfac.get_n() > 0) tdimp = 2;
  else                           tdimp = 3;

  int ibins = -1;

  // Create the point, set info for localization 
  int iseed, iref;
  ego obj = NULL;
  double algnd[3];
  if(tdimp == 1){
    int iedge = cav.lcedg[0];
    METRIS_ASSERT(iedge >= 0);
    cav.ipins = msh.newpoitopo(1,iedge);
    ibins = msh.newbpotopo(cav.ipins,1,iedge);
    iseed = iedge;
    iref = msh.edg2ref[iedge];
    if(msh.CAD()) obj  = msh.CAD.cad2edg[iref];
  }else if(tdimp == 2){
    int iface = cav.lcfac[0];
    METRIS_ASSERT(iface >= 0);
    cav.ipins = msh.newpoitopo(2,iface);
    iseed = iface;
    iref  = msh.fac2ref[iface];
    if(msh.isboundary_faces()){
      ibins = msh.newbpotopo(cav.ipins,2,iface);
      if(msh.CAD()) obj = msh.CAD.cad2fac[iref];
    }
  }else{
    cav.ipins = msh.newpoitopo(3,ientt);
    iseed = ientt;
    iref = msh.tet2ref[ientt];
  }
  if(msh.CAD()) METRIS_ASSERT(obj != NULL 
                    || tdimp == 2 && !msh.isboundary_faces() || tdimp == 3);

  CPRINTF1(" - create ipins %d tdim = %d seed %d ref %d\n",cav.ipins,tdimp,iseed,iref);

  for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins,ii) = coop[ii];


  // Evaluate ipins on CAD, also get algnd for interpMetBack 
  int ip[2] = {ent2poi(ientt,lnoed[iedl][0]),
               ent2poi(ientt,lnoed[iedl][1])};
  if(ibins >= 0 && msh.CAD()){
    int ib[2];
    // Correct ibs : attach to ref or edge/face as needed
    for(int ii = 0; ii < 2; ii++){
      ib[ii] = msh.poi2ebp(ip[ii],tdimp,iseed,iref);
      METRIS_ASSERT(ib[ii] >= 0);
    }

    for(int ii = 0; ii < 2; ii++) msh.bpo2rbi(ibins,ii) = 
        bar1*msh.bpo2rbi[ib[0]][ii] + (1.0 - bar1)*msh.bpo2rbi[ib[1]][ii];

    double result[18];
    METRIS_ASSERT(obj != NULL);
    ierro = EG_evaluate(obj, msh.bpo2rbi[ibins], result);
    if(ierro != 0){
      ierro = INS2D_ERR_EGEVALUATE; 
      goto cleanup;
    }
    if(DOPRINTS2()){
      CPRINTF2("EG_evaluate orig = ");
      dblAr1(msh.idim,msh.coord[cav.ipins]).print();
      CPRINTF2("new = ");
      dblAr1(msh.idim,result).print();
    }
    for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins,ii) = result[ii];

    if(tdimp == 1){
      for(int ii = 0; ii < msh.idim; ii++) algnd[ii] = result[3+ii];
    }else{
      vecprod(&result[3], &result[6], algnd);
    }

    //// If the point moved away from the initial cavity, increase it. 
    //ierro = increase_cavity_validity(msh,cav,ithrd1);

    //if(DOPRINTS2()){
    //  writeMeshCavity("insert_cavity0.meshb", msh, cav);
    //  CPRINTF2("increase_cavity_validity after EGADS  ierro %d ipins = %d \n",ierro,cav.ipins);
    //}
    //if(ierro != 0){
    //  ierro = INS2D_ERR_INCCAV2D2;
    //  goto cleanup;
    //} 
  }else if(ibins >= 0 && !msh.CAD()){ 
    METRIS_ASSERT(tdimp <= 2);
    // No reevaluation, but initialize algnd to edge tangent 
    CPRINTF1(" - discrete algnd initialization tdimp %d \n",tdimp);
    if(tdimp == 1){
      double dum[2];
      double bary[2] = {bar1, 1.0 - bar1};
      // To compute at higher degree, copy more vertices into ip
      //CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
        eval1<2,1>(msh.coord, ip,
                   msh.getBasis(), DifVar::Bary, DifVar::None,
                   bary, dum, algnd, NULL);
      //}}CT_FOR1(ideg);
    }else if(msh.idim == 3){
      getnorfacP1(ent2poi[ientt],msh.coord,algnd);
    }
  }

  try{

  if(!msh.param->ins_lazy_interp){
    ierro = msh.interpMetBack(cav.ipins, tdimp, iseed, iref, algnd);
    if(ierro != 0){
      //printf("debug interpMetBack error %d\n",ierro);
      //wait();
      ierro = INS2D_ERR_INTERPMETBACK;
      goto cleanup;
    }
  }

  }catch(const MetrisExcept& e){

    printf("Exception in interpMetBack from insertEdge\n");
    printf(" insertion dim was %d\n",tdimp);
    printf(" using iseed iref %d %d initial seed ientt %d tdim %d\n",
           iseed,iref,ientt,tdim);
    printf("ipins = %d\n",cav.ipins);
    MPRINTF("-- cavity ncedg = %d ncfac = %d nctet = %d ipins = %d \n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n(),cav.ipins);
    MPRINTF("   npoin %d nedge %d nface %d nelem %d\n",msh.npoin,msh.nedge,msh.nface,msh.nelem);
    if(cav.lcedg.get_n() > 0){
      MPRINTF(" - Edge cavity: ");
      cav.lcedg.print();
    }
    if(cav.lcfac.get_n() > 0){
      MPRINTF(" - Face cavity: ");
      cav.lcfac.print();
    }
    if(cav.lctet.get_n() > 0){
      MPRINTF(" - Tetra cavity: ");
      cav.lctet.print();
    }

    for(int tdimn = 1; tdimn <= 3; tdimn++){
      intAr1 &lcent = cav.lcent(tdimn);
      int ncent = lcent.get_n();
      if(ncent <= 0) continue;
      intAr2 &ent2poi = msh.ent2poi(tdimn);

      if(tdimn == 1){
        MPRINTF(" - Edge cavity: \n");
      }else if(tdimn == 2){
        MPRINTF(" - Face cavity: \n");
      }else{
        MPRINTF(" - Tetra cavity: \n");
      }
      int nnode = msh.nnode(tdimn);
      for(int ientt : lcent){
        MPRINTF("%d : ",ientt);
        for(int ii = 0; ii < nnode; ii++){
          printf(" %d ",ent2poi(ientt,ii));
        }
        printf("\n");
      }
    }
    writeMeshCavity("inscavity1",msh,cav);
    throw(e);
  }

  ncavcorr = 0;
  ierro = 0;
  do{

    int nced0 = cav.lcedg.get_n();
    int ncfa0 = cav.lcfac.get_n();
    int ncte0 = cav.lctet.get_n();
   
    if(DOPRINTS2()) writeMeshCavity("insert_cavity0."+std::to_string(ncavcorr), 
                                  msh,cav);
    CPRINTF1(" - initial cavity nedge %d nface %d nelem %d\n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());

    nprem = increase_cavity_lenedg(msh,cav,opts,cav.ipins,ithrd1,ithrd2);
    if(nprem < 0){
      ierro = INS2D_ERR_SHORTEDG;
      CPRINTF1(" - with ops.allow_remove_points == false, short edge would be created\n");
      goto fixpoint;
    }
    CPRINTF1(" - +len cavity size %d nprem = %d\n", cav.lcfac.get_n(),nprem); 
    if(DOPRINTS2()){
      writeMeshCavity("insert_cavity1."+std::to_string(ncavcorr), 
                                  msh,cav);
    }


    ierro = increase_cavity(msh, cav, true, ithrd1, ithrd2);
    if(ierro != 0){
      CPRINTF1(" - +cav error %d\n",ierro);
      ierro = INS2D_ERR_INCCAV2D;
    }
    CPRINTF1(" - +cav nedge %d nface %d nelem %d\n",
             cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n());
    if(DOPRINTS2()){
      writeMeshCavity("insert_cavity2."+std::to_string(ncavcorr), 
                                  msh,cav);
    }
   
    if(ierro <= 0) break;

    fixpoint:

    if(ibins >= 0){
      ierro = INS2D_ERR_BDRYNOCORR; 
      CPRINTF1(" - Cannot correct boundary point in insertEdge\n");
      if(msh.param->interactive){
        printf("## WAIT HERE INS2D_ERR_BDRYNOCORR\n");
        wait();
      }
      goto cleanup;
    }

    METRIS_ASSERT(tdim == msh.idim);
    // Try relocate ipins to cavity barycenter using only lowest dim and CAD if
    // available. 

    cav.lcedg.set_n(nced0);
    cav.lcfac.set_n(ncfa0);
    cav.lctet.set_n(ncte0);

    int ierr2 = aux_movePointCav(msh, cav, tdimp, iseed, iref, algnd);

    if(ierr2 != 0){
      ierro = INS2D_ERR_INTERPMETBACK;
      CPRINTF1(" - Failed to move point in insertEdge\n");
      break;
    }

    if(DOPRINTS2()) writeMeshCavity("insert_cavity3."+std::to_string(ncavcorr)+".meshb", 
                                    msh,cav);

  }while(ierro > 0 && ncavcorr++ < mcavcorr);
  if(ierro > 0) goto cleanup;


  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd1);
  }}CT_FOR1(ideg);

  if(ierro > 0) lerro[ierro-1]++;
  #ifndef NDEBUG
  if(ierro == CAV_ERR_NOBPO){
    printf("## CAV_ERR_NOBPO error \n");
    wait();
  }
  #endif

  if(ierro != 0) ierro = INS2D_ERR_CAVITYOPERATOR;

  if(info.done){

    if(msh.param->ins_lazy_interp){
      ierro = msh.interpMetBack(cav.ipins, tdimp, iseed, iref, algnd);
      if(ierro != 0){
        printf("debug interpMetBack error %d\n",ierro);
        wait();
        //ierro = INS2D_ERR_INTERPMETBACK;
        //goto cleanup;
      }
      METRIS_ENFORCE(ierro == 0);
    }

    CPRINTF1("-- END insertEdge ipins = %d  \n",cav.ipins);
    #ifndef NDEBUG
      if(DOPRINTS2()) writeMesh("debug_insert1.meshb",msh);
    #endif
    return -1; // Return did op
  }

  cleanup:
  msh.killpoint(cav.ipins);

  return ierro;
}



template int insertEdge<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                         int tdim, int ientt, int iedl, double *coop, double bar1, 
                         MshCavity &cav, CavWrkArrs &work, 
                         intAr1 &lerro, int ithrd1, int ithrd2);
template int insertEdge<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                         int tdim, int ientt, int iedl, double *coop, double bar1, 
                         MshCavity &cav, CavWrkArrs &work, 
                         intAr1 &lerro, int ithrd1, int ithrd2);






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

  if(!msh.param->ins_lazy_interp){
    // reinterp metric. This is always interior case, no need for ref of bdry dir
    ierro = msh.interpMetBack(cav.ipins,tdimp,iseed,iref,algnd);
    if(ierro != 0){
      CPRINTF1(" - interpMetBack failed ierro = %d \n",ierro);
      ierro = INS2D_ERR_INTERPMETBACK;
    }
  }
  
  return ierro;
}

















} // end namespace
