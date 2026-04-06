//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../cavity/msh_cavity.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../io_libmeshb.hxx"
#include "../utils/aux_pp_inc.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"
#include "../msh_checktopo.hxx"


namespace Metris{


void MshCavity::print(const MeshBase &msh, int iforce) const{
  GETVDEPTH(msh.param);

  if(!DOPRINTS1() && iforce <= 0) return;

  if(DOPRINTS1() || iforce >= 1){

    MPRINTF(" - cavity ipins {} pdim {} ncedg {} ncfac {} nctet {}\n",ipins,
            ipins >= 0 ? msh.getpoitdim(ipins) : -1,
            lcedg.get_n(),lcfac.get_n(),lctet.get_n());

    if(this->lcedg.get_n() > 0){
      MPRINTF(" - Edge cavity: {}\n",this->lcedg);
    }
    if(this->lcfac.get_n() > 0){
      MPRINTF(" - Face cavity: {}\n",this->lcfac);
    }
    if(this->lctet.get_n() > 0){
      MPRINTF(" - Tetra cavity: {}\n",this->lctet);
    }
  }

  if(DOPRINTS2() || iforce >= 2){
    for(int tdimn = 1; tdimn <= msh.get_tdim(); tdimn++){
      const intAr1 &lcent = this->lcent(tdimn);
      int ncent = lcent.get_n();
      if(ncent <= 0) continue;
      const intAr2 &ent2poi = msh.ent2poi(tdimn);

      if(tdimn == 1){
        MPRINTF(" - Edge cavity: \n");
      }else if(tdimn == 2){
        MPRINTF(" - Face cavity: \n");
      }else{
        MPRINTF(" - Tetra cavity: \n");
      }
      int nnode = msh.nnode(tdimn);
      for(int ientt : lcent){
        MPRINTF("  {} : {}\n",ientt,intAr1(nnode,ent2poi[ientt]));
      }
    }
  }
}

// -------------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------------- //
// For prints: "level 3" routine. level 0 = adapMesh level 1 = msh_collapse.. level 2 = low_collapse.. level 3 = cavity
// 2 spaces per level
template <class MFT, int ideg>
int cavity_operator(Mesh<MFT> &msh ,
                    MshCavity  &cav,
                    CavOprOpt  &opts,
                    CavWrkArrs &work,
                    CavOprInfo &info,
                    int ithread){

  METRIS_ASSERT_MSG(cav.inewp == 0 || cav.inewp == 1, "Caller must set cav.inewp to 0 if new point, 1 otherwise.");

  try{
  INCVDEPTH(msh.param);
  info.done = false;

  //static int nwarnprt = 0;
  //if(nwarnprt++ < 10){
  //  printf("## WARNING REMOVE THIS FILE WRITE\n");
  //  if(DOPRINTS3()) writeMesh("cavinifull",msh);
  //}

  METRIS_ENFORCE_MSG(opts.max_increase_cav_geo <= 1,"Implement cavity correction")

  CPRINTF1("-- START cavity_operator ncedg = {} ncfac = {} nctet = {} ipins = {} \n",
           cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n(),cav.ipins);
  CPRINTF1("   npoin {} nedge {} nface {} nelem {}\n",msh.npoin,msh.nedge,msh.nface,msh.nelem);

  cav.print(msh);

  if(DOPRINTS2()) writeMeshCavity("cavity0",msh,cav);
  if(DOPRINTS3()) writeMesh("cavinifull",msh,true);

	if(cav.ipins < 0 || cav.ipins >= msh.npoin)
		METRIS_THROW_MSG("ipins out of bounds\n");

	int ierro = CAV_NOERR;

	// In lbad, store: [2*i + 0] = entity rank in cavity
	//                 [2*i + 1] = entity type: 0 = edg, 1 = fac, 2 = tet
  intAr2 &lbad = work.lbad;

	int nbpo0 = msh.nbpoi,
	    npoi0 = msh.npoin,
	    nedg0 = msh.nedge,
	    nfac0 = msh.nface,
	    nele0 = msh.nelem;

  double qmax;

  // reset tag arrays in case is about to overflow
  if (msh.tag[ithread] >= INT_MAX - 1000000) {
    CPRINTF1("   tag value about to overflow... reseting tag arrays\n");
    msh.reset_tags();
  }

  cav.maxtag = msh.tag[ithread];

  /*  --------------   Correctness checks -------------------- */
	if((cav.lcedg.get_n() > 0 && msh.isboundary_edges())
  || (cav.lcfac.get_n() > 0 && msh.isboundary_faces()) ){
		int ibpoi = msh.poi2bpo[cav.ipins];
    if(ibpoi < 0){
      CPRINTF1("## ERROR CAV_ERR_NOBPO\n");
      ierro = CAV_ERR_NOBPO;
      goto cleanup;
    }
    int ityp = msh.bpo2ibi(ibpoi,1);
    if(cav.lcedg.get_n() > 0 && msh.isboundary_edges() && ityp > 1){
      CPRINTF1("## ERROR CAV_ERR_TDIMN\n");
      ierro = CAV_ERR_TDIMN;
      goto cleanup;
    }
	}

  // This only implements norempts for now.
	ierro = check_cavity_topo(msh, cav, opts, ithread);
	if(ierro > 0) goto cleanup;

  // The minimum value we need to set to (note right ++ is fine as we need only >=)
  cav.maxtag = MAX(cav.maxtag,++msh.tag[ithread]);

   /*  -------------- Generate final cavity --------------------
   		 For typent in (line|face|tetra) do
         Generate typent cavity + new typent-1 elements boundary
         Reconnect bdry to ipins
   */

	ierro = reconnect_lincav<MFT, ideg>(msh, cav, opts, ithread);
	if(ierro > 0) goto cleanup;
  CPRINTF1("-- reconnect_lincav done nedg0 = {} nedge = {} npoi0 = {} npoin = {}\n",
           nedg0,msh.nedge,npoi0,msh.npoin);

	ierro = reconnect_faccav<MFT, ideg>(msh, cav, opts, work, nedg0, &qmax, ithread);
  if(msh.get_tdim() == 2) info.qmax_end = qmax;
	if(ierro > 0) goto cleanup;
  CPRINTF1("-- reconnect_faccav done nfac0 = {} nface = {} \n",nfac0,msh.nface);

	ierro = reconnect_tetcav<MFT, ideg>(msh, cav, opts, info, nfac0, &qmax, ithread);
  if(ierro > 0) goto cleanup;
  if(msh.get_tdim() == 3) info.qmax_end = qmax;
  CPRINTF1("-- reconnect_tetcav done nele0 = {} nelem = {} \n",nele0,msh.nelem);

  // Update ibpois in case surface for new points (HO and ipins if new and dim < 2)
  try{
    ierro = update_bpois_newp(msh, cav, work, npoi0, nfac0, ithread);
  }catch(const MetrisExcept& e){
    writeMesh("cavenderr",msh,true,nedg0,nfac0,nele0);
    throw(e);
  }

	/*  -------------- Validity correction and quality checks ------------------
	  If this doesn't pass, don't invest on expensive optimization
	*/

	ierro = correct_cavity<MFT,ideg>(msh,cav,opts,npoi0,nedg0,nfac0,nele0,
                                   lbad,work,ithread);
  if(lbad.get_n() > 0) ierro = CAV_ERR_CORRECTCAV;
  if(ierro > 0) goto cleanup;

  if(opts.dryrun){
    if((opts.qmax_suf < 0 || qmax > opts.qmax_suf) &&
       (opts.qmax_iff < 0 || qmax > opts.qmax_iff)){
      ierro = CAV_ERR_DRYFAIL1;
      goto cleanup;
    }
  }
  if(opts.qmax_nec > 0 && qmax > opts.qmax_nec){
    CPRINTF1(" # specified qmax_nec = {} and got qmax = {} -> reject\n",
             opts.qmax_nec,qmax);
    ierro = CAV_ERR_DRYFAIL2;
    goto cleanup;
  }


	ierro = update_cavity<MFT,ideg>(msh, cav, work, npoi0, nedg0, nfac0, nele0, ithread);
	if(ierro > 0) goto cleanup;

  info.done = true;


  CPRINTF1("-- Cavity successful exit\n");
  //if(msh.nface >= 13326){
  //  printf("## DEBUG WAIT AT 13326 creation\n");
  //  writeMesh("13326",msh);
  //  wait();
  //}
	goto finish;


  //-------- cleanup (error case)
	cleanup:
  msh.tag[ithread] = cav.maxtag;
  CPRINTF1("-- Cavity error ierro = {} \n",ierro);
  if(DOPRINTS2()){
    // Write out the cavity.
    writeMesh("cavenderr",msh,true,nedg0,nfac0,nele0);
  }
	//METRIS_THROW_MSG(
  //  "Get rid of bpoi entries of existing points? Do these exist? Check ierro = "<<ierro);
  msh.tag[ithread]++;
  if(msh.isboundary_faces()){
    for(int iface = nfac0; iface < msh.nface; iface++){
      msh.fac2tag(ithread,iface) = msh.tag[ithread];
    }
  }
  if(msh.isboundary_edges()){
    for(int iedge = nedg0; iedge < msh.nedge; iedge++){
      msh.edg2tag(ithread,iedge) = msh.tag[ithread];
    }
  }
  for(int tdim = 1; tdim <= msh.get_tdim(); tdim++){
    int nent0 = tdim == 1 ? nedg0 : tdim == 2 ? nfac0 : nele0;
    for(int ientt = nent0; ientt < msh.nentt(tdim); ientt++){
      msh.ent2poi(tdim)(ientt,0) = -1;
    }
  }
  for(int ibpoi = nbpo0; ibpoi < msh.nbpoi; ibpoi++){
    INCVDEPTH(msh.param);
    int ip = msh.bpo2ibi(ibpoi,0);
    if( ip < 0 ) continue;
    if(msh.poi2tag(ithread,ip) >= msh.tag[ithread]) continue;
    CPRINTF1(" - remove ipoin {} bpois of new entities\n",ip);
    msh.poi2tag(ithread,ip) = msh.tag[ithread];
    msh.rembpotag(ip,ithread);
  }
	msh.set_nbpoi(nbpo0);
	msh.set_npoin(npoi0);
	msh.set_nedge(nedg0);
	msh.set_nface(nfac0);
	msh.set_nelem(nele0);

  cav.maxtag = MAX(cav.maxtag, msh.tag[ithread]);


  finish:

  msh.tag[ithread] = cav.maxtag;

  // cav.ipins may not be correctly in the mesh
  if(msh.param->dbgfull && ierro == 0) check_topo(msh,ithread);

  //static int nwarnprt = 0;
  //if(nwarnprt++ < 10) printf("## DEBUG REMOVE PRINT in msh_cavity\n");
  //if(msh.npoin == 26240 && ierro == 0){
  //  printf("## EXCEPTIONAL CHECK TOPO CALL\n");
  //  check_topo(msh,ithread);
  //}


  //static int nafter600 = 0;
  //if(msh.nface >= 600){
  //  nafter600++;
  //  PRINTF("## DEBUG 600 faces + {} op, write mesh\n",nafter600);
  //  writeMesh("iface600_" + std::to_string(nafter600),msh);
  //  if(nafter600 == 7) exit(0);
  //  //wait();
  //}

	return ierro;
  }catch(const MetrisExcept& e){
    printf("## EXCEPTION IN CAVITY \n");
    cav.print(msh,2);
    writeMeshCavity("cavity_except",msh,cav);
    throw(e);
  }
}
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int cavity_operator<MetricFieldAnalytical, n >(Mesh<MetricFieldAnalytical>&,MshCavity&,CavOprOpt&,CavWrkArrs&,CavOprInfo&,int);\
template int cavity_operator<MetricFieldFE        , n >(Mesh<MetricFieldFE        >&,MshCavity&,CavOprOpt&,CavWrkArrs&,CavOprInfo&,int);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


}// End namespace
