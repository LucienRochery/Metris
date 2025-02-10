//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../cavity/msh_cavity.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_pp_inc.hxx"
#include "../mprintf.hxx"


namespace Metris{

// The boundary here is a set of edges. They will be reconnected to 
// ipins. 
// Note: if the initial cavity includes edges, these have been reconnected
// and the result is in lnewed
// Triangles, likewise. 
template <class MFT, int ideg>
int reconnect_tetcav([[maybe_unused]]Mesh<MFT> &msh, 
                     [[maybe_unused]] const MshCavity& cav, 
                     [[maybe_unused]] CavOprOpt &opts, 
                     [[maybe_unused]] double *qmax, 
                     [[maybe_unused]] int ithread){
  METRIS_THROW_MSG(TODOExcept(), "Implement reconnect_tetcav");
	return 0;
}


// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int reconnect_tetcav<MetricFieldAnalytical, n >(Mesh<MetricFieldAnalytical> &msh, const MshCavity& cav, CavOprOpt &opts, double *qmax, int ithread);\
template int reconnect_tetcav<MetricFieldFE        , n >(Mesh<MetricFieldFE        > &msh, const MshCavity& cav, CavOprOpt &opts, double *qmax, int ithread); 
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

// -------------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------------- //
// For prints: "level 3" routine. level 0 = adapMesh level 1 = msh_collapse.. level 2 = low_collapse.. level 3 = cavity 
// 2 spaces per level 
template <class MFT, int ideg>
int cavity_operator(Mesh<MFT> &msh , 
                   MshCavity  &cav,
                   CavOprOpt  &opts  ,
                    CavWrkArrs &work  ,
                    CavOprInfo &info  ,
                   int ithread){
  GETVDEPTH(msh);
  info.done = false;

  METRIS_ENFORCE_MSG(opts.max_increase_cav_geo <= 1, 
                     "Implement cavity correction")

  CPRINTF1("-- cavity_operator start ncedg = %d ncfac = %d nctet = %d ipins = %d \n",
    cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n(),cav.ipins);
  
  if(DOPRINTS1()){
    if(cav.lcedg.get_n() > 0){
      CPRINTF1(" - Edge cavity: ");
      cav.lcedg.print(cav.lcedg.get_n());
    }
    if(cav.lcfac.get_n() > 0){
      CPRINTF1(" - Face cavity: ");
      cav.lcfac.print(cav.lcfac.get_n());
    }
    if(cav.lctet.get_n() > 0){
      CPRINTF1(" - Tetra cavity: ");
      cav.lctet.print(cav.lctet.get_n());
    }
  }

  if(DOPRINTS2()){
    for(int tdimn = 1; tdimn <= 3; tdimn++){
      intAr1 &lcent = cav.lcent(tdimn);
      int ncent = lcent.get_n();
      if(ncent <= 0) continue;
      intAr2 &ent2poi = msh.ent2poi(tdimn);

      if(tdimn == 1){
        CPRINTF2(" - Edge cavity: \n");
      }else if(tdimn == 2){
        CPRINTF2(" - Face cavity: \n");
      }else{
        CPRINTF2(" - Tetra cavity: \n");
      }
      int nnode = msh.nnode(tdimn);
      for(int ientt : lcent){
        CPRINTF2("%d : ",ientt);
        for(int ii = 0; ii < nnode; ii++){
          printf(" %d ",ent2poi(ientt,ii));
        }
        printf("\n");
      }
    }
    writeMeshCavity("cavity0",msh,cav);
  }

	if(cav.ipins < 0 || cav.ipins >= msh.npoin) 
		METRIS_THROW_MSG(WArgExcept(),"ipins out of bounds\n");

	// The outcav must only store elements that will be put in the mesh
	// It will not include everything coming from gen_lcav, gen_faccav, gen_tetcav.
	// Those are stored separately and are only used to generate new elements. 
	//MshCavity outcav(cav.lctet.size(),cav.lcfac.size(),cav.lcedg.size());
	int ierro = CAV_NOERR;
	//int mbad = 100, nbad = 0;
	//RoutineWorkMemory<int> iwrk(msh.iwrkmem);
	//intAr2 lbad(mbad,2,iwrk.allocate(2*mbad));

  //int nbad; 
  intAr2 &lbad = work.lbad;

	int iinva = 0;
	int niter_incr = 0;
	// In lbad, store: [2*i + 0] = entity rank in cavity
	//                 [2*i + 1] = entity type: 0 = edg, 1 = fac, 2 = tet


	int nbpo0 = msh.nbpoi,
	    npoi0 = msh.npoin,
	    nedg0 = msh.nedge,
	    nfac0 = msh.nface,
	    nele0 = msh.nelem;

  double qmax;

  /*  --------------   Correctness checks -------------------- */
	if(cav.lcedg.get_n() > 0 && msh.isboundary_edges() || cav.lcfac.get_n() > 0 && msh.isboundary_faces() ){
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

    // Only if the point is pure surface need the normal be given
    // Otherwise we'll compute on the fly for each surface. 
    // This begs refactoring at some point 
    if(cav.lcfac.get_n() > 0 && msh.idim >= 3 && ityp == 2){
      METRIS_ASSERT(cav.nrmal != NULL);
    }
    // This can happen legitimately. (but the op is not valid)
		//if(ibpoi < 0) {
    //  METRIS_THROW_MSG(WArgExcept(),
		//	"When (re)inserting a boundary point, provide an ibpoi with (u,v)s\n ipins = "
    //  <<cav.ipins<<" ibpoi = "<<ibpoi);
    //}
	}


  // This only implements norempts for now. 
	ierro = check_cavity_topo(msh, cav, opts, ithread);
	if(ierro > 0) goto cleanup;


	// Goto constraint, decl needs to be above
	iinva = 0;
	niter_incr = 0;
	do{
   /*  -------------- Generate final cavity -------------------- 
   		------- For typent in (line|face|tetra) do 
       -------   Generate typent cavity + new typent-1 elements boundary
       -------   Reconnect bdry to ipins
   */

			ierro = reconnect_lincav<MFT, ideg>(msh, cav, opts, ithread);
			if(ierro > 0) goto cleanup;

      CPRINTF1("-- reconnect_lincav done nedg0 = %d nedge = %d npoi0 = %d npoin = %d\n",
               nedg0,msh.nedge,npoi0,msh.npoin);

      //printf("## DEBUG: costly call to check_topo after lincav\n");
      //check_topo(msh,nbpo0, npoi0, nedg0, nfac0, nele0);

			ierro = reconnect_faccav<MFT, ideg>(msh, cav, opts, work, nedg0, &qmax, ithread);

      {
        bool check_qua = (opts.qmax_nec > 0 && msh.get_tdim() == 2)
                      || (opts.qmax_suf > 0 && msh.get_tdim() == 2)
                      || (opts.qmax_iff > 0 && msh.get_tdim() == 2);
        if(check_qua) CPRINTF1(" - after reconnect_faccav qmax = %f \n",qmax);
      }
      if(msh.get_tdim() == 2) info.qmax_end = qmax;

      //printf("## DEBUG: costly call to check_topo after lincav\n");
      //check_topo(msh,nbpo0, npoi0, nedg0, nfac0, nele0);

			if(ierro > 0) goto cleanup;

      CPRINTF1("-- reconnect_faccav done nfac0 = %d nface = %d \n",nfac0,msh.nface);


			ierro = reconnect_tetcav<MFT, ideg>(msh, cav, opts, &qmax, ithread);
      if(ierro > 0) goto cleanup; 
	
	/*  -------------- Fast validity correction --------------------
	  If this doesn't pass, don't invest on expensive optimization 
		Instead, increase cavity and restart.
	*/


		ierro = correct_cavity_fast<MFT,ideg>(msh,cav,opts,npoi0,nedg0,
                                                      nfac0,nele0,lbad,work,ithread);
    if(ierro > 0) goto cleanup;

    int nbad = lbad.get_n();

		if(nbad > 0){
      if(DOPRINTS1()){
        CPRINTF1(" - debug nbad = %d \n",nbad);
        lbad.print(nbad);
        CPRINTF1("## Reject cavity\n");
        CPRINTF1(" - debug write mesh:\n");
        writeMesh("debug_failed_correctcav.meshb",msh, false, nedg0, nfac0, nele0);
      }
      ierro = CAV_ERR_INCORRECTIBLE;
      goto cleanup;
		}

	}while(iinva > 0 && niter_incr < opts.max_increase_cav_geo);

	//ierro = make_cavity_unit(msh,outcav,opts);
	//if(ierro > 0){
	//	ierro += 10*(ierr0);
	//	goto cleanup;
	//}
	//ierr0++;

  if(opts.dryrun){
    if((opts.qmax_suf < 0 || qmax > opts.qmax_suf) && 
       (opts.qmax_iff < 0 || qmax > opts.qmax_iff)){
      ierro = CAV_ERR_DRYFAIL1;
      goto cleanup;
    }
  }
  if(opts.qmax_nec > 0 && qmax < opts.qmax_nec){
    ierro = CAV_ERR_DRYFAIL2;
    goto cleanup;
  }


	ierro = update_cavity<MFT,ideg>(msh, cav, work, npoi0, nedg0, nfac0, nele0, ithread);
	if(ierro > 0) goto cleanup;

  info.done = true;


  CPRINTF1("-- Cavity successful exit\n");
	return 0;


  //-------- cleanup (error case)
	cleanup:
  CPRINTF1("-- Cavity error ierro = %d \n",ierro);
  if(DOPRINTS2()){
    // Write out the cavity. 
    writeMesh("cavenderr",msh,true,nedg0,nfac0,nele0);
  }
	//METRIS_THROW_MSG(TODOExcept(), 
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
  for(int ibpoi = nbpo0; ibpoi < msh.nbpoi; ibpoi++){
    int ip = msh.bpo2ibi(ibpoi,0);
    if( ip < 0 ) continue;
    if(msh.poi2tag(ithread,ip) >= msh.tag[ithread]) continue;
    msh.poi2tag(ithread,ip) = msh.tag[ithread];
    msh.rembpotag(ip,ithread);
  }
	msh.set_nbpoi(nbpo0);
	msh.set_npoin(npoi0);
	msh.set_nedge(nedg0);
	msh.set_nface(nfac0);
	msh.set_nelem(nele0);
	return ierro;
}
// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int cavity_operator<MetricFieldAnalytical, n >(Mesh<MetricFieldAnalytical>&,MshCavity&,CavOprOpt&,CavWrkArrs&,CavOprInfo&,int);\
template int cavity_operator<MetricFieldFE        , n >(Mesh<MetricFieldFE        >&,MshCavity&,CavOprOpt&,CavWrkArrs&,CavOprInfo&,int);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


}// End namespace
