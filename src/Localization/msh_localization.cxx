//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "msh_localization.hxx"
#include "low_localization.hxx"

#include "../msh_lag2bez.hxx"
#include "../aux_utils.hxx"
#include "../low_topo.hxx"
#include "../low_geo.hxx"
#include "../aux_timer.hxx"
//#include "msh_metric.hxx"
#include "../io_libmeshb.hxx"
#include "../Boundary/low_projsurf.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

namespace Metris{


// Very rudimentary routine to be enhanced with a kd-tree in the future 
// Localize points, interpolate metric and store seed 
// Needs two slots in tag arrays 
// if ipoi0 > 0, it should be the number of P1 points !
template<class MetricFieldType, int bdeg>
void interpFrontBack(Mesh<MetricFieldType> &msh, MeshBack &bak, int ipoi0){
	if(bak.getBasis() == FEBasis::Lagrange && bak.curdeg > 1) 
			METRIS_THROW_MSG(WArgExcept(), "Back should be in Bézier format!");

  METRIS_ENFORCE_MSG(msh.idim == msh.get_tdim(), "Mesh is surface or line in plane.");

  METRIS_ENFORCE(METRIS_MAXTAGS >= 2);

  msh.setBasis(FEBasis::Lagrange);

  int ierro;

	bak.met.setSpace(MetSpace::Log);
	bak.setBasis(FEBasis::Bezier);

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



	msh.tag[0]++;

	int gdim  = msh.idim;

  // Will increase in size if needed 
  intAr1 lerro(100);
  lerro.set_n(0); 

  double t0 = get_wall_time();
  for(int ipoin = ipoi0; ipoin < msh.npoin; ipoin++){
    if(msh.poi2ent(ipoin,0) < 0) continue;

    int pdim = msh.getpoitdim(ipoin);
    // Corner not an element, start at 1 if lower. 
    pdim = MAX(pdim, 1);

    // Localize in elements of same dimension 
    const int tdim = pdim;

    if(pdim == 0) continue;


    int iref = -1;
    int iseed = -1;
    double algnd[3];

    // Corner is simply a copy from corresponding corned in backmesh
    if(pdim == 0){
      printf("## WATCH OUT THIS IS A HACK Assuming %d front = %d back\n",
        ipoin,ipoin);
      int nnmet = (msh.idim * (msh.idim + 1)) / 2;
      for(int ii = 0; ii < nnmet; ii++) msh.met(ipoin,ii) = bak.met(ipoin,ii);
      continue;
    }

    if(pdim < msh.idim){
      if(tdim == 2){
        METRIS_THROW_MSG(TODOExcept(), 
                                "Implement get face dir interpFrontBack in 3D")
      }
      int ibpoi = msh.poi2bpo[ipoin];
      METRIS_ASSERT(ibpoi >= 0);
      for(;ibpoi != -1; ibpoi = msh.bpo2ibi(ibpoi,3)){
        int itype = msh.bpo2ibi(ibpoi,1);
        if(itype == tdim){
          iseed = msh.bpo2ibi(ibpoi,2);
          iref  = msh.edg2ref[iseed];

          if(msh.CAD()){
            // Also compute the tangent at
            double result[18];
            ego obj = msh.CAD.cad2edg[iref];
            ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
            METRIS_ENFORCE(ierro == 0);
            for(int ii = 0; ii < msh.idim; ii++) algnd[ii] = result[3+ii];
          }else{
            int ipoi1 = msh.edg2poi(iseed,0);
            int ipoi2 = msh.edg2poi(iseed,1);
            for(int ii = 0; ii < msh.idim; ii++){
              algnd[ii] = msh.coord(ipoi1,ii) - msh.coord(ipoi2,ii);
            }
          }
        }
      }
      METRIS_ASSERT(iref != -1);
      METRIS_ASSERT(iseed != -1);
    }else{
      // Get seed, ref
      iseed = getpoifac(msh,ipoin);
      METRIS_ASSERT(iseed != -1);
      iref = msh.fac2ref[iseed];
      METRIS_ASSERT(iref != -1);
    }


    //if(msh.CAD()){
      ierro = msh.interpMetBack(ipoin, tdim, iseed, iref, algnd);
    //}else{
    //  ierro = msh.interpMetBack(ipoin, tdim, iseed, iref, NULL);
    //}
    if(ierro != 0){
      printf("## interpMetBack failed ierro = %d\n",ierro);
      printf("Failed to find ipoin %d iseed %d tdim %d iref %d pdim %d \n",
             ipoin,iseed,tdim,iref,pdim);
      printf("ibpoi %d \n",msh.poi2bpo[ipoin]);
      for(int ii = 0; ii < 3; ii++){
        printf("iseed vertex %d \n",msh.fac2poi(iseed,ii));
        printf("poi2bak = %d \n",msh.poi2bak(msh.fac2poi(iseed,ii),1));
      }
      writeMesh("interpMetDebug", msh);
      #ifndef NDEBUG
      msh.param->iverb = 20;
      printf("Try to localize coop = ");
      dblAr1(msh.idim,msh.coord[ipoin]).print();

      int ipdbg = msh.bak->newpoitopo(-1,-1);
      msh.bak->template newbpotopo<0>(ipdbg,ipdbg);
      for(int ii = 0; ii < msh.idim; ii++) 
        msh.bak->coord[ipdbg][ii] = msh.coord(ipoin,ii);
      writeMesh("interpMetDebug.back",*(msh.bak));
      ierro = msh.interpMetBack(ipoin, tdim, iseed, iref, algnd);
      #endif

      METRIS_THROW(TopoExcept());
    }


  } // for int ipoin
  double t1 = get_wall_time();
  printf("-- Interp Back -> Front time %f pt/s %d nerror %d \n",t1-t0,
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
        for(int ii = 0 ;ii < nnode; ii++){
          int ipoi2 = ent2poi(iebal,ii);
          // These are the points that have failed 
          if(msh.poi2tag(0,ipoi2) >= msh.tag[0]) continue;

          int ieleg = msh.poi2bak(ipoi2,tdim-1);
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
            msh.poi2bak(ipoin,tdim-1) = ieleg;
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




template <int gdim, int tdim, int ideg>
int locMesh(MeshBase &msh, int *ientt, 
           const double* coop, int pdim, const double *uvsrf,
           int iref, const double* algnd_, 
	         double* coopr, double* bary, 
	         double tol, int ithrd, bool iexpensive){

  METRIS_ASSERT(pdim > 0);
  METRIS_ASSERT(pdim == msh.get_tdim() || uvsrf != NULL);

  int ierro = 0;
  const int iverb = msh.param->iverb;

  // Instead of barycentrics, consider scalar product of displacement 
  // with opposite edge or face and take minimum as the best nei crit. 
  const bool dir_nei_criterion = true;

  //printf("Debug set iptr in locMesh = 4 inp guess = %d \n",*ientt);
  //iverb = 4;

  static_assert(gdim == 1 || gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim && tdim > 0);

  //if constexpr(tdim < gdim) METRIS_ASSERT(algnd_ != NULL);

  //constexpr int nnode = msh.nnode(tdim);
  constexpr auto entnpps = ENTNPPS(tdim);
  constexpr int nnode = entnpps[ideg];

  double algnd[gdim];
  if(algnd_ != NULL){
    double algnd_norm = getnrml2<gdim>(algnd_);
    if(algnd_norm < 1e-32) METRIS_THROW_MSG(GeomExcept(), "Singular algnd");
    for(int ii = 0; ii < gdim; ii++) algnd[ii] = algnd_[ii] / sqrt(algnd_norm);
  }

  int nentt = msh.nentt(tdim);
        intAr2 &ent2tag = msh.ent2tag(tdim);
  const intAr2 &ent2ent = msh.ent2ent(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const intAr1 &ent2ref = msh.ent2ref(tdim);
  if(iref >= 0) METRIS_ASSERT_MSG(ent2ref[*ientt] == iref,
       "Provided ref "<<iref<<" is not seed ref = "<<ent2ref[*ientt]);


  //double tolcur = MAX(1.0e-2,tol+1.0e-16);
  //if constexpr(ideg == 1) tolcur = tol;
  //// With fast Newton this does not matter, in fact probably faster
  const double tolcur = tol;

	if(*ientt < 0 || *ientt >= nentt){
		printf("## locMeshVol inva ini guess %d, use 1\n",*ientt);
    //if(msh.param->dbgfull){
      printf("## WAIT HERE\n");
      wait();
      METRIS_THROW(GeomExcept());
    //}
		*ientt = 0;
	}

  // Only use BB crit if interior expected. 
	if(tdim == gdim && locMeshQuick<gdim>(msh,coop)) return LOC_ERR_OUTBB;


  if(iverb >= 3) printf("-- START locMesh gdim %d tdim %d ideg = %d guess %d "
                       "search coor0 = ",gdim,tdim,ideg,*ientt);
  if(iverb >= 3) dblAr1(gdim,coop).print();


	if constexpr(ideg > 1){

    ierro = locMesh<gdim,tdim,1>(msh,ientt,coop,pdim,uvsrf,iref,algnd_,
                                 coopr,bary,tol,ithrd);

    //if(ierro != 0){
    //  if(iverb >= 3) printf("##Failed P1 localization. Attempting P%d\n",ideg);
    //  return LOC_ERR_FAILP1;
    //}

    if(iverb >= 3){
      double dist = geterrl2<gdim>(coopr,coop);
      printf(" -> P1 loc done ientt = %d bary = ",*ientt);
      dblAr1(tdim+1,bary).print();
      printf(" dist = %15.7e\n",dist);
    }

    //if(gdim == tdim){
    //  double tol1 = getepsent<gdim>(msh, tdim, *ientt);
    //  ierro = inveval<gdim,ideg>(msh,*ientt,coop,coopr,bary,tolcur*tol1);
    //}else if(tdim == 2){
    //  METRIS_THROW_MSG(TODOExcept(), "Implement projptfac in low_projsurf")
    //}else{
    //  ierro = projptedg<gdim,ideg>(msh, coop, *ientt, bary, coopr);
    //  ierro = 0;
    //  if(bary[0] < -Constants::baryTol || bary[0] > 1 + Constants::baryTol)
    //    ierro = 1;
    //  // redundant but who knows with floating point
    //  if(bary[1] < -Constants::baryTol || bary[1] > 1 + Constants::baryTol)
    //    ierro = 1;
    //  //METRIS_THROW_MSG(TODOExcept(), "Implement Pk projptedg");
    //  //METRIS_THROW_MSG(TODOExcept(), "Implement Pk algnd handling");
    //}
    //if(ierro == 0){
    //  return 0;
    //}else{
    //  printf(" - Failed Pk localization in P1 element got ierro = %d \n",ierro);
    //}
    
	}

  //intAr1 lnext(10);
  intAr1 &iwrkarr = gdim == 2 ? msh.fac2iwk : msh.tet2iwk;
  intAr1 lnext(iwrkarr.size(),&iwrkarr[0]);
  lnext.set_n(0);
  lnext.stack(*ientt);

  // For direction crit in tdim == gdim case (select best neighbour)
  double edg1[gdim], edg2[gdim];
  // For half-space of neighbour determination in tdim < gdim case 
  double nrm1[gdim], nrm2[gdim];

  int niter = 0;
  // This is to manage the recursive calls to lower dim and use 
  // dim-specific tags. All adjusted in the end. 
  int maxtag = msh.tag[ithrd];

  int ifnd;
  // Only for HO elements: localize with successively smaller tolerances 
  //do{

    ifnd = 0;
	  msh.tag[ithrd]++;
    maxtag = MAX(maxtag,msh.tag[ithrd]);
	  double tol1;
    int ntry = 0;
	  while(!ifnd){
      ntry++;
      if(ntry > msh.nentt(tdim)){
        printf("ntry = %d nentt %d \n",ntry, msh.nentt(tdim));
        METRIS_THROW_MSG(AlgoExcept(),"TRIED ALL ELEMENTS !");
      }
      while(lnext.n1_ > 0){
        *ientt = lnext.pop(); 
  	  	METRIS_ASSERT(*ientt >= 0 && *ientt < nentt);
        METRIS_ASSERT(!isdeadent(*ientt,ent2poi));
  	  	ent2tag(ithrd,*ientt) = msh.tag[ithrd];
  	  	niter++;

  	  	tol1 = getepsent<gdim>(msh, gdim, *ientt);

        if constexpr(ideg > 1){
          for(int ii = 0; ii < gdim + 1 ; ii++) bary[ii] = 1.0 / (gdim + 1);
        }

        #ifndef NDEBUG
        if(iverb >= 3){
          //printf("   - try in ientt = %d with tol %f \n",*ientt,tol1*tolcur);
          //printf("Vertices : ");
          //intAr1(nnode,ent2poi[*ientt]).print();
          //for(int ii = 0; ii < nnode; ii++){
          //  printf("%d : ",ent2poi(*ientt,ii));
          //  dblAr1(gdim,msh.coord[ent2poi(*ientt,ii)]).print();
          //}
        }
        #endif

        if(gdim == tdim){

          ierro = inveval<gdim,ideg>(msh,*ientt,coop,coopr,bary,tolcur*tol1);

        }else if(tdim == 2){

          METRIS_THROW_MSG(TODOExcept(), "Implement projptfac in low_projsurf")

        }else{

          if(pdim == 1){ 
            // If point is line then use t coordinate
            int ipoi1 = ent2poi(*ientt,0);
            int ipoi2 = ent2poi(*ientt,1);
            int ibpo1 = msh.poi2bpo[ipoi1];
            int ibpo2 = msh.poi2bpo[ipoi2];
            METRIS_ASSERT(ibpo1 >= 0 && ibpo2 >= 0);

            // Get the correct t for each point if corner
            if(msh.bpo2ibi(ibpo1,1) < 1) 
              ibpo1 = getent2bpo(msh, ibpo1, *ientt, tdim);
            if(msh.bpo2ibi(ibpo2,1) < 1) 
              ibpo2 = getent2bpo(msh, ibpo2, *ientt, tdim);

            if(iverb >= 3) printf(" ipoi1 = %d ibpo1 = %d ipoi2 = %d ibpo2 = %d\n",
              ipoi1,ibpo1,ipoi2,ibpo2);

            if(ibpo1 < 0 || ibpo2 < 0){
              ibpo1 = msh.poi2bpo[ipoi1];
              printf(" ipoi1 = %d dump all ibpo1 start at %d \n",ipoi1,ibpo1);
              for(ibpo1 =  msh.poi2bpo[ipoi1]; ibpo1 >= 0; ibpo1 =  msh.bpo2ibi(ibpo1,3)){
                printf(" %d :",ibpo1);
                intAr1(nibi,msh.bpo2ibi[ibpo1]).print();
              }
              ibpo2 = msh.poi2bpo[ipoi2];
              printf(" ipoi2 = %d dump all ibpo2 start at %d \n",ipoi2,ibpo2);
              for(ibpo2 =  msh.poi2bpo[ipoi2]; ibpo2 >= 0; ibpo2 =  msh.bpo2ibi(ibpo2,3)){
                printf(" %d :",ibpo2);
                intAr1(nibi,msh.bpo2ibi[ibpo2]).print();
              }
              METRIS_THROW(TopoExcept());
            }

            double t1 = msh.bpo2rbi(ibpo1,0);
            double t2 = msh.bpo2rbi(ibpo2,0);
            double tp = *uvsrf;

            if(abs(t2-t1) < 1.0e-16) METRIS_THROW_MSG(GeomExcept(),
                                      "t coordinates too close");

            ierro = 0;
            bary[0] = (t2 - tp) / (t2 - t1);
            bary[1] = (tp - t1) / (t2 - t1);

            if(iverb >= 3) printf("tdim 1/pdim 1: bary = %15.7e %15.7e\n",
                                  bary[0],bary[1]);

            for(int ii = 0; ii < 2; ii++){
              if( bary[ii] >   - Constants::baryTol 
              &&  bary[ii] < 1 + Constants::baryTol ) continue;
              ierro = 1;
            }

            // If the t coord is in the element, then do a proper projection 
            // to get the real bary. 
            if(ierro == 0){

              ierro = projptedg<gdim,ideg>(msh, coop, *ientt, bary, coopr);
              //if constexpr (ideg > 1)
              //  METRIS_THROW_MSG(TODOExcept(), "Implement Pk projptedg");
              if(iverb >= 3) printf(" - found t in %d new bary = " 
                "%15.7e %15.7e ierro = %d\n",*ientt,bary[0],bary[1],ierro);

              if(ierro == 0){
                bool okbar2 = true;
                for(int ii = 0; ii < 2; ii++){
                  if( bary[ii] >   - Constants::baryTol 
                  &&  bary[ii] < 1 + Constants::baryTol ) continue;
                  okbar2 = false;
                }
                if(!okbar2){
                  printf("## T FITS BUT BARY IS WRONG AFTER PROJ! \n");
                  printf("bary = %15.7e %15.7e\n",bary[0],bary[1]);
                  METRIS_THROW(GeomExcept());
                }
              }else{
                METRIS_THROW(GeomExcept());
              }

            }

          }else{

            ierro = projptedg<gdim,ideg>(msh, coop, *ientt, bary, coopr);
            //if constexpr (ideg > 1)
            //  METRIS_THROW_MSG(TODOExcept(), "Implement Pk projptedg");

            if(ierro == 0 && algnd_ != NULL){
              double dum[gdim], tanedg[gdim];
              eval1<gdim,ideg>(msh.coord, ent2poi[*ientt],
                               msh.getBasis(), DifVar::Bary, DifVar::None,
                               bary, dum, tanedg, NULL);
              double tanedg_norm = getnrml2<gdim>(tanedg);
              if(tanedg_norm < 1e-32) METRIS_THROW_MSG(GeomExcept(), 
                                                       "Singular algnd");
              for(int ii = 0; ii < gdim; ii++) 
                tanedg[ii] = tanedg[ii] / sqrt(tanedg_norm);

              double dtprd = getprdl2<gdim>(tanedg, algnd);
              double dev = 1 - abs(dtprd);

              double maxdev = msh.get_geodev(1);
              // If the mesh is a MeshBack, we have geodev for each edge. 
              if(msh.meshClass() == MeshClass::MeshBack){
                maxdev = ((MeshBack &) msh).edg2dev[*ientt];
              }

              if(iverb >= 3) printf("bdry 1: dtprd %15.7e dev %15.7e <?= %15.7e algnd = %f %f" 
                " tanedg = %f %f \n",dtprd,dev,maxdev,algnd[0],algnd[1],tanedg[0],tanedg[1]);
              if(dev > maxdev) ierro = 2;
            }

          }// if(pdim == 1) // else

        }// if (gdim == tdim) // else


        if(ierro == 0){
          if(iverb >= 3){
            printf("  - END niter = %d ierro %d ientt %d tdim %d bary ",niter,ierro,*ientt,tdim);
            dblAr1(gdim+1,bary).print();
          } 
          ifnd = 1;
          break;
        }


        if(iverb >= 3) printf("   - not in %d got bary = ",*ientt);
        if(iverb >= 3) dblAr1(tdim + 1,bary).print();

        // Initially, we were using minimum barycentric coordinate as the criterion.
        // This is ok for isotropic elements. But highly anisotropic means 
        // a closer to 0 bary can in fact 
        double bmin = 1.0e30;
  	  	int    imin = -1;
        // imax is probably a better strategy...
        double bmax = -1.0e30;
        int    imax = -1;

        // If using direction crit compute P1 centroid 
        if(dir_nei_criterion && tdim > 1){
          double coom[gdim];
          for(int jj = 0; jj < gdim; jj++){
            coom[jj] = 0.0;
            for(int ii = 0; ii < tdim + 1; ii++){
              coom[jj] += msh.coord(ent2poi(*ientt,ii), jj) / (tdim + 1);
            }
            edg2[jj] = coom[jj] - coopr[jj];
          }
          double nrm = getnrml2<gdim>(edg2);
          nrm = sqrt(nrm);
          if(nrm < 1.0e-16) METRIS_THROW_MSG(GeomExcept(),"Zero length displ "<< nrm);
          for(int jj = 0; jj < gdim; jj++) edg2[jj] /= nrm;
        }

        if(iexpensive){
          int idx[4] = {0, 1, 2, 3};
          //sortupto8_dec(bary,idx,tdim+1);
          sortupto8_dec<double,tdim+1>(bary,idx);
          // From worst to best 
          for(int ii = 0 ; ii < tdim + 1 ; ii++){
            int i = idx[ii]; 
            //if(bary[i] > tolcur*tol1) continue;
            //if(bary[i] > 1) continue;
            int ienei = ent2ent[*ientt][i];

            if(iverb >= 3) printf("Test neighbour %d = %d \n",i,ienei);
            if(ienei < 0) continue;
            if(iverb >= 3) printf("ienei = %d tettag = %d tag = %d \n",ienei,
                                           ent2tag(ithrd,ienei),msh.tag[ithrd]);
            if(iref >= 0 && ent2ref[ienei] != iref) continue;
            if(ent2tag(ithrd,ienei) >= msh.tag[ithrd] ) continue;
            lnext.stack(ienei);
            if(bary[i] < bmin){
              bmin = bary[i];
              imin = i;
            }
          }
        }else{
          if constexpr(tdim < gdim){
            // For edges, just use the tangent 
            if constexpr(tdim == 1){
              for(int ii = 0; ii < gdim; ii++) 
                nrm1[ii] = msh.coord(ent2poi(*ientt,1),ii) 
                         - msh.coord(ent2poi(*ientt,0),ii);
            }else{
              METRIS_THROW_MSG(TODOExcept(), "tdim < gdim case tdim != 1");
            }
          }
          for(int ii = 0 ; ii < tdim + 1 ; ii++){
            int ienei = ent2ent(*ientt,ii);
            if(iverb >= 3) printf("\n    - check ienei %d bary %15.7e ",ienei,bary[ii]);
            if(ienei < 0) continue;
            if(iverb >= 3) printf(" iref %d =? %d ",ent2ref[ienei], iref);
            if(iref >= 0 && ent2ref[ienei] != iref) continue;
            if(iverb >= 3) printf(" nei tag? %d ",ent2tag(ithrd,ienei) >= msh.tag[ithrd]);
            if(ent2tag(ithrd,ienei) >= msh.tag[ithrd] ) continue;

            // if sg = 1, apply normal computation
            // if -1, opposite sign
            // sg = 0 if "more or less orthogonal", in which case
            // check the neighbour in doubt.
            int sg = 1;
            if(tdim < gdim && tdim < pdim){
              // In this case, the ii-th neighbour is not guaranteed 
              // to be in the half-space bary[ii] < 0
              // We need to check the scalar product of normal with current
              if constexpr(tdim == 1){
                for(int ii = 0; ii < gdim; ii++) 
                  nrm2[ii] = msh.coord(ent2poi(ienei,1),ii) 
                           - msh.coord(ent2poi(ienei,0),ii);
              }else{
                METRIS_THROW_MSG(TODOExcept(),
                               "Surface half-space determination neighbour loc")
              }
              double dtprd = getprdl2<gdim>(nrm1,nrm2);
              if(abs(dtprd) < Constants::baryTol)  sg = 0;
              else if(dtprd < -Constants::baryTol) sg = 0;
              //else if(dtprd < -Constants::baryTol) sg = -1;
              // We can have a pertinent "fold back" (i.e. negative sg per the
              // previous law), but the barycentric is also negative. 
              // Since sg = 0 lead to never skip, simply put sg = 0 if not > 0
              if(iverb >= 3) printf(" dtprd = %f sg = %d ",dtprd,sg);
            }

            if(sg != 0 && sg*bary[ii] > -Constants::baryTol)continue;

            if(dir_nei_criterion && tdim > 1){

              if(tdim == 2){
                int ipoi1 = msh.fac2poi(*ientt, lnoed2[ii][0]);
                int ipoi2 = msh.fac2poi(*ientt, lnoed2[ii][1]);
                for(int jj = 0; jj < gdim; jj++) 
                  edg1[jj] = msh.coord(ipoi1,jj) - msh.coord(ipoi2,jj);
                double nrm = getnrml2<gdim>(edg1);
                nrm = sqrt(nrm);
                if(nrm < 1.0e-16) METRIS_THROW_MSG(GeomExcept(),"Zero length edge "<< nrm);
                double dtprd = getprdl2<gdim>(edg1,edg2);

                // We want to keep the one that has least scalar product 
                // i.e. edge is least aligned with displacement. 
                // Call 1 - that deviation and put that in bmax. 

                double dev = 1 - abs(dtprd) / nrm;
                if(iverb >= 3) printf(" dev = %15.7e ",dev);

                if(dev >= bmax){
                  bmax = dev;
                  imax = ii;
                }


              }else if(tdim == 3){
                METRIS_THROW_MSG(TODOExcept(),"Implement getnorface and use that in locMesh");
              }

            }else{

              if(sg == 1 && bary[ii] < bmin
              || bmin > 1.0e29){
                if(sg == 1) bmin = bary[ii];
                else        bmin = 1.0e30;
                imin = ii;
              }
              if(sg == 1 && bary[ii] > bmax
              || bmax < -1.0e29){
                if(sg == 1) bmax = bary[ii];
                else        bmax = -1.0e30;
                imax = ii;
              }

            }
          }
          if(iverb >= 3) printf("\n");
        }

        //METRIS_ASSERT_MSG(imin != -1,"NO ELIGIBLE NEXT ELEMENT ideg = " << ideg)
        if(imax == -1){
        //if(imin == -1){
          if(iverb >= 3) printf("END no candidates\n");
          ierro = LOC_ERR_ALLPOS;

          // If this dimension is higher than 1, try lower dimension
          // This is akin to a projection but better 
          // As there are cases where we legitimately haven't found 
          // the element the projected falls on. 
          
          if constexpr(tdim > 1){
            int ientf = -1;
            for(int ii = 0; ii < tdim + 1; ii++){
              if(bary[ii] < Constants::baryTol){
                // Look at case with 2 negative
                METRIS_ENFORCE(ientf == -1);
                if constexpr (tdim == 2){
                  ientf = msh.fac2edg(*ientt, ii);
                }else{
                  ientf = msh.tet2fac(*ientt, ii);
                }
              }
            }


            double coopf[gdim], barf[tdim];
            ierro = 0;
            if(ientf >= 0){
              if(iverb >= 3) printf(" - restart loc dim %d from %d \n",tdim-1,ientf);
              // We could decrement then increment after but only in ideg = 1
              // More generally, we keep a max tag and will set in the end 
              int tag0 = msh.tag[ithrd];
              int ierr2 = locMesh<gdim, tdim-1, ideg>(msh, &ientf, coop, pdim, NULL,
                                   -1, NULL, coopf, barf, tol, ithrd, iexpensive);
              maxtag = MAX(maxtag, msh.tag[ithrd]);
              msh.tag[ithrd] = tag0;
              if(ierr2 > 0 && ierr2 != LOC_ERR_ALLPOS) ierro = LOC_ERR_PROJ;
              if(iverb >= 3) printf("  - lower dim tried but failed ierro %d \n",ierr2);
            }else{
              if(iverb >= 3) printf("  - lower dim not pertinent \n");
              ierro = LOC_ERR_PROJ;
            }

            if(ierro == LOC_ERR_PROJ){
              goto cleanup;
            }

            for(int ii = 0; ii < gdim; ii++) coopr[ii] = coopf[ii];

            int ient2;
            if constexpr(tdim == 2){
              ient2 = msh.edg2fac[ientf];
            }else{
              ient2 = msh.fac2tet[ientf][0];
            }

            if(iverb >= 3) printf("  - ientt = %d -> %d \n",*ientt,ient2);
            *ientt = ient2;

            // If this entity has been seen before, call it quits. Get bary:
            if(ent2tag(ithrd,ient2) >= msh.tag[ithrd] || ierro == 0){

              if constexpr(tdim == 2){
                int ipoi1 = msh.edg2poi(ientf, 0);
                int ipoi2 = msh.edg2poi(ientf, 1);
                int iedl  = getedgfac(msh, *ientt, ipoi1, ipoi2);
                METRIS_ASSERT(iedl >= 0);
                bary[iedl] = 0;
                if(ipoi1 == msh.fac2poi(*ientt,lnoed2[iedl][0])){
                  bary[lnoed2[iedl][0]] = barf[0];
                  bary[lnoed2[iedl][1]] = barf[1];
                }else{
                  bary[lnoed2[iedl][0]] = barf[1];
                  bary[lnoed2[iedl][1]] = barf[0];
                }
              }else{
                METRIS_THROW_MSG(TODOExcept(), "Get bary from facet");
              }

              // Maybe an error, maybe not, but not a standard run. 
              // Caller needs to check the distance and make a decision. 
              ierro = LOC_WARN_PROJ;
              goto cleanup;

            }else{

              if(iverb >= 3) printf("  - detach tdim = %d -> %d from ientt %d \n",
                                    tdim-1,tdim,ient2);

              lnext.stack(*ientt);
              continue;

            }

          } 

          //// Project 
          //double sum = 0;
          //for(int ii = 0; ii < tdim + 1; ii++){
          //  if(bary[ii] < Constants::baryTol) bary[ii] = 0;
          //  if(bary[ii] > 1 - Constants::baryTol) bary[ii] = 1;
          //  sum += bary[ii];
          //}
          //for(int ii = 0; ii < tdim + 1; ii++) bary[ii] /= sum;
          //evalf(msh.coord, ent2poi[*ientt],  
          //      msh.getBasis(), DifVar::None, DifVar::None, 
          //      bary, coopr, NULL, NULL);

          //if(iverb >= 3){
          //  double dist = sqrt(geterrl2<gdim>(coopr,coor0));
          //  printf("  - project dist = %15.7e \n",dist);
          //}
          ierro = LOC_ERR_ALLPOS;
          goto cleanup;
        } // end if imax == -1

        // In other case, already stacked 
        if(!iexpensive){
          if(iverb >= 3) printf("   - imax = %d bmax = %f\n",imax,bmax);
          //lnext.stack(ent2ent[*ientt][imin]);
          lnext.stack(ent2ent[*ientt][imax]);
        }
  	  	//*ientt = ent2ent[*ientt][imin];
        
        //for(int ii = 0; ii < gdim + 1 ; ii++) bary[ii] = 1.0 / (gdim + 1);
      }
	  }//end while(!ifnd)
    //if constexpr(ideg == 1) break;
    //tolcur /= 10.0;
    //if(iverb >= 3) printf("  - restart loop with tol = %f > %f\n",tolcur,tol);
  //}while(tolcur >= tol); 

  cleanup:
  msh.tag[ithrd] = maxtag;
  return ierro;
}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int locMesh<1, 1, n >(MeshBase &msh, int *ientt, const double* coop, \
                               int pdim, const double* uvsrf,\
                               int iref, const double* algnd,\
                               double* coopr, double* bary, double tol,\
                               int ithrd, bool iexpensive);\
template int locMesh<2, 1, n >(MeshBase &msh, int *ientt, const double* coop, \
                               int pdim, const double* uvsrf,\
                               int iref, const double* algnd,\
                               double* coopr, double* bary, double tol,\
                               int ithrd, bool iexpensive);\
template int locMesh<2, 2, n >(MeshBase &msh, int *ientt, const double* coop, \
                               int pdim, const double* uvsrf,\
                               int iref, const double* algnd,\
                               double* coopr, double* bary, double tol,\
                               int ithrd, bool iexpensive);\
template int locMesh<3, 1, n >(MeshBase &msh, int *ientt, const double* coop, \
                               int pdim, const double* uvsrf,\
                               int iref, const double* algnd,\
	                            double* coopr, double* bary, double tol,\
                              int ithrd, bool iexpensive);\
template int locMesh<3, 2, n >(MeshBase &msh, int *ientt, const double* coop, \
                               int pdim, const double* uvsrf,\
                               int iref, const double* algnd,\
                               double* coopr, double* bary, double tol,\
                               int ithrd, bool iexpensive);\
template int locMesh<3, 3, n >(MeshBase &msh, int *ientt, const double* coop, \
                               int pdim, const double* uvsrf,\
                               int iref, const double* algnd,\
                               double* coopr, double* bary, double tol,\
                               int ithrd, bool iexpensive);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


//#include <src/msh_localization.ixx>


} // End namespace

