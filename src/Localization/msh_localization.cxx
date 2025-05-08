//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "msh_localization.hxx"
#include "low_localization.hxx"

#include "../msh_lag2bez.hxx"
#include "../utils/aux_misc.hxx"
#include "../linalg/det.hxx"
#include "../low_topo.hxx"
#include "../low_normal.hxx"
#include "../low_geo.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/aux_timer.hxx"
#include "../utils/mprintf.hxx"
//#include "msh_metric.hxx"
#include "../io_libmeshb.hxx"
#include "../Boundary/low_projsurf.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

namespace Metris{




template <int gdim, int tdim, int ideg>
int locMesh(MeshBase &msh, int *ientt, 
           const double* coop, int pdim, const double *uvsrf,
           int iref, const double* algnd_, 
	         double* coopr, double* bary, 
	         double tol, int ithrd, bool iexpensive){
  GETVDEPTH(msh.param);

  METRIS_ASSERT(pdim > 0);
  METRIS_ASSERT(pdim == msh.get_tdim() || uvsrf != NULL);

  int ierro = 0;

  // Instead of barycentrics, consider scalar product of displacement 
  // with opposite edge or face and take minimum as the best nei crit. 
  static int nwarnprt = 0;
  const bool dir_nei_criterion = tdim != 3;
  if(!dir_nei_criterion && nwarnprt++ < 10 && msh.param->iverb != 0) 
    MPRINTF("## WARNING dir_nei_criterion disabled\n");


  //printf("Debug set iptr in locMesh = 4 inp guess = %d \n",*ientt);
  //iverb = 4;

  static_assert(gdim == 1 || gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim && tdim > 0);

  //if constexpr(tdim < gdim) METRIS_ASSERT(algnd_ != NULL);

  //constexpr int nnode = msh.nnode(tdim);
  //constexpr auto entnpps = ENTNPPS(tdim);
  //constexpr int nnode = entnpps[ideg];

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
    #ifndef NDEBUG
		MPRINTF("## locMeshVol inva ini guess %d, use 1\n",*ientt);
    //if(msh.param->dbgfull){
      printf("## WAIT HERE\n");
      wait();
      METRIS_THROW(GeomExcept());
    //}
    #endif
		*ientt = 0;
	}

  // Only use BB crit if interior expected. 
	if(tdim == gdim && locMeshQuick<gdim>(msh,coop)) return LOC_ERR_OUTBB;


  if(DOPRINTS1()){
    CPRINTF1("-- START locMesh gdim %d tdim %d ideg = %d guess %d "
                         "search coor0 = ",gdim,tdim,ideg,*ientt);
    dblAr1(gdim,coop).print();
  }


	if constexpr(ideg > 1){

    ierro = locMesh<gdim,tdim,1>(msh,ientt,coop,pdim,uvsrf,iref,algnd_,
                                 coopr,bary,tol,ithrd);

    //if(ierro != 0){
    //  if(iverb >= 3) printf("##Failed P1 localization. Attempting P%d\n",ideg);
    //  return LOC_ERR_FAILP1;
    //}

    if(DOPRINTS1()){
      double dist = geterrl2<gdim>(coopr,coop);
      CPRINTF1(" -> P1 loc done ientt = %d bary = ",*ientt);
      dblAr1(tdim+1,bary).print();
      CPRINTF1(" - dist = %15.7e\n",dist);
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
  intAr1 &lnext = msh.iwork; 
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
    int ntry = 0;
	  while(!ifnd){
      INCVDEPTH(msh.param);
      ntry++;
      if(ntry > msh.nentt(tdim)){
        MPRINTF("ntry = %d nentt %d \n",ntry, msh.nentt(tdim));
        METRIS_THROW_MSG(AlgoExcept(),"TRIED ALL ELEMENTS !");
      }
      while(lnext.n1_ > 0){
        *ientt = lnext.pop(); 
      METRIS_ASSERT(*ientt >= 0 && *ientt < nentt);
      METRIS_ASSERT(!isdeadent(*ientt,ent2poi));
      ent2tag(ithrd,*ientt) = msh.tag[ithrd];
      niter++;


        if constexpr(ideg > 1){
          for(int ii = 0; ii < gdim + 1 ; ii++) bary[ii] = 1.0 / (gdim + 1);
        }


        if(gdim == tdim){

          {
          INCVDEPTH(msh.param);
          ierro = inveval<gdim,ideg>(msh,*ientt,coop,coopr,bary,tolcur);
          }
          CPRINTF2(" - called inveval dim %d deg %d ientt %d tol %e got ierro %d bary = ",
                   gdim,ideg,*ientt,tolcur,ierro);
          if(DOPRINTS2()){
            dblAr1(gdim+1,bary).print();
            CPRINTF2(" got coopr ");
            dblAr1(gdim,coopr).print();
          }

        }else if(tdim == 2){

          // Unlike tdim 1, we're not doing (u,v)-space localization

          ierro = projptfac<ideg>(msh, coop, *ientt, bary, coopr);

          if(ierro == 0 && algnd_ != NULL){
            double dum[gdim], jmat[2][gdim];
            eval2<gdim,ideg>(msh.coord, ent2poi[*ientt],
                             msh.getBasis(), DifVar::Bary, DifVar::None,
                             bary, dum, jmat[0], NULL);
            double norfac[3];
            vecprod(jmat[0], jmat[1], norfac);
            if(normalize_vec<3>(norfac))METRIS_THROW_MSG(GeomExcept(),"Singular norfac");

            double dtprd = getprdl2<gdim>(norfac, algnd);
            double dev = 1 - abs(dtprd);

            double maxdev = msh.get_geodev(2);
            // If the mesh is a MeshBack, we have geodev for each edge. 
            if(msh.meshClass() == MeshClass::MeshBack){
              maxdev = ((MeshBack &) msh).fac2dev[*ientt];
            }

            CPRINTF1(" - bdry 2: dtprd %15.7e dev %15.7e <?= %15.7e algnd = %f %f %f" 
              " norfac = %f %f %f\n",dtprd,dev,maxdev,algnd[0],algnd[1],algnd[2],
               norfac[0],norfac[1],norfac[2]);
            if(dev > maxdev) ierro = 2;
          }
        }else{ // tdim == 1

          if(pdim == 1 && msh.CAD()){
            // If point is line then use t coordinate
            int ipoi1 = ent2poi(*ientt,0);
            int ipoi2 = ent2poi(*ientt,1);
            int ibpo1 = msh.poi2ebp(ipoi1,1,*ientt,-1);
            int ibpo2 = msh.poi2ebp(ipoi2,1,*ientt,-1);
            METRIS_ASSERT(ibpo1 >= 0 && ibpo2 >= 0);

            CPRINTF1(" ipoi1 = %d ibpo1 = %d ipoi2 = %d ibpo2 = %d\n",
                     ipoi1,ibpo1,ipoi2,ibpo2);

            if(ibpo1 < 0 || ibpo2 < 0){
              ibpo1 = msh.poi2bpo[ipoi1];
              MPRINTF(" - ipoi1 = %d dump all ibpo1 start at %d \n",ipoi1,ibpo1);
              for(ibpo1 =  msh.poi2bpo[ipoi1]; ibpo1 >= 0; ibpo1 =  msh.bpo2ibi(ibpo1,3)){
                MPRINTF(" %d :",ibpo1);
                intAr1(nibi,msh.bpo2ibi[ibpo1]).print();
              }
              ibpo2 = msh.poi2bpo[ipoi2];
              MPRINTF(" ipoi2 = %d dump all ibpo2 start at %d \n",ipoi2,ibpo2);
              for(ibpo2 =  msh.poi2bpo[ipoi2]; ibpo2 >= 0; ibpo2 =  msh.bpo2ibi(ibpo2,3)){
                MPRINTF(" %d :",ibpo2);
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

            CPRINTF1(" - t1 = %f t2 = %f tp = %f bary = %f %f \n",
                     t1,t2,tp, bary[0],bary[1]);

            for(int ii = 0; ii < 2; ii++){
              if( bary[ii] >   - Constants::baryTol 
              &&  bary[ii] < 1 + Constants::baryTol ) continue;
              ierro = 1;
            }

            // If the t coord is in the element, then do a proper projection 
            // to get the real bary. 

            // In some cases, there is no "real" bary on the edge. 
            if(ierro == 0){


              //ierro = projptedg<gdim,ideg>(msh, coop, *ientt, bary, coopr);
              //if constexpr (ideg > 1)
              //  METRIS_THROW_MSG(TODOExcept(), "Implement Pk projptedg");
              CPRINTF1(" - found t in %d w/ t bary = %15.7e %15.7e ierro = %d\n",
                       *ientt,bary[0],bary[1],ierro);

              // Compute coopr
              eval1<gdim,ideg>(msh.coord,ent2poi[*ientt],msh.getBasis(),DifVar::None,
                               DifVar::None,bary,coopr,NULL,NULL);


              //if(ierro == 0){
              //  bool okbar2 = true;
              //  for(int ii = 0; ii < 2; ii++){
              //    if( bary[ii] >   - Constants::baryTol 
              //    &&  bary[ii] < 1 + Constants::baryTol ) continue;
              //    okbar2 = false;
              //  }
              //  if(!okbar2){
              //    printf("## T FITS BUT BARY IS WRONG AFTER PROJ! \n");
              //    printf("bary = %15.7e %15.7e\n",bary[0],bary[1]);
              //    METRIS_THROW(GeomExcept());
              //  }
              //}else{
              //  printf("ierro = %d \n",ierro);
              //  printf("bary = ");
              //  dblAr1(gdim+1,bary).print();
              //  METRIS_THROW(GeomExcept());
              //}

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

              CPRINTF1(" - bdry 1: dtprd %15.7e dev %15.7e <?= %15.7e algnd = %f %f" 
                " tanedg = %f %f \n",dtprd,dev,maxdev,algnd[0],algnd[1],tanedg[0],tanedg[1]);
              if(dev > maxdev) ierro = 2;
            }

          }// if(pdim == 1) // else

        }// if (gdim == tdim) // else


        if(ierro == 0){
          if(DOPRINTS1()){
            CPRINTF1("  - END niter = %d ierro %d ientt %d tdim %d bary ",niter,ierro,*ientt,tdim);
            dblAr1(gdim+1,bary).print();
          } 
          ifnd = 1;
          break;
        }


        CPRINTF1(" - not in %d got bary = ",*ientt);
        if(DOPRINTS1()) dblAr1(tdim + 1,bary).print();

        // Initially, we were using minimum barycentric coordinate as the criterion.
        // This is ok for isotropic elements. But highly anisotropic means 
        // a closer to 0 bary can in fact 
        //double bmin = 1.0e30;
        //int    imin = -1;
        // imax is probably a better strategy...
        double bmax = -1.0e30;
        int    imax = -1;

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

            CPRINTF1(" - Test neighbour %d = %d \n",i,ienei);
            if(ienei < 0) continue;
            CPRINTF1(" - ienei = %d tettag = %d tag = %d \n",ienei,
                                           ent2tag(ithrd,ienei),msh.tag[ithrd]);
            if(iref >= 0 && ent2ref[ienei] != iref) continue;
            if(ent2tag(ithrd,ienei) >= msh.tag[ithrd] ) continue;
            lnext.stack(ienei);
            //if(bary[i] < bmin){
            //  bmin = bary[i];
            //  //imin = i;
            //}
          }
        }else{
          if constexpr(tdim < gdim){
            // For edges, just use the tangent 
            if constexpr(tdim == 1){
              for(int ii = 0; ii < gdim; ii++) 
                nrm1[ii] = msh.coord(ent2poi(*ientt,1),ii) 
                         - msh.coord(ent2poi(*ientt,0),ii);
              if(DOPRINTS2()){
                CPRINTF2(" - using nrm1 = ");
                dblAr1(gdim,nrm1).print();
              }
            }else{
              getnorfacP1(msh.fac2poi[*ientt], msh.coord, nrm1);
              METRIS_ENFORCE(!normalize_vec<3>(nrm1));
            }
          }
          for(int ii = 0 ; ii < tdim + 1 ; ii++){
            int ienei = ent2ent(*ientt,ii);
            if(DOPRINTS1() && ii > 0) MPRINTF("\n");
            CPRINTF1(" - check ienei %d bary %15.7e ",ienei,bary[ii]);
            if(ienei < 0) continue;
            CPRINTF1(" iref %d =? %d ",ent2ref[ienei], iref);
            if(iref >= 0 && ent2ref[ienei] != iref) continue;
            CPRINTF1(" nei tag? %d ",ent2tag(ithrd,ienei) >= msh.tag[ithrd]);
            if(ent2tag(ithrd,ienei) >= msh.tag[ithrd] ) continue;
            if(DOPRINTS1()) printf("\n");

            // if sg = 1, apply normal computation
            // if -1, opposite sign
            // sg = 0 if "more or less orthogonal", in which case
            // check the neighbour in doubt.
            int sg = 1;
            if(tdim < gdim){ //  && tdim < pdim
              // In this case, the ii-th neighbour is not guaranteed 
              // to be in the half-space bary[ii] < 0
              // We need to check the scalar product of normal with current
              if constexpr(tdim == 1){
                for(int ii = 0; ii < gdim; ii++) 
                  nrm2[ii] = msh.coord(ent2poi(ienei,0),ii) 
                           - msh.coord(ent2poi(ienei,1),ii);
                if(DOPRINTS2()){
                  CPRINTF2(" - using nrm2 = ");
                  dblAr1(gdim,nrm2).print();
                }
              }else{
                getnorfacP1(msh.fac2poi[*ientt], msh.coord, nrm2);
                METRIS_ENFORCE(!normalize_vec<3>(nrm2));
              }
              double dtprd = getprdl2<gdim>(nrm1,nrm2);
              // lines are not oriented
              if(tdim == 1){
                // - If this is first neighbour (ii == 0)
                //   - If ent2poi(ientt,1-ii) == ent2poi(ienei,0): +
                //   -                                         1 : -
                // - If second neighbour (ii == 1)
                //   - If ent2poi(ientt,1-ii) == ent2poi(ienei,0): -
                //   -                                         1 : +
                // Hence if ent2poi(ientt,1-ii) == ent2poi(ienei,ii)  : +
                //       if ent2poi(ientt,1-ii) == ent2poi(ienei,1-ii): -
                if(ent2poi(ienei,ii) == ent2poi(*ientt,1-ii)){
                  dtprd *= -1;
                }
              }
              if(abs(dtprd) < Constants::baryTol)  sg = 0;
              else if(dtprd < -Constants::baryTol) sg = -1;
              //else if(dtprd < -Constants::baryTol) sg = -1;
              // We can have a pertinent "fold back" (i.e. negative sg per the
              // previous law), but the barycentric is also negative. 
              // Since sg = 0 lead to never skip, simply put sg = 0 if not > 0
              CPRINTF1(" - dtprd = %f sg = %d ",dtprd,sg);
            }

            if(sg != 0 && sg*bary[ii] > -Constants::baryTol) continue;

            if(dir_nei_criterion && tdim > 1){


              // If using direction crit compute P1 centroid 
              double coom[gdim] = {}; // value-init to 0

              for(int ii = 0; ii < tdim + 1; ii++){
                for(int jj = 0; jj < gdim; jj++){
                  coom[jj] += msh.coord(ent2poi(*ientt,ii), jj) / (tdim + 1);
                }
              }

              for(int jj = 0; jj < gdim; jj++)
                edg2[jj] = coom[jj] - coopr[jj];


              if(DOPRINTS1()){
                CPRINTF1(" - dbg edg2 nrm %e = ",getnrml2<gdim>(edg2));
                dblAr1(gdim,edg2).print();
              }
              METRIS_ENFORCE(!normalize_vec<gdim>(edg2));


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
                CPRINTF1(" dev = %15.7e ",dev);

                if(dev >= bmax){
                  bmax = dev;
                  imax = ii;
                }


              }else if(tdim == 3){

                METRIS_THROW_MSG(TODOExcept(),"Implement getnorface and use that in locMesh");
              }

            }else{

              //if(sg == 1 && bary[ii] < bmin
              //|| bmin > 1.0e29){
              //  if(sg == 1) bmin = bary[ii];
              //  else        bmin = 1.0e30;
              //  //imin = ii;
              //}
              if(sg != 0 && sg*bary[ii] > bmax
              || bmax < -1.0e29){
                if(sg == 0) bmax = -1.0e30;
                else        bmax = sg*bary[ii];
                imax = ii;
              }

            }
          }
          CPRINTF1("\n");
        }

        //METRIS_ASSERT_MSG(imin != -1,"NO ELIGIBLE NEXT ELEMENT ideg = " << ideg)
        if(imax == -1){
        //if(imin == -1){
          CPRINTF1("-- END no candidates\n");
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
              CPRINTF1(" - restart loc dim %d from %d \n",tdim-1,ientf);
              // We could decrement then increment after but only in ideg = 1
              // More generally, we keep a max tag and will set in the end 
              int tag0 = msh.tag[ithrd];
              int ierr2 = locMesh<gdim, tdim-1, ideg>(msh, &ientf, coop, pdim, NULL,
                                   -1, NULL, coopf, barf, tol, ithrd, iexpensive);
              maxtag = MAX(maxtag, msh.tag[ithrd]);
              msh.tag[ithrd] = tag0;
              if(ierr2 > 0 && ierr2 != LOC_ERR_ALLPOS) ierro = LOC_ERR_PROJ;
              CPRINTF1(" # lower dim tried but failed ierro %d \n",ierr2);
            }else{
              CPRINTF1(" # lower dim not pertinent \n");
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

            CPRINTF1(" - ientt = %d -> %d \n",*ientt,ient2);
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
                METRIS_THROW_MSG(TODOExcept(), "Get bary from facet in case tdim = "<<tdim);
              }

              // Maybe an error, maybe not, but not a standard run. 
              // Caller needs to check the distance and make a decision. 
              ierro = LOC_WARN_PROJ;
              goto cleanup;

            }else{

              CPRINTF1(" - detach tdim = %d -> %d from ientt %d \n", tdim-1,tdim,ient2);

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
          CPRINTF1(" - imax = %d bmax = %f\n",imax,bmax);
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

