//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "DIRECT.hxx"
#include "../Mesh/MeshBase.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"


#include <sstream>
#include <fstream>

namespace Metris{


/*
This one optimizes within n tetrahedra at once, weighing their 
respective levels and costfun evaluations together.
Thought for an edge shell but can deal with any set of tetras.

The particular elements are not given to this algo, instead their 
number is given. This means that nblob reference tetrahedra are 
split which corresponds to as many sets of barycentric coordinates
for a list of tetras known by the caller. Only indexes are stored in
leval(2,:).

args.ent2pol stores: 1-4 vertices, 5 level, 6 the global host tetra

-- Can probably get rid of leval and store nfac0 = nfacl at 
beginning: since iface is replaced when split by center
triangle, the new ones are all at nfac0 ++   -> DONE

-- Can save space, getting rid of peval:
   - have user store values directly in args.fuelt(ieval + i)
   - have user compute triangle barycenter themselves 
      (same computations in the end except for optimum storing
           -> must recompute then)

To use:
Initialize iflag = 0
START:
  call directtri_drivertype()
  if(iflag < 0) GOTO END
  do i = 1,neval
     peval(4) = costfun(peval(1),peval(2),peval(3))
  enddo
GOTO START
END: optimum at barmin(1:2) with value fmin
*/
//-- Convex hull information
// - args.lhull,args.rhull: args.lhull points to triangle, args.rhull stores min f value
// - args.fuelt: evaluated function at face barycenter
// - peval: 1-3 coordinates of bary, 4 value
// - leval[ii]; which element this relates to
void DIBLOB(DIBLOB_args &args,
            int idim, int nblob, 
            intAr2 &leval, dblAr2 &peval, dblAr1 &feval,
            int *ifmin, double *barmin ,double *fmin){

  if(idim == 3){
    METRIS_THROW_MSG(TODOExcept(), "Improve DIBLOB dim 3 splits and run test_DIRECT");
  }

  METRIS_ASSERT(idim == 2 || idim == 3);

  GETVDEPTH(args.param);

  METRIS_ASSERT(args.miter > 0);


  if(args.iflag == 0){
    CPRINTF1("-- START DIBLOB\n");
    args.ent2pol.allocate(10*nblob,idim+1+2); // over-allocate a little
    args.ent2pol.set_n(nblob);
    args.fminhist.allocate(10*nblob);
    args.fminhist.set_n(nblob);
    args.fminhist.fill(0);

    args.coorl.allocate(40*nblob,idim); // over-allocate a little
    args.coorl.set_n(4*nblob);

    args.iflag = 1;
    args.niter = 0;

    for(int ii = 0; ii < nblob; ii++){
      int npoi0 = 4*ii;

      args.ent2pol(ii, 0) = npoi0;
      args.ent2pol(ii, 1) = npoi0+1;
      args.ent2pol(ii, 2) = npoi0+2;

      if(idim == 2){
        args.coorl(npoi0 + 0, 0) = 1.0;
        args.coorl(npoi0 + 0, 1) = 0.0;

        args.coorl(npoi0 + 1, 0) = 0.0;
        args.coorl(npoi0 + 1, 1) = 1.0;

        args.coorl(npoi0 + 2, 0) = 0.0;
        args.coorl(npoi0 + 2, 1) = 0.0;
      }else{
        args.coorl(npoi0 + 0, 0) = 0.0;
        args.coorl(npoi0 + 0, 1) = 0.0;
        args.coorl(npoi0 + 0, 2) = 0.0;

        args.coorl(npoi0 + 1, 0) = 1.0;
        args.coorl(npoi0 + 1, 1) = 0.0;
        args.coorl(npoi0 + 1, 2) = 0.0;

        args.coorl(npoi0 + 2, 0) = 0.0;
        args.coorl(npoi0 + 2, 1) = 1.0;
        args.coorl(npoi0 + 2, 2) = 0.0;

        args.coorl(npoi0 + 3, 0) = 0.0;
        args.coorl(npoi0 + 3, 1) = 0.0;
        args.coorl(npoi0 + 3, 2) = 1.0;

        args.ent2pol(ii, 3) = npoi0+3;
      }

      args.ent2pol(ii ,idim + 1) = 0;
      args.ent2pol(ii ,idim + 2) = ii;
    }// for ii

    // First round, evaluate on all the input elements
    leval.allocate(nblob,2);
    leval.set_n(nblob);
    peval.allocate(nblob,idim+1);
    peval.set_n(nblob);
    feval.allocate(nblob);
    feval.set_n(nblob);
    for(int ii = 0; ii < nblob; ii++){
      // -- This correspondence only holds at the first iter
      leval(ii,0) = ii;
      leval(ii,1) = ii;
      double sum = 0;
      for(int jj = 0; jj < idim; jj++){
        peval(ii,jj) = args.coorl(args.ent2pol(ii,0),jj)
                     + args.coorl(args.ent2pol(ii,1),jj)
                     + args.coorl(args.ent2pol(ii,2),jj);
        if(idim == 3) peval(ii,jj) += args.coorl(args.ent2pol(ii,3),jj);
        peval(ii,jj) /= (idim + 1);
        sum += peval(ii,jj);
      }
      peval(ii,idim) = 1 - sum; // for convenience
    }

    args.fmin_pre = 1.0e30;
    *fmin = 1.0e30;
    *ifmin= -1;

    args.lhull.allocate(args.miter);
    args.rhull.allocate(args.miter);

    return;
  }// if iflag == 0

  (args.niter)++;
  CPRINTF1("-- ENTER DIBLOB niter = {}/{} \n",args.niter,args.miter);
  if(args.niter >= args.miter){
    CPRINTF1("# DIBLOB exceeded number of iterations\n");
    args.iflag = -3;
    return;
  }

  int neval = leval.get_n();
  // This just increases the size, we already have elements allocated and set
  // up to old number of elements.
  args.fuelt.allocate(args.ent2pol.get_n());
  args.fuelt.set_n(args.ent2pol.get_n());
  // Update the fuelt values with the new evaluations. 
  int ielmin = -1;
  double fvmin = 1.0e30;
  for(int ii = 0; ii < neval; ii++){
    INCVDEPTH(args.param);
    int ietes = leval(ii,0);
    CPRINTF1(" - eval {} elt (loc) {} value {}\n",ii,ietes,feval[ii]);
    args.fuelt[ietes] = feval[ii];
    if(feval[ii] <= fvmin) fvmin = feval[ii];
    if(feval[ii] <= *fmin){
      //CPRINTF1(" -> fmin update\n");
      ielmin = ietes;
      args.fmin_pre = *fmin;
      *fmin  = feval[ii];
      *ifmin = leval(ii,1);
      CPRINTF1(" -> update fmin = {} ifmin = {} from eval {}\n",*fmin, *ifmin, ii);

      double sum = (barmin[0] = peval(ii,0));
      sum += (barmin[1] = peval(ii,1));
      if(idim == 3) sum += (barmin[2] = peval(ii,2));
      barmin[idim] = 1 - sum; // convenience
    }
  }// for ii < neval

  if(ielmin >= 0){
    args.fminhist[ielmin]++;
    if(args.fminhist[ielmin] >= args.nloc_switch && args.nloc_switch > 0 && !args.iloc_mode){
      CPRINTF1(" - iter {} element {} and parents have provided fmin {} iters in a row -> switch to local mode\n", 
        args.niter, ielmin, args.fminhist[ielmin]);
      args.iloc_mode = true;
    }
  }else if(args.iloc_mode){
    // In local mode, we expect a fmin update every iteration.
    // In the future, we can simply handle this error by reverting to non-local
    //printf("## Debug min eval {} fmin {} \n",fvmin, *fmin);
    //METRIS_ASSERT(args.iloc_mode == false);
    CPRINTF1(" - iter {} fmin not updated in local mode, switch back to global\n",args.niter);
    args.iloc_mode = false;
  }

  CPRINTF1(" - DIBLOB iter {} decrease {} vals {} {} \n", 
    args.niter, *fmin/args.fscale, *fmin, args.fscale);

  if( abs(*fmin) <= args.ftol*abs(args.fscale) ){
    CPRINTF1("-- END DIBLOB: ftol reached\n");
    args.iflag = -2;
    return;
  }

  if( (args.fmin_pre - *fmin) <= args.dftol*args.fscale ){
    CPRINTF1("-- END DIBLOB: dftol reached, decrease = {} < {}\n",
             (args.fmin_pre - *fmin) / args.fscale, args.dftol);
    args.iflag = -4;
    return;
  }
  

  leval.set_n(0);
  peval.set_n(0);
  feval.set_n(0);

  // In local mode, we only need to split the fmin providing element. 
  if(args.iloc_mode){
    METRIS_ASSERT(ielmin >= 0);
    int ieglo = args.ent2pol(ielmin,idim+2);
    int ilev  = args.ent2pol(ielmin,idim+1);
    int nele0 = args.ent2pol.get_n();
    if(idim == 2){
      aux_DIRECT_splittri(args,ielmin,ieglo,ilev);
    }else{
      aux_DIRECT_splittet(args,ielmin,ieglo,ilev);
    }

    aux_DIRECT_newevals(args, idim, nele0, ilev, ieglo, leval, peval, feval);

    CPRINTF1(" - ilev {} ielmin {} split\n",ilev, ielmin);
    return;
  }



  args.lhull.set_n(args.niter);
  args.rhull.set_n(args.niter);
  for(int ii = 0; ii < args.niter; ii++){
    args.lhull[ii] = -1;
    args.rhull[ii] = 1.0e30;
  }

  int nelel = args.ent2pol.get_n();
  int minlv = args.niter-1;
  int maxlv = -1;// Should be args.niter - 1 unless we change the splitting logic
  int nhull = 0;
  // Initialize the hull to the minimum evaluation at each level. 
  for(int ielem = 0; ielem < nelel; ielem++){
    if(ielem != ielmin) args.fminhist[ielem] = 0;

    INCVDEPTH(args.param);
    int ilev = args.ent2pol(ielem,idim+1);
    METRIS_ASSERT(ilev >= 0);
    if(args.fuelt[ielem] <= args.rhull[ilev]){
      if(args.lhull[ilev] < 0) minlv = MIN(minlv, ilev);
      if(args.lhull[ilev] < 0) maxlv = MAX(maxlv, ilev);
      if(args.lhull[ilev] < 0) nhull++;
      args.rhull[ilev] = args.fuelt[ielem];
      args.lhull[ilev] = ielem;
    }
    CPRINTF1(" - check element {} ilev = {} f = {} flev {}\n",ielem,ilev,
             args.fuelt[ielem],args.rhull[ilev]);
  }

  // No longer true with switch to local mode and back
  //METRIS_ASSERT_MSG(maxlv == args.niter - 1,"maxlv = "<<maxlv<<" args.niter = "<<args.niter);

  //// cull minlv
  //minlv = MAX(minlv, maxlv-10);


  if(DOPRINTS2()){
    CPRINTF2(" - {} potential hull points\n",args.niter);
    for(int ilev = minlv; ilev <= maxlv; ilev++){
      INCVDEPTH(args.param)
      if(args.lhull[ilev] < 0) continue;
      CPRINTF2("   - ilev {} ielem {} value {}\n",ilev,args.lhull[ilev],args.rhull[ilev]);
    }
  }

  if(DOPRINTS3()){
    // Dump to python file to plot
    std::vector<std::ostringstream> strlev(maxlv-minlv+1);
    intAr1 istrlev(maxlv-minlv+1);
    istrlev.fill(-1);

    for(int ilev = minlv; ilev <= maxlv; ilev++){
      strlev[ilev-minlv] << "vals" << ilev << "=[";
    }

    for(int ielem = 0; ielem < nelel; ielem++){
      int ilev = args.ent2pol(ielem,idim+1);
      if(ilev < minlv) continue;
      if(istrlev[ilev-minlv] == -1){
        istrlev[ilev-minlv] = 0;
      }else{
        strlev[ilev-minlv] << ",";
      }
      strlev[ilev-minlv] << args.fuelt[ielem];
    }

    for(int ilev = minlv; ilev <= maxlv; ilev++){
      strlev[ilev-minlv] << "]\n";
    }

    std::ostringstream str;

    str<<"import matplotlib.pyplot as plt\n";

    for(int ilev = minlv; ilev <= maxlv; ilev++){
      str << strlev[ilev-minlv].str() << "\n";
    }

    //str << "plt.figure(1)\n";

    //for(int ilev = minlv; ilev <= maxlv; ilev++){
    //  str << "plt.plot(["<<ilev<<"] * len(vals"<<ilev<<"), vals"<<ilev<<", 'o')\n";
    //}

    //str << "plt.show()\n";

    str << "plt.figure(1)\n";
    for(int ilev = minlv; ilev <= maxlv; ilev++){
      str << "plt.plot(["<<ilev<<"] * len(vals"<<ilev<<"), vals"<<ilev<<", 'o')\n";
    }
    str << "plt.yscale('log')\n";
    str << "plt.show()\n";


    std::ofstream f;
    std::string fname = args.param->outmPrefix + "DIRECT_graph" 
                        + std::to_string(args.niter) + ".py";
    f.open(fname.c_str(), std::ios::out);
    f << str.str();
    f.close();
  }


  // -- Scan potential points for the convex hull from left to right
  // Compare to segment given by (ilmin,flmin) -- (ilmax,flmax)
  // -> above: reject
  // -> under: accept and update ilmin,flmin

  CPRINTF1(" - scanning through potential hull pts with {} <= level <= {}\n",minlv,maxlv);



  if(args.niter > 4)
    CPRINTF1("## DEBUG ielem 14 ilev {} \n",args.ent2pol(14,idim+1));

  // Begin by dealing with the endpoints: these are always in the hull.
  // This loop can handle one (if only one total) or two end points
  for(int ilev = minlv; ilev <= maxlv; ilev += maxlv - minlv){
    INCVDEPTH(args.param);

    int ielem = args.lhull[ilev];
    METRIS_ASSERT(args.ent2pol(ielem,idim+1) == ilev);
    int ieglo = args.ent2pol(ielem,idim+2);
    if(ielem < 0){
      MPRINTF("## VERY STRANGE !\n");
      MPRINTF("niter = {} \n",args.niter);
      MPRINTF("hull is: {}\n",args.lhull);
      METRIS_THROW(TopoExcept())
    }


    int nele0 = args.ent2pol.get_n();
    if(idim == 2){
      aux_DIRECT_splittri(args,ielem,ieglo,ilev);
    }else{
      aux_DIRECT_splittet(args,ielem,ieglo,ilev);
    }

    aux_DIRECT_newevals(args, idim, nele0, ilev, ieglo, leval, peval, feval);

    CPRINTF1(" - ilev {} ielem {} split\n",ilev, ielem);

    if(maxlv == minlv) break;
  }

  CPRINTF1(" - split ends neval = {} \n",leval.get_n());

  if(nhull <= 2) return;

  double xlmin = 1.0;
  double xlmax = args.niter;
  double flmin = args.rhull[minlv];
  double flmax = args.rhull[maxlv];

  CPRINTF1(" - finishing convex hull\n");
  int nsplit = 0;
  for(int ilev = minlv + 1; ilev < maxlv; ilev++){
    INCVDEPTH(args.param);
    int ielem = args.lhull[ilev];
    if(ielem < 0) continue;

    int ieglo = args.ent2pol(ielem,idim+2);

    double xlcur = ilev;
    double flcur = args.rhull[ilev];

    // -- Compute area of triangle Xmin Xcur Xmax 
    // if trigonometric, then Xcur belongs to convex hull
    double artri = (xlmin - xlmax)*(flcur-flmax)
                 - (flmin - flmax)*(xlcur-xlmax);

    if(artri < -1.0e16) continue;

    nsplit++;

    CPRINTF1(" - ilev {} ielem {} split\n",ilev,ielem);

    // -- In this case, add the triangle for splitting and update 
    // left hull point

    xlmin = xlcur;
    flmin = flcur;

    int nele0 = args.ent2pol.get_n();
    if(idim == 2){
      aux_DIRECT_splittri(args,ielem,ieglo,ilev);
    }else{
      aux_DIRECT_splittet(args,ielem,ieglo,ilev);
    }

    aux_DIRECT_newevals(args, idim, nele0, ilev, ieglo, leval, peval, feval);

  }// for ilev
  CPRINTF1(" - split ends neval = {} ; total nsplit {} / {}\n",leval.get_n(),nsplit,nhull);

  return;
}


void aux_DIRECT_splittri(DIBLOB_args &args, int ielem,int ieglo,int ilev){

  int npoi0 = args.coorl.get_n();

  for(int ii = 0; ii < 3; ii++){
    int ip1 = args.ent2pol(ielem,lnoed2[ii][0]);
    int ip2 = args.ent2pol(ielem,lnoed2[ii][1]);
    args.coorl.inc_n();
    for(int kk = 0; kk < 2; kk++)
      args.coorl(npoi0+ii,kk) = (args.coorl(ip1,kk) + args.coorl(ip2,kk))/2.0;
  }

  int ip1 = args.ent2pol(ielem,0);
  int ip2 = args.ent2pol(ielem,1);
  int ip3 = args.ent2pol(ielem,2);


  args.ent2pol(ielem,0) = npoi0 + 0;
  args.ent2pol(ielem,1) = npoi0 + 1;
  args.ent2pol(ielem,2) = npoi0 + 2;
  args.ent2pol(ielem,3) = ilev  + 1;

  int nele0 = args.ent2pol.get_n();

  args.ent2pol.inc_n();
  args.ent2pol(nele0+0,0) = ip1;
  args.ent2pol(nele0+0,1) = npoi0 + 2;
  args.ent2pol(nele0+0,2) = npoi0 + 1;
  args.ent2pol(nele0+0,3) = ilev  + 1;
  args.ent2pol(nele0+0,4) = ieglo;
  args.fminhist.stack(args.fminhist[ielem]);

  args.ent2pol.inc_n();
  args.ent2pol(nele0+1,0) = ip2;
  args.ent2pol(nele0+1,1) = npoi0 + 0;
  args.ent2pol(nele0+1,2) = npoi0 + 2;
  args.ent2pol(nele0+1,3) = ilev  + 1;
  args.ent2pol(nele0+1,4) = ieglo;
  args.fminhist.stack(args.fminhist[ielem]);

  args.ent2pol.inc_n();
  args.ent2pol(nele0+2,0) = ip3;
  args.ent2pol(nele0+2,1) = npoi0 + 1;
  args.ent2pol(nele0+2,2) = npoi0 + 0;
  args.ent2pol(nele0+2,3) = ilev  + 1;
  args.ent2pol(nele0+2,4) = ieglo;
  args.fminhist.stack(args.fminhist[ielem]);

  return;
}



void aux_DIRECT_splittet(DIBLOB_args &args, int ielem,int ieglo,int ilev){

  int ityp = 2;

  // New split: dumb but covers the whole tet. 
  // 1 tet shares the centroid, then 4 tets that overlap with this one. 

  if(ityp == 1){

    int npoi0 = args.coorl.get_n();

    // Face points for the middle tetra
    for(int ii = 0; ii < 4; ii++){
      int ip1 = args.ent2pol(ielem,lnofa3[ii][0]);
      int ip2 = args.ent2pol(ielem,lnofa3[ii][1]);
      int ip3 = args.ent2pol(ielem,lnofa3[ii][2]);
      args.coorl.inc_n();
      for(int kk = 0; kk < 3; kk++)
        args.coorl(npoi0+ii,kk) = (args.coorl(ip1,kk) + args.coorl(ip2,kk) + args.coorl(ip3,kk))/3.0;
    }

    // Middle point for the other tets
    int ipmid = args.coorl.get_n();
    int ip1 = args.ent2pol(ielem,0);
    int ip2 = args.ent2pol(ielem,1);
    int ip3 = args.ent2pol(ielem,2);
    int ip4 = args.ent2pol(ielem,3);

    args.ent2pol(ielem,0) = npoi0 + 1;
    args.ent2pol(ielem,1) = npoi0 + 0;
    args.ent2pol(ielem,2) = npoi0 + 2;
    args.ent2pol(ielem,3) = npoi0 + 3;
    args.ent2pol(ielem,4) = ilev  + 1;


    npoi0 = args.coorl.get_n();
    args.coorl.inc_n();
    for(int kk = 0; kk < 3; kk++)
      args.coorl(npoi0,kk) = ( args.coorl(ip1,kk) + args.coorl(ip2,kk) 
                             + args.coorl(ip3,kk) + args.coorl(ip4,kk))/4.0;


    int nele0 = args.ent2pol.get_n();

    args.ent2pol.inc_n();
    args.ent2pol(nele0+0,0) = ipmid;
    args.ent2pol(nele0+0,1) = ip2;
    args.ent2pol(nele0+0,2) = ip3;
    args.ent2pol(nele0+0,3) = ip4;
    args.ent2pol(nele0+0,4) = ilev + 1;
    args.ent2pol(nele0+0,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+1,0) = ip1;
    args.ent2pol(nele0+1,1) = ipmid;
    args.ent2pol(nele0+1,2) = ip3;
    args.ent2pol(nele0+1,3) = ip4;
    args.ent2pol(nele0+1,4) = ilev + 1;
    args.ent2pol(nele0+1,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+2,0) = ip1;
    args.ent2pol(nele0+2,1) = ip2;
    args.ent2pol(nele0+2,2) = ipmid;
    args.ent2pol(nele0+2,3) = ip4;
    args.ent2pol(nele0+2,4) = ilev + 1;
    args.ent2pol(nele0+2,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+3,0) = ip1;
    args.ent2pol(nele0+3,1) = ip2;
    args.ent2pol(nele0+3,2) = ip3;
    args.ent2pol(nele0+3,3) = ipmid;
    args.ent2pol(nele0+3,4) = ilev + 1;
    args.ent2pol(nele0+3,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

  }else{// if ityp == 1

    // Same as previous but 1 + 12 instead of 1 + 4
    // Should not go to flat as quickly. 

    int npoi0 = args.coorl.get_n();

    // Face points for the middle tetra
    int ipf[4];
    for(int ii = 0; ii < 4; ii++){
      int ip1 = args.ent2pol(ielem,lnofa3[ii][0]);
      int ip2 = args.ent2pol(ielem,lnofa3[ii][1]);
      int ip3 = args.ent2pol(ielem,lnofa3[ii][2]);
      ipf[ii] = args.coorl.get_n();
      args.coorl.inc_n();
      for(int kk = 0; kk < 3; kk++)
        args.coorl(npoi0+ii,kk) = (args.coorl(ip1,kk) + args.coorl(ip2,kk) + args.coorl(ip3,kk))/3.0;
    }

    // Middle point for the other tets
    int ipmid = args.coorl.get_n();
    int ip1 = args.ent2pol(ielem,0);
    int ip2 = args.ent2pol(ielem,1);
    int ip3 = args.ent2pol(ielem,2);
    int ip4 = args.ent2pol(ielem,3);

    args.ent2pol(ielem,0) = ipf[1];
    args.ent2pol(ielem,1) = ipf[0];
    args.ent2pol(ielem,2) = ipf[2];
    args.ent2pol(ielem,3) = ipf[3];
    args.ent2pol(ielem,4) = ilev  + 1;


    npoi0 = args.coorl.get_n();
    args.coorl.inc_n();
    for(int kk = 0; kk < 3; kk++)
      args.coorl(npoi0,kk) = ( args.coorl(ip1,kk) + args.coorl(ip2,kk) 
                             + args.coorl(ip3,kk) + args.coorl(ip4,kk))/4.0;


    int nele0 = args.ent2pol.get_n();

    args.ent2pol.inc_n();
    args.ent2pol(nele0+0,0) = ipmid;
    args.ent2pol(nele0+0,1) = ipf[0];
    args.ent2pol(nele0+0,2) = ip3;
    args.ent2pol(nele0+0,3) = ip4;
    args.ent2pol(nele0+0,4) = ilev + 1;
    args.ent2pol(nele0+0,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+1,0) = ipmid;
    args.ent2pol(nele0+1,1) = ip2;
    args.ent2pol(nele0+1,2) = ipf[0];
    args.ent2pol(nele0+1,3) = ip4;
    args.ent2pol(nele0+1,4) = ilev + 1;
    args.ent2pol(nele0+1,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+2,0) = ipmid;
    args.ent2pol(nele0+2,1) = ip2;
    args.ent2pol(nele0+2,2) = ip3;
    args.ent2pol(nele0+2,3) = ipf[0];
    args.ent2pol(nele0+2,4) = ilev + 1;
    args.ent2pol(nele0+2,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);



    args.ent2pol.inc_n();
    args.ent2pol(nele0+3,0) = ipf[1];
    args.ent2pol(nele0+3,1) = ipmid;
    args.ent2pol(nele0+3,2) = ip3;
    args.ent2pol(nele0+3,3) = ip4;
    args.ent2pol(nele0+3,4) = ilev + 1;
    args.ent2pol(nele0+3,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+4,0) = ip1;
    args.ent2pol(nele0+4,1) = ipmid;
    args.ent2pol(nele0+4,2) = ipf[1];
    args.ent2pol(nele0+4,3) = ip4;
    args.ent2pol(nele0+4,4) = ilev + 1;
    args.ent2pol(nele0+4,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+5,0) = ip1;
    args.ent2pol(nele0+5,1) = ipmid;
    args.ent2pol(nele0+5,2) = ip3;
    args.ent2pol(nele0+5,3) = ipf[1];
    args.ent2pol(nele0+5,4) = ilev + 1;
    args.ent2pol(nele0+5,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);




    args.ent2pol.inc_n();
    args.ent2pol(nele0+6,0) = ipf[2];
    args.ent2pol(nele0+6,1) = ip2;
    args.ent2pol(nele0+6,2) = ipmid;
    args.ent2pol(nele0+6,3) = ip4;
    args.ent2pol(nele0+6,4) = ilev + 1;
    args.ent2pol(nele0+6,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+7,0) = ip1;
    args.ent2pol(nele0+7,1) = ipf[2];
    args.ent2pol(nele0+7,2) = ipmid;
    args.ent2pol(nele0+7,3) = ip4;
    args.ent2pol(nele0+7,4) = ilev + 1;
    args.ent2pol(nele0+7,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+8,0) = ip1;
    args.ent2pol(nele0+8,1) = ip2;
    args.ent2pol(nele0+8,2) = ipmid;
    args.ent2pol(nele0+8,3) = ipf[2];
    args.ent2pol(nele0+8,4) = ilev + 1;
    args.ent2pol(nele0+8,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);



    args.ent2pol.inc_n();
    args.ent2pol(nele0+9,0) = ipf[3];
    args.ent2pol(nele0+9,1) = ip2;
    args.ent2pol(nele0+9,2) = ip3;
    args.ent2pol(nele0+9,3) = ipmid;
    args.ent2pol(nele0+9,4) = ilev + 1;
    args.ent2pol(nele0+9,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+10,0) = ip1;
    args.ent2pol(nele0+10,1) = ipf[3];
    args.ent2pol(nele0+10,2) = ip3;
    args.ent2pol(nele0+10,3) = ipmid;
    args.ent2pol(nele0+10,4) = ilev + 1;
    args.ent2pol(nele0+10,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

    args.ent2pol.inc_n();
    args.ent2pol(nele0+11,0) = ip1;
    args.ent2pol(nele0+11,1) = ip2;
    args.ent2pol(nele0+11,2) = ipf[3];
    args.ent2pol(nele0+11,3) = ipmid;
    args.ent2pol(nele0+11,4) = ilev + 1;
    args.ent2pol(nele0+11,5) = ieglo;
    args.fminhist.stack(args.fminhist[ielem]);

  }


  #if 0
  // Old split, wrong? Check 
  int npoi0 = args.coorl.get_n();

  for(int ii = 0; ii < 4; ii++){
    int ip1 = args.ent2pol(ielem,lnofa3[ii][0]);
    int ip2 = args.ent2pol(ielem,lnofa3[ii][1]);
    int ip3 = args.ent2pol(ielem,lnofa3[ii][2]);
    args.coorl.inc_n();
    for(int kk = 0; kk < 3; kk++)
      args.coorl(npoi0+ii,kk) = (args.coorl(ip1,kk) + args.coorl(ip2,kk) + args.coorl(ip3,kk))/3.0;
  }

  int ip1 = args.ent2pol(ielem,0);
  int ip2 = args.ent2pol(ielem,1);
  int ip3 = args.ent2pol(ielem,2);
  int ip4 = args.ent2pol(ielem,3);

  args.ent2pol(ielem,0) = npoi0 + 0;
  args.ent2pol(ielem,1) = npoi0 + 2;
  args.ent2pol(ielem,2) = npoi0 + 3;
  args.ent2pol(ielem,3) = npoi0 + 1;
  args.ent2pol(ielem,4) = ilev  + 1;

  int nele0 = args.ent2pol.get_n();

  args.ent2pol.inc_n();
  args.ent2pol(nele0+0,0) = ip3;
  args.ent2pol(nele0+0,1) = ip4;
  args.ent2pol(nele0+0,2) = npoi0 + 0;
  args.ent2pol(nele0+0,3) = npoi0 + 1;
  args.ent2pol(nele0+0,4) = ilev + 1;
  args.ent2pol(nele0+0,5) = ieglo;
  args.fminhist.inc_n();
  args.fminhist[nele0+0] = args.fminhist[ielem];

  args.ent2pol.inc_n();
  args.ent2pol(nele0+1,0) = ip1;
  args.ent2pol(nele0+1,1) = ip4;
  args.ent2pol(nele0+1,2) = npoi0 + 2;
  args.ent2pol(nele0+1,3) = npoi0 + 1;
  args.ent2pol(nele0+1,4) = ilev  + 1;
  args.ent2pol(nele0+1,5) = ieglo;
  args.fminhist.inc_n();
  args.fminhist[nele0+1] = args.fminhist[ielem];

  args.ent2pol.inc_n();
  args.ent2pol(nele0+2,0) = ip1;
  args.ent2pol(nele0+2,1) = ip3;
  args.ent2pol(nele0+2,2) = npoi0 + 3;
  args.ent2pol(nele0+2,3) = npoi0 + 1;
  args.ent2pol(nele0+2,4) = ilev  + 1;
  args.ent2pol(nele0+2,5) = ieglo;
  args.fminhist.inc_n();
  args.fminhist[nele0+2] = args.fminhist[ielem];

  args.ent2pol.inc_n();
  args.ent2pol(nele0+3,0) = ip3;
  args.ent2pol(nele0+3,1) = npoi0 + 0;
  args.ent2pol(nele0+3,2) = ip2;
  args.ent2pol(nele0+3,3) = npoi0 + 3;
  args.ent2pol(nele0+3,4) = ilev  + 1;
  args.ent2pol(nele0+3,5) = ieglo;
  args.fminhist.inc_n();
  args.fminhist[nele0+3] = args.fminhist[ielem];

  args.ent2pol.inc_n();
  args.ent2pol(nele0+4,0) = ip1;
  args.ent2pol(nele0+4,1) = ip2;
  args.ent2pol(nele0+4,2) = npoi0 + 3;
  args.ent2pol(nele0+4,3) = npoi0 + 2;
  args.ent2pol(nele0+4,4) = ilev  + 1;
  args.ent2pol(nele0+4,5) = ieglo;
  args.fminhist.inc_n();
  args.fminhist[nele0+4] = args.fminhist[ielem];

  args.ent2pol.inc_n();
  args.ent2pol(nele0+5,0) = ip4;
  args.ent2pol(nele0+5,1) = npoi0 + 2;
  args.ent2pol(nele0+5,2) = ip2;
  args.ent2pol(nele0+5,3) = npoi0 + 0;
  args.ent2pol(nele0+5,4) = ilev  + 1;
  args.ent2pol(nele0+5,5) = ieglo;
  args.fminhist.inc_n();
  args.fminhist[nele0+5] = args.fminhist[ielem];

  return;
  #endif
}



void aux_DIRECT_newevals(DIBLOB_args &args, int idim, int nele0, int ilev, int ieglo,
                         intAr2 &leval, dblAr2 &peval, dblAr1 &feval){
  GETVDEPTH(args.param);

  for(int ielem = nele0; ielem < args.ent2pol.get_n(); ielem++){
    int ieval = peval.get_n();
    CPRINTF1("   - ask new eval iele loc {} glo {} lev {} split {} nodes {}\n", 
             ielem, ieglo,ilev,ielem,intAr1(idim+1,args.ent2pol[ielem]));
    peval.inc_n();
    leval.inc_n();
    feval.inc_n();
    double sum = 0;
    for(int jj = 0; jj < idim; jj++){
      peval(ieval,jj) = args.coorl(args.ent2pol(ielem, 0),jj) 
                      + args.coorl(args.ent2pol(ielem, 1),jj) 
                      + args.coorl(args.ent2pol(ielem, 2),jj);
      if(idim == 3) peval(ieval,jj) += args.coorl(args.ent2pol(ielem, 3),jj);
      peval(ieval,jj) /= (idim + 1);
      sum += peval(ieval,jj);
    }
    peval(ieval,idim) = 1 - sum; // for convenience
    leval(ieval,0) = ielem;
    leval(ieval,1) = ieglo;
  }
}








}// end namespace
