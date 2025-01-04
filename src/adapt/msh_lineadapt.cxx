//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_lineadapt.hxx"
#include "msh_lineforce.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../low_lenedg.hxx"
#include "../msh_structs.hxx"
#include "../low_geo.hxx"
#include "../aux_utils.hxx"
#include "../msh_checktopo.hxx"
#include "../aux_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../adapt/low_increasecav.hxx"
#include "../cavity/msh_cavity.hxx"


namespace Metris{

// Adapt lines "frontally" using CAD
template<class MFT>
void adaptGeoLines(Mesh<MFT> &msh, int ithrd1, int ithrd2){
  int iverb = msh.param->iverb;
  if(!msh.CAD()) return;

  if(iverb >= 1) printf("-- Start adaptGeoLines.\n");

  if(msh.param->dbgfull)  check_topo(msh);

  //printf("## DEBUG corner point 2 print all bpois\n");
  //for(int ibpoi = msh.poi2bpo[2]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
  //  printf("%d = %f : ",ibpoi,msh.bpo2rbi[ibpoi][0]);
  //  intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
  //}

  MetSpace ispac0 = msh.met.getSpace();
  msh.met.setSpace(MetSpace::Exp);

  METRIS_ASSERT(msh.CAD.ncaded > 0); 

  //const double tola = 1.0e-14;
  const double edgtol = 0.1;
  const int medgit = 100;
  // Absolute tolerance for whether a point is on the geometry or not 
  // (temporary)
  const double geotol = msh.param->geo_abstoledg;
  // Multiply and divide target length by this factor for admissible new
  // edge length interval 
  const double lentolfac = msh.param->geo_lentolfac;
  if(lentolfac < 1.0){
    METRIS_THROW_MSG(WArgExcept(),"lentolfac < 1! = "<<lentolfac);
  }
  // Maximum number of iterations for length bisection 
  const int miter_bisection = 1000;

  RoutineWorkMemory iwrk(msh.iwrkmem);

  // From each ref, point onto a corner that has an edge, of this ref, inciding
  intAr1 ref2cor(msh.CAD.ncaded,iwrk.allocate(msh.CAD.ncaded));

  for(int iref = 0; iref < msh.CAD.ncaded; iref++) ref2cor[iref] = -1;

  int nseen = 0;
  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;

    int iref = msh.edg2ref[iedge];
    if(ref2cor[iref] >= 0) continue;

    bool icnt = false;
    for(int ii = 0; ii < 2; ii++){
      int ip = msh.edg2poi(iedge,ii);
      int ib = msh.poi2bpo[ip];
      if(msh.bpo2ibi[ib][1] == 0){
        nseen++;
        ref2cor[iref] = ip;
        icnt = true;
        break;
      }
    }

    if(nseen == msh.CAD.ncaded) break;

    if(icnt) continue;
  }

  METRIS_ASSERT_MSG(nseen == msh.CAD.ncaded, "nseen = "<<nseen<<" ncaded = "
    <<msh.CAD.ncaded);


  //printf("## DEBUG point 1 ibpoi & stuff\n");
  //for(int ibpoi = msh.poi2bpo[1]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
  //  printf("%d t = %f : ",ibpoi,msh.bpo2rbi(ibpoi,0));
  //  intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
  //}
  //wait();


  // Start by forcing all points on the geometry 
  reinsertLines<MFT>(msh,ithrd1,ithrd2);
  if(iverb >= 2) writeMesh("lineforce",msh);

  dblAr1 crv_lens(msh.CAD.ncaded);
  getCADCurveLengths(msh, (lentolfac - 1.0), crv_lens);


  int mcfac = 100, mcedg = 10;
  MshCavity cav(0,mcfac,mcedg);

  CavOprOpt opts;
  CavOprInfo info;
  CavWrkArrs work;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = true;
  opts.allow_remove_points = true;
  opts.dryrun = false;
  opts.iverb  = iverb;
  opts.geodev1 = 1.0; // lax


  int CADedgtag = 1;
  std::map<ego,int> CADedg2tag, CADedg2ref;
  for(int iref = 0; iref < msh.CAD.ncaded; iref++){
    ego edg = msh.CAD.cad2edg[iref];
    CADedg2ref[edg] = iref;
    CADedg2tag[edg] = 0;
  }

  int nEGrro = 0;
  for(int iloop = 0; iloop < msh.CAD.ncadlp; iloop++){

    // Make a map of edge ego orientations 
    ego loop = msh.CAD.cad2lop[iloop], *lchild, geom; 
    int oclass,mtype,nchild,*senses;
    int ierro = EG_getTopology(loop,&geom,&oclass,&mtype,NULL,
                           &nchild,&lchild,&senses);
    if(ierro != 0){
      print_EGADS_error("EG_getTopology (LOOP)",ierro);
      METRIS_THROW(TopoExcept());
    }
    METRIS_ENFORCE_MSG(nchild == msh.CAD.ncaded || msh.CAD.ncadlp > 1," nchild = "<<nchild);
    std::map<ego,int> edgorient;
    for(int ii = 0; ii < nchild; ii++){
      ego edg = lchild[ii];
      edgorient[edg] = senses[ii];
    }


    double result[18],  sz[2], len, nrmal[3];
    int nloop;
    // Loop over CAD edges and remesh each one 
    for(int iCADed = 0; iCADed < nchild; iCADed++){
      ego obj = lchild[iCADed];
      int iref = CADedg2ref[obj];
      if(CADedg2tag[obj] >= CADedgtag) continue;

      if(iverb >= 1) printf("\n   - Loop %d adapt line %d / %d \n", 
                            iloop, iref+1, msh.CAD.ncaded);

      // Initially we want length one. As we discretize and converge to the curve
      // length, we can pinpoint the optimal spacing. 
      //double tarlen = 1.0;
      double tarlen = crv_lens[iref] / round(crv_lens[iref]);

      // STEP 0: check if edge is already correctly meshed. This also kickstarts
      // tarlen
      bool doskip = true;

      // icor0: seed corner for this CAD edge. 
      int icor0 = ref2cor[iref];
      int ibcr0 = msh.poi2bpo[icor0];
      METRIS_ASSERT(ibcr0 >= 0);
      ibcr0 = getref2bpo(msh,ibcr0,iref,1);
      int iedc0 = msh.bpo2ibi[ibcr0][2];
      METRIS_ASSERT(iedc0 >= 0);

      int iedge = iedc0;
      int ipoi0 = icor0;
      // Compute curve length, as well as min / max edges length. 
      // Also evaluate edge endpoints on CAD edge and compute distance 
      double min_len = 1.0e30, max_len = -1.0, avg_len = 0.0;
      //double crv_len = 0;
      int npavg = 0;
      while(true){

        bool ifin = false;
        int iver = getveredg<1>(iedge,msh.edg2poi,ipoi0);
        // Next point:
        ipoi0 = msh.edg2poi(iedge,1-iver);
        // Evaluate 
        int ibpo0 = msh.poi2bpo[ipoi0];
        if(msh.bpo2ibi(ibpo0,1) == 0){
          ifin = true;
          // If corner, get the bpo for this edge. 
          ibpo0 = getent2bpo(msh,ibpo0,iedge,1);
        }
        ierro = EG_evaluate(obj, msh.bpo2rbi[ibpo0], result);
        double dst; 
        if(msh.idim == 2){
          dst = geterrl2<2>(msh.coord[ipoi0], result);
        }else{
          dst = geterrl2<3>(msh.coord[ipoi0], result);
        }

        // If distance is too large, points are not on the geometry, get on with
        // the normal procedure. 
        if(dst > geotol*geotol){
          if(iverb >= 2) printf("   - Point %d not on geometry, dist = %15.7e > "
                                      "%15.7e = tol\n", ipoi0, sqrt(dst), geotol);
          doskip = false;
          //break;
        }

        if(msh.idim == 2){
          len = getlenedg_geosz<MFT,2,1>(msh,msh.edg2poi[iedge],sz);
        }else{
          len = getlenedg_geosz<MFT,3,1>(msh,msh.edg2poi[iedge],sz);
        }

        if(iverb >= 2) printf("     - initial edge %d (%d,%d) length %15.7e\n",
          iedge,msh.edg2poi(iedge,0),msh.edg2poi(iedge,1),len);

        iedge = msh.edg2edg(iedge, iver);

        //crv_len += len;
        avg_len += len;
        npavg++;

        min_len = MIN(min_len, len);
        max_len = MAX(max_len, len);

        if(ifin) break;
      }

      // If, so far, no reason for skipping has been found, check with tarlen and
      // edge lengths. 
      if(doskip){
        avg_len /= npavg; 
        //tarlen = crv_len / round(crv_len);
        double lentolabs = tarlen * (lentolfac - 1.0); 
        // Check that the min and max of the edge lengths are within the tolerance
        // to tarlen. 
        doskip = min_len >= tarlen - lentolabs  
              && max_len <= tarlen + lentolabs;
        //doskip = min_len >= tarlen / lentolfac  
        //      && max_len <= tarlen * lentolfac;
        if(iverb >= 2) printf("   - Points on geometry, now check length: min = %f" 
          " avg = %f max = %f new tarlen = %f -> doskip = %d crv_len = %f  \n", 
          min_len, avg_len, max_len, tarlen, doskip, crv_lens[iref]);
      }

      // Move on to next reference. 
      if(doskip) continue;

      // Get CAD parameter range
      double range[2];
      int iperi;
      ierro = EG_getRange(obj,range,&iperi);
      if(ierro != 0){
        printf("  ## EG_getRange failed %d \n",ierro);
        METRIS_ASSERT(ierro == 0);
        continue;
      }
      if(range[0] > range[1]){
        double tmp = range[0];
        range[0] = range[1];
        range[1] = tmp;
      }

      if(iverb >= 2) printf("   - range %f , %f \n",range[0],range[1]);


      for(int nedgit = 0; nedgit < medgit; nedgit++){

        // We'll adjust the target length depending on the lengths effectively obtained.
        // This is to avoid propagating error until the last point. 
        double adjusted_tarlen = tarlen;
        double lentolabs = tarlen * (lentolfac - 1.0); 

        // This stuff has to enter the loop because the ibcr0 changes, also iedc0
        ibcr0 = msh.poi2bpo[icor0];
        // Get the corresponding ib on the correct tdimn 1 ref 
        ibcr0 = getref2bpo(msh,ibcr0,iref,1);

        METRIS_ASSERT(ibcr0 >= 0);

        // Seed edge 
        iedc0 = msh.bpo2ibi[ibcr0][2];
        METRIS_ASSERT(iedc0 >= 0);

        if(iverb >= 2) printf("   - iter %d seed corner %d (t = %f) edge = %d (%d, %d) tarlen = %f\n",
          nedgit,icor0,msh.bpo2rbi[ibcr0][0],iedc0,msh.edg2poi(iedc0,0),msh.edg2poi(iedc0,1),tarlen);

        METRIS_ASSERT_MSG(msh.edg2poi(iedc0,0) == icor0 || msh.edg2poi(iedc0,1) == icor0,
          "Corner not in seed edge ! ib = "<<ibcr0<<" entries = "
          <<msh.bpo2ibi(ibcr0,0)<<" " <<msh.bpo2ibi(ibcr0,1)<<" "
          <<msh.bpo2ibi(ibcr0,2)<<" " <<msh.bpo2ibi(ibcr0,3)<<"\n");


        // Which end is the corner?
        int irnge = -1;
        double drnge = abs(range[1] - range[0]);
        for(int ii = 0; ii < 2; ii++){
          if(abs(range[ii] - msh.bpo2rbi[ibcr0][0]) < 1.0e-6 * drnge) irnge = ii;
        }
        METRIS_ENFORCE_MSG(irnge != -1,"## CORNERS IN MESH HAVE WRONG CAD EDGE "
        "PARAMETRIC COORDINATES !\n" <<" icor = "<<icor0<<" range = "<<range
        [0]<<" - "<<range[1]<<" this t = "<<msh.bpo2rbi[ibcr0][0]);


        if(iverb >= 2) printf("   - CAD edge %d inner iteration %d tarlen = %f \n",
                              iref,nedgit,adjusted_tarlen);

        // Error on existing edge points, used to determine whether re-iteration
        // is needed
        double erredg = 0;
        int    npterr = 0;

        // Departing from ipoi0, with edge see iedg0, generate insertions /
        // collapses until reaching other end of CAD edge. 
        int ipoi0 = icor0; 
        int ibpo0 = ibcr0;
        int iedg0 = iedc0;

        //crv_len = 0;
        bool   continue_ref = true;
        int    npins        = 0;
        while(continue_ref){

          // Initialize / reset cavity
          cav.lcfac.set_n(0);
          cav.lcedg.set_n(0);
          cav.ipins = -1;
          cav.nrmal = NULL;

          // Tag is shared with aux_walk_line
          msh.tag[ithrd1]++;

          // Walk from edge to edge until the target length is reached. 
          bool ifin, inewp;
          int iedge, iface, ibpnw;
          // End points of edge to be considered (computed the len of)
          int edg2pol[2];
          ierro = aux_walk_line(msh, cav, obj, 
                                ipoi0, iedg0, adjusted_tarlen,
                                edg2pol, &ibpnw, &iedge, &iface, &len, sz, &ifin,
                                &erredg, &npterr, &nEGrro, ithrd1);
          if(ierro != 0) goto cleanup1;

          // ---- An insertion or reinsertion is going to be carried out. 
          // ifin just means last seen point is a corner. But perhaps there 
          // is in fact room for a new point. 
          inewp = len > adjusted_tarlen*sqrt(2) || !ifin;


          if(inewp){ // New point

            npins++;
            // Bisection: create new point, put into ipins, return final len and
            // update adjusted_tarlen. 
            ierro = gen_newp_line(msh,cav,obj,ibpo0,ibpnw,iedge,iface,edg2pol,sz,
                                  miter_bisection,tarlen,lentolabs,
                                  &adjusted_tarlen,&len,&nEGrro);
            if(ierro != 0) goto cleanup1;

          }else{ // Both len <= tar and ifin, then reinsert corner with this cavity. 

            cav.ipins = edg2pol[1];
            continue_ref = false;
            if(iverb >= 2) printf("   - Terminating edge insertions for iref = %d."
                            "Inserted %d points. Last length %f error %15.7e  \n",
                            iref,npins,len, len-tarlen);
            if(iverb >= 2){
              std::string fname = "debug_lineadapt" + std::to_string(iref) + ".meshb";
              writeMesh(fname,msh);
            }

          }

          if(iverb >= 4)
            writeMeshCavity("debug_lineadap0_cav",msh,cav);

          // Proceed to insertion We have our ipins, edge cavity also. Now extend
          // triangle cavity from edg2fac seeds
          nloop = 0;
          do{
            ierro = increase_cavity2D(msh,msh.fac2ref[iface],msh.coord[cav.ipins],
                                      opts,cav,ithrd1);
            nloop++;
            if(nloop > 10) METRIS_THROW_MSG(TopoExcept(), "Too many cavity increases !");
          }while(ierro == 2);

          increase_cavity_Delaunay(msh, cav, cav.ipins, ithrd1);

          if(iverb >= 4)
            writeMeshCavity("debug_lineadap1_cav",msh,cav);


          if(iverb >= 3){
            printf(" insert ipins = %d (%f,%f) list lcedg, lcfac:\n",cav.ipins,
              msh.coord(cav.ipins,0), msh.coord(cav.ipins,1));
            cav.lcedg.print();
            cav.lcfac.print();
          }


          // If the point is a corner, we probably added triangles against edges
          // of the other corner refs. We need to add them to the cavity 
          // Note: OR NOT !
          if(!inewp){
            
            #if 0
            METRIS_ASSERT(ifin);
            int ib = msh.poi2bpo[cav.ipins];
            METRIS_ASSERT(msh.bpo2ibi[ib][1] == 0);
            ib = msh.bpo2ibi[ib][3];
            METRIS_ASSERT(msh.bpo2ibi[ib][1] == 1);
            // Tag all the refs
            do{
              if(msh.bpo2ibi[ib][1] == 1){
                int itmp = msh.bpo2ibi[ib][2];
                int iref2 = msh.edg2ref[itmp];
                if(iref2 != iref) msh.ced2tag[ithrd1][iref2] = msh.tag[ithrd1];
              }
              ib = msh.bpo2ibi[ib][3];
            }while(ib >= 0);


            for(int ifacl = 0; ifacl < cav.lcfac.get_n(); ifacl++){
              int ifac2 = cav.lcfac[ifacl];
              for(int ii = 0; ii < 3; ii++){
                int ip1 = msh.fac2poi(ifac2,lnoed2[ii][0]);
                int ip2 = msh.fac2poi(ifac2,lnoed2[ii][1]);
                int itmp = getedgglo(msh,ip1,ip2);
                if(itmp < 0) continue;
                int iref2 = msh.edg2ref[itmp];
                if(msh.ced2tag[ithrd1][iref2] >= msh.tag[ithrd1]){
                  cav.lcedg.stack(itmp); 
                  if(iverb >= 3) printf("   - Add edge %d to cavity (ref %d) \n",
                                                                      itmp,iref2);
                }
              }
            }
            #endif

          }
          

          if(msh.idim >= 3 && iface >= 0){
            getnorfacP1(msh.fac2poi[iface],msh.coord,nrmal);
            cav.nrmal = nrmal;
          }

          // If insertion fails, try inserting a Steiner point 
          // First insertion is regular. 
          for(int isteiner = 0; isteiner <= 1; isteiner++){

            if(iverb >= 3){
              printf("   - Starting insert ipins = %d cav ncedg = %d ncfac = %d \n",    
                                     cav.ipins,cav.lcedg.get_n(),cav.lcfac.get_n());
            }

            if(iverb >= 4){
              int *refold = (int *) malloc(msh.nface*sizeof(int));
              METRIS_ENFORCE(refold != NULL);
              for(int ii = 0; ii < msh.nface; ii++){
                refold[ii] = msh.fac2ref[ii];
                msh.fac2ref[ii] = 2;
              }

              for(int ii = 0; ii < cav.lcfac.get_n(); ii++){
                msh.fac2ref[cav.lcfac[ii]] = 3;
              }
              // Add a corner at ipoin 
              int ipoin = msh.newpoitopo(-1,-1);
              int ibpoi = msh.template newbpotopo<0>(ipoin,ipoin);
              for(int ii = 0; ii < msh.idim; ii++) 
                msh.coord[ipoin][ii] = msh.coord[cav.ipins][ii] ;
              writeMesh("debug_lineadap0.meshb",msh);
              for(int ii = 0; ii < nibi; ii++) msh.bpo2ibi[ibpoi][ii]  = -1;
              msh.set_npoin(msh.npoin-1);
              msh.set_nbpoi(msh.nbpoi-1);
              //printf("Wait here \n");
              //wait();
              for(int ii = 0; ii < msh.nface; ii++){
                msh.fac2ref[ii] = refold[ii];
              }
              free(refold);

              writeMeshCavity("debug_lineadap2_cav",msh,cav);
            }

            CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
              ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd1);
            }}CT_FOR1(ideg);

            if(ierro == 0 && isteiner == 1 && iverb >= 3){
              printf("   - Steiner reinsertion successful\n");
              //wait();
            }
            if(ierro == 0) break;


            // Insertion failed. If this is the Steiner try, throw error. 
            if(isteiner > 0){
              if(iverb >= 3){
                #ifndef NDEBUG
                printf("## Steiner insertion cavity error %d wait \n",ierro);
                if(msh.param->dbgfull) wait();
                #endif
              }
              goto cleanup1;
            }else if(iverb >= 3){
              printf("   - going to try Steiner insertion \n");
              #ifndef NDEBUG
              if(msh.param->dbgfull) wait();
              #endif
            }

            double norpoi[3];
            ierro = getnorpoiCAD(msh, cav.ipins, edgorient, norpoi);

            if(iverb >= 3){
              printf("   - got CAD = ");
              dblAr1(msh.idim, norpoi).print();
            }

            //msh.set_npoin(msh.npoin-1);
            //msh.set_nbpoi(msh.nbpoi-1);

            #if 0
            // Otherwise, let's try inserting a point along the normal to the 
            // point. 
            METRIS_ASSERT(msh.bpo2ibi[ibpnw][1] == 1);
            ierro = EG_evaluate(obj, msh.bpo2rbi[ibpnw], result);
            METRIS_ASSERT(ierro == EGADS_SUCCESS); 
            if(ierro != 0){
              if(iverb >= 3) printf("## EG_evaluate error %d \n",ierro);
              wait();
              goto cleanup1;
            }

            double *du = &result[3];
            METRIS_ENFORCE(msh.idim == 2);

            // ingoing normal 
            double dv[2]; 
            dv[0] = -du[1];
            dv[1] = du[0];

            if(iverb >= 4){
              for(int ii = 0; ii < 18 ;ii++){
                printf("resu %d = %f \n",ii,result[ii]);
              }
              printf("iref %d ipend %d \n",iref,edg2pol[0]);
            }
            #endif

            // Get an approximate edge length. We don't really care to be very 
            // close to 1 here. 
            double len2; 
            if(msh.met.getSpace() == MetSpace::Log){
              if(msh.idim == 2){
                len2 = getlenedg_log<2>(norpoi, msh.met[cav.ipins]);
              }else{
                len2 = getlenedg_log<3>(norpoi, msh.met[cav.ipins]);
              }
            }else{
              if(msh.idim == 2){
                len2 = getlenedg<2>(norpoi, msh.met[cav.ipins]);
              }else{
                len2 = getlenedg<3>(norpoi, msh.met[cav.ipins]);
              }
            }

            //int isens = edgorient[obj];
            if(iverb >= 3) 
              printf("   - Steiner nor = %f %f len2 = %f orient %d \n",
                                  norpoi[0],norpoi[1],len2,edgorient[obj]);

            int ipins_old = cav.ipins;
            cav.ipins = msh.newpoitopo(msh.get_tdim(), -1);

            int ncedg = cav.lcedg.get_n();
            cav.lcedg.set_n(0);

            bool stsuc = false;
            double step = 1.0;
            for(int iststp = 0; iststp < 10; iststp++){

              for(int ii = 0; ii < msh.idim; ii++){
                msh.coord(cav.ipins,ii) = msh.coord(ipins_old, ii)
                                        - step *  norpoi[ii] / len2;  // + isens *
              }

              if(iverb >= 3){
                printf("   - Steiner insertion attempt step %f coop = ",step);
                dblAr1(msh.idim,msh.coord[cav.ipins]).print();
              }

              // Interpolate in mesh interior -> msh.get_tdim()
              ierro = msh.interpMetBack(cav.ipins,msh.get_tdim(),iface,-1,NULL);

              if(ierro != 0){
                if(iverb >= 3) printf("## Steiner attempt interpmet error %d \n",ierro);
                // If this fails, relax down to 0 for a few iterations. Probably
                // overshot. 
                // Note: not any more since we're relaxing up for visibility 
                break;
              }


              // Remove the edges from the cavity. This is very easily done by 
              // setting lcedg.Set_n = 0 and restoring.


              ierro = increase_cavity2D(msh,msh.fac2ref[iface],msh.coord[cav.ipins],
                                        opts,cav,ithrd1);

              if(ierro != 0){
                step *= 1.5;
                if(iverb >= 3) printf("## Steiner attempt increase_cavity error %d \n",ierro);
                if(iverb >= 3) writeMeshCavity("debug_lineadap0_cav_Steiner"+
                                               std::to_string(iststp),msh,cav);
                continue;
              }


              CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
                ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd1);
              }}CT_FOR1(ideg);

              if(ierro != 0){
                step *= 1.5;
                if(iverb >= 3) printf("## Steiner attempt cavity error %d \n",ierro);
                continue; 
              }

              stsuc = true;
              // Restore the edges. Get rid of the triangles (they're dead)
              cav.lcedg.set_n(ncedg);
              cav.lcfac.set_n(0);
              cav.ipins = ipins_old;
              // Now get the faces from the edges 
              msh.tag[ithrd2]++;
              for(int iedgc : cav.lcedg){
                int ifac1 = msh.edg2fac[iedgc];
                METRIS_ASSERT(ifac1 >= 0);
                if(msh.fac2tag(ithrd2,ifac1) >= msh.tag[ithrd2]) continue;
                msh.fac2tag(ithrd2,ifac1) = msh.tag[ithrd2];
                cav.lcfac.stack(ifac1);
              }

              if(iverb >= 3) printf("   - Steiner insertion succeeded \n");
              if(iverb >= 4) writeMesh("debug_Steiner.meshb",msh);
              break;
            }

            if(!stsuc){
              if(iverb >= 3){
                printf("## Steiner attempt failed error %d \n",ierro);
                #ifndef NDEBUG
                wait();
                #endif
              }
              METRIS_ASSERT_MSG(cav.ipins == msh.npoin - 1, "(1) cav.ipins == msh.npoin - 1");
              msh.set_npoin(msh.npoin-1);
              // This ipins is not an nbpoi
              cav.ipins = ipins_old;
              goto cleanup1;
            }


            #if 0
            double metl[6];
            const int nnmet = (msh.idim*(msh.idim+1))/2;
            for(int ii = 0; ii < nnmet; ii++) metl[ii] = msh.met(cav.ipins,ii);
            if(msh.met.getSpace() == MetSpace::Log){
              double lmet[6];
              for(int ii = 0; ii < nnmet; ii++) lmet[ii] = metl[ii];
              if(msh.idim == 2){
                getexpmet_cpy<2>(lmet, metl);
              }else{
                getexpmet_cpy<3>(lmet, metl);
              }
            }

            double len2;
            if(msh.idim==2){
              len2 = getlenedgsq<2>(dv,metl);
            }else{
              len2 = getlenedgsq<3>(dv,metl);
            }
            #endif






          }


          if(msh.param->dbgfull) check_topo(msh);
          if(iverb >= 4){
            writeMesh("debug_lineadap1.meshb",msh);
          }

          if(!continue_ref) break;

          // -- Prepare next seeds.
          //  - start from ipins 
          ipoi0 = cav.ipins;
          ibpo0 = msh.poi2bpo[ipoi0];

          // Should be a line point, if corner we exited already (!continue_ref)
          METRIS_ASSERT_MSG(msh.bpo2ibi[ibpo0][1] == 1, 
            "corner despite continue_ref? itype = "<<msh.bpo2ibi[ibpo0][1]<<
            " ibpo0 = "<<ibpo0<<" ipoin = "<<msh.bpo2ibi(ibpo0,0));

          if(iverb >= 3) printf("   - Start from ipoi0 = %d ibpo0 = %d \n",
            ipoi0,ibpo0);

          iedg0 = -1;
          //  - edge seed is whichever of 2 last goes way from icor0
          for(int ii = 0; ii < 2; ii++){
            iedge = msh.nedge - 2 + ii;
            if(iverb >= 3) printf("    - check iedge %d = (%d,%d)\n",
              iedge,msh.edg2poi(iedge,0),msh.edg2poi(iedge,1));
            int ip; 
            if(msh.edg2poi(iedge,0) == cav.ipins){
              ip = msh.edg2poi(iedge,1);
            }else if(msh.edg2poi(iedge,1) == cav.ipins){
              ip = msh.edg2poi(iedge,0);
            }else{
              METRIS_THROW_MSG(TopoExcept(),"Check again edge indices after insertion");
            }
            int ib = msh.poi2bpo[ip];
            METRIS_ASSERT(ib >= 0);

            if(msh.bpo2ibi[ib][1] == 0){
              ib = getent2bpo(msh,ib,iedge,1);
            }
            METRIS_ASSERT(ib >= 0);

            double t = msh.bpo2rbi[ib][0];

            if(iverb >= 3) printf("    - Candidate point %d t = %f orig = %f \n",
                                   ip,t,msh.bpo2rbi(ibpo0,0));
            // Corner is range[irnge]
            // If this point is away from initial corner respective to ipins
            if(irnge == 0 && t > msh.bpo2rbi(ibpo0,0)
            || irnge == 1 && t < msh.bpo2rbi(ibpo0,0)){
              // Found our next seed: iedge
              iedg0 = iedge;
              break;
            }
          } 

          if(iedg0 < 0){
            printf("## iedg0 = %d \n",iedg0);
            printf("iedge = %d \n",iedge);
            printf("edg2 poi = ");
            intAr1(2,msh.edg2poi[iedge]).print();
            writeMesh("debug_iedg0", msh);
          }
          METRIS_ASSERT(iedg0 >= 0);


          continue;

          cleanup1: 
          if(cav.ipins >= 0 && !ifin){
            // This removes any debug points that were put here.
            //msh.set_nbpoi(msh.poi2bpo[cav.ipins]-1);
            //msh.set_npoin(cav.ipins-1);
            break;
          }

        } // End loop operations

        if(iverb >= 3){
          std::string fname = "geolines_ref" + std::to_string(iref) 
                            + "_iter" + std::to_string(nedgit);
          writeMesh(fname,msh);
        }

        // Compute new target length. If this varies, there'll be a restart no matter what. 
        // double tarle0 = tarlen;
        // Mesh in excess
        //tarlen = crv_len / round(crv_len);  

        double tarle0 = tarlen;

        // Now that the curve length is known well enough (from getCADCurveLengths)
        // the only source of tarlen correction is that the mesh is coarse
        // In that case, there might be an error on the last edge length. 
        // The tarlen is corrected to even that error out. 
        int nedcrv = crv_lens[iref] / tarlen;
        if(iverb >= 2) printf("   - End len %f correct tarlen %f -> %f \n",
                              len,tarlen,tarlen+(len-tarlen) / nedcrv);
        // Damping -> in the future, we need to be less aimless here 
        tarlen += (len-tarlen) / nedcrv / 10; 

        //break;

        //if(abs(tarlen - tarle0) > abs(1.0 - lentolfac)){
        if(abs(tarlen - len) > abs(1.0 - lentolfac)){
          if(iverb >= 1) printf("   - iter %d correct tarlen %f -> %f, err = %f" 
                        " > tol = %f \n",nedgit,tarle0,tarlen,len-tarlen,abs(1.0 - lentolfac));
          continue;
        }else if(iverb >= 1){
          if(iverb >= 1) printf("   - iter %d end line: |err| = %f < tol = %f\n",nedgit,
                                abs(len-tarlen),abs(1.0 - lentolfac));

        }

        // Restart by default, not if no pt inexact
        if(npterr == 0) break;
        erredg /= npterr;
        if(erredg < edgtol) break;


      } // End repeat same CAD edge; negit 


    } // End CAD edges

  } // End for iloop


  msh.met.setSpace(ispac0);
  msh.cleanup();

}

template void adaptGeoLines<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, int ithrd1, int ithrd2);
template void adaptGeoLines<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, int ithrd1, int ithrd2);


// Binary refinement until tolerance is met after split. 
// ref2cor should, for each edge reference, give its two end points. 
template<class MFT>
void getCADCurveLengths(Mesh<MFT> &msh, double tol, dblAr1 &crv_len){

  const int nref = msh.CAD.ncaded;
  const int iverb = msh.param->iverb; 
  crv_len.set_n(nref);

  MetSpace ispac0 = msh.met.getSpace();
  msh.met.setSpace(MetSpace::Log);

  if(iverb >= 1) printf(" - START getCADCurveLengths\n");

  double result[18];
  const int nnmet = (msh.idim*(msh.idim + 1)) / 2;

  // add two dummy points 
  int ipon[2]; 
  ipon[0] = msh.newpoitopo(-1,-1);
  ipon[1] = msh.newpoitopo(-1,-1);
  int edg2pol[2] = {ipon[0], ipon[1]};
  double sz[2];

  intAr1 ref2ned(nref);
  ref2ned.set_n(nref);
  ref2ned.fill(0);
  for(int iedge = 0; iedge < msh.nedge; iedge++){
    int iref = msh.edg2ref[iedge];
    ref2ned[iref] ++;
  }


  //for(int ii = 0; ii < nref; ii++) crv_len0[ii] = crv_len[ii];

  // Initialize crv_len using initial mesh. This is simply for normalization 
  // purposes. 
  // We also will need a copy, so the normalization doesn't depend on the 
  // tolerance or otherwise execution. 
  dblAr1 crv_len0(nref);
  crv_len0.set_n(nref);
  crv_len0.fill(0);
  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;
    int iref = msh.edg2ref[iedge]; 

    double len;
    if(msh.idim == 2){
      len = getlenedg_geosz<MFT,2,1>(msh,msh.edg2poi[iedge],sz);
    }else{
      len = getlenedg_geosz<MFT,3,1>(msh,msh.edg2poi[iedge],sz);
    }

    crv_len0[iref] += len;
  }

  if(iverb >= 1){
    for(int iref = 0; iref < nref; iref++)
      printf("   - line %d len0 = %f\n",iref, crv_len0[iref]);
  }


  // Refine each edge until the successive lengths are close enough 
  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;
    int iref = msh.edg2ref[iedge]; 
    ego obj = msh.CAD.cad2edg[iref];

    int ipoi1 = msh.edg2poi(iedge,0);
    int ipoi2 = msh.edg2poi(iedge,1);

    int ibpo1 = msh.poi2bpo[ipoi1];
    if(msh.bpo2ibi(ibpo1,1) == 0) ibpo1 = getent2bpo(msh,ibpo1,iedge,1);
    int ibpo2 = msh.poi2bpo[ipoi2];
    if(msh.bpo2ibi(ibpo2,1) == 0) ibpo2 = getent2bpo(msh,ibpo2,iedge,1);

    double t0_1 = msh.bpo2rbi(ibpo1, 0);
    double t0_2 = msh.bpo2rbi(ibpo2, 0); 

    if(iverb >= 3){
      printf("  - iedge %d ib1 %d ib2 %d t0_1 %f t0_2 %f \n",
                          iedge,ibpo1,ibpo2,t0_1,t0_2);
      printf("ipoi1 = %d ipoi2 = %d \n",ipoi1,ipoi2);
      printf("ib1 full bpo:\n");
      for(int ibpoi = msh.poi2bpo[ipoi1];ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
        printf(" ibpoi = %d t = %10.3e : ",ibpoi,msh.bpo2rbi(ibpoi,0));
        intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
      }
      printf("ib2 full bpo:\n");
      for(int ibpoi = msh.poi2bpo[ipoi2];ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
        printf(" ibpoi = %d t = %10.3e : ",ibpoi,msh.bpo2rbi(ibpoi,0));
        intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
      }
    }

    double du[3];
    for(int ii = 0; ii < msh.idim; ii++) du[ii] = msh.coord(ipoi2,ii)
                                                - msh.coord(ipoi1,ii);
    double dunrm = msh.idim == 2 ? getnrml2<2>(du) : getnrml2<3>(du);
    METRIS_ASSERT(dunrm >= 1.0e-32);
    double x0 = msh.idim == 2 ? getprdl2<2>(msh.coord[ipoi1], du)
                              : getprdl2<3>(msh.coord[ipoi1], du);
    double x1 = msh.idim == 2 ? getprdl2<2>(msh.coord[ipoi2], du)
                              : getprdl2<3>(msh.coord[ipoi2], du);

    METRIS_ASSERT(abs(x1 - x0) > 1.0e-15) 

    int ndiv = 2;
    double edg_len0 = 0; 
    // if EG_eval proves costly, we can cache the results from the previous iteration. 
    // (as we're doubling samples each iteration) -> it's also a simple idiv%2 == 0 check. 
    while(true){

      int iwhich = 0; // which one is the previous in coop
      for(int ii = 0; ii < msh.idim; ii++) 
        msh.coord(ipon[iwhich],ii) = msh.coord(ipoi1,ii);
      for(int ii = 0; ii < nnmet   ; ii++) 
        msh.met(  ipon[iwhich],ii) = msh.met(ipoi1,ii);

      if(iverb >= 3){
        printf("       - init ipon[%d] at coord:",iwhich);
        dblAr1(msh.idim,msh.coord[ipoi1]).print();
        printf("       - met:");
        dblAr1(nnmet,msh.met[ipoi1]).print();
      }

      double tprev = t0_1;
      // Likewise, we can cache back localization results to be used sequentially
      // -> as we refine, a miss means a new seed to insert between prev and next.
      //int ibseed   = msh.poi2bak(ipoi1,msh.idim-1); 
      double edg_len1 = 0;
      for(int idiv = 0; idiv < ndiv; idiv++){
        double tnext = tprev + (t0_2 - t0_1)/ndiv; 

        if(idiv == ndiv - 1){
          for(int ii = 0; ii < msh.idim; ii++) 
            msh.coord(ipon[1 - iwhich],ii) = msh.coord(ipoi2,ii);
          for(int ii = 0; ii < nnmet   ; ii++) 
            msh.met(  ipon[1 - iwhich],ii) = msh.met(ipoi2,ii);
        }else{
          int ierro = EG_evaluate(obj, &tnext, result);
          METRIS_ENFORCE(ierro == EGADS_SUCCESS);


          if(iverb >= 4){
            printf("       - idiv %d / %d t = %f coord = ", idiv, ndiv, tnext);
            dblAr1(msh.idim, result).print();
          }

          for(int ii = 0; ii < msh.idim; ii++) 
            msh.coord(ipon[1 - iwhich],ii) = result[ii];

          // Project point on edge 
          //double dv[3];
          //for(int ii = 0; ii < msh.idim; ii++) dv[ii] = msh.coord(ipoi2,ii)
          //                                            - msh.coord(ipon[1 - iwhich],ii);
          double dtprd = 0;                                                    
          if(msh.idim == 2){
            dtprd = getprdl2<2>(msh.coord[ipon[1 - iwhich]],du);
          }else{
            dtprd = getprdl2<3>(msh.coord[ipon[1 - iwhich]],du);
          }

          dtprd = (dtprd - x0) / (x1 - x0);



          for(int ii = 0; ii < nnmet; ii++) 
             msh.met[ipon[1 - iwhich]][ii]
                      =        dtprd  * msh.met(ipoi1,ii) 
                      + (1.0 - dtprd) * msh.met(ipoi2,ii);

          if(iverb >= 3){
            printf("       - idiv %d / %d dtprd %f param = %f coord: ",
                   idiv,ndiv,dtprd,tnext);
            dblAr1(msh.idim,msh.coord[ipon[1 - iwhich]]).print();
            printf("       - met:");
            dblAr1(nnmet,msh.met[ipon[1 - iwhich]]).print();
          }


          #ifndef NDEBUG
          if(dtprd < -10){
            msh.template newbpotopo<0>(ipon[1-iwhich],-1);
            msh.template newbpotopo<0>(ipon[  iwhich],-1);

            printf("## FAILURE idiv = %d / %d dtprd = %f iedge = %d\n",
              idiv, ndiv, dtprd, iedge);
            printf("x0 = %15.7e x1 = %15.7e \n",x0,x1);

            printf("       - ipoi1 %d at coord:",ipoi1);
            dblAr1(msh.idim,msh.coord[ipoi1]).print();
            printf("       - ipoi2 %d at coord:",ipoi2);
            dblAr1(msh.idim,msh.coord[ipoi2]).print();
            printf("       - ipoi1 met:");
            dblAr1(nnmet,msh.met[ipoi1]).print();
            printf("       - ipoi2 met:");
            dblAr1(nnmet,msh.met[ipoi2]).print();

            printf("       - ipon %d at coord:",ipon[1 - iwhich]);
            dblAr1(msh.idim,msh.coord[ipon[1 - iwhich]]).print();
            printf(" -> comp at t = %15.7e iref = %d \n",tnext,iref);


            printf(" ib1 bpo2ibi: ");
            intAr1(nibi, msh.bpo2ibi[ibpo1]).print();
            printf(" ib2 bpo2ibi: ");
            intAr1(nibi, msh.bpo2ibi[ibpo2]).print();

            printf(" parametric coords ipoi1 = %15.7e ipoi2 = %15.7e\n",
                    t0_1,t0_2);

            printf(" full bpo list for ipoi1:\n");
            int ib1 = msh.poi2bpo[ipoi1];
            for(;ib1 >= 0; ib1 = msh.bpo2ibi(ib1,3)){
              printf(" %d : ",ib1);
              intAr1(nibi, msh.bpo2ibi[ib1]).print();
              printf("   r = %15.7e \n", msh.bpo2rbi(ib1,0));
            }
            printf(" full bpo list for ipoi2:\n");
            int ib2 = msh.poi2bpo[ipoi2];
            for(;ib2 >= 0; ib2 = msh.bpo2ibi(ib2,3)){
              printf(" %d : ",ib2);
              intAr1(nibi, msh.bpo2ibi[ib2]).print();
              printf("   r = %15.7e \n", msh.bpo2rbi(ib2,0));
            }


            ierro = EG_evaluate(obj, &t0_1, result);
            printf("       - ipoi1 reeval at coord:");
            dblAr1(msh.idim,result).print();

            ierro = EG_evaluate(obj, &t0_2, result);
            printf("       - ipoi2 reeval at coord:");
            dblAr1(msh.idim,result).print();

            writeMesh("debug_CADCurveLengths.meshb", msh);
            METRIS_THROW(GeomExcept());
          }
          #endif

          //double lmet[6];
          //for(int ii = 0; ii < nnmet; ii++) 
          //   lmet[ii] =        dtprd  * msh.met(ipoi1,ii) 
          //            + (1.0 - dtprd) * msh.met(ipoi2,ii);
          //if(msh.idim == 2){
          //  getexpmet_cpy<2>(lmet, msh.met[ipon[1 - iwhich]]);
          //}else{
          //}

          // No need! The front is the back still (ran first thing adaptMesh)
          // Localize new point in back mesh
          //int isdnew; 
          //ierro = msh.interpMetBack(msh.coord[ipon[1 - iwhich]],ibseed,&isdnew,
          //                          msh.met[ipon[1 - iwhich]]);
          //METRIS_ENFORCE(ierro == 0);
          //ibseed = isdnew; 
        }


        double len;
        if(msh.idim == 2){
          len = getlenedg_geosz<MFT,2,1>(msh,edg2pol,sz);
        }else{
          len = getlenedg_geosz<MFT,3,1>(msh,edg2pol,sz);
        }

        if(iverb >= 3) printf("       - len = %15.7e\n",len);

        #ifndef NDEBUG
        if(std::isnan(len)){
          printf("NaN length ! \n");
          writeMesh("debug_NaNlen",msh);
          msh.param->iverb = 5;
          if(msh.idim == 2){
            len = getlenedg_geosz<MFT,2,1>(msh,edg2pol,sz);
          }else{
            len = getlenedg_geosz<MFT,3,1>(msh,edg2pol,sz);
          }
          METRIS_THROW(GeomExcept());
        }
        #endif

        edg_len1 += len;


        iwhich = 1 - iwhich;
        tprev  = tnext; 
      }
      if(iverb >= 3) printf("     - iref = %d ndiv = %d length = %f\n",
                                                            iref,ndiv,edg_len1);
      // Worst case, two successive are both over/undershooting. Then the actual 
      // error is twice that between the two. 
      // Normalize by number of edges in reference, ensuring the total curve 
      // length tolerance is met. 
      // Also normalize by initial curve length estimate : DO NOT normalize by 
      // this edge's length ! otherwise this'll blow up as mesh is refined. 
      if(abs(edg_len1 - edg_len0) * ref2ned[iref] < tol * crv_len0[iref] / 2){
        if(iverb >= 3) printf("     -> converged\n");
        break;
      }
      ndiv *= 2;
      METRIS_ENFORCE_MSG(ndiv < 1e6, "Infinite loop? Increase geotol. Last len ="
        << edg_len0 << " new = "<<edg_len1<< " ref2ned = "<<ref2ned[iref]
        << " tol = "<<tol<<" crv_len0 = "<<crv_len0[iref]);
      edg_len0 = edg_len1;
      METRIS_ENFORCE_MSG(ndiv > 0, "Integer overflow");
    } // while(true)

    crv_len[iref] += edg_len0;

  } // for int iedge 


  if(iverb >= 1){
    for(int iref = 0; iref < nref; iref++)
      printf("   - END line %d/%d len = %f\n",iref, nref, crv_len[iref]);
  }


  msh.met.setSpace(ispac0);
}


template void getCADCurveLengths<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, double tol, dblAr1 &crv_len);
template void getCADCurveLengths<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, double tol, dblAr1 &crv_len);





/* --------------------------------------------
AUXILIARY FUNCTIONS 
   -------------------------------------------- */

template<class MFT>
static int aux_walk_line(Mesh<MFT> &msh, MshCavity& cav, ego obj, 
                         int ipoistart, int iedgstart, double adjusted_tarlen, 
                         int *edg2pol, int *ibpnw, int *iedgend, int *ifacend, 
                         double *lenend, double *szend, bool *ifin, 
                         double *erredg, int *npterr, int *nEGrro, 
                         int ithrd1){

  const int iverb = msh.param->iverb;
  if(iverb >= 3) 
    printf("    - aux_walk_line start ipoistart %d iedgstart %d = %d, %d\n",
    ipoistart,iedgstart,msh.edg2poi(iedgstart,0),msh.edg2poi(iedgstart,1));

  int iface, ibpoi, iedge = iedgstart;
  edg2pol[0] = ipoistart;
  double coor0[3], result[18];
  int ipoiprev = ipoistart; 
  while(true){

    // Current edge will be in cavity no matter what
    cav.lcedg.stack(iedge); 



    // Add adjacent face. This will seed the cavity extension later on
    iface = msh.edg2fac[iedge];
    METRIS_ASSERT(iface >= 0 || msh.nface == 0);
    if(iface >= 0 && msh.fac2tag(ithrd1,iface) < msh.tag[ithrd1]){
      msh.fac2tag(ithrd1,iface) = msh.tag[ithrd1];
      int iedl = getedgfac(msh,iface,msh.edg2poi(iedge,0),msh.edg2poi(iedge,1));
      METRIS_ASSERT(iedl >= 0);
      cav.lcfac.stack(iface); 
    }



    // Get next point for length computation 
    int iipo2 = -1;
    for(int ii = 0; ii < 2; ii++){
      if(msh.edg2poi(iedge,ii) == ipoiprev) iipo2 = ii;
    }
    METRIS_ASSERT(iipo2 != -1);
    edg2pol[1] = msh.edg2poi[iedge][1-iipo2];

    if(iverb >= 3) 
    printf("      - iedge %d = %d, %d iface = %d nxtp = %d \n",
           iedge,msh.edg2poi(iedge,0),msh.edg2poi(iedge,1), iface, iipo2);

    // Evaluate this point on the CAD edge using its param coord 
    ibpoi = msh.poi2bpo[edg2pol[1]];
    *ifin = false;
    // If it's a corner, this is the end of the line. 
    if(msh.bpo2ibi[ibpoi][1] < 1){
      *ifin = true;
      // Still need the length to update tarlen next outer iteration. 
      ibpoi = getent2bpo(msh,ibpoi,iedge,1);
    }
    METRIS_ASSERT(ibpoi >= 0);
    METRIS_ASSERT(msh.bpo2ibi[ibpoi][1] == 1);


    if(!(*ifin)){ // Corner does not need evaluating, and end is iff corner
      int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
      METRIS_ASSERT(ierro == 0);
      if(ierro != 0){ // This failed 
        if(iverb >= 2) printf("   - EG_evaluate error %d \n",ierro);
        (*nEGrro) ++;
        return 1;
      }

      // Accumulate error (distance between eval'd and original). Also
      // virtually update coordinates (revert to prior copy) for
      // length computation. 
      for(int ii = 0; ii < msh.idim; ii++){
        coor0[ii] = msh.coord[edg2pol[1]][ii];
        *erredg  += abs(msh.coord[edg2pol[1]][ii] - result[ii]);
        msh.coord[edg2pol[1]][ii] = result[ii];
      }
      (*npterr)++;
    }

    if(msh.idim == 2){
      *lenend = getlenedg_geosz<MFT,2,1>(msh,edg2pol,szend);
    }else{
      *lenend = getlenedg_geosz<MFT,3,1>(msh,edg2pol,szend);
    }
    if(iverb >= 2) printf("      - len next %f\n",*lenend);

    //// Accumulate curve length for tarlen update 
    //crv_len += len;

    // Only if the point was previously moved should it be reverted
    if(!(*ifin)){
      for(int ii = 0; ii < msh.idim; ii++){
        msh.coord[edg2pol[1]][ii] = coor0[ii];
      }
    }

    // Stop walking if nowhere to go, or if target length met.
    if(*lenend > adjusted_tarlen || (*ifin)) break;


    ipoiprev = msh.edg2poi[iedge][1-iipo2]; // New previous next point.
    if(iverb >= 2) printf("      - iedge %d -> %d len %f\n",iedge,msh.edg2edg(iedge,iipo2),*lenend);
    iedge = msh.edg2edg(iedge,iipo2);    // New next edge 

  } // End while(true)


  METRIS_ASSERT(iface >= 0);
  *ifacend = iface;

  METRIS_ASSERT(iedge >= 0);
  *iedgend = iedge; 

  METRIS_ASSERT(ibpoi >= 0);
  *ibpnw = ibpoi; 

  return 0;
}



template<class MFT>
static int gen_newp_line(Mesh<MFT> &msh, MshCavity& cav, ego obj,
                         int ibpo0, int ibpnw,
                         int iedgseed, int ifacseed, const int *edg2pol, 
                         double *sz, int miter_bisection,
                         double tarlen, double lentolabs, 
                         double *adjusted_tarlen,  // in/out
                         double *lenend, int *nEGrro){

  const double tola = 1.0e-14;

  
  const int iverb = msh.param->iverb;
  const int nnmet = (msh.idim*(msh.idim+1))/2;

  cav.ipins = msh.newpoitopo(1,iedgseed);
  int ibins = msh.template newbpotopo<1>(cav.ipins,iedgseed);


  // Get optimal location. First guess, t. Then bisection 
  double a = sz[0] / sz[1], t;
  if(abs(a-1.0) < tola){ // Constant size 
    t = 1.0 - *adjusted_tarlen / sz[0];
  }else{
    double loga = log(a);
    t = log(1.0 + loga * *adjusted_tarlen / (2.0 * sz[1])) / loga;
  }

  int ierro;
  double result[18];
  int edg2po2[2];
  edg2po2[0] = edg2pol[0];
  edg2po2[1] = cav.ipins;
  double t0 = 0, t1 = 1;
  bool fndt = false; 

  double dlen_best = 1.0e30, adjl_best, t_best;
  double coop_best[3], met_best[6];

  double len;
  if(iverb >= 2){
    printf("   - start bisect from pt %d extrem u %f %f tarlen = %f t ini %f\n",
            edg2pol[0], msh.bpo2rbi(ibpo0,0), msh.bpo2rbi(ibpnw,0), 
            *adjusted_tarlen, t);
    printf("   - ibpo0 = %d t = %f : ",ibpo0, msh.bpo2rbi(ibpo0,0));
    intAr1(nibi, msh.bpo2ibi[ibpo0]).print();
    printf("   - From full list:\n");
    for(int ii = msh.poi2bpo[msh.bpo2ibi(ibpo0,0)]; ii >= 0; ii = msh.bpo2ibi(ii,3)){
      printf("%d = %f : ",ii,msh.bpo2rbi[ii][0]);
      intAr1(nibi,msh.bpo2ibi[ii]).print();
    }
    printf("   - ibpon = %d t = %f: ",ibpnw, msh.bpo2rbi(ibpnw,0));
    intAr1(nibi, msh.bpo2ibi[ibpnw]).print();
    printf("   - From full list:\n");
    for(int ii = msh.poi2bpo[msh.bpo2ibi(ibpnw,0)]; ii >= 0; ii = msh.bpo2ibi(ii,3)){
      printf("%d = %f : ",ii,msh.bpo2rbi[ii][0]);
      intAr1(nibi,msh.bpo2ibi[ii]).print();
    }
    #ifndef NDEBUG
    if(iverb >= 3) writeMesh("bisec0", msh);
    #endif
    printf("   - interpMetBack starts from iedgseed = %d (%d,%d)\n",iedgseed,
      msh.edg2poi(iedgseed,0),msh.edg2poi(iedgseed,1));
  }
  for(int itfnd = 0; itfnd < miter_bisection; itfnd++){
    msh.bpo2rbi(ibins,0) =        t * msh.bpo2rbi(ibpo0,0) 
                         + (1.0 - t)* msh.bpo2rbi(ibpnw,0);

    ierro = EG_evaluate(obj, msh.bpo2rbi[ibins], result);
    METRIS_ASSERT(ierro == 0);
    if(ierro != 0){ // This failed 
      if(iverb >= 2) printf("   - EG_evaluate error %d \n",ierro);
      nEGrro ++;
      return 1;
    }
    for(int ii = 0; ii < msh.idim; ii++) msh.coord[cav.ipins][ii] = result[ii];

    // -- Interpolate metric at insertion point, only if new
    //ierro = msh.interpMetBack(msh.coord[cav.ipins],ifacseed,
    ierro = msh.interpMetBack(cav.ipins,1,iedgseed,
                              msh.edg2ref[iedgseed],&result[3]);
    #ifndef NDEBUG
      if(iverb >= 3){
        printf("   - computed metric ");
        dblAr1(nnmet,msh.met[cav.ipins]).print();
      }
      if(ierro != 0){
        msh.template newbpotopo<0>(cav.ipins);
        writeMesh("debug_interpMetBack"+std::to_string(cav.ipins),msh);
      }
      METRIS_ASSERT_MSG(ierro == 0, "interpMetBack failed");
    #endif

    if(ierro != 0){
      if(iverb >= 3){
        printf("Wait error interpMetBack %d \n",ierro);
        #ifndef NDEBUG
        wait();
        #endif
      }
      return 1;
    }

    if(msh.idim == 2){
      len = getlenedg_geosz<MFT,2,1>(msh,edg2po2,sz);
    }else{
      len = getlenedg_geosz<MFT,3,1>(msh,edg2po2,sz);
    }

    if(len > *adjusted_tarlen) t0 = t;
    else                       t1 = t;
    

    if(iverb >= 2) printf("     - tried t = %f param %f got len = %f coop = %f %f new t0 %f t1 %f \n",
                          t,msh.bpo2rbi(ibins,0),len,result[0],result[1],t0,t1);
    t = (t0 + t1) / 2;

    if(std::isnan(len)){
      printf("## NaN len !!\n");
      printf("edg2po2 = ");
      intAr1(2,edg2po2).print();
      printf("metrics at points:\n");
      for(int ii = 0; ii < 2; ii++){
        int ipdbg = edg2po2[ii];
        if(ipdbg < 0) continue;
        printf("%d : ",ipdbg);
        dblAr1(nnmet,msh.met[ipdbg]).print();
      }
      exit(1);
      wait();
    }

    if(abs(len - *adjusted_tarlen) < dlen_best){
      dlen_best = abs(len - *adjusted_tarlen);
      adjl_best = tarlen + (*adjusted_tarlen - len); 
      for(int ii = 0; ii < msh.idim; ii++) coop_best[ii] = msh.coord(cav.ipins,ii);
      t_best = msh.bpo2rbi(ibins,0);
      for(int ii = 0; ii < nnmet; ii++) met_best[ii] = msh.met(cav.ipins,ii);
    }

    if(len < *adjusted_tarlen + lentolabs 
    && len > *adjusted_tarlen - lentolabs){
      if(iverb >= 2) printf("   - end bisect: len = target %f +- %f, new "
        "tarlen = %f error %15.7e (adj) nedge %d \n", *adjusted_tarlen, lentolabs, 
        tarlen + (*adjusted_tarlen - len), *adjusted_tarlen - len, msh.nedge); 
      *adjusted_tarlen = tarlen + (*adjusted_tarlen - len); 
      fndt = true;
      break;
    }
  } // bisection 

  *lenend = len;

  if(msh.param->dbgfull){
    METRIS_ASSERT(fndt); 
  }else if(!fndt){
    if(iverb >= 2) printf("   - end bisect !fndt: len = target %f +- %f, new "
                 "tarlen = %f error %15.7e (adj) nedge %d \n", *adjusted_tarlen, 
                 lentolabs, adjl_best, dlen_best, msh.nedge); 
    *adjusted_tarlen = adjl_best;
    for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins,ii) = coop_best[ii];
    msh.bpo2rbi(ibins,0) = t_best;
    for(int ii = 0; ii < nnmet; ii++) msh.met(cav.ipins,ii) = met_best[ii];
  }
  return 0;
}

} // end namespace
