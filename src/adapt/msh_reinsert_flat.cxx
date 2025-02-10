//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "msh_reinsert_flat.hxx"
#include "../low_geo.hxx"
#include "../mprintf.hxx"
#include "../msh_structs.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../io_libmeshb.hxx"
#include "../adapt/low_collapse.hxx"
//#include "../Boundary/low_projsurf.hxx"


namespace Metris{

// Reinsert vertices that create almost flat elements. 
template<class MFT, int idim, int ideg>
int reinsertFlat(Mesh<MFT> &msh, bool allow_collapse, int ithread){
  GETVDEPTH(msh);
  constexpr int gdim = idim;
  constexpr int tdim = idim; 
  constexpr int nnmet = (gdim * (gdim + 1)) / 2;

  if(tdim != 2) METRIS_THROW_MSG(TODOExcept(), "Implement reinsertFlat dim 3");

  // For now, make it an option later 
  const double hgttol = 1.0e-8;
  const int miter = 10;
  //const double prjtol = 1.0e-10;
  //// Interior ones are swaps: avoid 
  // Note: they are not always
  //bool ibdryonly = true;

  //// Boundary insertion of interior point requires collapse, unless we can 
  //// reversibly newbpo on an existing point. This sounds hazardous. 
  //// The idea is particularly badly inspired as a change of implementation 
  //// of newpotopo could break it!
  //// And adding the entry manually is similarly not satisfactory. 
  //// Better to update info in smooth than to do this. 
  //if(ibdryonly && !allow_collapse) return 0;

  intAr2 &ent2poi = msh.ent2poi(tdim);
  //intAr2 &ent2ent = msh.ent2ent(tdim);
  //intAr2 &ent2tag = msh.ent2tag(tdim);

  //writeMesh("debug00", msh);

  MetSpace ispac0 = msh.met.getSpace();
  msh.met.setSpace(MetSpace::Log);

  // Can have up to one edge in the cavity
  int mcfac = 100, mcedg = 1; 
  int nwork = mcfac + mcedg;
  RoutineWorkMemory ipool(msh.iwrkmem);
  int *iwork = ipool.allocate(nwork);
  MshCavity cav(0,mcfac,mcedg,nwork,iwork);
  CavOprOpt  opts;
  CavWrkArrs work;
  opts.allow_topological_correction = false;
  opts.skip_topo_checks = false;
  opts.allow_remove_points = true;
  //opts.dryrun   = true;
  opts.dryrun   = false;

  const int merror = CAV_ERR_NERROR;
  intAr1 lerror(merror);
  lerror.set_n(merror);
  lerror.fill(0); 
  int nerro_tot = 0;
  int noper_tot = 0;

  double metl[nnmet];
  double height[tdim + 1];
  double bar0[tdim + 1];
  for(int ii = 0; ii < tdim + 1; ii++) bar0[ii] = 1.0 / (tdim + 1.0);

  int nent0 = 0;
  int nent1 = msh.nentt(tdim);
  for(int niter = 0; niter < miter; niter++){
    CPRINTF1(" - niter %d / %d nent0 -> nent1 %d -> %d\n",niter,miter,nent0,nent1);
    int noper  = 0;
    int nerro = 0;
    for(int ientt = nent0; ientt < nent1; ientt++){
      INCVDEPTH(msh);
      if(isdeadent(ientt,ent2poi)) continue;

      //if(ibdryonly){
      //  bool onbdry = false;
      //  for(int ifa = 0; ifa < tdim + 1; ifa++){
      //    if(ent2ent(ientt,ifa) >= 0) continue;
      //    onbdry = true;
      //  }

      //  if(!onbdry) continue;
      //}

      msh.met.getMetBary(AsDeg::P1,DifVar::None, 
                         MetSpace::Exp, 
                         ent2poi[ientt], tdim, bar0, metl, NULL);

      getheightentP1_aniso<gdim>(ent2poi[ientt], msh.coord, metl, height);

      if(DOPRINTS1()){
        MPRINTF(" - ientt %d heights ",ientt);
        dblAr1(tdim + 1, height).print();
      }


      for(int ifa = 0; ifa < tdim + 1; ifa++){
        int ipoin = ent2poi(ientt,ifa);


        if(height[ifa] >= hgttol) continue;
        // Low height element: insert itself and neighbour in cavity, as well 
        // as all points in ball of ipoin 

        CPRINTF1(" - Collapse point %d \n", ipoin);
        // Just collapse the point. 
        collversurf(msh,ientt,ifa,msh.param->iverb,lerror);
        break;

        #if 0
        int ifnei = ent2ent(ientt,ifa);
        if(ifnei >= 0 && ibdryonly) continue;


        if(iverb >= 1) printf("   - Reinsert pt %d face %d (%d,%d,%d)" 
                              "ied %d height = %15.7e\n",ipoin, ientt,
                               msh.fac2poi(ientt,0), msh.fac2poi(ientt,1),
                               msh.fac2poi(ientt,2), ifa, height[ifa]);
 
        cav.reset();
        if(ifnei >= 0){
          cav.ipins = ipoin;  
        }else{
          int ibins; 
          int iedge, iface;
          if(tdim == 2){
            iedge = msh.facedg2glo(ientt, ifa);
            cav.lcedg.stack(iedge);
            cav.ipins = msh.newpoitopo(1, iedge);
            ibins = msh.newbpotopo(cav.ipins,1,iedge);
          }else if(tdim == 3){
            iface = msh.tetfac2glo(ientt, ifa);
            cav.ipins = msh.newpoitopo(2, iface);
            ibins = msh.newbpotopo(cav.ipins,2,iface);
            METRIS_THROW_MSG(TODOExcept(), "Implement projpoifac (get (u,v))");
          }
          for(int ii = 0; ii < gdim; ii++) 
            msh.coord(cav.ipins,ii) = msh.coord(ipoin,ii);
          for(int ii = 0; ii < nnmet; ii++) 
            msh.met(cav.ipins,ii) = msh.met(ipoin,ii);
          if(tdim == 2){
            double bary[2], coopr[3];
            ierro = projptedgCAD<idim>(msh, msh.coord[cav.ipins], prjtol, 
                                       iedge, msh.bpo2rbi[ibins], bary, coopr);
            if(ierro > 0){
              if(iverb >= 2) 
                printf(" ## msh_reinsert_flat projptedg ierro %d \n",ierro);
              nerro++;
              continue;
            }
            double dist = geterrl2<gdim>(coopr, msh.coord[cav.ipins]);
            if(dist > 1.0e-12){

              //#ifndef NDEBUG 
              printf("Large dist ? %15.7e \n",dist);
              for(int ii = 0; ii < msh.idim; ii++){
                printf(" %15.7e ",coopr[ii]);
              }
              printf("\n");
              for(int ii = 0; ii < msh.idim; ii++){
                printf(" %15.7e ",msh.coord[cav.ipins][ii]);
              }
              printf("\n");

              ipoin = msh.newpoitopo(-1,-1);
              int ibpoi = msh.newbpotopo(ipoin,0,ipoin);
              for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipoin,ii) = coopr[ii];


              ipoin = msh.newpoitopo(-1,-1);
              ibpoi = msh.newbpotopo(ipoin,0,ipoin);
              for(int ii = 0; ii < msh.idim; ii++) 
                msh.coord(ipoin,ii) = msh.coord(cav.ipins,ii);
              writeMesh("dbg_projptedg.meshb",msh);

              METRIS_THROW_MSG(GeomExcept(),"LArge dist "<<dist);
              //#else
              //if(iverb >= 3) printf("  - Large dist %f continue \n",dist);
              //continue;
              //#endif

            }

          }else if(tdim == 3){
            METRIS_THROW_MSG(TODOExcept(), "Implement use projpoifac");
          }
        }

        if constexpr(idim == 2){
          bool imani; 
          int iopen;
          ball2(msh, ipoin, ientt, cav.lcfac, cav.lcedg, &iopen, &imani, ithread);
          if(ifnei >= 0){
            cav.lcfac.stack(ifnei);
          }else{
            int ifedg = msh.facedg2glo(ientt, ifa);
            cav.lcedg.stack(ifedg);
          }
        }else{
          METRIS_THROW(TODOExcept());
        }

        ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithread);
        if(ierro > 0){
          nerro++;
          lerror[ierro-1]++;
        }else{
          noper++;
          // Triangle ceased to exist
          break; 
        }

        #endif

      } // for ifa 
    } // for ientt 

    noper_tot += noper;
    nerro_tot += nerro;

    if(noper == 0) break;

    nent0 = nent1;
    nent1 = msh.nentt(tdim);
    CPRINTF1(" - iter %d noper = %d nerro = %d \n", niter, noper, nerro);

  }// for(int niter)

  CPRINTF1(" - End noper = %d nerro = %d \n", noper_tot, nerro_tot);
  if(DOPRINTS2() && nerro_tot > 0){
    CPRINTF1(" - ierro list:\n");
    for(int ii = 0; ii < merror; ii++){
      if(lerror[ii] == 0) continue;
      CPRINTF1("   - ierro = %d : %d \n",ii+1,lerror[ii]);
    }
  }


  msh.met.setSpace(ispac0);

  return noper_tot;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template int reinsertFlat<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                                       bool,int ithread);\
template int reinsertFlat<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                                       bool,int ithread);\
template int reinsertFlat<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                                       bool,int ithread);\
template int reinsertFlat<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                                       bool,int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


}//namespace Metris
