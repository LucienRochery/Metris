//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../adapt/low_collapse.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../low_geo.hxx"
#include "../mprintf.hxx"
#include "../msh_structs.hxx"
#include "../io_libmeshb.hxx"
#include "../msh_checktopo.hxx"
#include "../adapt/low_increasecav.hxx"

namespace Metris{

/*
Collapse a vertex in short edge 
Cavity is passed in to reuse allocations
*/
template<class MFT>
int colledgsurf(Mesh<MFT>& msh, int iface, int iedl, double qmax_suf, 
                MshCavity &cav, CavWrkArrs &work, 
                intAr1 &lerro, int ithrd1, int ithrd2){
  if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(), "Implement + tet nelem = "<<msh.nelem)

  GETVDEPTH(msh);


  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  //METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  //METRIS_ASSERT(ithrd1 != ithrd3);
  //METRIS_ASSERT(ithrd2 != ithrd3);

 
  CavOprOpt  opts;
  CavOprInfo info;
  opts.allow_topological_correction = true; // To fetch missing edges
  opts.skip_topo_checks = false;
  opts.allow_remove_points = true;
  //opts.dryrun   = true;
  opts.dryrun   = false;
  opts.qmax_suf = qmax_suf;
  cav.lcedg.set_n(0);
  cav.lcfac.set_n(0);
  cav.lctet.set_n(0);

  int ierro = 0;

  int ip1 = msh.fac2poi(iface,lnoed2[iedl][0]);
  int ip2 = msh.fac2poi(iface,lnoed2[iedl][1]);


  CPRINTF1("-- START colledgsurf iface = %d iedl = %d = (%d,%d) \n",iface,iedl,ip1,ip2);
  if(DOPRINTS2()) writeMesh("debug_collapse0.meshb",msh);

  // Track best insertion dry run to run it at the end
  double qmabest = 1.0e30;
  int ivebest = -1; // p collapse vertex of 
  int ipibest = -1; // p insert 

  double nrmal[3];

  int tdim[2], tdimc;
  for(int ii = 0; ii < 2 ; ii ++){
    int ipoin = msh.fac2poi(iface,lnoed2[iedl][ii]);
    tdim[ii] = msh.getpoitdim(ipoin);
    if(DOPRINTS2() && tdim[ii] == 0){
      CPRINTF2("Topo dim 0 point %d \n",ipoin);
      int ib = msh.poi2bpo[ipoin];
      CPRINTF2("poi2bpo = %d bpo2ibi = ",ib);
      intAr1(nibi,msh.bpo2ibi[ib]).print();
    }
  }
  CPRINTF1(" - topo dims %d %d \n",tdim[0],tdim[1]);

  tdimc = MAX(tdim[0], tdim[1]);

  bool imani;
  int  iopen;
  for(int iver = 0; iver < 2; iver++){
    // Collapse this one 
    int ipcol = msh.fac2poi(iface,lnoed2[iedl][iver]);

    // Collapse only among the highest-dimensional points
    if(tdim[iver] < tdimc) continue;

    int ibcol = msh.poi2bpo[ipcol];
    if(ibcol >= 0){
      if(msh.bpo2ibi(ibcol,1) == 0) continue; // Skip corners
    }

    ierro = ball2(msh,ipcol,iface,
                  cav.lcfac,cav.lcedg,
                  &iopen,&imani,ithrd1);
    if(!imani) METRIS_THROW_MSG(TODOExcept(), "Non manifold not handled");
    METRIS_ASSERT(ierro == 0);

    if(msh.idim == 3){
      // We'll not use the ipins normal but rather this one as we already have the ball. 
      // Just treat it as P1
      getnorballref<1>(msh,cav.lcfac,-1,nrmal);
      cav.nrmal = nrmal;
    }


    CPRINTF1(" - try collapse poi = %d ball nface = %d nedge = %d \n",
                                ipcol,cav.lcfac.get_n(),cav.lcedg.get_n());

    // Try the cavity call with different ipins in neighbours of ipcol
    msh.tag[ithrd1]++;
    // We'll stack on top of the ball, so we need to be able to prune to nbalf
    // with each attempt, as well as restrict search of ipins to ball 
    int nbalf = cav.lcfac.get_n();
    int nbale = cav.lcedg.get_n();
    for(int iicfc = 0; iicfc < nbalf; iicfc++){
      int icfac = cav.lcfac[iicfc]; 
      METRIS_ASSERT(icfac >= 0 && icfac < msh.nface);
      METRIS_ASSERT(!isdeadent(icfac,msh.fac2poi));

      // Doesn't change but easy to get it here 

      for(int ive2 = 0; ive2 < 3; ive2 ++){
        int ipins = msh.fac2poi(icfac,ive2);
        if(ipins == ipcol) continue;
        if(msh.poi2tag(ithrd1,ipins) >= msh.tag[ithrd1]) continue;
        msh.poi2tag(ithrd1,ipins) = msh.tag[ithrd1];

        //if(msh.isboundary_edges() && cav.lcedg.get_n() > 0 
        //|| msh.isboundary_faces() && cav.lcfac.get_n() > 0){
        //  int ibpoi = msh.poi2bpo[ipins];
        //  if(ibpoi < 0) continue;
        //}

        cav.ipins = ipins;
        cav.lctet.set_n(0);
        cav.lcfac.set_n(nbalf); // Revert to simple ball 
        cav.lcedg.set_n(nbale); // Revert 
        CPRINTF1(" - try reinsert point %d tag = %d vs %d \n",
                             ipins,msh.poi2tag(ithrd1,ipins),msh.tag[ithrd1]);


        if(DOPRINTS2()) writeMeshCavity("collapse_cavity0.meshb", msh, cav, ithrd2);

        // Increase cavity with Delaunay criterion
        ierro = increase_cavity_Delaunay(msh, cav, ithrd2);
        if(ierro != 0) continue;

        //int nprem = increase_cavity_lenedg(msh,cav,ipins,ithrd2,ithrd3);
        ierro = increase_cavity2D(msh,cav,ithrd2);
        if(ierro > 0) continue;

        if(DOPRINTS2()) writeMeshCavity("collapse_cavity1.meshb", msh, cav, ithrd2);


        CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
          ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
        }}CT_FOR1(ideg);

        if(ierro > 0){
          lerro[ierro-1]++;
        }
        //printf("Debug wait;\n");
        // If operation was done, out
        if(info.done){
          CPRINTF1("-- END colledgsurf successful using ipcol = %d ipins = %d \n",ipcol,ipins);
          if(DOPRINTS2()) writeMesh("debug_collapse1.meshb",msh);
          msh.poi2ent(ipcol,0) = -1;
          msh.poi2ent(ipcol,1) = -1;
        }
        if(msh.param->dbgfull){
          check_topo(msh, msh.nbpoi, msh.npoin, msh.nedge, msh.nface, msh.nelem, ithrd2); 
        }
        if(info.done) return 0;

        CPRINTF1(" - return qmax = %f \n",info.qmax_end);
        if(info.qmax_end < qmabest && ierro == 0){
          CPRINTF1(" - new best quality !\n");
          qmabest = info.qmax_end;
          ivebest = iver;
          ipibest = ipins;
        }
      }
    }
  }

  if(ivebest == -1){
    // Failed to find a reinsertion candidate. 
    // Let's not rely on exceptions as failure here is not invalid behaviour 
    // and could happen quite frequently 
    return 1;
  }

  cav.ipins = ipibest;
  opts.dryrun = false;

  if(ivebest == 0){ // Otherwise ball hasn't changed !
    int ipcol = msh.fac2poi(iface,lnoed2[iedl][ivebest]);
    ierro = ball2(msh,ipcol,iface,
                  cav.lcfac,
                  cav.lcedg,
                  &iopen,&imani,ithrd1);
    // Increase cavity with Delaunay criterion
    ierro = increase_cavity_Delaunay(msh, cav, ithrd2);
    METRIS_ASSERT(ierro == 0);

    //int nprem = increase_cavity_lenedg(msh,cav,ipins,ithrd2,ithrd3);
    ierro = increase_cavity2D(msh,cav,ithrd2);
    METRIS_ASSERT(ierro == 0);
  }

  if(DOPRINTS2()) writeMeshCavity("collapse_cavity0.meshb", msh, cav, ithrd2);


  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
  }}CT_FOR1(ideg);

  METRIS_ASSERT(ierro == 0);
  METRIS_ASSERT(abs(info.qmax_end - qmabest) < 1.0e-6);
  METRIS_ASSERT(info.done);

  return 0;
}

template int colledgsurf<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                                          int iface, int iedl, double qmax_suf, 
                                          MshCavity &cav, CavWrkArrs &work, 
                                          intAr1 &lerro, int ithrd1, int ithrd2);
template int colledgsurf<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                                          int iface, int iedl, double qmax_suf, 
                                          MshCavity &cav, CavWrkArrs &work, 
                                          intAr1 &lerro, int ithrd1, int ithrd2);


template<class MFT>
int collversurf(Mesh<MFT>& msh, int iface, int iver, double qmax_suf, 
                MshCavity &cav, CavWrkArrs &work, 
                intAr1 &lerro, int ithrd1, int ithrd2){
  if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(), "Implement + tet nelem = "<<msh.nelem)
  GETVDEPTH(msh);

  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  //METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  //METRIS_ASSERT(ithrd1 != ithrd3);
  //METRIS_ASSERT(ithrd2 != ithrd3);

  CavOprOpt  opts;
  CavOprInfo info;
  opts.allow_topological_correction = true; // To fetch missing edges
  opts.skip_topo_checks = false;
  opts.allow_remove_points = true;
  //opts.dryrun   = true;
  opts.dryrun   = false;
  opts.qmax_suf = qmax_suf;

  int ierro = 0;

  // Track best insertion dry run to run it at the end
  double qmabest = 1.0e30;
  int ivebest = -1; // p collapse vertex of 
  //int ipibest = -1; // p insert 

  double nrmal[3];

  bool imani;
  int  iopen;
  // Collapse this one 
  int ipcol = msh.fac2poi(iface,iver);

  CPRINTF1("-- START colledgsurf iface = %d iver = %d ipoin = %d \n",iface,iver,ipcol);
  if(DOPRINTS2()) writeMesh("debug_collapse0.meshb",msh);


  int ibcol = msh.poi2bpo[ipcol];
  if(ibcol >= 0){
    if(msh.bpo2ibi(ibcol,1) == 0) return 1;
  }

  ierro = ball2(msh,ipcol,iface,
                cav.lcfac,cav.lcedg,
                &iopen,&imani,ithrd1);
  if(!imani) METRIS_THROW_MSG(TODOExcept(), "Non manifold not handled");
  METRIS_ASSERT(ierro == 0);

  if(msh.idim == 3){
    // We'll not use the ipins normal but rather this one as we already have the ball. 
    // Just treat it as P1
    getnorballref<1>(msh,cav.lcfac,-1,nrmal);
    cav.nrmal = nrmal;
  }


  CPRINTF1(" - try collapse poi = %d ball nface = %d nedge = %d \n",
                              ipcol,cav.lcfac.get_n(),cav.lcedg.get_n());

  // Try the cavity call with different ipins in neighbours of ipcol
  msh.tag[ithrd1]++;
  // We'll stack on top of the ball, so we need to be able to prune to nbalf
  // with each attempt, as well as restrict search of ipins to ball 
  int nbalf = cav.lcfac.get_n();
  int nbale = cav.lcedg.get_n();
  for(int iicfc = 0; iicfc < nbalf; iicfc++){
    int icfac = cav.lcfac[iicfc]; 
    METRIS_ASSERT(icfac >= 0 && icfac < msh.nface);
    METRIS_ASSERT(!isdeadent(icfac,msh.fac2poi));

    // Doesn't change but easy to get it here 

    for(int ive2 = 0; ive2 < 3; ive2 ++){
      int ipins = msh.fac2poi(icfac,ive2);
      if(ipins == ipcol) continue;
      if(msh.poi2tag(ithrd1,ipins) >= msh.tag[ithrd1]) continue;
      msh.poi2tag(ithrd1,ipins) = msh.tag[ithrd1];

      //if(msh.isboundary_edges() && cav.lcedg.get_n() > 0 
      //|| msh.isboundary_faces() && cav.lcfac.get_n() > 0){
      //  int ibpoi = msh.poi2bpo[ipins];
      //  if(ibpoi < 0) continue;
      //}

      cav.ipins = ipins;
      cav.lctet.set_n(0);
      cav.lcfac.set_n(nbalf); // Revert to simple ball 
      cav.lcedg.set_n(nbale); // Revert 
      CPRINTF1(" - try reinsert point %d tag = %d vs %d \n",
                           ipins,msh.poi2tag(ithrd1,ipins),msh.tag[ithrd1]);


      if(DOPRINTS2()) writeMeshCavity("collapse_cavity0.meshb", msh, cav, ithrd2);

      // Increase cavity with Delaunay criterion
      ierro = increase_cavity_Delaunay(msh, cav, ithrd2);
      if(ierro > 0) continue;

      //int nprem = increase_cavity_lenedg(msh,cav,ipins,ithrd2,ithrd3);
      ierro = increase_cavity2D(msh,cav,ithrd2);
      if(ierro > 0) continue;

      if(DOPRINTS2()) writeMeshCavity("collapse_cavity1.meshb", msh, cav, ithrd2);


      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
        ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
      }}CT_FOR1(ideg);
      if(ierro > 0){
        lerro[ierro-1]++;
      }
      //printf("Debug wait;\n");
      // If operation was done, out
      if(info.done){
        CPRINTF1("-- END colledgsurf successful using ipcol = %d ipins = %d \n",ipcol,ipins);
        if(DOPRINTS2()) writeMesh("debug_collapse1.meshb",msh);
        msh.poi2ent(ipcol,0) = -1;
        msh.poi2ent(ipcol,1) = -1;
      }
      if(msh.param->dbgfull){
        check_topo(msh, msh.nbpoi, msh.npoin, msh.nedge, msh.nface, msh.nelem, ithrd2); 
      }
      if(info.done) return 0;

      CPRINTF1(" - return qmax = %f \n",info.qmax_end);
      if(info.qmax_end < qmabest && ierro == 0){
        CPRINTF1(" - new best quality !\n");
        qmabest = info.qmax_end;
        ivebest = iver;
      }
    }
  }

  if(ivebest == -1){
    // Failed to find a reinsertion candidate. 
    // Let's not rely on exceptions as failure here is not invalid behaviour 
    // and could happen quite frequently 
    return 1;
  }

  if(ivebest == 0){ // Otherwise ball hasn't changed !
    ierro = ball2(msh,ipcol,iface,
                  cav.lcfac,
                  cav.lcedg,
                  &iopen,&imani,ithrd1);
    // Increase cavity with Delaunay criterion
    ierro = increase_cavity_Delaunay(msh, cav, ithrd2);
    METRIS_ASSERT(ierro == 0);
      
    //int nprem = increase_cavity_lenedg(msh,cav,ipins,ithrd2,ithrd3);
    ierro = increase_cavity2D(msh,cav,ithrd2);
    METRIS_ASSERT(ierro == 0);
  }

  if(DOPRINTS2()) writeMeshCavity("collapse_cavity0.meshb", msh, cav, ithrd2);

  opts.dryrun = false;

  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
  }}CT_FOR1(ideg);

  METRIS_ASSERT(ierro == 0);
  METRIS_ASSERT(abs(info.qmax_end - qmabest) < 1.0e-6);
  METRIS_ASSERT(info.done);

  return 0;
}



template int collversurf<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                   int iface, int iver, double qmax_suf, 
                   MshCavity &cav, CavWrkArrs &work, 
                   intAr1 &lerro, int ithrd1, int ithrd2);
template int collversurf<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                   int iface, int iver, double qmax_suf, 
                   MshCavity &cav, CavWrkArrs &work, 
                   intAr1 &lerro, int ithrd1, int ithrd2);


} // end namespace

