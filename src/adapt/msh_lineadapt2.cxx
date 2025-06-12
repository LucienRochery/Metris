//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_lineadapt.hxx"
#include "msh_lineforce.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../low_lenedg.hxx"
#include "../msh_structs.hxx"
#include "../low_geo.hxx"
#include "../utils/aux_timer.hxx"
#include "../utils/aux_misc.hxx"
#include "../msh_checktopo.hxx"
#include "../aux_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../adapt/low_increasecav.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../utils/mprintf.hxx"


namespace Metris{

// Adapt lines "frontally" using CAD
// In this rewrite, we first generate points in t space, and only then insert 
// them in bulk. 
template<class MFT>
void adaptGeoLines2(Mesh<MFT> &msh){
  GETVDEPTH(msh.param);
  int ithrd1 = 0;
  int ithrd2 = 1;
  int ithrd3 = 2;

  if(!msh.CAD()) return;

  CPRINTF1("-- Start adaptGeoLines.\n");

  if(msh.param->dbgfull) check_topo(msh, ithrd1);

  MetSpace ispac0 = msh.met.getSpace();
  msh.met.setSpace(MetSpace::Exp);

  METRIS_ASSERT(msh.CAD.ncaded > 0); 

  // Multiply and divide target length by this factor for admissible new
  // edge length interval 
  const double lentolfac = msh.param->geo_lentolfac;
  if(lentolfac < 1.0){
    METRIS_THROW_MSG(WArgExcept(),"lentolfac < 1! = "<<lentolfac);
  }

  // From each ref, point onto a corner that has an edge, of this ref, inciding
  intAr2 ref2cor(msh.CAD.ncaded, 2);
  dblAr2 ref2rng(msh.CAD.ncaded, 2);
  ref2cor.fill(-1);

  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;

    int iref = msh.edg2ref[iedge];
    if(ref2cor(iref,0) >= 0 && ref2cor(iref,1) >= 0) continue;

    for(int ii = 0; ii < 2; ii++){
      int ip = msh.edg2poi(iedge,ii);
      int ib = msh.poi2bpo[ip];
      if(msh.bpo2ibi(ib,1) == 0){
        ib = msh.poi2ebp(ip,1,iedge,-1);
        METRIS_ASSERT(ib >= 0);
        double tcorn = msh.bpo2rbi(ib,0);
        int iupd = 0;
        if(ref2cor(iref,0) >= 0) iupd = 1;
        ref2cor(iref,iupd) = ip;
        ref2rng(iref,iupd) = tcorn;
        break;
      }
    }
  }

  // Order by ascending param
  for(int iref = 0; iref < msh.CAD.ncaded; iref++){
    // There can be missing edges (degenerate)
    if(ref2cor(iref,0) < 0) continue;
    METRIS_ASSERT(ref2cor(iref,0) >= 0 && ref2cor(iref,1) >= 0);

    double t1 = ref2rng(iref,0);
    double t2 = ref2rng(iref,1);

    if(t1 >= t2){
      ref2rng(iref,0) = t2;
      ref2rng(iref,1) = t1;
      int itmp = ref2cor(iref,0);
      ref2cor(iref,0) = ref2cor(iref,1);
      ref2cor(iref,1) = itmp;
    }
  }



  dblAr1 crv_lens(msh.CAD.ncaded);
  {
  double t0 = get_wall_time();
  getCADCurveLengths(msh, (lentolfac - 1.0), crv_lens);
  double t1 = get_wall_time();
  CPRINTF1(" - getCADCurveLengths time %f \n",t1-t0);
  }


  int mcfac = 100, mcedg = 10;
  MshCavity cav(0,mcfac,mcedg);

  CavOprOpt opts;
  CavWrkArrs work;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = false;
  opts.allow_remove_points = true;
  opts.dryrun = false;
  opts.geodev1 = 1.0; // lax

  int CADedgtag = 1;
  std::map<ego,int> CADedg2tag, CADedg2ref;
  for(int iref = 0; iref < msh.CAD.ncaded; iref++){
    ego edg = msh.CAD.cad2edg[iref];
    CADedg2ref[edg] = iref;
    CADedg2tag[edg] = 0;
  }

  dblAr1 lnewt(100);
  intAr1 &ledge = msh.iwork;
  ledge.set_n(0);

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

    // Loop over CAD edges and remesh each one 
    for(int iCADed = 0; iCADed < nchild; iCADed++){
      INCVDEPTH(msh.param);

      ego obj = lchild[iCADed];
      if(CADedg2tag[obj] >= CADedgtag) continue;
      CADedg2tag[obj] = CADedgtag;

      int iref = CADedg2ref[obj];
      int icor0 = ref2cor(iref,0);
      if(icor0 < 0){
        CPRINTF1(" - Loop %d line %d / %d is degenerate -> skip\n",
                 iloop, iref+1, msh.CAD.ncaded);
        continue;
      }

      CPRINTF1(" - Loop %d adapt line %d / %d \n", iloop, iref+1, msh.CAD.ncaded);

      //double range[2];
      //ego CADed = msh.CAD.cad2edg[iref];
      //int iperi;
      //int ierro = EG_getRange(CADed,range,&iperi);
      //if(ierro != 0){
      //  printf("  ## EG_getRange failed %d \n",ierro);
      //  METRIS_ASSERT(ierro == 0);
      //  return;
      //}
      //if(range[0] > range[1]){
      //  double tmp = range[0];
      //  range[0] = range[1];
      //  range[1] = tmp;
      //}

      // Generate points on the curve 
      genPointsCurve<MFT>(msh, iref, icor0, crv_lens[iref], ref2rng[iref], lnewt, ledge);

      if(lnewt.get_n() == 0){
        CPRINTF1(" # Warning line %d with len %f -> no points to insert\n",
                  iref, crv_lens[iref]);
        continue;
      }
      // Insert them 
      insPointsCurve<MFT>(msh, iref, ref2rng[iref], ref2cor[iref], lnewt, ledge, ithrd1, ithrd2, ithrd3);

      if(DOPRINTS2()) writeMesh("lineadapt" + std::to_string(iref),msh);
    } // for iCADed
  } // for iloop

  msh.met.setSpace(ispac0);
  msh.cleanup();
}


template void adaptGeoLines2<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh);
template void adaptGeoLines2<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh);





} //namespace Metris