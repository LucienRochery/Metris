//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include <bunit/common_setup.hxx>

#include <boost/timer/progress_display.hpp>

#include "../src/ho_constants.hxx"
#include "../src/io_libmeshb.hxx"
#include "../src/utils/aux_misc.hxx"
#include "../src/quality/low_metqua.hxx"
#include "../SANS/Surreal/SurrealS.h"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
#include "../src/low_geo/misc.hxx"
#include "../src/Adaptation/low_cavqual.cxx"
#include <boost/hana.hpp> 

using namespace Metris;

typedef MetricFieldAnalytical MFT;

// -- Test metqua3_xi_d 
// Constant metric fields should yield reliable derivatives in all cases 
// In non constant metric fields, derivatives only defined for DoFs in back element
// interiors... 

int get_nrempt(MeshBase &msh, const MshCavity &cav, int ithrd1){
  int tdim = cav.get_tdim();
  const intAr1 &lcent = cav.lcent(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const intAr2 &ent2ent = msh.ent2ent(tdim);
  intAr2 &ent2tag = msh.ent2tag(tdim);
  
  msh.tag[ithrd1]++;

  for(int ientt : lcent) ent2tag(ithrd1, ientt) = msh.tag[ithrd1];
  for(int ientt : lcent){
    for(int ifa = 0; ifa < tdim + 1; ifa++){
      int ienei = ent2ent(ientt, ifa);
      if(ienei >= 0 && ent2tag(ithrd1,ienei) == msh.tag[ithrd1]) continue;
      for(int iver = 0; iver < tdim + 1; iver++){
        if(iver == ifa) continue;
        int ipoin = ent2poi(ientt,iver);
        msh.poi2tag(ithrd1, ipoin) = msh.tag[ithrd1];
      }
    }
  }
  // Now count the internal points. Boundary points are tagged.
  int nrempt = 0;
  for(int ientt : lcent){
    for(int iver = 0; iver < tdim + 1; iver++){
      int ipoin = ent2poi(ientt,iver);
      if(msh.poi2tag(ithrd1, ipoin) == msh.tag[ithrd1]) continue;
      msh.poi2tag(ithrd1, ipoin) = msh.tag[ithrd1];
      // An internal point, and one not seen yet.
      nrempt++;
    }
  }
  return nrempt;
}

BOOST_AUTO_TEST_CASE(test_eval3) 
{

  // bool is whether straight
  std::vector<std::string> meshes = {
     METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
    ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.500"
    }; 


  CavOprOpt opts;
  CavOprInfo info;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = false;
  opts.dryrun = false;
  //opts.allow_remove_points = icollapse; 
  opts.allow_remove_points = true; 
  opts.allow_remove_points_superdim = true; // For boundary
  opts.qmax_nec = -1;
  opts.qmax_suf = -1;
  opts.qmax_iff = -1;

  const int ithrd1 = 0;
  const int ithrd2 = 1;

  for(std::string mesh : meshes)
  { 

    cargHandler arg("-in " + mesh + "  -anamet 1 -verb 0 -vdepth 0 -prefix tmp/");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g); 
    run.degElevate();


    std::cout<<"------------------------------------------------\n";
    std::cout<<"Mesh "<<mesh<<"\n";

    int tdim = msh.get_tdim();
    double pi = 3.141592653589793238462643383279502884;
    const double opt_dens = msh.get_tdim() == 2 ? pi / 4 : 0.54;
    const intAr2 &ent2poi = msh.ent2poi(tdim);
    const int nedgl = (tdim*(tdim+1))/2;
    const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                              tdim == 2 ? lnoed2[0] : lnoed3[0]);

    const int medge = tdim == 2 ? msh.nface : msh.nelem;
    HshTab_I2I ledge;
    ledge.reserve(medge);

    double t0s = get_cpu_time();
    int nedge_tot = 0;
    for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
      INCVDEPTH(msh.param);
      if(isdeadent(ientt,ent2poi)) continue;
      for(int ied = 0; ied < nedgl; ied++){
        INCVDEPTH(msh.param);

        // Check edge already seen
        int ip1 = ent2poi(ientt, lnoed(ied,0));
        int ip2 = ent2poi(ientt, lnoed(ied,1));
        auto key = stup2(ip1,ip2);
        if(ledge.find(key) != ledge.end()) continue;

        nedge_tot++;

        double sz[2];
        double len = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,ientt,tdim,ied,sz) :
                                     getlenedg_geosz<MFT,3,1>(msh,ientt,tdim,ied,sz);
        if(len < sqrt(2)) continue;
        ledge[key] = ientt;

        // New long edge, add to stack
        CPRINTF1(" + long edge %d %d seed %d len = %e\n",ip1,ip2,ientt,len);

      }// for ied
    }// for ientt
    double t1s = get_cpu_time();
    printf(" - init time %f nlong = %d\n",t1s-t0s,(int)ledge.size());

    //msh.param->iverb = 5;
    //msh.param->ivdepth = 10;

    double dens0, dens1;
    double qumin0, qumin1, qumax0, qumax1, quavg0, quavg1;
    int iedge = 0;
    for(auto edge_it : ledge){
      INCVDEPTH(msh.param);
      int ientt = edge_it.second;
      METRIS_ASSERT(ientt >= 0);

      std::tuple<int,int> key = edge_it.first;
      int ip1 = std::get<0>(key);
      int ip2 = std::get<1>(key);

      MshCavity cav(100,100,1);
      int iopen;
      shell(msh,ip1,ip2,tdim,ientt,cav.lcedg,cav.lcfac,cav.lctet,&iopen);
      if(cav.get_tdim() < tdim) continue;
      cav.inewp = 0; // Start at no point, to get initial density
      collrejcav_dens(msh,cav,&dens0,&dens1,ithrd1,ithrd2);
      
      MPRINTF(" - Init cav edge %d face %d tetra %d dens0 %f dens1 %f optimal %f\n",
              cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n(),dens0,dens1,opt_dens);

      cav.ipins = msh.newpoitopo(tdim,ientt);
      printf("Debug ipins = %d \n",cav.ipins);
      cav.inewp = 1; // Now flag as new point for future density computations
      for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins,ii) = (msh.coord(ip1,ii) + msh.coord(ip2,ii))/2.0;
      int ierro = msh.interpMetBack(cav.ipins,tdim,ientt,msh.ent2ref(tdim)[ientt],NULL);
      if(ierro != 0){
        printf("## INTERPMETBACK ERROR %d \n",ierro);
        continue;
      }

      double coor0[3];
      for(int ii = 0; ii < msh.idim; ii++) coor0[ii] = msh.coord(cav.ipins,ii);
      double met0[6];
      const int nnmet = (msh.idim*(msh.idim+1))/2;
      for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(cav.ipins,ii);

      collrejcav_dens(msh,cav,&dens0,&dens1,ithrd1,ithrd2);
      MPRINTF(" - With +1 point, dens0 %f dens1 %f optimal %f\n",dens0,dens1,opt_dens);

      getquacav(msh,cav,&qumin0,&qumin1,&qumax0,&qumax1,&quavg0,&quavg1,ithrd1);
      MPRINTF(" - Init quacav: min0 %f min1 %f max0 %f max1 %f avg0 %f avg1 %f\n",
              qumin0,qumin1,qumax0,qumax1,quavg0,quavg1);

      iedge++;

      writeMeshCavity("cav0_" + std::to_string(iedge),msh,cav);

      int ncedg0 = cav.lcedg.get_n();
      int ncfac0 = cav.lcfac.get_n();
      int nctet0 = cav.lctet.get_n();


      static std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp;
      double lenqua_short_max = 1;
      ierro = setCavityInsertion(msh,cav,opts,10,lenqua_short_max,nocomp,ithrd1,ithrd2);
      getquacav(msh,cav,&qumin0,&qumin1,&qumax0,&qumax1,&quavg0,&quavg1,ithrd1);
      printf(" - setCavityInsertion  returned %d ncent %d qumin %f -> %f max %f -> %f avg %f -> %f\n",ierro,
             cav.lcent(tdim).get_n(),qumin0,qumin1,qumax0,qumax1,quavg0,quavg1); 
      writeMeshCavity("cav1_" + std::to_string(iedge),msh,cav);


      for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins,ii) = coor0[ii];
      for(int ii = 0; ii < nnmet; ii++) msh.met(cav.ipins,ii) = met0[ii];
      cav.lcedg.set_n(ncedg0);
      cav.lcfac.set_n(ncfac0);
      cav.lctet.set_n(nctet0);
      ierro = setCavityInsertion2(msh,cav,ientt,10,ithrd1,ithrd2);
      getquacav(msh,cav,&qumin0,&qumin1,&qumax0,&qumax1,&quavg0,&quavg1,ithrd1);
      printf(" - setCavityInsertion2 returned %d ncent %d qumin %f -> %f max %f -> %f avg %f -> %f\n",ierro,
             cav.lcent(tdim).get_n(),qumin0,qumin1,qumax0,qumax1,quavg0,quavg1);
      writeMeshCavity("cav2_" + std::to_string(iedge),msh,cav);

             

      wait();
      



      continue;

      ncedg0 = cav.lcedg.get_n();
      ncfac0 = cav.lcfac.get_n();
      nctet0 = cav.lctet.get_n();

      for(int niter = 0; niter < 5; niter++){

        cav.lcedg.set_n(ncedg0);
        cav.lcfac.set_n(ncfac0);
        cav.lctet.set_n(nctet0);

        MPRINTF(" - long edge %d %d seed %d\n",ip1,ip2,ientt);
        ierro = movePointCavLen<MFT>(msh,cav,tdim,ientt,msh.ent2ref(tdim)[ientt],5,0);
        MPRINTF(" - ierro = %d\n",ierro);

        ierro = increase_cavity_Delaunay(msh, cav, -1, ithrd1);
        METRIS_ENFORCE(ierro == 0);
        ierro = increase_cavity(msh, cav, false, ithrd1, ithrd2);
        METRIS_ENFORCE(ierro == 0);
        collrejcav_dens(msh,cav,&dens0,&dens1,ithrd1,ithrd2);
        const int ncent = cav.lcent(tdim).get_n();

        int nrempt = get_nrempt(msh, cav, ithrd1);

        getquacav(msh,cav,&qumin0,&qumin1,&qumax0,&qumax1,&quavg0,&quavg1,ithrd1);
        MPRINTF(" - iter %d + del, quacav: min0 %f min1 %f max0 %f max1 %f avg0 %f avg1 %f\n",
                niter,qumin0,qumin1,qumax0,qumax1,quavg0,quavg1);         
        MPRINTF(" - iter %d + del, dens0 %f dens1 %f nentt %d nremp %d\n",niter,dens0,dens1,ncent,nrempt);
      }

      


      writeMeshCavity("cav1_" + std::to_string(iedge),msh,cav);

      // Do one step adding short points and check quality afterwards (from prints)

      int nprem = increase_cavity_lenedg(msh,cav,opts,cav.ipins,ithrd1,ithrd2);
      MPRINTF(" - +remp nedge %d nface %d nelem %d nprem = %d\n",
              cav.lcedg.get_n(),cav.lcfac.get_n(),cav.lctet.get_n(),nprem);

      writeMeshCavity("cav2_" + std::to_string(iedge),msh,cav);

      collrejcav_dens(msh,cav,&dens0,&dens1,ithrd1,ithrd2);
      MPRINTF(" - With + %d nprem, dens0 %f dens1 %f optimal %f\n",nprem,dens0,dens1,opt_dens);

      ierro = movePointCavLen<MFT>(msh,cav,tdim,ientt,msh.ent2ref(tdim)[ientt],10,0);
      MPRINTF(" - ierro = %d\n",ierro);

      writeMeshCavity("cav3_" + std::to_string(iedge),msh,cav);
      wait();

    }


  }// fortest cases

}// end boost test case











