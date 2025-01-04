//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_lineforce.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "low_increasecav.hxx"
#include "../cavity/msh_cavity.hxx"

#include "../low_geo.hxx"
#include "../low_topo.hxx"
#include "../io_libmeshb.hxx"


namespace Metris{

// Evaluate points on lines and forcibly reinsert them in the mesh, breaking 
// everything on the way 
template<class MFT>
void reinsertLines(Mesh<MFT> &msh, int ithrd1, int ithrd2){
  MshCavity cav(0,100,2);
  cav.lcedg.set_n(2);
  // There could be one single face for 2 edges even not corner 
  cav.lcfac.set_n(0);

  CavOprInfo info;
  CavWrkArrs work;
  CavOprOpt opts;
  opts.dryrun = false;
  opts.allow_remove_points = true;


  // Absolute tolerance for whether a point is on the geometry or not 
  // (temporary)
  const double geotol = msh.param->geo_abstoledg;
  const int iverb = msh.param->iverb;
  const double vtol = msh.param->vtol;

  double result[18];


  msh.tag[ithrd1]++;
  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;

    int iref = msh.edg2ref[iedge];
    ego obj = msh.CAD.cad2edg[iref];

    // In case two insertions are needed -> edge is killed on first insert
    int ip[2] = {msh.edg2poi(iedge,0), msh.edg2poi(iedge,1)};

    for(int iver = 0; iver < 2; iver++){
      int ipoin = ip[iver];
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;


      int ibpoi = msh.poi2bpo[ipoin];
      // Don't move corners, assumed in position 
      if(msh.bpo2ibi(ibpoi,1) == 0) continue;
      METRIS_ASSERT(msh.bpo2ibi(ibpoi,1) == 1);

      METRIS_ASSERT_MSG(msh.edg2ref[msh.bpo2ibi(ibpoi,2)] == iref,
       "ipoin="<<ipoin<<" ibpoi = "<<ibpoi<<" bpo2ibi[2] = "
        <<msh.bpo2ibi(ibpoi,2)<< 
        " iref = "<<iref<<" edg2ref = "<<msh.edg2ref[msh.bpo2ibi(ibpoi,2)]);


      // Proceed to reinsertion? First test distance. 
      int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
      if(ierro != 0){
        printf("  ## EG_getRange failed %d \n",ierro);
        METRIS_THROW_MSG(GeomExcept()," EG_Evaluate failed")  
      }

      double dst; 
      if(msh.idim == 2){
        dst = geterrl2<2>(msh.coord[ipoin], result);
      }else{
        dst = geterrl2<3>(msh.coord[ipoin], result);
      }

      if(dst < geotol*geotol) continue;

      // Not on geometry, reinsert 
      if(iverb >= 2) printf("   - Point %d not on geometry ref %d, dist = %15.7e > "
                                  "%15.7e = tol\n", ipoin, iref,sqrt(dst), geotol);
      if(iverb >= 3){
        writeMesh("dbg_geometry_pt"+std::to_string(ipoin),msh);
      }

      cav.lcfac.set_n(0);
      cav.lcedg.set_n(0);
      cav.ipins = ipoin;

      double coor0[3];
      for(int ii = 0; ii < msh.idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
      for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipoin,ii) = result[ii];

      // Start from ball2 so we can check if a reinsertion is even required. 
      // (no-op)
      // Edge could have been killed if ii == 1
      int iedg1 = msh.poi2ent(ipoin,0);
      METRIS_ASSERT(msh.poi2ent(ipoin,1) == 1);
      int iface = msh.edg2fac[iedg1];
      bool imani;
      int  iopen;
      ball2(msh,ipoin,iface,cav.lcfac,cav.lcedg,&iopen,&imani,ithrd2);
      METRIS_ASSERT(iopen);
      METRIS_ASSERT(cav.lcedg.get_n() == 2);

      bool ireins = false;
      for(int icfac : cav.lcfac){
        // Dummy normal 
        bool iflat;
        if(msh.idim == 2) getmeasentP1<2,2>(msh.fac2poi[icfac], msh.coord, vtol, &result[4], &iflat);
        else              getmeasentP1<3,2>(msh.fac2poi[icfac], msh.coord, vtol, &result[4], &iflat);

        if(iflat){
          ireins = true;
          break; 
        }
      }

      if(!ireins) continue;

      ierro = increase_cavity2D(msh,msh.fac2ref[iface],msh.coord[cav.ipins],opts,cav,ithrd2);
      if(ierro != 0){
        std::cout<<"increase_cavity2D failed "<< ierro<<"\n";
        printf("Trying to reinsert ipoin %d iedge = %d iedg1 %d iface %d iref %d \n",
          cav.ipins, iedge, iedg1, iface, iref);
        writeMesh("increase_fail", msh);
        writeMeshCavity("increase_cav2D_fail", msh, cav);
        METRIS_THROW(GeomExcept());
      }


      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
        ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
      }}CT_FOR1(ideg);

      METRIS_ENFORCE_MSG(ierro == 0, "Cavity failed ierro = "<<ierro);


    }
  }

}

template
void reinsertLines(Mesh<MetricFieldAnalytical> &msh, int ithrd1, int ithrd2);
template
void reinsertLines(Mesh<MetricFieldFE        > &msh, int ithrd1, int ithrd2);


}// end Namespace
