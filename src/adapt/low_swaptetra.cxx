//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_swap.hxx"
#include "msh_swap.hxx" // for swapOptions

#include "../Mesh/Mesh.hxx"

#include "../linalg/det.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../utils/aux_misc.hxx"
#include "../low_lenedg.hxx"
#include "../utils/mprintf.hxx"
#include "../low_normal.hxx"
#include "../io_libmeshb.hxx"
#include "../quality/low_metqua.hxx"

namespace Metris{


// tet swaps: try face and edge swaps. 
// if lazy, take first quality improving operation. 
// Tet swaps can only be quality, not length, based (for now and perhaps ever).
template<class MFT, int ideg>
int swaptetra(Mesh<MFT>& msh, int itetr, swapOptions opt, 
             MshCavity &cav, CavWrkArrs &work, 
             double *qnrm0_, double *qnrm1_, 
             int ithread){
  INCVDEPTH(msh.param);
  constexpr int tdim = 3;
  constexpr int gdim = 3;
  constexpr AsDeg asdmet = AsDeg::P1;

  const bool ilazy = true;

  METRIS_ENFORCE_MSG(ilazy == true, "Non-lazy swaptetra not implemented");

  const int qpnorm = opt.swap_norm;
  METRIS_ENFORCE_MSG(qpnorm == 0, 
           "Only infinity norm of element quality is supported in tetra swaps");
  const int qpower = msh.param->opt_power;

  CavOprOpt opts;
  opts.allow_topological_correction = false;
  opts.skip_topo_checks = true;
  opts.allow_remove_points = false;
  opts.cache_tetra_quality = true;

  // Precompute initial tetra quality and store in cavity hash table.
  cav.reset(); // this resets the quality hash table
  double quael = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,itetr,qpower,
                                       qpnorm,1.0);
  auto key = stup4(msh.tet2poi(itetr,0),msh.tet2poi(itetr,1),
                   msh.tet2poi(itetr,2),msh.tet2poi(itetr,3));
  cav.qtetr[key] = quael;

  if(ilazy){
    // In case lazy, simply set a quality threshold. 
    opts.dryrun = false;
    opts.qmax_nec = quael*0.99; // Only accept if new cavity improves on this tet.
  }else{
    opts.dryrun = true;
    opts.qmax_suf = quael*0.5; // Dryrun but accept if final quality considerably better.
  }


  for(int ifa = 0; ifa < 4; ifa++){
    int ierro = aux_swaptetface<MFT,ideg>(msh, itetr, ifa, quael, cav, opts, work, 
                                          qnrm0_, qnrm1_, ithread);
    if(ierro < 0){
      CPRINTF1("- Accepted tet %d swap face %d quality %e -> %e\n",
               itetr, ifa, *qnrm0_, *qnrm1_);
      return 0;
    }
  }

  return 1;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template int swaptetra<MetricFieldAnalytical,n>(Mesh<MetricFieldAnalytical>& msh, \
                                    int itetr, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1,int ithread);\
template int swaptetra<MetricFieldFE        ,n>(Mesh<MetricFieldFE        >& msh, \
                                    int itetr, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1,int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


// Returns > 0 if error, 0 if no error and no operation, -1 if operation done.
template<class MFT, int ideg>
int aux_swaptetface(Mesh<MFT>& msh, int itetr, int ifacl, double quae1,
                    MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,
                    double *qnrm0_, double *qnrm1_,
                    int ithread){
  GETVDEPTH(msh.param);
  constexpr int tdim = 3;
  constexpr int gdim = 3;
  constexpr AsDeg asdmet = AsDeg::P1;

  const int qpnorm = msh.param->opt_pnorm;
  const int qpower = msh.param->opt_power;

  METRIS_ASSERT_MSG(qpnorm == 0, "Implement proper quality norms (number of elts varies, weighting)")

  double &qnrm0 = *qnrm0_;
  double &qnrm1 = *qnrm1_;

  int itet2 = msh.tet2tet(itetr, ifacl);
  if(itet2 < 0) return 1;

  int idom1 = msh.tet2ref[itetr];
  if(idom1 != msh.tet2ref[itet2]){
    CPRINTF1("# Tetra face swap error: crosses refs %d != %d\n",idom1, msh.tet2ref[itet2]);
    return 1;
  }

  // Compute second tetrahedron quality for computing qnrm0. 
  // We cache this quality for the following cases:
  // - The caller had done dryruns and is now effecting the operation -> skip computation
  // - The caller will eventually call an edge-based swap, its qnrm0 will involve
  // this element if the edge is on the face.
  auto key = stup4(msh.tet2poi(itet2,0),msh.tet2poi(itet2,1),
                   msh.tet2poi(itet2,2),msh.tet2poi(itet2,3));
  auto tt = cav.qtetr.find(key);
  double quae2;
  if(tt != cav.qtetr.end()){
    quae2 = tt->second;
    CPRINTF2(" - found cached quality for neighbour tet %d: %e\n",itet2,quae2);
  }else{
    quae2 = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,itet2,
                                  qpower,qpnorm,1.0);
    cav.qtetr[key] = quae2;
  }

  if(qpnorm == 0){
    qnrm0 = MAX(quae1,quae2);
  }else if (qpnorm > 0){
    qnrm0 = pow(pow(quae1,qpnorm) + pow(quae2,qpnorm), 1.0/qpnorm);
  }


  int ifa2 = -1;
  for(int itmp = 0; itmp < 4; itmp++){
    if(msh.tet2tet(itet2,itmp) != itetr) continue;
    ifa2 = itmp;
    break;
  }
  METRIS_ASSERT(ifa2 >= 0);

  // Get point opposite face in second tetra. 
  int ipopp = msh.tet2poi(itet2,ifa2);

  // Do not call reset()! We need the quality hash table to persist between calls.
  cav.lcedg.set_n(0);
  cav.lcfac.set_n(0);
  cav.lctet.set_n(0);

  cav.ipins = ipopp;

  cav.lctet.stack(itetr);
  cav.lctet.stack(itet2);


  CavOprInfo info;
  int ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithread);
  qnrm1 = info.qcav3;
  CPRINTF1("- aux_swaptetface called cavity, ierro = %d info.done = %d qnrm1 = %e\n",
           ierro,info.done,qnrm1);

  if(info.done) return -1;

  return ierro;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template int aux_swaptetface<MetricFieldAnalytical,n>(Mesh<MetricFieldAnalytical>& msh, \
                        int itetr, int ifacl, double quae1,\
                        MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,\
                        double *qnrm0_, double *qnrm1_,\
                        int ithread);\
template int aux_swaptetface<MetricFieldFE,n>(Mesh<MetricFieldFE>& msh, \
                        int itetr, int ifacl, double quae1,\
                        MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,\
                        double *qnrm0_, double *qnrm1_,\
                        int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



} // end namespace