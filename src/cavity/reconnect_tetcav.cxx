//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../cavity/msh_cavity.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../aux_topo.hxx"
#include "../ho_constants.hxx"
#include "../low_geo/measure.hxx"
#include "../quality/low_metqua.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/aux_misc.hxx"
#include "../io_libmeshb.hxx"

namespace Metris{

int aux_check_tetbdry(MeshBase &msh, MshCavity& cav, int ielen, int ifa0, int *cfatag, int ithread);

template <class MFT, int ideg>
int reconnect_tetcav(Mesh<MFT> &msh,
                     MshCavity& cav,
                     CavOprOpt  &opts,
                     CavOprInfo &info,
                     int nfac0, double *qmax, int ithread){

  const int nctet = cav.lctet.get_n();
  if(nctet <= 0) return 0;

  bool check_qua = (opts.qmax_nec > 0 && msh.get_tdim() == 3)
                || (opts.qmax_suf > 0 && msh.get_tdim() == 3)
                || (opts.qmax_iff > 0 && msh.get_tdim() == 3);

  //// When would we not though?
  //bool check_val = opts.fast_reject
  //              || opts.max_increase_cav_geo <= 0
  //              || check_qua;

  const int qpnorm = msh.param->opt_pnorm;
  *qmax = -1;
  info.qcav3 = qpnorm == 0 ? -1 : 0;
  info.objective_numerator_end = 0.;
  info.objective_target_weight_end = 0.;
  info.objective_element_count_end = 0;
  int nqual_cav3 = 0;
  double targetWeightCav3 = 0.;

  // We need to increment this independently per element, avoid conflict with
  // ithread
  int cfatag = msh.tag[ithread];

  GETVDEPTH(msh.param);
  CPRINTF1(" - START reconnect_tetcav check_qua = {} qpnorm {} qpower {}\n",
            check_qua,qpnorm,msh.param->opt_power);

  const int nele0 = msh.nelem;

  // In the future, handle connex components here. For now, unique ref.
  // If the cavity has two+ tet connex components, then ipins is necessarily
  // a tdim 2 or 1 point. It is on the interface beetween the connex cpts.
  // In that case, the faces nfac0 through nface do not generate tetrahedra.
  // Hence the only case where tets are created with new faces is if point
  // is interior. In that case it is quite simple, the tetrahedra are oriented
  // same as face with ipins in first position, and we are always the first fac2tet.

  for(int ielem : cav.lctet){
    msh.tet2tag(ithread, ielem) = msh.tag[ithread];
  }

  // Hash faces containing ipins for neighbour updates.
  // Can be boundary faces, or new interior faces.
  // As they all share ipins, it is enough to hash the remaining two vertices.
  // Such faces hit at most 2 elements; we need only keep track of one, as
  // when we are the second, we update info and can trash the entry.
  HshTab_I2I2 facHsh;
  for(int iface = nfac0; iface < msh.nface; iface++){
    int iver = msh.template getverfac<1>(iface, cav.ipins);
    METRIS_ASSERT(iver >= 0);
    int ip1 = msh.fac2poi(iface,lnoed2[iver][0]);
    int ip2 = msh.fac2poi(iface,lnoed2[iver][1]);
    auto key = stup2(ip1,ip2);
    facHsh[key] = {2, iface};
  }

  // "Hash" edges containing ipins for HO updates
  // Loose arrays we'll point into using poi2tag on the
  // non-ipins vertex.
  // Helper struct that resets poi2tag on destruction.
  struct edgHashHO_t{
    edgHashHO_t(MeshBase &msh_, const MshCavity& cav, int nfac0, int ithread_) :
      msh(msh_), ithread(ithread_), ipins(cav.ipins),
      // Compilers complain about std::move here but it is necessary for
      // how work arrays work: only one live at a time, and update the tracking with
      // the move constructor.
      // Note: we could allow copy constructor for untracked work arrays, but let's not
      // complicate things.
      edgHshDim(std::move(msh.get_iwork(cav.lctet.get_n()*(msh.curdeg > 1)))),
      edgHshElt(std::move(msh.get_iwork(cav.lctet.get_n()*(msh.curdeg > 1)))),
      edgHshEdg(std::move(msh.get_iwork(cav.lctet.get_n()*(msh.curdeg > 1)))){
        edgHshElt.set_n(0);
        edgHshDim.set_n(0);
        edgHshEdg.set_n(0);
        // Initialize with new triangles; new edges from edges have been
        // now used by triangles, so this is enough.
        int nhash = 0;
        for(int iface = nfac0; iface < msh.nface; iface++){
          for(int ied = 0; ied < 3; ied++){
            int ipoin = -1;
            if(ipins == msh.fac2poi(iface,lnoed2[ied][0])){
              ipoin = msh.fac2poi(iface,lnoed2[ied][1]);
            }else if(ipins == msh.fac2poi(iface,lnoed2[ied][1])){
              ipoin = msh.fac2poi(iface,lnoed2[ied][0]);
            }else continue;
            if(msh.poi2tag(ithread,ipoin) >= msh.tag[ithread]) continue;
            msh.poi2tag(ithread,ipoin) = msh.tag[ithread] + nhash;
            nhash++;
            edgHshDim.stack(2);
            edgHshElt.stack(iface);
            edgHshEdg.stack(ied);
          }
        }
      };
    ~edgHashHO_t(){
      int nhash = edgHshDim.get_n();
      for(int ihash = 0; ihash < nhash; ihash++){
        int tdim  = edgHshDim[ihash];
        int ientt = edgHshElt[ihash];
        int ied   = edgHshEdg[ihash];
        int ipoi1 = msh.ent2poi(tdim)(ientt, lnoed[tdim-1][ied][0]);
        int ipoi2 = msh.ent2poi(tdim)(ientt, lnoed[tdim-1][ied][1]);
        // cav.ipins also, doesn't matter
        msh.poi2tag(ithread,ipoi1) = msh.tag[ithread];
        msh.poi2tag(ithread,ipoi2) = msh.tag[ithread];
      }
    }
    // tuple return values can be conveniently taken as
    // auto [a, b ,c] = func();
    std::tuple<int,int,int> operator[](int ihash) const{
      return std::make_tuple(edgHshDim[ihash], edgHshElt[ihash], edgHshEdg[ihash]);
    }

    int ihash(int ipoi1, int ipoi2){
      int ipoin = ipoi1 == ipins ? ipoi2 : ipoi1;
      return msh.poi2tag(ithread,ipoin) - msh.tag[ithread];
    }
    void add(int tdim, int ientt, int ied){
      edgHshDim.stack(tdim);
      edgHshElt.stack(ientt);
      edgHshEdg.stack(ied);
      int ipoi1 = msh.ent2poi(tdim)(ientt, lnoed[tdim-1][ied][0]);
      int ipoi2 = msh.ent2poi(tdim)(ientt, lnoed[tdim-1][ied][1]);
      int ipoin = ipoi1 == ipins ? ipoi2 : ipoi1;
      msh.poi2tag(ithread,ipoin) = msh.tag[ithread] + edgHshDim.get_n() - 1;
    }

    MeshBase& msh;
    const int ithread;
    const int ipins;
    intWrkAr1 edgHshDim; // Host element tdim
    intWrkAr1 edgHshElt; // Host element
    intWrkAr1 edgHshEdg; // Local edge
  };
  edgHashHO_t edgHashHO(msh, cav, nfac0, ithread);

  // Hash external edges for control point copies.
  // A tetra can share no face (but an edge) with an external element.
  // If we don't hash these edges, we're reliant on order of control point copies
  // and it can take arbitrarily many steps to go from an tetra that copied
  // the ctrl point to another.
  HshTab_I2I2 extEdgHsh;
  if constexpr (ideg > 1){
    // Any external edge is necessarily in an external face.
    // Hence there is at least one cavity element that sees it as
    // an edge within one of its external faces.
    for(int ielem : cav.lctet){
      INCVDEPTH(msh.param);
      for(int ifa = 0; ifa < 4; ifa++){
        int ienei = msh.tet2tet(ielem,ifa);
        CPRINTF3(" - extEdgHsh check iecav {} ifa {} : {} {} {}\n",
                 ielem, ifa, msh.tet2poi(ielem,lnofa3[ifa][0]),
                             msh.tet2poi(ielem,lnofa3[ifa][1]),
                             msh.tet2poi(ielem,lnofa3[ifa][2]));
        if(ienei >= 0 && msh.tet2tag(ithread,ienei) >= msh.tag[ithread]) continue;
        for(int ii = 0; ii < 3; ii++){
          int ied = ledfa3[ifa][ii];
          int ipoi1 = msh.tet2poi(ielem,lnoed3[ied][0]);
          int ipoi2 = msh.tet2poi(ielem,lnoed3[ied][1]);
          auto key = stup2(ipoi1, ipoi2);
          // emplace adds if key doesn't exist
          extEdgHsh.emplace(key, std::make_pair(ielem, ied));
          CPRINTF2(" - extEdgHsh add edge {} {} ielem {} ied {} = {} {}\n",
                   ipoi1,ipoi2,ielem,ied,ipoi1,ipoi2);
        }
      }
    }

  }


  // To check if all points are on the boundary.
  const int pdim_ipins = msh.getpoitdim(cav.ipins);

  // Next loop over cavity boundary faces to generate new tets. If there are
  // new faces, they do not generate tetrahedra, as those would have ipins twice.
  for(int iele0 : cav.lctet){
    INCVDEPTH(msh.param);

    for(int ifa0 = 0; ifa0 < 4; ifa0++){

      int ienei = msh.tet2tet(iele0,ifa0);

      if(ienei >= 0 && msh.tet2tag(ithread,ienei) >= msh.tag[ithread]){
        CPRINTF1(" - skip tetra creation from {} face {} due to cavity neighbour {} \n",
                 iele0, ifa0, ienei);
        continue;
      }

      // Also skip if face contains ipins already.
      bool iskip = false;
      for(int ii = 0; ii < 3; ii++){
        if(msh.tet2poi(iele0, lnofa3[ifa0][ii]) != cav.ipins) continue;
        iskip = true;
        break;
      }
      if(iskip){
        CPRINTF1(" - skip tetra creation from {} face {} which contains ipins\n",
                 iele0, ifa0);
        continue;
      }

      // Additionally, skip any faces that coincide with triangles which were
      // in the cavity. These have disappeared from the mesh.
      {// namespace
        auto key = stup3(msh.tet2poi(iele0, lnofa3[ifa0][0]),
                         msh.tet2poi(iele0, lnofa3[ifa0][1]),
                         msh.tet2poi(iele0, lnofa3[ifa0][2]));
        auto tt = msh.facHshTab.find(key);
        if(tt != msh.facHshTab.end()){
          int iface = tt->second;
          if(msh.fac2tag(ithread,iface) >= msh.tag[ithread]){
            CPRINTF1(" - skip tetra creation from {} face {} is glo face {} in cavity\n",
                     iele0, ifa0 ,iface);
            continue;
          }
        }
      }

      // Cavity boundary face, create tet.
      int ielen = msh.nelem;
      msh.set_nelem(msh.nelem+1);

      msh.tet2tet(ielen, 0) = -1;
      msh.tet2tet(ielen, 1) = -1;
      msh.tet2tet(ielen, 2) = -1;
      msh.tet2tet(ielen, 3) = -1;
      msh.tet2tet(ielen, ifa0) = ienei;
      msh.tet2poi(ielen, ifa0) = cav.ipins;
      msh.tet2poi(ielen, lnofa3[ifa0][0]) = msh.tet2poi(iele0, lnofa3[ifa0][0]);
      msh.tet2poi(ielen, lnofa3[ifa0][1]) = msh.tet2poi(iele0, lnofa3[ifa0][1]);
      msh.tet2poi(ielen, lnofa3[ifa0][2]) = msh.tet2poi(iele0, lnofa3[ifa0][2]);

      if constexpr(ideg > 1){
        int nnode = getnnod3(ideg);
        for(int ii = 4; ii < nnode; ii++){
          msh.tet2poi(ielen,ii) = -1;
        }
      }

      CPRINTF1(" - new tetra = {} from iele0 = {} ifa = {} vertices: {} {} {} {}\n",
               ielen,iele0,ifa0,msh.tet2poi(ielen,0),msh.tet2poi(ielen,1),
               msh.tet2poi(ielen,2),msh.tet2poi(ielen,3));

      // High-order nodes are still unset here, so only the provisional P1
      // skeleton can be checked. correct_cavity assigns those nodes and then
      // applies the full-dimensional classifier to the completed tetrahedron.
      double meas0;
      if(!isvalideltP1<3,3>(msh, ielen,NULL,&meas0)){
        CPRINTF1(" - iflat ! return ip1 ip2 ip3 ip4 = {} {} {} {} meas = {:15.7e} \n",
                 msh.tet2poi(ielen,0),msh.tet2poi(ielen,1),msh.tet2poi(ielen,2),
                 msh.tet2poi(ielen,3),meas0);
        return CAV_ERR_FLATTET;
      }
      CPRINTF2(" - new tetra volume {} \n",meas0);


      // If the tet has all vertices on boundary, check faces orientations.
      if(pdim_ipins < 3){
        int ierro = aux_check_tetbdry(msh, cav, ielen, ifa0, &cfatag, ithread);
        if(ierro != 0) return ierro;
      }

      #ifdef STEPDISTANCE
      constexpr QuaFun iquaf = QuaFun::StepDistance;
      #else
      constexpr QuaFun iquaf = QuaFun::SizeShape;
      #endif

      if(check_qua){
        // Regardless of degree, verify underlying P1 element is decent enough
        double quael;
        if(opts.cache_tetra_quality){
          auto key = stup4(msh.tet2poi(ielen,0),msh.tet2poi(ielen,1),
                           msh.tet2poi(ielen,2),msh.tet2poi(ielen,3));
          auto tt = cav.qtetr.find(key);
          if(tt != cav.qtetr.end()){
            quael = tt->second;
            CPRINTF2(" - found cached quality\n");
          }else{
            #ifdef TESTQUALITYALGO
            CPRINTF1("Using distance-based qual in cavity operator\n");
            quael = metqua<MFT,3,3,iquaf>(msh,AsDeg::P1,AsDeg::P1,ielen,1.0);
            #else
            quael = metqua<MFT,3,3>(msh,AsDeg::P1,AsDeg::P1,ielen,1.0);
            #endif
            cav.qtetr[key] = quael;
          }
        }else{
          #ifdef TESTQUALITYALGO
          CPRINTF1("Using dsitance-based qual in cavity operator\n")
          quael = metqua<MFT,3,3,iquaf>(msh,AsDeg::P1,AsDeg::P1,ielen,1.0);
          #else
          quael = metqua<MFT,3,3>(msh,AsDeg::P1,AsDeg::P1,ielen,1.0);
          #endif
        }
        CPRINTF1(" - new tetra {} = {} {} {} {} from {} conf error = {} \n",
          ielen,
          msh.tet2poi(ielen,0), msh.tet2poi(ielen,1), msh.tet2poi(ielen,2),
          msh.tet2poi(ielen,3) ,iele0,quael);
        *qmax = MAX(quael, *qmax);
        #ifdef STEPDISTANCE
        if(msh.param->step_distance_cavity_target_average){
          const double objectiveWeight =
              step_distance_element_target_weight<MFT,3,3>(
                  msh,AsDeg::P1,ielen);
          info.objective_numerator_end += objectiveWeight*quael;
          info.objective_target_weight_end += objectiveWeight;
        }else{
          info.objective_numerator_end += quael;
        }
        #else
        info.objective_numerator_end += quael;
        #endif
        info.objective_element_count_end++;
        if(qpnorm == 0){
          info.qcav3 = MAX(info.qcav3, quael);
        }else{
          const double poweredQuality = pow(abs(quael),qpnorm);
          #ifdef STEPDISTANCE
          if(msh.param->step_distance_cavity_target_average){
            const double targetWeight =
                step_distance_element_target_weight<MFT,3,3>(
                    msh,AsDeg::P1,ielen);
            info.qcav3 += step_distance_region_contribution(
                poweredQuality,targetWeight,true);
            targetWeightCav3 += targetWeight;
          }else{
            info.qcav3 += poweredQuality;
          }
          #else
          info.qcav3 += poweredQuality;
          #endif
        }
        nqual_cav3++;
        if(quael > opts.qmax_nec && opts.qmax_nec > 0.0){
          CPRINTF1(" # quael = {} > {} = qmax_nec -> reject\n",quael, opts.qmax_nec);
          return CAV_ERR_QMAXNEC; // Run rejected
        }
        if(quael > opts.qmax_iff && opts.qmax_iff > 0.0) return CAV_ERR_QMAXIFF; // Run rejected
        if(quael < 0) return CAV_ERR_QFACNEG;  // Run rejected
      }

      msh.tet2ref[ielen] = msh.tet2ref[iele0];


      if constexpr (ideg > 1){

        // Copy ctrl pts from mother element
        // Capping at ideg-1 and preventing ii == jj == 0 skips the vertices.
        for(int ii = 0; ii <= ideg-1; ii++){
          for(int jj = (int)(ii == 0); jj <= ideg-1; jj++){
            int kk = ideg - ii - jj;

            int inode = ifa0 == 0 ? mul2nod(0,ii,jj,kk)
                      : ifa0 == 1 ? mul2nod(kk,0,ii,jj)
                      : ifa0 == 2 ? mul2nod(jj,kk,0,ii)
                                  : mul2nod(ii,jj,kk,0);

            CPRINTF1(" - new tet node {} copy {} from spawning face {} of {}\n",
                    inode,msh.tet2poi(iele0, inode), ifa0, iele0);

            msh.tet2poi(ielen,inode) = msh.tet2poi(iele0, inode);
          }
        }

        // Copy ctrl pts from external edges
        int idx1[4], idx2[4];
        // Edges in face opposite ipins are not the only candidates
        bool taged[6] = {};
        for(int ied = 0; ied < 6; ied++){
          int ipoi1 = msh.tet2poi(ielen,lnoed3[ied][0]);
          int ipoi2 = msh.tet2poi(ielen,lnoed3[ied][1]);
          auto key = stup2(ipoi1, ipoi2);
          auto tt = extEdgHsh.find(key);
          CPRINTF2(" - new tet {} check ied {} in extEdgHsh : found = {}\n",ielen,ied,tt != extEdgHsh.end());
          if(tt == extEdgHsh.end()) continue;
          taged[ied] = true;
          const auto [iecav, ied1] = tt->second;
          // Copy ied1 of iecav into ied of ielen
          CPRINTF1(" - new tet {} edge {} copy ext edge {} of {}\n",
                   ielen, ied, ied1, iecav);
          for(int ii = 0; ii < 4; ii++) idx1[ii] = 0;
          for(int ii = 0; ii < 4; ii++) idx2[ii] = 0;
          bool isori = true;
          if(ipoi1 == msh.tet2poi(iecav, lnoed3[ied1][1])) isori = false;
          for(int ii = 1; ii <= ideg-1; ii++){
            idx1[lnoed3[ied][0]] = ii;
            idx1[lnoed3[ied][1]] = ideg - ii;
            idx2[lnoed3[ied1][0]] = isori ? ii        : ideg - ii;
            idx2[lnoed3[ied1][1]] = isori ? ideg - ii : ii;
            int inod1 = mul2nod(3,idx1);
            int inod2 = mul2nod(3,idx2);
            CPRINTF2(" - copy ctrl pt {}{}{}{} inode {} = {} from {}{}{}{} inode {}\n",
                     idx1[0],idx1[1],idx1[2],idx1[3], inod1,
                     msh.tet2poi(iecav, inod2),
                     idx2[0],idx2[1],idx2[2],idx2[3], inod2);
            msh.tet2poi(ielen, inod1) = msh.tet2poi(iecav, inod2);
          }
        }

        // Copy control points from edges containing ipins stored in edgHashHO
        for(int ie0 = 0; ie0 < 3; ie0++){
          // edge impinging on ipins
          int ied = ledno3[ifa0][ie0];
          if(taged[ied]) continue;
          int ipoi1 = msh.tet2poi(ielen,lnoed3[ied][0]);
          int ipoi2 = msh.tet2poi(ielen,lnoed3[ied][1]);
          int ihash = edgHashHO.ihash(ipoi1, ipoi2);
          CPRINTF1(" - new tet {} check ied {} : {} {}, ihash = {}\n",
                  ielen, ied, ipoi1, ipoi2, ihash);

          if(ihash >= 0){
            // These are ints
            const auto [tdim, iele1, ied1] = edgHashHO[ihash];
            const intAr2& ent2poi = msh.ent2poi(tdim);

            CPRINTF1(" - copy from tdim {} elt {} edge {} : {} {}\n",
                    tdim, iele1, ied1,
                    ent2poi(iele1, lnoed[tdim-1][ied1][0]),
                    ent2poi(iele1, lnoed[tdim-1][ied1][1]));
            METRIS_ASSERT(tdim == 2 || tdim == 3);
            METRIS_ASSERT(iele1 >= 0 && iele1 < msh.nentt(tdim));
            METRIS_ASSERT(ied1 >= 0 && ied1 < (tdim == 2 ? 3 : 6));

            // Copy control points
            int idx1[4] = {}, idx2[4] = {};
            bool isori = true;
            if(ipoi1 == ent2poi(iele1, lnoed[tdim-1][ied1][1])) isori = false;
            for(int ii = 1; ii <= ideg-1; ii++){
              idx1[lnoed3[ied][0]] = ii;
              idx1[lnoed3[ied][1]] = ideg - ii;
              idx2[lnoed[tdim-1][ied1][0]] = isori ? ii        : ideg - ii;
              idx2[lnoed[tdim-1][ied1][1]] = isori ? ideg - ii : ii;
              int inod1 = mul2nod(3,idx1);
              int inod2 = mul2nod(tdim,idx2);
              CPRINTF3("idx2 = {} {} {} {} inod2 = {}\n",idx2[0],idx2[1],idx2[2],idx2[3],inod2);
              CPRINTF3(" - copy ctrl pt iele1 node {} = {} to node {}\n",
                      inod2,ent2poi(iele1, inod2), inod1);
              msh.tet2poi(ielen, inod1) = ent2poi(iele1, inod2);
            }
          }else{ // if ihash < 0
            // Create new HO nodes
            int idx[4] = {};
            int ipoi1 = msh.tet2poi(ielen,lnoed3[ied][0]);
            int ipoi2 = msh.tet2poi(ielen,lnoed3[ied][1]);
            for(int ii = 1; ii <= ideg-1; ii++){
              idx[lnoed3[ied][0]] = ii;
              idx[lnoed3[ied][1]] = ideg - ii;
              double ci = ((double)ii)/ideg;
              double cj = ((double)(ideg - ii))/ideg;
              int inode = mul2nod(3,idx);
              int ipnew = msh.newpoitopo(PointType::CtrlPt, 3, ielen);
              CPRINTF3(" - new tet {} edge {} ctrl pt inode {} : {}\n",ielen,ied,inode,ipnew);
              for(int kk = 0; kk < msh.idim; kk++){
                msh.coord(ipnew,kk) = ci*msh.coord(ipoi1,kk)
                                    + cj*msh.coord(ipoi2,kk);
              }
              msh.tet2poi(ielen,inode) = ipnew;
            }
            edgHashHO.add(3, ielen, ied);

          }// if ihash >= 0
        }
      }

      // Loop over the three faces impinging on ipins:
      // - update neighbors
      // - create/copy HO nodes
      for(int ied = 0; ied < 3; ied++){
        // ifa is the node we want to avoid (ipins)
        // lnoed[ied][0] and lnoed[ied][1] will span all combinations of
        // 2 among 0,1,2
        // to this we add 1, and go over ifa + 1, ifa + 2 and ifa + 3
        // modulo 4, which spans all pairs of vertices not ifa.
        int ip1 = msh.tet2poi(ielen, (ifa0 + lnoed2[ied][0]+1)%4);
        int ip2 = msh.tet2poi(ielen, (ifa0 + lnoed2[ied][1]+1)%4);

        // This makes the face opposite ifa + ied.
        int ifan = (ifa0 + ied + 1)%4;

        CPRINTF1(" - with ifa0 {} have ifan {} \n",ifa0, ifan);

        auto key = stup2(ip1, ip2);
        auto tt = facHsh.find(key);
        if(tt != facHsh.end()){
          const auto [tdimf, ielef] = tt->second;
          METRIS_ASSERT(tdimf == 2 || tdimf == 3);
          METRIS_ASSERT(ielef >= 0 && ielef < msh.nentt(tdimf));

          CPRINTF1(" - found internal or new boundary face tdim {} entt {} \n",
                   tdimf,ielef);

          // Copy the nodes from the entity. Could be dimension 2 or 3.
          // Also update neighbours
          if(tdimf == 2){

            CPRINTF1(" - bdry face nodes {} {} {} \n",msh.fac2poi(ielef,0)
                     ,msh.fac2poi(ielef,1),msh.fac2poi(ielef,2));
            CPRINTF1("   match with {} {} {} \n",msh.tet2poi(ielen,lnofa3[ifan][0])
              ,msh.tet2poi(ielen,lnofa3[ifan][1]),msh.tet2poi(ielen,lnofa3[ifan][2]));

            if constexpr (ideg > 1) cpy_glofac2tetfac<ideg>(msh, ielef, ielen, ifan);

            // This is a new face. Update the fac2tet
            METRIS_ASSERT(ielef >= nfac0);
            if(msh.fac2tet(ielef,0) < 0){
              // Nothing here yet
              msh.fac2tet(ielef,0) = ielen;
            }else if(msh.fac2tet(ielef,0) >= nele0){
              // 0 points to a new one, put in 1
              msh.fac2tet(ielef,1) = ielen;
            }else{
              METRIS_THROW_MSG( "New face points to dead element")
            }

          }else{
            METRIS_ASSERT(ielef >= nele0);
            // Copy from a tetrahedron. For this we need to know the ifa for
            // the tet.
            int ifaf = getfactet(msh, ielef, ip1, ip2, cav.ipins);
            METRIS_ASSERT(ifaf >= 0);

            if constexpr (ideg > 1) cpy_tetfac2tetfac<ideg>(msh, ielef, ifaf, ielen, ifan);

            // Also update the neighbour. This is a cavity interior face and both are
            // new elements.
            msh.tet2tet(ielen,ifan) = ielef;
            msh.tet2tet(ielef,ifaf) = ielen;
          }

        }else{// tt != facHsh.end()
          // Create control points and add this interior face to hash table.
          // Loop over face indices all but vertices and edges
          int idx[4] = {};
          for(int ii = 1; ii <= ideg - 2; ii++){
            for(int jj = 1; jj <= ideg - 1 - ii; jj++){
              int kk = ideg - ii - jj;
              idx[lnofa3[ifan][0]] = ii;
              idx[lnofa3[ifan][1]] = jj;
              idx[lnofa3[ifan][2]] = kk;
              int inod3 = mul2nod(3,idx);
              if(msh.tet2poi(ielen,inod3) >= 0){
                CPRINTF1(" - new tet {} face {} ctrl pt {} already exists: {}\n",
                         ielen,ifan,inod3,msh.tet2poi(ielen,inod3));
                continue;
              }
              int ipnew = msh.newpoitopo(PointType::CtrlPt, 3, ielen);
              msh.tet2poi(ielen,inod3) = ipnew;
              for(int mm = 0; mm < 3; mm++){
                msh.coord(ipnew,mm) = ((double)ii)/ideg*msh.coord(msh.tet2poi(ielen,lnofa3[ifan][0]),mm)
                                    + ((double)jj)/ideg*msh.coord(msh.tet2poi(ielen,lnofa3[ifan][1]),mm)
                                    + ((double)kk)/ideg*msh.coord(msh.tet2poi(ielen,lnofa3[ifan][2]),mm);
              }
              CPRINTF3(" - new tet {} face {} ctrl pt {} coord {} {} {} = {}xcoord({}) + {}xcoord({}) + {}xcoord({})\n",ielen,ifan,ipnew,
                       msh.coord(ipnew,0),msh.coord(ipnew,1),msh.coord(ipnew,2),
                      ((double)ii)/ideg, msh.tet2poi(ielen,lnofa3[ifan][0]),
                      ((double)jj)/ideg, msh.tet2poi(ielen,lnofa3[ifan][1]),
                      ((double)kk)/ideg, msh.tet2poi(ielen,lnofa3[ifan][2]) );
            }
          }

          // Additionally, add the new face to hash table.
          facHsh[key] = {3, ielen};

        }// tt != facHsh.end()

      }// for ied


      // Create interior control points. These are all those of the form (1,1,1,1) + (i,j,k,l)
      for(int ii = 0; ii < ideg - 4; ii++){
        for(int jj = 0; jj < ideg - 4; jj++){
          for(int kk = 0; kk < ideg - 4; kk++){
            int ll = 1 - ii - jj - kk;
            int inod3 = mul2nod(1+ii,1+jj,1+kk,1+ll);
            int ipnew = msh.newpoitopo(PointType::CtrlPt,3, ielen);
            msh.tet2poi(ielen,inod3) = ipnew;

            for(int mm = 0; mm < msh.idim; mm++){
              msh.coord(ipnew,mm) = (1.0+ii)/ideg*msh.coord(msh.tet2poi(ielen,0),mm)
                                  + (1.0+jj)/ideg*msh.coord(msh.tet2poi(ielen,1),mm)
                                  + (1.0+kk)/ideg*msh.coord(msh.tet2poi(ielen,2),mm)
                                  + (1.0+ll)/ideg*msh.coord(msh.tet2poi(ielen,3),mm);
            }
          }
        }
      }// for ii

      for(int ii = 0; ii < METRIS_MAXTAGS; ii++) msh.tet2tag(ii,ielen) = 0;


    }// for iver
  }// for iele0

  if(check_qua && qpnorm > 0){
    #ifdef STEPDISTANCE
    if(nqual_cav3 > 0){
      info.qcav3 = step_distance_region_objective(
          info.qcav3,targetWeightCav3,
          msh.param->step_distance_cavity_target_average);
    }
    #endif
    info.qcav3 = pow(info.qcav3, 1.0 / qpnorm);
    CPRINTF1(" - Final tetra cavity quality = {:.3f}\n",info.qcav3);
  }

  return 0;
}


// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int reconnect_tetcav<MetricFieldAnalytical, n >(Mesh<MetricFieldAnalytical> &msh,\
 MshCavity& cav, CavOprOpt &opts, CavOprInfo &info, int nfac0, double *qmax, int ithread);\
template int reconnect_tetcav<MetricFieldFE        , n >(Mesh<MetricFieldFE        > &msh,\
 MshCavity& cav, CavOprOpt &opts, CavOprInfo &info, int nfac0, double *qmax, int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




int aux_check_tetbdry(MeshBase &msh, MshCavity& cav, int ielen, int ifa0, int *cfatag, int ithread){
  GETVDEPTH(msh.param);

  for(int ii = 0; ii < 3; ii++){
    int pdim = msh.getpoitdim(msh.tet2poi(ielen, lnofa3[ifa0][ii]));
    if(pdim < 3) continue;
    CPRINTF1(" face point {} has dim 3 -> skip\n",msh.tet2poi(ielen, lnofa3[ifa0][ii]));
    return 0;
  }

  // The previous test was more lenient: it allowed tetrahedra to have
  // all four vertices on the boundary, provided that the tet faces
  // were reverse oriented to the boundary faces they supported.
  // This meant the tet was inside the domain.
  // One issue is a tetrahedron could be created with only one (perhaps even none)
  // face on the boundary, and later be "uncovered".
  // New criterion is more conservative but does not depend on future topology
  // we simply check if no four vertices are on the same face reference.

  // tag[ithread] has been incremented for other uses already, and no
  // conflict with cfa2tag:
  (*cfatag)++;
  cav.maxtag = MAX(cav.maxtag,*cfatag + 4);
  for(int iver = 0; iver < 4; iver++){
    //INCVDEPTH(msh.param);
    int ipoin = msh.tet2poi(ielen, iver);
    for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
      int tdimb = msh.bpo2ibi(ibpoi,1);
      if(tdimb != 2) continue;
      int iface = msh.bpo2ibi(ibpoi,2);
      int ireff = msh.fac2ref[iface];
      int nseenm1 = msh.cfa2tag(ithread,ireff) - *cfatag;
      CPRINTF1(" - check ver {} ipoin {} face {} ref {} already seen {}x\n",
                iver,ipoin,iface,ireff,nseenm1+1);
      // Same ref can be seen twice by same vertex, meaning it has already
      // been updated.
      METRIS_ASSERT_MSG(nseenm1 <= iver, "nseenm1 = {} iver = {}", nseenm1, iver);

      // We must be careful to skip refs that have been skipped by at
      // least one prior point, because our method of counting is
      // setting tag to tag + iver. Otherwise, if only the last point sees
      // a face ref, we'd think all did.
      if(nseenm1 + 1 < iver && iver != 0){
        CPRINTF1(" - inactive ref -> skip\n");
        continue;
      }

      // If has already been seen 3 times before,
      if(nseenm1 + 1 == 3 && iver == 3){
        CPRINTF1(" # error 4 points on same face ref\n");
        return CAV_ERR_BDRYTET2;
      }
      // An edge point can see the same face ref several times. We mustn't
      // increment.
      msh.cfa2tag(ithread,ireff) = *cfatag + iver;
    }
  }
  *cfatag += 4;

  #if 0
  // All vertices are on the boundary. Go over faces and compare orientation
  bool nogood = false;
  for(int ifal = 0; ifal < 4; ifal++){
    int iface = msh.tet2fac(ielen, ifal);
    CPRINTF1(" - tet face {} -> glo face {}\n",ifal, iface);
    if(iface < 0) continue;
    if(msh.fac2tet(iface,0) >= 0 && msh.fac2tet(iface,1) >= 0){
      if(msh.param->dbgfull)
        METRIS_THROW_MSG("TODO: Handle internal surface case here");
      CPRINTF1(" # REJECT: tet face {} matches glo face {} which is internal\n",
        ifal, iface);
      return CAV_ERR_BDRYTET2;
    }
    int ip1 = msh.tet2poi(ielen, lnofa3[ifal][0]);
    int ip2 = msh.tet2poi(ielen, lnofa3[ifal][1]);
    int ip3 = msh.tet2poi(ielen, lnofa3[ifal][2]);
    int jp1 = msh.fac2poi(iface, 0);
    int jp2 = msh.fac2poi(iface, 1);
    int jp3 = msh.fac2poi(iface, 2);
    CPRINTF1(" - glo face {} = {} {} {}: check opposite orientation\n",
            iface, jp1, jp2, jp3);
    if(ip1 == jp1){
      if(ip2 == jp2) nogood = true;
    }else if(ip1 == jp2){
      if(ip2 == jp3) nogood = true;
    }else if(ip1 == jp3){
      if(ip2 == jp1) nogood = true;
    }else{
      METRIS_THROW_MSG( "Tet local <-> global face mismatch");
    }
    if(nogood) break;
  }
  if(nogood){
    CPRINTF1(" # REJECT: tet face matches glo face orientation\n");
    return CAV_ERR_BDRYTET;
  }
#endif
  return 0;
}


}
