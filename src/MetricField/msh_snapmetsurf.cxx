//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "msh_snapmetsurf.hxx"
#include "../Mesh/MeshMetric.hxx"
#include "../low_geo.hxx"
#include "../utils/mprintf.hxx"
#include "../low_normal.hxx"

#include <unordered_set>

namespace Metris{


// If surface normal direction close enough to metric direction closest to it,
// snap the metric. 
// Repeat for edge tangents. 
template<class MetricFieldType>
void snapMetSurf(MeshMetric<MetricFieldType> &msh,
                 int ithread){

  double surfdir[3];
  const int gdim = msh.idim;

  // Store keys (ipoin, iref) to keep track of operations done (lower dim)
  std::unordered_set<std::pair<int, int>, tup2_hash::hash> poireftag[2];

  intAr1 ntryref1(msh.CAD.ncaded), ntryref2(msh.CAD.ncadfa);
  intAr1 ncorref1(msh.CAD.ncaded), ncorref2(msh.CAD.ncadfa);
  dblAr1 maxerrminref1(msh.CAD.ncaded), maxerrminref2(msh.CAD.ncaded);
  dblAr1 minerrminref1(msh.CAD.ncaded), minerrminref2(msh.CAD.ncaded);
  ntryref1.fill(0);
  ntryref2.fill(0);
  ncorref1.fill(0);
  ncorref2.fill(0);
  maxerrminref1.fill(-1);
  maxerrminref2.fill(-1);
  minerrminref1.fill( 1);
  minerrminref2.fill( 1);

  for(int tdim = gdim-1; tdim >= 1; tdim--){
    INCVDEPTH(msh.param);
    int nentt = msh.nentt(tdim);
    const intAr2& ent2poi = msh.ent2poi(tdim);
    const intAr1& ent2ref = msh.ent2ref(tdim);
    intAr1& ntryref = tdim == 1 ? ntryref1 : ntryref2;
    intAr1& ncorref = tdim == 1 ? ncorref1 : ncorref2;
    dblAr1& maxerrminref = tdim == 1 ? maxerrminref1 : maxerrminref2;
    dblAr1& minerrminref = tdim == 1 ? minerrminref1 : minerrminref2;

    msh.tag[ithread]++;
    int ncorr = 0, ntry = 0, nerro = 0;

    double maxerrmin = -1, minerrmin = 1;

    for(int ientt = 0; ientt < nentt; ientt++){
      INCVDEPTH(msh.param);
      int iref = ent2ref[ientt];

      for(int iver = 0; iver < tdim + 1; iver++){
        int ipoin = ent2poi(ientt,iver);

        // Normals need only be set once for surface interior, once per 
        // (ipoin,iref) pair otherwise. 
        int pdim = msh.getpoitdim(ipoin);
        if(pdim == tdim){
          if(msh.poi2tag(ithread,ipoin) >= msh.tag[ithread]) continue;
          msh.poi2tag(ithread,ipoin) = msh.tag[ithread];
        }else{
          std::pair<int, int> key = {ipoin, iref};
          if(poireftag[tdim-1].find(key) != poireftag[tdim-1].end()) continue;
          poireftag[tdim-1].insert(key);
        }
        ntry++;
        ntryref[iref]++;

        if(tdim == 2){
          getnorpoiref(msh, ipoin, iref, surfdir);
        }else{
          gettanpoiref(msh, ipoin, iref, surfdir);
        }

        double eigval[3], eigvec[9];
        if(gdim == 2) geteigsym<2>(msh.met[ipoin], eigval, eigvec);
        else          geteigsym<3>(msh.met[ipoin], eigval, eigvec);

        int imin = -1;
        double errmin = 1.0e30;
        //printf("Debug surfdir = ");
        //dblAr1(gdim, surfdir).print();
        for(int ii = 0; ii < gdim; ii++){
          double dtprd = gdim == 2 ? getprdl2<2>(&eigvec[gdim*ii], surfdir)
                                   : getprdl2<3>(&eigvec[gdim*ii], surfdir);
          double err = MAX(1 - abs(dtprd),0);
          if(err < errmin){
            errmin = err;
            imin = ii;
          }
        }

        CPRINTF1(" - ipoin %d iref %d tdim %d errmin %e\n",ipoin,iref,tdim,errmin);
        maxerrmin = MAX(maxerrmin, errmin);
        minerrmin = MIN(minerrmin, errmin);

        maxerrminref[iref] = MAX(maxerrminref[iref], errmin);
        minerrminref[iref] = MIN(minerrminref[iref], errmin);

        if(errmin >= msh.param->met_snap_tol*2) continue;

        ncorr++;
        ncorref[iref]++;

        // To snap, we simply replace the closest eigenvector with surfdir
        // Then we apply Gram-Schmidt, which gives us the closest vectors in L2
        // norm to the original. 

        for(int ii = 0; ii < gdim; ii++) eigvec[gdim*imin + ii] = surfdir[ii];

        int ivec1 = (imin + 1)%gdim;
        int ivec2 = (ivec1+ 1)%gdim;

        double dtprd = gdim == 2 ? getprdl2<2>(&eigvec[gdim*ivec1], surfdir)
                                 : getprdl2<3>(&eigvec[gdim*ivec1], surfdir);
        for(int ii = 0; ii < gdim; ii++) 
          eigvec[gdim*ivec1 + ii] -= dtprd*surfdir[ii];
        int ierro = gdim == 2 ? normalize_vec<2>(&eigvec[gdim*ivec1])
                              : normalize_vec<3>(&eigvec[gdim*ivec1]);
        METRIS_ASSERT(ierro == 0);
        if(ierro != 0){
          nerro++;
          continue;
        }


        double dtpr1 = gdim == 2 ? getprdl2<2>(&eigvec[gdim*ivec2], surfdir)
                                 : getprdl2<3>(&eigvec[gdim*ivec2], surfdir);
        double dtpr2 = gdim == 2 ? getprdl2<2>(&eigvec[gdim*ivec2], &eigvec[gdim*ivec1])
                                 : getprdl2<3>(&eigvec[gdim*ivec2], &eigvec[gdim*ivec1]);
        for(int ii = 0; ii < gdim; ii++) 
          eigvec[gdim*ivec2 + ii] -= dtpr1*surfdir[ii]
                                   + dtpr2*eigvec[gdim*ivec1 + ii];
        ierro = gdim == 2 ? normalize_vec<2>(&eigvec[gdim*ivec2])
                          : normalize_vec<3>(&eigvec[gdim*ivec2]);
        METRIS_ASSERT(ierro == 0);
        if(ierro != 0){
          nerro++;
          continue;
        }

        if(gdim == 2) eig2met<2>(eigval, eigvec, msh.met[ipoin]);
        else          eig2met<3>(eigval, eigvec, msh.met[ipoin]);

      }// for iver
    }// for ientt
    CPRINTF1("-- DONE tdim = %d snapped %d/%d (%4.1f%%) metrics to surface\n"
             ,tdim,ncorr,ntry,ncorr/(double)ntry*100);
    CPRINTF1("        errmin min = %5.2e max = %5.2e\n",minerrmin,maxerrmin)
    if(nerro > 0) CPRINTF1("        %d/%d (%4.1f%%) errors\n",
                           nerro,ntry,nerro/(double)ntry*100);

    if(DOPRINTS2()){
      int nref = tdim == 1 ? msh.CAD.ncaded : msh.CAD.ncadfa;
      for(int iref = 0; iref < nref; iref++){
        if(ncorref[iref] == 0) continue;
        CPRINTF2("   - iref %d snapped %d/%d (%4.1f%%) err min = %9.2e max = %9.2e\n",iref,
               ncorref[iref],ntryref[iref],ncorref[iref]/(double)ntryref[iref]*100,
               minerrminref[iref],maxerrminref[iref]);
      }
    }
  }// for tdim

}

template
void snapMetSurf<MetricFieldAnalytical>(MeshMetric<MetricFieldAnalytical> &msh,
                                        int ithread);
template
void snapMetSurf<MetricFieldFE        >(MeshMetric<MetricFieldFE        > &msh,
                                        int ithread);

}//namespace
