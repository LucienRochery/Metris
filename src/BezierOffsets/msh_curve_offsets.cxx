//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_curve_offsets.hxx"
#include "low_gaps.hxx"
#include "msh_bez2gaps.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../Mesh/Mesh.hxx"

#include "../quality/msh_metqua.hxx"
#include "../ho_constants.hxx"
#include "../LPopt/msh_maxccoef.hxx"

#include "../io_libmeshb.hxx"
#include "../aux_exceptions.hxx"
#include "../aux_histogram.hxx"
#include "../aux_topo.hxx"
#include "../aux_timer.hxx"

namespace Metris{
	 


#define BEZGAPS_LP

template <class MFT, int gdim, int ideg>
int curveMeshOffsets(Mesh<MFT> &msh, bool icorr){
  constexpr int tdim = gdim;
  static_assert(gdim == 2 || gdim == 3);

  double t0 = get_wall_time();


	if(ideg > 2) METRIS_THROW_MSG(TODOExcept(), "Implement ideg > 2 curveMeshOffsets.");

  if constexpr(ideg <= 1) return 0;

	FEBasis ibas0 = msh.getBasis();
	msh.setBasis(FEBasis::Bezier);
  MetSpace ispac0 = msh.met.getSpace();
	msh.met.setSpace(MetSpace::Log);

  int iverb = msh.param->iverb;


  constexpr int nnmet = (gdim*(gdim+1))/2;
  const intAr2 &ent2poi = gdim == 2 ? msh.fac2poi : msh.tet2poi;
  const int nentt = gdim == 2 ? msh.nface : msh.nelem;


	double offset[3];
	int ithread = 0; // This can easily be parallelized.

  #ifdef BEZGAPS_LP
    intAr1 idx_point;
    dblAr2 coor0;
    int npopt = 0;

    if(icorr){
      idx_point.allocate(msh.npoin);
      idx_point.set_n(msh.npoin);
      idx_point.fill(-1);
      coor0.allocate(msh.npoin,gdim);
      coor0.set_n(msh.npoin);


      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        for(int ii = 0; ii < gdim; ii++) coor0(ipoin,ii) = msh.coord(ipoin,ii) ;
      }
    }

  #endif

  bez2gaps<ideg>(msh);
    

	msh.tag[ithread]++;
	for(int ientt = 0; ientt < nentt; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;

		// Not a typo: number of edges is also number of metric components
		for(int ie = 0; ie < nnmet; ie++){
			int inode = gdim + 1 + ie;
			int ipoih = ent2poi(ientt,inode);
			if(msh.poi2bpo[ipoih] >= 0) continue;
			if(msh.poi2tag(ithread,ipoih) >= msh.tag[ithread]) continue;
			msh.poi2tag(ithread,ipoih) = msh.tag[ithread];
      getBezOffsetsEdge<MFT, gdim, ideg>(msh,gdim,ent2poi[ientt], ie,offset);
	   for(int ii = 0; ii < gdim; ii++) msh.coord(ipoih,ii) = offset[ii];
      #ifdef BEZGAPS_LP
        if(icorr){
          idx_point[ipoih] = npopt;
          npopt++;
        }
      #endif
		}
	}

  gaps2bez<ideg>(msh);



  if(iverb >= 2){
    msh.setBasis(FEBasis::Lagrange);
    writeMesh("crv0.meshb",msh);
    double qmin, qmax, qavg;
    bool iinva;
    dblAr1 lquae, dum = {0.1, 0.9};
    if(iverb >= 1){
      if(msh.idim == 2){
        getmetquamesh<MFT,2,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
      }else{
        getmetquamesh<MFT,3,AsDeg::Pk>(msh,&iinva,&qmin,&qmax,&qavg,&lquae);
      }
      print_histogram(msh,lquae,IntrpTyp::Geometric,dum,"q","Element quality");
    }
    msh.setBasis(FEBasis::Bezier);
  } 


  if(!icorr) return 0;


  #ifndef BEZGAPS_LP
  // Pre-correction: control points should not leave their shells. 
  msh.tag[ithread]++;
  for(int ientt = 0; ientt < nentt; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;

    for(int ie = 0; ie < nnmet; ie++){
      int inode = gdim + 1 + ie;
      int ipoih = ent2poi(ientt,inode);
      if(msh.poi2bpo[ipoih] >= 0) continue;
      if(msh.poi2tag(ithread,ipoih) >= msh.tag[ithread]) continue;
      msh.poi2tag(ithread,ipoih) = msh.tag[ithread];

      double bary[tdim+1];
      if constexpr(tdim == 2){

        // Start with current element: could be the limiting one. 
        // This is only if bary opposite vertex > 0
        inventP1<gdim>(msh.fac2poi[ientt],msh.coord,msh.coord[ipoih],bary);
        // Set + tol as we can leave some breathing room 
        if(bary[ie] >= Constants::baryTol){
          // Next, check if the other two are positive; if so, then nothing to do.
          int ie1 = (ie  + 1)%(tdim + 1);
          int ie2 = (ie1 + 1)%(tdim + 1);
          if(bary[ie1] < Constants::baryTol || bary[ie2] < Constants::baryTol){
            // "truncate" barycentrics 
            if(bary[ie1] < Constants::baryTol){
              bary[ie1] = Constants::baryTol;
            }
            if(bary[ie2] < Constants::baryTol){
              bary[ie2] = Constants::baryTol;
            }
            double sum = 0;
            for(int ii = 0; ii < tdim + 1; ii++) sum += bary[ii];
            if(sum <= 1.0e-16){
              METRIS_THROW_MSG(GeomExcept(), "Zero or negative bary sum");
            }
            for(int ii = 0; ii < tdim + 1; ii++) bary[ii] /= sum;

            double eval[gdim];
            eval2<gdim,1>(msh.coord, 
                          msh.fac2poi[ientt],
                          msh.getBasis(), DifVar::None, DifVar::None,
                          bary, eval, 
                          NULL, NULL);
            for(int ii = 0; ii < gdim; ii++) msh.coord(ipoih,ii) = eval[ii];
          }
        }else{
          // Check if the neighbour is limiting 
          int ient2 = msh.fac2fac(ientt,ie);
          if(ient2 < 0){
            METRIS_THROW_MSG(TODOExcept(),"Boundary or non-manifold case in msh_curve_offsets");
          }else{
            inventP1<gdim>(msh.fac2poi[ient2],msh.coord,msh.coord[ipoih],bary);
            int je = -1;
            for(int ii = 0; ii < tdim + 1; ii++){
              if(msh.fac2fac(ient2,ii) == ientt){
                je = ii;
                break;
              }
            }
            METRIS_ASSERT(je >= 0);
            // Set + tol as we can leave some breathing room 
            if(bary[ie] >= Constants::baryTol){
              // Next, check if the other two are positive; if so, then nothing to do.
              int ie1 = (je  + 1)%(tdim + 1);
              int ie2 = (ie1 + 1)%(tdim + 1);
              if(bary[ie1] < Constants::baryTol || bary[ie2] < Constants::baryTol){
                // "truncate" barycentrics 
                if(bary[ie1] < Constants::baryTol){
                  bary[ie1] = Constants::baryTol;
                }
                if(bary[ie2] < Constants::baryTol){
                  bary[ie2] = Constants::baryTol;
                }
                double sum = 0;
                for(int ii = 0; ii < tdim + 1; ii++) sum += bary[ii];
                if(sum <= 1.0e-16){
                  METRIS_THROW_MSG(GeomExcept(), "Zero or negative bary sum");
                }
                for(int ii = 0; ii < tdim + 1; ii++) bary[ii] /= sum;

                double eval[gdim];
                eval2<gdim,1>(msh.coord, 
                              msh.fac2poi[ient2],
                              msh.getBasis(), DifVar::None, DifVar::None,
                              bary, eval, 
                              NULL, NULL);
                for(int ii = 0; ii < gdim; ii++) msh.coord(ipoih,ii) = eval[ii];
              }
            }
          }
        }



      }else if (tdim == 3){
        METRIS_THROW_MSG(TODOExcept(),"Shell2 and correction in msh_curve_offsets");
      }

    }
  }


  if(iverb >= 2){
    msh.setBasis(FEBasis::Lagrange);
    writeMesh("crv1.meshb",msh);
  }
  double t1 = get_wall_time();
  printf("Curving time %f \n",t1-t0);

  maximizeCcoef<ideg,2,2>(msh, OptDoF::HO, LPMethod::IPM, LPLib::alglib);

  #else
    double t1 = get_wall_time();
    printf("Curving time %f \n",t1-t0);

    dblAr2 pos_ctrlp(npopt,gdim);
    pos_ctrlp.set_n(npopt);
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      int irank = idx_point[ipoin];
      if(irank < 0) continue;
      for(int ii = 0; ii < gdim; ii++) pos_ctrlp(irank,ii) = msh.coord(ipoin,ii);
      for(int ii = 0; ii < gdim; ii++) msh.coord(ipoin,ii) = coor0(ipoin,ii);
    }
    writeMesh("crv2.meshb",msh);

    bezGapsLP<gdim,tdim,ideg>(msh,idx_point,pos_ctrlp,LPMethod::IPM, LPLib::alglib);

  #endif

  double t2 = get_wall_time();
  printf("Correction time %f \n",t2-t1);

  msh.setBasis(ibas0);
  msh.met.setSpace(ispac0);

	return 0;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template int curveMeshOffsets<MetricFieldAnalytical,2,n>(\
       Mesh<MetricFieldAnalytical> &msh, bool icorr);\
template int curveMeshOffsets<MetricFieldAnalytical,3,n>(\
       Mesh<MetricFieldAnalytical> &msh, bool icorr);\
template int curveMeshOffsets<MetricFieldFE        ,2,n>(\
       Mesh<MetricFieldFE        > &msh, bool icorr);\
template int curveMeshOffsets<MetricFieldFE        ,3,n>(\
       Mesh<MetricFieldFE        > &msh, bool icorr);
#define BOOST_PP_LOCAL_LIMITS     (2, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


//#include <src/quality/opt_metqua.ixx>

} // End namespace
