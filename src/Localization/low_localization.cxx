//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "low_localization.hxx"

//#include <nlopt.h>

#include "../low_geo.hxx"
#include "../opt_generic.hxx"
#include "../ho_constants.hxx"
#include "../low_eval.hxx"

#include "../linalg/eigen.hxx"
#include "../linalg/matprods.hxx"
#include "../linalg/invmat.hxx"
#include "../linalg/sym3idx.hxx"
#include "../Mesh/MeshBase.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../aux_pp_inc.hxx"


#include <random>

namespace Metris{

// Return 0 if converged and in element
//        1 if converged but outside
//        2 if unconverged
template<int gdim, int ideg> 
int inveval(const MeshBase &msh,
            const int* ent2pol,
            const dblAr2 &coord,
            const double* coor0, 
            double* __restrict__ coopr, double* __restrict__ bary,
            double tol){

  static_assert(gdim == 1 || gdim == 2 || gdim == 3);
  constexpr int tdim = gdim;
  constexpr auto evalf = tdim == 1 ? eval1<gdim,ideg> : 
                         tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;
  const int iverb = msh.param->iverb;

  inventP1<gdim>(ent2pol, coord, coor0, bary);

  if constexpr(ideg == 1){

    evalf(coord,ent2pol,msh.getBasis(),DifVar::None,
          DifVar::None,bary,coopr,NULL,NULL);

    goto checkbar;

  }else{

    newton_drivertype_args<gdim> args;
    args.stpmin = 1.0e-12;
    args.ratnew = 0.5; // LS step decrease factor 
    args.iprt = msh.param->iverb - 1;
    args.maxit = 500;
    int iflag = 0, ihess;
    constexpr int nhess = (gdim*(gdim + 1))/2;
    double xcur[tdim+1] = {0}, fcur, gcur[gdim], hcur[nhess];
    double jmat[tdim][gdim],hmat[nhess][gdim];
    double bar0[tdim+1];



    //restart:

    for(int ii = 0; ii < tdim + 1; ii++) bar0[ii] = bary[ii];
    if(iverb >= 3) printf("    - init bary0 = ");
    if(iverb >= 3) dblAr1(tdim+1,bar0).print();

    double dist;
    DifVar d2flag;
    while(true){
      int ierro = optim_newton_drivertype<gdim>(args, &xcur[1], &fcur, gcur, hcur, &iflag, &ihess);

      if(ierro > 0){
        if(iverb >= 1) printf("    ## optim_newton_drivertype error %d\n",ierro);
        return 1 + ierro;
      }
      if(iflag <= 0) {
        if(iverb >= 3) printf("     - iflag = 0 termination\n");
        break;
      }

      xcur[0] = 0;
      for(int ii = 1; ii < tdim + 1; ii++) xcur[0] -= xcur[ii];

      // Turns out not necessary
      #if 0
      // Renormalize xcur in order not to exit the element
      double fac = 1;
      bool dofac = false;
      for(int ii = 0; ii < tdim + 1; ii++){
        if(abs(xcur[ii]) < 1.0e-12) continue;

        if(bar0[ii] + xcur[ii] > 1.0 - Constants::baryTol){
          dofac = true;
          fac = MIN(fac, abs( (1 - bar0[ii]) / xcur[ii] ) );
          if(iverb >= 3) printf("    - debug tmp 1 = %f before, after %f %f \n", 
              abs( (1 - bar0[ii]) / xcur[ii] ) , bar0[ii], bar0[ii] + xcur[ii]);
        }

        if(bar0[ii] + xcur[ii] < Constants::baryTol){
          dofac = true;
          fac = MIN(fac, abs( -bar0[ii] / xcur[ii] ) );
          if(iverb >= 3) printf("    - debug tmp 2 = %f before, after %f %f w cor %f \n", 
              abs( -bar0[ii] / xcur[ii] ) , bar0[ii], bar0[ii] + xcur[ii]
              , bar0[ii] + abs( -bar0[ii] / xcur[ii] )*xcur[ii]);
        }

        //if(iverb >= 3) printf("    - debug tmp 2 = %f before, after %f %f \n", tmp, bar0[ii], bar0[ii] + xcur[ii]);
      }

      if(dofac){
        if(iverb >= 3) printf("  - applying factor %f to descent ",fac);
        if(iverb >= 3) dblAr1(tdim+1,xcur).print();
        for(int ii = 0; ii < tdim + 1; ii++) xcur[ii] *= fac;
      }
      #endif

      for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = bar0[ii] + xcur[ii];


      d2flag = ihess > 0 ? DifVar::Bary : DifVar::None;
      evalf(coord,ent2pol,msh.getBasis(),DifVar::Bary,
            d2flag,bary,coopr,jmat[0],hmat[0]);
      dist = geterrl2<gdim>(coopr,coor0);


      if(iverb >= 3){
        printf("    - niter %d dist = %15.7e  dbar = ",args.niter,dist);
        dblAr1(tdim+1, xcur).print();
        printf("    - bary = ");
        dblAr1(tdim+1,bary).print();
      }

      if(dist < tol*tol){
        if(iverb >= 3) printf("    - dist < tol termination\n");
        goto checkbar;
      }

      // f = ||FK(xi) - X0||^2 / 2
      // dif = (J_{i:}, (FK(xi) - X0))
      // dijf = (H_{ij:}, (FK(xi) - X0)) + (J_{i:}, J_{j:})
      fcur = dist / 2;
      for(int ii = 0; ii < gdim; ii++) 
        gcur[ii] = getprdl2<gdim>(jmat[ii], coopr)  
                 - getprdl2<gdim>(jmat[ii], coor0);
      if(ihess > 0){
        for(int ii = 0; ii < gdim; ii++){
          for(int jj = ii; jj < gdim; jj++){
            hcur[sym3idx(ii,jj)] = getprdl2<gdim>(hmat[sym3idx(ii,jj)], coopr)
                                 - getprdl2<gdim>(hmat[sym3idx(ii,jj)], coor0)
                                 + getprdl2<gdim>(jmat[ii],jmat[jj]);
          }
        }
      }// if ihess

    }

    //// Never have we passed the small dist test. Try preconditionning and restart

    //constexpr auto evalp1 = tdim == 1 ? eval1<gdim,1> : 
    //                        tdim == 2 ? eval2<gdim,1> : eval3<gdim,1>;
    //evalf(coord,ent2pol,msh.getBasis(),DifVar::Bary,
    //      DifVar::None,bary,coopr,jmat[0],NULL);
    //// Apply mapping phi : X -> J^{-1} (X - X0)
    //// In particular sought point is set to 0

  }

  return 2;

  checkbar:

  bool iout = false;
  for(int ii = 0 ;ii < gdim+1; ii++){
    if(bary[ii] >   - Constants::baryTol 
    && bary[ii] < 1 + Constants::baryTol) continue;
    iout = true;
  }

  return iout;
}
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int inveval<1, n>(const MeshBase&, const int*, const dblAr2 &,const double*, \
                        double* __restrict__, double* __restrict__, double);\
template int inveval<2, n>(const MeshBase&, const int*, const dblAr2 &,const double*, \
                        double* __restrict__, double* __restrict__, double);\
template int inveval<3, n>(const MeshBase&, const int*, const dblAr2 &,const double*, \
                        double* __restrict__, double* __restrict__, double);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()





template<int gdim, int ideg> 
int inveval(MeshBase &msh, int ientt, const double* coor0, 
            double* __restrict__ coopr, double* __restrict__ bary,
            double tol0){
  constexpr int tdim = gdim;

  const int iverb = msh.param->iverb;

  double tol = tol0;
  if constexpr(ideg > 1){
    double eps = getepsent<gdim>(msh,gdim,ientt);
    tol *= eps;
  }

  if(iverb >= 3)
    printf("    -- Start inveval_linesearch ientt %d adj tol fac %15.7e\n",ientt,tol/tol0);
  

  const intAr2 &ent2poi = msh.ent2poi(tdim);
  int ierro = inveval<gdim,ideg>(msh,ent2poi[ientt],msh.coord,coor0,coopr,bary,tol);
  if(ierro < 2) return ierro;
  // Only if ideg > 1 should we be here 
  METRIS_ASSERT(ideg > 1);
  if(iverb >= 3) printf("    - inveval0 unconverged, precondition\n");


  double jmatP1[tdim][gdim];
  constexpr int nnode = tdim == 1 ? edgnpps[ideg] :
                        tdim == 2 ? facnpps[ideg] : tetnpps[ideg];
  int ent2pol[nnode];
  for(int ii = 0; ii < nnode; ii++) ent2pol[ii] = ii;

  msh.poi2rwk.allocate(gdim*nnode);
  dblAr2 coorl(nnode,gdim, &msh.poi2rwk[0]);
  coorl.set_n(nnode);


  constexpr auto evalp1 = tdim == 1 ? eval1<gdim,1> : 
                          tdim == 2 ? eval2<gdim,1> : eval3<gdim,1>;
  evalp1(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::Bary,
         DifVar::None,bary,coopr,jmatP1[0],NULL);
  invmat<gdim>(jmatP1[0]);
  double dum[gdim];
  for(int ii = 0; ii < nnode; ii++){
    int ipoin = ent2poi(ientt,ii);
    for(int jj = 0; jj < gdim; jj++) dum[jj] = msh.coord(ipoin,jj)
                                             - coor0[jj];
    tmatXvec<gdim>(jmatP1[0],dum,coorl[ii]);
  }

  for(int ii = 0; ii < gdim ;ii++) dum[ii] = 0;

  ierro = inveval<gdim,ideg>(msh,ent2pol,coorl,dum,coopr,bary,tol0);
  

  // Bary is correct but coopr needs reevaluating
  constexpr auto evalf = tdim == 1 ? eval1<gdim,ideg> : 
                         tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;
  evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
        DifVar::None,bary,coopr,NULL,NULL);

  double dist = geterrl2<gdim>(coopr,coor0);
  if(iverb >= 3) printf("     - precond call ret ierro %d dist = %15.7e\n",ierro,dist);

  if(dist >= tol*tol){
    if(iverb >= 3) printf("    ## still unconverged !!\n");
    return 2;
  }

  return ierro;
}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int inveval<1, n>(MeshBase&, int, const double*, double* __restrict__, double* __restrict__, double);\
template int inveval<2, n>(MeshBase&, int, const double*, double* __restrict__, double* __restrict__, double);\
template int inveval<3, n>(MeshBase&, int, const double*, double* __restrict__, double* __restrict__, double);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



template<int gdim>
int locMeshQuick(MeshBase &msh, const double *coor0){
	for(int ii = 0; ii < gdim; ++ii){
		if(coor0[ii] > msh.bb[ii][1] || coor0[ii] < msh.bb[ii][0]){
			return 1;
		}
	}
	return 0;
}

template int locMeshQuick<1>(MeshBase &msh, const double *coor0);
template int locMeshQuick<2>(MeshBase &msh, const double *coor0);
template int locMeshQuick<3>(MeshBase &msh, const double *coor0);


} // End namespace

