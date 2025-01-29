//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "low_localization.hxx"

//#include <nlopt.h>

#include "../low_geo.hxx"
#include "../aux_utils.hxx"
#include "../opt_generic.hxx"
#include "../ho_constants.hxx"
#include "../low_eval.hxx"

#include "../linalg/eigen.hxx"
#include "../linalg/matprods.hxx"
#include "../linalg/invmat.hxx"
#include "../linalg/symidx.hxx"
#include "../Mesh/MeshBase.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../aux_pp_inc.hxx"

#include <nlopt.hpp>

#include <random>

namespace Metris{

// Return 0 if converged and in element
//        1 if converged but outside
//        2 if unconverged
template<int gdim, int ideg> 
int inveval_badNewton0(const MeshBase &msh,
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
    args.wlfc1 = 1e-4;
    args.wlfc2 = 0.9;
    int iflag = 0, ihess;
    constexpr int nhess = (gdim*(gdim + 1))/2;
    //double dbar[tdim+1] = {0}, fcur, gcur[gdim], hcur[nhess];
    double fcur, gcur[gdim], hcur[nhess];
    double jmat[tdim][gdim],hmat[nhess][gdim];
    double bar0[tdim+1];

    double buf[5*gdim];
    dblAr1 buf_(5*gdim,buf);
    truncated_newton_work workTNCG(gdim,buf_);


    //restart:

    for(int ii = 0; ii < tdim + 1; ii++) bar0[ii] = bary[ii];
    if(iverb >= 3) printf("    - init bary0 = ");
    if(iverb >= 3) dblAr1(tdim+1,bar0).print();

    double dist;
    DifVar d2flag;
    // Function is ||F_K(xi^0 + dbar) - X||^2 / 2
    while(true){
      //int ierro = optim_newton_drivertype<gdim>(args, &bary[1], &fcur, gcur, hcur, &iflag, &ihess);
      
      args.wlfc1 = 0.4;
      int ierro = optim_newton_drivertype_TNCG<gdim>(args, workTNCG, &bary[1], &fcur, gcur, hcur, &iflag, &ihess);

      if(ierro > 0){
        if(iverb >= 3) printf("    ## optim_newton_drivertype error %d\n",ierro);
        return 1 + ierro;
      }
      if(iflag <= 0) {
        if(iverb >= 3) printf("     - iflag = 0 termination\n");
        break;
      }

      // To preserve bary sum, dbar must sum to 0
      //dbar[0] = 0;
      //for(int ii = 1; ii < tdim + 1; ii++) dbar[0] -= dbar[ii];
      bary[0] = 1;
      for(int ii = 1; ii < tdim + 1; ii++) bary[0] -= bary[ii];

      // Turns out not necessary
      #if 0
      // Renormalize dbar in order not to exit the element
      double fac = 1;
      bool dofac = false;
      for(int ii = 0; ii < tdim + 1; ii++){
        if(abs(dbar[ii]) < 1.0e-12) continue;

        if(bar0[ii] + dbar[ii] > 1.0 - Constants::baryTol){
          dofac = true;
          fac = MIN(fac, abs( (1 - bar0[ii]) / dbar[ii] ) );
          if(iverb >= 3) printf("    - debug tmp 1 = %f before, after %f %f \n", 
              abs( (1 - bar0[ii]) / dbar[ii] ) , bar0[ii], bar0[ii] + dbar[ii]);
        }

        if(bar0[ii] + dbar[ii] < Constants::baryTol){
          dofac = true;
          fac = MIN(fac, abs( -bar0[ii] / dbar[ii] ) );
          if(iverb >= 3) printf("    - debug tmp 2 = %f before, after %f %f w cor %f \n", 
              abs( -bar0[ii] / dbar[ii] ) , bar0[ii], bar0[ii] + dbar[ii]
              , bar0[ii] + abs( -bar0[ii] / dbar[ii] )*dbar[ii]);
        }

        //if(iverb >= 3) printf("    - debug tmp 2 = %f before, after %f %f \n", tmp, bar0[ii], bar0[ii] + dbar[ii]);
      }

      if(dofac){
        if(iverb >= 3) printf("  - applying factor %f to descent ",fac);
        if(iverb >= 3) dblAr1(tdim+1,dbar).print();
        for(int ii = 0; ii < tdim + 1; ii++) dbar[ii] *= fac;
      }
      #endif

      //for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = bar0[ii] + dbar[ii];


      d2flag = ihess > 0 ? DifVar::Bary : DifVar::None;
      evalf(coord,ent2pol,msh.getBasis(),DifVar::Bary,
            d2flag,bary,coopr,jmat[0],hmat[0]);
      dist = geterrl2<gdim>(coopr,coor0);


      if(iverb >= 3){
        //printf("    - niter %d dist = %15.7e  dbar = ",args.niter,dist);
        //dblAr1(tdim+1, bary).print();
        //printf("    - bary = ");
        //dblAr1(tdim+1,bary).print();
        printf("    - niter %d dist = %15.7e  bary = ",args.niter,dist);
        dblAr1(tdim+1, bary).print();
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
            hcur[sym2idx(ii,jj)] = getprdl2<gdim>(hmat[sym2idx(ii,jj)], coopr)
                                 - getprdl2<gdim>(hmat[sym2idx(ii,jj)], coor0)
                                 + getprdl2<gdim>(jmat[ii],jmat[jj]);
          }
        }
      }// if ihess

    }

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
template int inveval_badNewton0<1, n>(const MeshBase&, const int*, const dblAr2 &,const double*, \
                        double* __restrict__, double* __restrict__, double);\
template int inveval_badNewton0<2, n>(const MeshBase&, const int*, const dblAr2 &,const double*, \
                        double* __restrict__, double* __restrict__, double);\
template int inveval_badNewton0<3, n>(const MeshBase&, const int*, const dblAr2 &,const double*, \
                        double* __restrict__, double* __restrict__, double);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()







template<int gdim, int ideg> 
int inveval_badNewton(MeshBase &msh, int ientt, 
                      const double* coor0, 
                      double* __restrict__ coopr, 
                      double* __restrict__ bary,
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
  int ierro = inveval_badNewton0<gdim,ideg>(msh,ent2poi[ientt],msh.coord,coor0,coopr,bary,tol);
  if(ierro < 2) return ierro;
  // Only if ideg > 1 should we be here 
  METRIS_ASSERT(ideg > 1);
  if(iverb >= 3) printf("    - inveval0 unconverged, precondition\n");

  //printf("DEBUG WAIT HERE\n");
  //wait();


  double jmatP1[tdim][gdim];
  constexpr int nnode = tdim == 1 ? edgnpps[ideg] :
                        tdim == 2 ? facnpps[ideg] : tetnpps[ideg];
  int ent2pol[nnode];
  for(int ii = 0; ii < nnode; ii++) ent2pol[ii] = ii;

  // Just make a static buffer. This will probably blow up at large degrees
  if constexpr(ideg > 3) METRIS_THROW_MSG(TODOExcept(), "WATCH OUT FOR LARGE STATIC BUFFER");
  double buf[gdim*nnode];
  dblAr2 coorl(nnode,gdim,buf);
  coorl.set_n(nnode);


  constexpr auto evalp1 = tdim == 1 ? eval1<gdim,1> : 
                          tdim == 2 ? eval2<gdim,1> : eval3<gdim,1>;
  evalp1(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::Bary,
         DifVar::None,bary,coopr,jmatP1[0],NULL);
  ierro = invmat<gdim>(jmatP1[0]);
  if(ierro != 0) return 2;

  double dum[gdim];
  for(int ii = 0; ii < nnode; ii++){
    int ipoin = ent2poi(ientt,ii);
    for(int jj = 0; jj < gdim; jj++) dum[jj] = msh.coord(ipoin,jj)
                                             - coor0[jj];
    tmatXvec<gdim>(jmatP1[0],dum,coorl[ii]);
  }

  for(int ii = 0; ii < gdim ;ii++) dum[ii] = 0;

  ierro = inveval_badNewton0<gdim,ideg>(msh,ent2pol,coorl,dum,coopr,bary,tol0);

  // Bary is correct but coopr needs reevaluating
  constexpr auto evalf = tdim == 1 ? eval1<gdim,ideg> : 
                         tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;
  evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
        DifVar::None,bary,coopr,NULL,NULL);

  double dist = geterrl2<gdim>(coopr,coor0);
  if(iverb >= 3) printf("     - precond call ret ierro %d dist = %15.7e\n",ierro,dist);


  //printf("DEBUG WAIT HERE after precond\n");
  //wait();

  if(dist >= tol*tol){
    if(msh.param->dbgfull && iverb >= 3) printf("    ## still unconverged !!\n");
    return 2;
  }

  return ierro;
}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int inveval_badNewton<1, n>(MeshBase&, int, const double*, double* __restrict__, double* __restrict__, double);\
template int inveval_badNewton<2, n>(MeshBase&, int, const double*, double* __restrict__, double* __restrict__, double);\
template int inveval_badNewton<3, n>(MeshBase&, int, const double*, double* __restrict__, double* __restrict__, double);
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






// Return 0 if converged and in element
//        1 if converged but outside
//        2 if unconverged
template<int gdim, int ideg> 
int inveval0(const MeshBase &msh,
                   const int* ent2pol,
                   const dblAr2 &coord,
                   const double* coor0, 
                   double* __restrict__ coopr, 
                   double* __restrict__ bary,
                   dblAr1 &work,
                   double tol){
/*nlopt::LD_TNEWTON_PRECOND_RESTART,
  nlopt::LD_TNEWTON_RESTART,
  nlopt::LD_TNEWTON */
  // NLOPT_LD_TNEWTON_PRECOND_RESTART.
  // Second, simplified versions NLOPT_LD_TNEWTON_PRECOND (same without restarting), 
  // NLOPT_LD_TNEWTON_RESTART (same without preconditioning), 
  // and NLOPT_LD_TNEWTON (same without restarting or preconditioning).
  constexpr int tdim = gdim;

  nlopt_algorithm algo = NLOPT_LD_TNEWTON_PRECOND_RESTART;

  invevalfun_data mydata(msh,ent2pol,coord,coor0,coopr);
  double fopt; 
  int nwork = luksan_pnet_worksize(gdim);
  work.allocate(nwork);
  work.set_n(nwork);
  double fstop = tol*tol/2;
  double ftol_rel = -1e30; 
  double ftol_abs = -1e30;
  double lb[gdim], ub[gdim];

  inventP1<gdim>(ent2pol, coord, coor0, bary);

  for(int ii = 0; ii < gdim; ii++){
    lb[ii] = -HUGE_VAL;
    ub[ii] =  HUGE_VAL;
  }
  int ierro = luksan_pnetS<gdim>(invevalfun_nlointf<gdim,ideg>, &mydata,
                                 lb, ub, /* bounds */
                                 &bary[1], /* in: initial guess, out: minimizer */
                                 &fopt,
                                 //int mf, /* subspace dimension (0 for default) */
                                 algo,
                                 work,
                                 fstop , ftol_rel, ftol_abs);

  bary[0] = 1;
  for(int ii = 0; ii < gdim; ii++){
    bary[0]   -= bary[ii+1];
  }

  constexpr auto evalf = tdim == 1 ? eval1<gdim,ideg> : 
                         tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;
  evalf(coord,ent2pol,msh.getBasis(),DifVar::None,
        DifVar::None,bary,coopr,NULL,NULL);
  double dist = geterrl2<gdim>(coopr,coor0);
  // Why is this mixed in with unconverged conditions?
  if(ierro == NLOPT_STOPVAL_REACHED
  || ierro == NLOPT_FTOL_REACHED
  || ierro == NLOPT_XTOL_REACHED) ierro = NLOPT_SUCCESS;

  if(ierro != NLOPT_SUCCESS){
    if(msh.param->iverb >= 3) 
      printf("     -luksan_pnet failed ierro %d \n",ierro);
    return 2;
  }


  if(dist >= tol*tol){
    if(msh.param->iverb >= 3) 
      printf("     -luksan_pnet did not converge dist %e >= %e reported %e \n",dist,tol*tol,fopt*2);
    return 2;
  }

  //if(dist >= tol*tol){
  //  if(msh.param->iverb >= 3) printf("dist > tol despite  failed ierro %d \n",ierro);
  //}

  bool iout = false;
  for(int ii = 0 ;ii < gdim+1; ii++){
    if(bary[ii] >   - Constants::baryTol 
    && bary[ii] < 1 + Constants::baryTol) continue;
    iout = true;
  }

  return iout;

}
#define BOOST_PP_LOCAL_MACRO(n)\
template int inveval0<2,n>(const MeshBase &msh,\
                           const int* ent2pol,\
                           const dblAr2 &coord,\
                           const double* coor0,\
                           double* __restrict__ coopr,\
                           double* __restrict__ bary,\
                           dblAr1 &work,\
                           double tol);\
template int inveval0<3,n>(const MeshBase &msh,\
                           const int* ent2pol,\
                           const dblAr2 &coord,\
                           const double* coor0,\
                           double* __restrict__ coopr,\
                           double* __restrict__ bary,\
                           dblAr1 &work,\
                           double tol);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



template<int gdim, int ideg> 
int inveval(MeshBase &msh, int ientt, 
                  const double* coor0, 
                  double* __restrict__ coopr, 
                  double* __restrict__ bary,
                  double tol0){
  constexpr int tdim = gdim;

  const int iverb = msh.param->iverb;

  double tol = tol0;
  if constexpr(ideg > 1){
    double eps = getepsent<gdim>(msh,gdim,ientt);
    tol *= eps;
    // we want to bound error on bary |d\xi|
    // all we know is error on x |dx|
    // use |d\xi| ~= |dx| |d\xi|/|dx|
    // and |d\xi| / |dx| = Frob(J_K^{-1}) (or any matrix norm)
    // Hence estimate err on xi |d\xi| ~= |dx| eps^{-1} ?
    // cdtion becomes  |dx|  < tol * eps
    // w eps = sqrt(sum l_i^2)

    // 
  }

  if(iverb >= 3)
    printf("    -- Start inveval ientt %d adj tol fac %15.7e\n",ientt,tol/tol0);
  

  const intAr2 &ent2poi = msh.ent2poi(tdim);
  dblAr1 &work = msh.rwork;
  int ierro = inveval0<gdim,ideg>(msh,ent2poi[ientt],msh.coord,coor0,coopr,bary,work,tol);
  if(ierro < 2) return ierro;
  // Only if ideg > 1 should we be here 
  METRIS_ASSERT(ideg > 1);
  if(iverb >= 3) printf("    - inveval0 unconverged, precondition\n");

  //printf("DEBUG WAIT HERE\n");
  //wait();


  double jmatP1[tdim][gdim];
  constexpr int nnode = tdim == 1 ? edgnpps[ideg] :
                        tdim == 2 ? facnpps[ideg] : tetnpps[ideg];
  int ent2pol[nnode];
  for(int ii = 0; ii < nnode; ii++) ent2pol[ii] = ii;

  // Just make a static buffer. This will probably blow up at large degrees
  if constexpr(ideg > 3) METRIS_THROW_MSG(TODOExcept(), "WATCH OUT FOR LARGE STATIC BUFFER");
  double buf[gdim*nnode];
  dblAr2 coorl(nnode,gdim,buf);
  coorl.set_n(nnode);


  constexpr auto evalp1 = tdim == 1 ? eval1<gdim,1> : 
                          tdim == 2 ? eval2<gdim,1> : eval3<gdim,1>;
  evalp1(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::Bary,
         DifVar::None,bary,coopr,jmatP1[0],NULL);
  ierro = invmat<gdim>(jmatP1[0]);
  if(ierro != 0) return 2;
  double dum[gdim];
  for(int ii = 0; ii < nnode; ii++){
    int ipoin = ent2poi(ientt,ii);
    for(int jj = 0; jj < gdim; jj++) dum[jj] = msh.coord(ipoin,jj)
                                             - coor0[jj];
    tmatXvec<gdim>(jmatP1[0],dum,coorl[ii]);
  }

  for(int ii = 0; ii < gdim ;ii++) dum[ii] = 0;

  ierro = inveval0<gdim,ideg>(msh,ent2pol,coorl,dum,coopr,bary,work,tol0/10);

  // Bary is correct but coopr needs reevaluating
  constexpr auto evalf = tdim == 1 ? eval1<gdim,ideg> : 
                         tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;
  evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
        DifVar::None,bary,coopr,NULL,NULL);

  double dist = geterrl2<gdim>(coopr,coor0);
  if(iverb >= 3) printf("     - precond call ret ierro %d dist = %15.7e\n",ierro,dist);


  //printf("DEBUG WAIT HERE after precond\n");
  //wait();

  if(dist >= tol*tol){
    if(iverb >= 3) printf("    ## still unconverged dist %e >= %e !!\n",dist,tol*tol);
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





// For debugging the derivatives
template<int gdim, int ideg> 
double invevalfun(const MeshBase &msh,
                  const int* ent2pol, 
                  const dblAr2 &coord, // this can differ from msh if precond
                  const double* coor0, 
                  const double* bary ,
                  int ihess,
                  double *coopr,
                  double *grad,
                  double *hess){
  constexpr int tdim = gdim;
  constexpr int nhess = (gdim*(gdim + 1))/2;
  double jmat[tdim][gdim],hmat[nhess][gdim];
  constexpr auto evalf = tdim == 1 ? eval1<gdim,ideg> : 
                         tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;

  DifVar d2flag = ihess > 0 ? DifVar::Bary : DifVar::None;
  DifVar d1flag = grad == NULL ? DifVar::None : DifVar::Bary;
  evalf(coord,ent2pol,msh.getBasis(),d1flag,
        d2flag,bary,coopr,jmat[0],hmat[0]);
  double dist = geterrl2<gdim>(coopr,coor0);


  // f = ||FK(xi) - X0||^2 / 2
  // dif = (J_{i:}, (FK(xi) - X0))
  // dijf = (H_{ij:}, (FK(xi) - X0)) + (J_{i:}, J_{j:})
  double fcur = dist / 2;
  if(grad != NULL){
    for(int ii = 0; ii < gdim; ii++) 
      grad[ii] = getprdl2<gdim>(jmat[ii], coopr)  
               - getprdl2<gdim>(jmat[ii], coor0);
  }
  if(ihess > 0){
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        hess[sym2idx(ii,jj)] = getprdl2<gdim>(hmat[sym2idx(ii,jj)], coopr)
                             - getprdl2<gdim>(hmat[sym2idx(ii,jj)], coor0)
                             + getprdl2<gdim>(jmat[ii],jmat[jj]);
      }
    }
  }// if ihess

  return fcur;
}
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double invevalfun<2,n>(const MeshBase &msh,\
                                const int* ent2pol,\
                                const dblAr2 &coord,\
                                const double* coor0, \
                                const double* bary ,\
                                int ihess,\
                                double *coopr,\
                                double *grad,\
                                double *hess);\
template double invevalfun<3,n>(const MeshBase &msh,\
                                const int* ent2pol,\
                                const dblAr2 &coord,\
                                const double* coor0, \
                                const double* bary ,\
                                int ihess,\
                                double *coopr,\
                                double *grad,\
                                double *hess);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

template<int gdim,int ideg>
double invevalfun_nlointf(unsigned int nvar, const double *x, 
                          double *grad, void *f_data){
  // The Fortran intfc calling this offsets the pointer?
  invevalfun_data *mydata = (invevalfun_data*)(f_data);
  double buf[gdim];
  double bary[gdim+1];
  bary[0] = 1;
  for(int ii = 0; ii < gdim; ii++){
    bary[1+ii] = x[ii];
    bary[0]   -= x[ii];
  } 
  double val = invevalfun<gdim,ideg>(mydata->msh,mydata->ent2pol,mydata->coord,
                            mydata->coor0,bary, 0, mydata->coopr, buf, NULL);
  //if constexpr (gdim == 2) printf(" x = %f %f f = %f g %f %f req %d \n",
  //                                 x[0],x[1],val,buf[0],buf[1],grad!=NULL);
  if(grad != NULL){
    for(int ii = 0; ii < gdim; ii++) grad[ii] = buf[ii];
  }
  return val;
}
#define BOOST_PP_LOCAL_MACRO(n)\
template double invevalfun_nlointf<2,n>(unsigned int nvar, \
                          const double *x, double *grad, void *f_data);\
template double invevalfun_nlointf<3,n>(unsigned int nvar, \
                          const double *x, double *grad, void *f_data);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




} // End namespace

