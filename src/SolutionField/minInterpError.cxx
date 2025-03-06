//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "SolutionField.hxx"
#include "minInterpError.hxx"
#include "interpError.hxx"

#include "../Mesh/Mesh.hxx"
#include "../ho_constants.hxx"
#include "../utils/aux_timer.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/aux_misc.hxx"
#include "../utils/bernstein_prod.hxx"
#include "../low_ccoef.hxx"
#include "../linalg/det.hxx"
#include "../opt_generic.hxx"
#include "../low_topo.hxx"
#include "../low_geo.hxx"

#include "codegen_lag2bez.hxx"

namespace Metris{

template <class MFT>
void minimizeInterpErrglo(Mesh<MFT> &msh, const SolutionFieldAnalytical &sol, 
                          int pdeg, int pnorm, int ithrd1, int ithrd2){
  METRIS_ASSERT(msh.idim == 2);
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
    CT_FOR0_INC(2,2,idim){if(idim == msh.idim){
      METRIS_ASSERT(idim == msh.get_tdim());
      printf("Debug pdeg = %d \n",pdeg);
      CT_FOR0_INC(1,METRIS_MAX_DEG,pdeg_){if(pdeg_ == pdeg){
        minimizeInterpErrglo0<MFT,idim,ideg,pdeg_>(msh,sol,pnorm,ithrd1,ithrd2);
      }}CT_FOR1(pdeg_);
    }}CT_FOR1(idim);
  }}CT_FOR1(ideg);
}


template
void minimizeInterpErrglo<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, 
                          const SolutionFieldAnalytical &sol, 
                          int pdeg, int pnorm, int ithrd1, int ithrd2);
template
void minimizeInterpErrglo<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, 
                          const SolutionFieldAnalytical &sol, 
                          int pdeg, int pnorm, int ithrd1, int ithrd2);



template<class MFT, int idim, int ideg, int pdeg>
void minimizeInterpErrglo0(Mesh<MFT> &msh, const SolutionFieldAnalytical &sol, 
                           int pnorm, int ithrd1, int ithrd2){
  GETVDEPTH(msh);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(pnorm == 1 || pnorm == 2);

  FEBasis ibas0 = msh.getBasis();
  msh.setBasis(FEBasis::Bezier);

  // parameters
  const int miter = 10;

  // boilerplate
  constexpr int tdim = idim;
  constexpr int gdim = idim;
  constexpr int nnode = getnnode(idim,ideg);
  constexpr int nedgl = (gdim*(gdim+1))/2;

  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const intAr2 &ent2ent = msh.ent2ent(tdim);
  const int nentt = msh.nentt(tdim);

  const int inod0 = tdim + 1; // only HO

  double t0 = get_wall_time();
  CPRINTF1("-- START minimizeInterpErrglo\n");

  // start
  intAr1 lball(10);
  intAr1 lnode(10); 
  double errGlo0, errGlo1;
  bool getErrGlo = DOPRINTS1();
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh);
    int nerro = 0, nsucc = 0;

    if(DOPRINTS1() || getErrGlo){
      double errGlo = pnorm == 1 ? interpErrGlo<idim,ideg,pdeg,1,true>(sol)
                                 : interpErrGlo<idim,ideg,pdeg,2,true>(sol);
      if(niter == 0) errGlo0 = errGlo;
      CPRINTF1("-- START iter %d/%d error %e ",niter,miter,errGlo);
      if(niter == 0 && DOPRINTS1()) printf("\n");
      else if(DOPRINTS1()) printf(" red %% %f\n",100*(errGlo0-errGlo)/errGlo0);
    }

    msh.tag[ithrd1]++;
    for(int ientt = 0; ientt < nentt; ientt++){
      if(isdeadent(ientt, ent2poi)) continue;

      INCVDEPTH(msh);
      for(int inode = inod0; inode < nnode; inode++){
        int ipoin = ent2poi(ientt,inode);
        if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];

        // skip bdry points
        if(msh.poi2bpo[ipoin] >= 0) continue;

        // If point is a node, use ball. Otherwise if edge and 3D, shell. Otherwise neighbour.
        if(inode < tdim + 1){
          if(tdim == 2){
            intAr1 dum;
            int iopen; 
            bool imani;
            ball2(msh, ipoin, ientt, lball, dum, &iopen, &imani, ithrd2);
            lnode.set_n(0);
            for(int ient2 : lball){
              int inod2 = msh.template getverfac<ideg>(ient2, ipoin);
              lnode.stack(inod2);
            }
            METRIS_ASSERT(imani && !iopen);
          }else{
            METRIS_THROW_MSG(TODOExcept(), "3D in minInterperr");
          }
        }else if(inode < tdim + 1 + nedgl * (ideg - 1)){
          // Edge node. 
          if(tdim == 2){
            lball.set_n(0);
            lnode.set_n(0);
            lball.stack(ientt);
            lnode.stack(inode);
            // Silencing a compiler warning, but we are never here if ideg == 1
            int iedgl = ideg == 1 ? -1 : inode - (tdim + 1) / (ideg - 1);
            int ient2 = ent2ent(ientt,iedgl);
            int inod2 = msh.template getverfac<ideg>(ient2, ipoin);
            METRIS_ASSERT(inod2 >= 0);
            lball.stack(ient2);
            lnode.stack(inod2);
          }else{
            METRIS_THROW_MSG(TODOExcept(), "3D in minInterperr");
          }
        }else{
          METRIS_THROW_MSG(TODOExcept(), "Interior 2D/3D or face in 3D");
        }

        if(DOPRINTS2()){
          CPRINTF2(" - ientt %d inode %d ipoin %d ball size %d\n",ientt,inode,ipoin,lball.get_n());
          CPRINTF2(" ball:");
          lball.print();
          CPRINTF2(" lnode:");
          lnode.print();
        }

        // We now have lball list of elements, lnode index of ipoin in these
        double errLp0, errLp1;
        int ierro = 
        minimizeInterpErrloc<MFT,idim,ideg,pdeg>(msh, sol, pnorm, lball, lnode, &errLp0, &errLp1);
        METRIS_ASSERT(errLp1 <= errLp0 || ierro);
        if(ierro != 0){
          CPRINTF1(" - returned ierro %d \n",ierro);
          nerro++;
          continue;
        }else{
          CPRINTF1(" - min interp err loc success %e -> %e, decrease %f%% \n",
                   errLp0, errLp1, (errLp1-errLp0)/errLp0*100);
          nsucc++;
        }

      }// for inode
    }// for ientt 
    CPRINTF1(" - END iter %d/%d nerro %d nsucc %d \n",niter,miter,nerro,nsucc);
  }// for niter

  if(DOPRINTS1()){
    double errGlo = pnorm == 1 ? interpErrGlo<idim,ideg,pdeg,1,true>(sol)
                               : interpErrGlo<idim,ideg,pdeg,2,true>(sol);
    errGlo1 = errGlo;
    double red = (errGlo0 - errGlo1) / errGlo0 * 100;
    double t1 = get_wall_time();
    CPRINTF1("-- END minimizeInterpErrglo error %e -> %e reduced%% %f \n",errGlo0,errGlo1,red);
    CPRINTF1("-- time = %f = %d elt/s\n",t1-t0,(int)(nentt/(t1-t0)));
  }

  msh.setBasis(ibas0);
}

#define BOOST_PP_LOCAL_MACRO(n)\
template \
void minimizeInterpErrglo0<MetricFieldAnalytical, 2, n, 1>(Mesh<MetricFieldAnalytical> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, int ithrd1, int ithrd2);\
template \
void minimizeInterpErrglo0<MetricFieldFE, 2, n, 1>(Mesh<MetricFieldFE> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, int ithrd1, int ithrd2);\
template \
void minimizeInterpErrglo0<MetricFieldAnalytical, 2, n, 2>(Mesh<MetricFieldAnalytical> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, int ithrd1, int ithrd2);\
template \
void minimizeInterpErrglo0<MetricFieldFE, 2, n, 2>(Mesh<MetricFieldFE> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

#if 0
template<int idim, int ideg>
double interpErrglo_constraint_nlopt(unsigned int ncstr, double *result, 
          unsigned int nvar, const double *xvar, double *grad, void *data)
{
  interpErrglo_constraint_data *d = (interpErrglo_constraint_data *) data;
  MeshBase<MFT> &msh = *(data->msh);
  const intAr1 &ldof = *(data->ldof);
  METRIS_ASSERT(nvar == ldof.get_n()*msh.idim);

  int idof = 0;
  for(int ipoin : ldof){
    for(int ii = 0; ii < msh.idim; ii++){
      msh.coord(ipoin,ii) = xvar[msh.idim*idof + ii];
    }
    idof++;
  }
  
  const int nentt = msh.nentt(idim);
  constexpr int jdeg = idim*(ideg - 1);
  constexpr int ncoef = getnnode(idim,jdeg);
  METRIS_ASSERT(ncstr == nentt*ncoef);

  for(int ientt = 0; ientt < nentt; ientt++){
    bool iinva;
    getsclccoef<idim,idim,ideg>(msh,ientt,NULL,&result[ncoef*ientt],&iinva);

    if(grad == NULL) continue;

    getccoef_dpoint<idim, ideg>(msh, ientt, inode, NULL, d_ccoef)

  }

  double a = d->a, b = d->b;
  if (grad) {
      grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
      grad[1] = -1.0;
  }
  return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}


template<class MFT, int idim, int ideg, int pdeg>
void minimizeInterpErrglo0_nlopt(Mesh<MFT> &msh, const SolutionFieldAnalytical &sol, 
                           int pnorm, int ithrd1, int ithrd2){
  GETVDEPTH(msh);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(pnorm == 1 || pnorm == 2);

  FEBasis ibas0 = msh.getBasis();
  msh.setBasis(FEBasis::Bezier);

  // parameters
  const int miter = 10;

  // boilerplate
  constexpr int tdim = idim;
  constexpr int gdim = idim;
  constexpr int nnode = getnnode(idim,ideg);
  constexpr int nedgl = (gdim*(gdim+1))/2;

  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const intAr2 &ent2ent = msh.ent2ent(tdim);
  const int nentt = msh.nentt(tdim);

  const int inod0 = tdim + 1; // only HO

  double t0 = get_wall_time();
  CPRINTF1("-- START minimizeInterpErrglo\n");

  // start
  intAr1 lball(10);
  intAr1 lnode(10); 
  double errGlo0, errGlo1;
  bool getErrGlo = DOPRINTS1();
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh);
    int nerro = 0, nsucc = 0;

    if(DOPRINTS1() || getErrGlo){
      double errGlo = pnorm == 1 ? interpErrGlo<idim,ideg,pdeg,1,true>(sol)
                                 : interpErrGlo<idim,ideg,pdeg,2,true>(sol);
      if(niter == 0) errGlo0 = errGlo;
      CPRINTF1("-- START iter %d/%d error %e ",niter,miter,errGlo);
      if(niter == 0 && DOPRINTS1()) printf("\n");
      else if(DOPRINTS1()) printf(" red %% %f\n",100*(errGlo0-errGlo)/errGlo0);
    }

    msh.tag[ithrd1]++;
    for(int ientt = 0; ientt < nentt; ientt++){
      if(isdeadent(ientt, ent2poi)) continue;

      INCVDEPTH(msh);
      for(int inode = inod0; inode < nnode; inode++){
        int ipoin = ent2poi(ientt,inode);
        if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];

        // skip bdry points
        if(msh.poi2bpo[ipoin] >= 0) continue;

        // If point is a node, use ball. Otherwise if edge and 3D, shell. Otherwise neighbour.
        if(inode < tdim + 1){
          if(tdim == 2){
            intAr1 dum;
            int iopen; 
            bool imani;
            ball2(msh, ipoin, ientt, lball, dum, &iopen, &imani, ithrd2);
            lnode.set_n(0);
            for(int ient2 : lball){
              int inod2 = msh.template getverfac<ideg>(ient2, ipoin);
              lnode.stack(inod2);
            }
            METRIS_ASSERT(imani && !iopen);
          }else{
            METRIS_THROW_MSG(TODOExcept(), "3D in minInterperr");
          }
        }else if(inode < tdim + 1 + nedgl * (ideg - 1)){
          // Edge node. 
          if(tdim == 2){
            lball.set_n(0);
            lnode.set_n(0);
            lball.stack(ientt);
            lnode.stack(inode);
            // Silencing a compiler warning, but we are never here if ideg == 1
            int iedgl = ideg == 1 ? -1 : inode - (tdim + 1) / (ideg - 1);
            int ient2 = ent2ent(ientt,iedgl);
            int inod2 = msh.template getverfac<ideg>(ient2, ipoin);
            METRIS_ASSERT(inod2 >= 0);
            lball.stack(ient2);
            lnode.stack(inod2);
          }else{
            METRIS_THROW_MSG(TODOExcept(), "3D in minInterperr");
          }
        }else{
          METRIS_THROW_MSG(TODOExcept(), "Interior 2D/3D or face in 3D");
        }

        if(DOPRINTS2()){
          CPRINTF2(" - ientt %d inode %d ipoin %d ball size %d\n",ientt,inode,ipoin,lball.get_n());
          CPRINTF2(" ball:");
          lball.print();
          CPRINTF2(" lnode:");
          lnode.print();
        }

        // We now have lball list of elements, lnode index of ipoin in these
        double errLp0, errLp1;
        int ierro = 
        minimizeInterpErrloc<MFT,idim,ideg,pdeg>(msh, sol, pnorm, lball, lnode, &errLp0, &errLp1);
        METRIS_ASSERT(errLp1 <= errLp0 || ierro);
        if(ierro != 0){
          CPRINTF1(" - returned ierro %d \n",ierro);
          nerro++;
          continue;
        }else{
          CPRINTF1(" - min interp err loc success %e -> %e, decrease %f%% \n",
                   errLp0, errLp1, (errLp1-errLp0)/errLp0*100);
          nsucc++;
        }

      }// for inode
    }// for ientt 
    CPRINTF1(" - END iter %d/%d nerro %d nsucc %d \n",niter,miter,nerro,nsucc);
  }// for niter

  if(DOPRINTS1()){
    double errGlo = pnorm == 1 ? interpErrGlo<idim,ideg,pdeg,1,true>(sol)
                               : interpErrGlo<idim,ideg,pdeg,2,true>(sol);
    errGlo1 = errGlo;
    double red = (errGlo0 - errGlo1) / errGlo0 * 100;
    double t1 = get_wall_time();
    CPRINTF1("-- END minimizeInterpErrglo error %e -> %e reduced%% %f \n",errGlo0,errGlo1,red);
    CPRINTF1("-- time = %f = %d elt/s\n",t1-t0,(int)(nentt/(t1-t0)));
  }

  msh.setBasis(ibas0);
}

#define BOOST_PP_LOCAL_MACRO(n)\
template \
void minimizeInterpErrglo0_nlopt<MetricFieldAnalytical, 2, n, 1>(Mesh<MetricFieldAnalytical> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, int ithrd1, int ithrd2);\
template \
void minimizeInterpErrglo0_nlopt<MetricFieldFE, 2, n, 1>(Mesh<MetricFieldFE> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, int ithrd1, int ithrd2);\
template \
void minimizeInterpErrglo0_nlopt<MetricFieldAnalytical, 2, n, 2>(Mesh<MetricFieldAnalytical> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, int ithrd1, int ithrd2);\
template \
void minimizeInterpErrglo0_nlopt<MetricFieldFE, 2, n, 2>(Mesh<MetricFieldFE> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

#endif






// Not thread safe
template<class MFT, int idim, int ideg, int pdeg>
int minimizeInterpErrloc(Mesh<MFT> &msh, const SolutionFieldAnalytical &sol, 
                         int pnorm, intAr1& lball, intAr1& lnode,
                         double *errLp0, double *errLp1){
  GETVDEPTH(msh);

  //const double qfac = 0.5; // How much we keep of the optimum
  //double qfac1 = qfac; // work

  const int miter = 100;

  static int nwarnprt = 0;
  constexpr int iexact = true;
  if(nwarnprt++ < 10 && iexact) printf("## WARNING iexact = true is expensive\n");

  constexpr int nhess = (idim*(idim+1))/2;
  static newton_drivertype_args<idim> args;
  args.maxit = miter;
  args.stpmin = 1.0e-8;
  //if(DOPRINTS2()) args.iprt = 4;
  int iflag = 0, ihess;
  double xcur[idim], fcur, gcur[idim], hess[nhess];
  for(int ii = 0; ii < idim; ii++) xcur[ii] = 0;
  const int nball = lball.get_n();

  const intAr2& ent2poi = msh.ent2poi(idim);
  int ipoin = ent2poi(lball[0], lnode[0]);
  METRIS_ASSERT(ipoin >= 0);

  constexpr int jdeg = idim*(ideg - 1);
  constexpr int ncoef = getnnode(idim,jdeg);
  bool iinva = false;
  double ccoef[ncoef];

  int ierro = 0;

  double coor0[idim];
  for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);

  for(int niter = 0; niter < miter; niter++){

    ierro = optim_newton_drivertype<idim>(args, xcur, &fcur, gcur, hess, &iflag, &ihess);
    
    CPRINTF1(" - newton ret ierro %d iflag %d xcur %f %f \n",ierro,iflag,xcur[0],xcur[1]);

    if(ierro > 0){
      CPRINTF1(" ## optim_newton_drivertype error %d\n",ierro);
      goto cleanup;
    }
    if(iflag <= 0){
      CPRINTF1(" - iflag = 0 termination\n");
      break;
    }

    fcur = 0;
    double dloc[idim], hloc[nhess];
    for(int ii = 0; ii < idim; ii++) 
      msh.coord(ipoin,ii) = coor0[ii] + xcur[ii];


    iinva = false;
    if constexpr (ideg == 1){
      for(int ientt : lball){
        getmeasentP1<idim,idim>(msh, ent2poi[ientt], NULL, &iinva);
        if(iinva) break;
      }
    }else{
      for(int ientt : lball){
        getsclccoef<idim,idim,ideg>(msh,ientt,NULL,ccoef,&iinva); 
        if(iinva) break;
      }
    }
    if(iinva){
      CPRINTF1(" - invalid config, ihess = %d \n",ihess);
      // Otherwise we got an invalid update. 
      METRIS_ASSERT(ihess == 0);
      if(ihess != 0){
        ierro = 2;
        goto cleanup;
      }
      fcur = 1.0e10;
      for(int ii = 0; ii < idim; ii++){
        gcur[ii] = xcur[ii] - args.xopt[ii];
      }
      continue;
    }


    for(int iball = 0; iball < nball; iball++){
      int ientt = lball[iball];
      int inode = lnode[iball];
      METRIS_ASSERT(ent2poi(ientt,inode) == ipoin);
      if(pnorm == 1){
        if(ihess == 0){
          fcur += interpErr<idim, ideg, pdeg, 1, iexact>(sol, ientt, inode, {dloc});
        }else{
          fcur += interpErr<idim, ideg, pdeg, 1, iexact>(sol, ientt, inode, {dloc, hloc});
        }
      }else{
        if(ihess == 0){
          fcur += interpErr<idim, ideg, pdeg, 2, iexact>(sol, ientt, inode, {dloc});
        }else{
          fcur += interpErr<idim, ideg, pdeg, 2, iexact>(sol, ientt, inode, {dloc, hloc});
        }
      }// if pnorm

      for(int ii = 0; ii < idim; ii++){
        gcur[ii] += dloc[ii];
      }

      if(ihess == 0) continue;

      for(int ii = 0; ii < nhess; ii++){
        hess[ii] += hloc[ii];
      }
    }// for iball

    if(niter == 0) *errLp0 = fcur;

    CPRINTF1(" - computed errLp ball = %e \n",fcur);
  }// for niter

  *errLp1 = args.fopt;
  for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii] + args.xopt[ii];
  iinva = false;
  #ifndef NDEBUG
  if constexpr (ideg == 1){
    for(int ientt : lball){
      getmeasentP1<idim,idim>(msh, ent2poi[ientt], NULL, &iinva);
      if(iinva) break;
    }
  }else{
    for(int ientt : lball){
      getsclccoef<idim,idim,ideg>(msh,ientt,NULL,ccoef,&iinva); 
      if(iinva) break;
    }
  }
  METRIS_ASSERT(!iinva);
  #endif

  //// "Back"track
  //for(int niter = 0; niter < 10; niter++){
  //  for(int ii = 0; ii < idim; ii++) 
  //    msh.coord(ipoin,ii) = (1-qfac)*coor0[ii] + qfac*args.xopt[ii];

  //  if constexpr (ideg == 1){
  //    for(int ientt : lball){
  //      getmeasentP1<idim,idim>(msh, ent2poi[ientt], NULL, &iinva);
  //      if(iinva) break;
  //    }
  //  }else{
  //    for(int ientt : lball){
  //      getsclccoef<idim,idim,ideg>(msh,ientt,NULL,ccoef,&iinva); 
  //      if(iinva) break;
  //    }
  //  }
  //  if(!iinva) break;

  //  qfac1 *= 1.1;
  //}

  //if(iinva){
  //  CPRINTF1("-- Gave up on relaxation")
  //  for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii] + args.xopt[ii];
  //}else{
  //  CPRINTF1("-- Final relaxation factor %f \n",qfac1)
  //}



  return 0;

  cleanup:
  for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];

  return ierro;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template \
int minimizeInterpErrloc<MetricFieldAnalytical, 2, n, 1>(Mesh<MetricFieldAnalytical> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, intAr1& lball, intAr1& lnode,\
                           double *errLp0, double *errLp1);\
template \
int minimizeInterpErrloc<MetricFieldFE, 2, n, 1>(Mesh<MetricFieldFE> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, intAr1& lball, intAr1& lnode,\
                           double *errLp0, double *errLp1);\
template \
int minimizeInterpErrloc<MetricFieldAnalytical, 2, n, 2>(Mesh<MetricFieldAnalytical> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, intAr1& lball, intAr1& lnode,\
                           double *errLp0, double *errLp1);\
template \
int minimizeInterpErrloc<MetricFieldFE, 2, n, 2>(Mesh<MetricFieldFE> &msh, \
                           const SolutionFieldAnalytical &sol, \
                           int pnorm, intAr1& lball, intAr1& lnode,\
                           double *errLp0, double *errLp1);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



} //namespace
