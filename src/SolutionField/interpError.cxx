//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "SolutionField.hxx"
#include "interpError.hxx"

#include "../Mesh/Mesh.hxx"
#include "../ho_constants.hxx"
#include "../utils/aux_misc.hxx"
#include "../utils/bernstein_prod.hxx"
#include "../low_ccoef.hxx"
#include "../linalg/det.hxx"

#include "codegen_lag2bez.hxx"

namespace Metris{


//template<class MFT>
//void minimizeInterpError(Mesh<MFT> &msh, SolutionFieldAnalytical &sol){
//  SolutionFieldAnalytical sol;
//}

template<int idim, int ideg, int pdeg, int pnorm, bool iexact>
double interpErrGlo(const SolutionFieldAnalytical &sol){
  const MeshBase& msh = *(sol.msh);
  int nentt = msh.nentt(idim);
  const intAr2& ent2poi = msh.ent2poi(idim);
  double errGlo = 0;
  for(int ientt = 0; ientt < nentt; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;
    errGlo += interpErr<idim,ideg,pdeg,pnorm,iexact>(sol,ientt,-1,{});
  }
  errGlo = pow(errGlo, 1.0 / (double) pnorm);
  return errGlo;
}
#define BOOST_PP_LOCAL_MACRO(n)\
template double interpErrGlo<2, 1, n, 1,false>(const SolutionFieldAnalytical &sol);\
template double interpErrGlo<2, 1, n, 2,false>(const SolutionFieldAnalytical &sol);\
template double interpErrGlo<2, 2, n, 1,false>(const SolutionFieldAnalytical &sol);\
template double interpErrGlo<2, 2, n, 2,false>(const SolutionFieldAnalytical &sol);\
template double interpErrGlo<2, 1, n, 1,true >(const SolutionFieldAnalytical &sol);\
template double interpErrGlo<2, 1, n, 2,true >(const SolutionFieldAnalytical &sol);\
template double interpErrGlo<2, 2, n, 1,true >(const SolutionFieldAnalytical &sol);\
template double interpErrGlo<2, 2, n, 2,true >(const SolutionFieldAnalytical &sol);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


// Evaluate interpolation error over an element. (tdim = msh.get_tdim())
// ideg is mesh degree maximum 2
// pdeg is interpolation degree maximum 2 
// pnorm is interp err pnorm maximum 2 
// Derivatives are derived in docs/Derivatives_interpolation_error.pdf
template<int idim, int ideg, int pdeg, int pnorm, bool iexact>
double interpErr(const SolutionFieldAnalytical &sol, int ielem,
                 int idof, std::initializer_list<double*> derr){
  static_assert(idim == 2); //change ordfac
  METRIS_ASSERT(idof < 0 || idof < getnnode(idim,ideg));
  METRIS_ASSERT(derr.size() == 0 || idof >= 0);

  const MeshBase &msh = *(sol.msh);
  METRIS_ASSERT(idim == msh.get_tdim());
  METRIS_ASSERT(ideg == msh.curdeg);

  constexpr auto ordent = ORDELT(idim);

  const int nderiv = derr.size();
  constexpr int nhess = (idim*(idim+1))/2;

  // Compute integral exactly using Bernstein basis
  // The Jacobian determinant in the Bernstein basis is easy -> ccoef
  // The interpolation error is computed in Lagrange then lag2bez'd
  // The power, if needed, is obtained using square_bernstein
  // Then the product of that to the Jacobian is computed using prod_bernstein

  constexpr int ideg_jdet = idim*(ideg - 1); // Jacobian det degree
  constexpr int ideg_intp = pdeg;  // interpolate degree
  // The function is assumed of degree pdeg + 1. This would then guarantee exact:
  //constexpr int ideg_err  = (pdeg + 1)*ideg; // Pointwise error
  constexpr int ideg_err  = !iexact ? MIN(pdeg+1,(pdeg + 1)*ideg)
                                    : (pdeg + 1)*ideg; // Pointwise error
  // But we can also truncate at pdeg + 1: (and should in the future, to test)
  //constexpr int ideg_err  = (pdeg + 1); 
  constexpr int ideg_errp = ideg_err*pnorm; // Pointwise error to the power pnorm
  constexpr int ideg_intgd= ideg_errp + ideg_jdet; // Integrand degree

  // Node counts. 
  constexpr int nnode_jdet = getnnode(idim,ideg_jdet);
  constexpr int nnode_intp = getnnode(idim,ideg_intp);
  constexpr int nnode_err  = getnnode(idim,ideg_err);
  constexpr int nnode_errp = getnnode(idim,ideg_errp);
  constexpr int nnode_intgd= getnnode(idim,ideg_intgd);
  constexpr int nnode_max  = nnode_intgd;

  double *d1err = NULL, *d2err = NULL;
  if(nderiv > 0){
    auto it = derr.begin();
    d1err = *it;
    for(int ii = 0; ii < idim; ii++) d1err[ii] = 0;
    if(nderiv > 1){
      it++;
      d2err = *it;
      for(int ii = 0; ii < nhess; ii++) d2err[ii] = 0;
    }
  }


  // Evaluation/lag2bez indices (identity)
  static intAr1 lfld(nnode_max);
  for(int ii = 0; ii < nnode_max; ii++) lfld[ii] = ii;

  // -- Compute the interpolate: f(F_K(\xi^\beta)) into rfld_intp, \xi^\beta intp nodes. 
  static dblAr2 rfld_intp(nnode_intp,1);
  static dblAr2 d1fld_intp(nnode_intp,idim );
  static dblAr2 d2fld_intp(nnode_intp,nhess);
  double bary[idim+1];
  for(int inode = 0; inode < nnode_intp; inode++){
    for(int ii = 0; ii < idim + 1; ii++) 
      bary[ii] = ordent[ideg_intp][inode][ii] / (double) ideg_intp;

    METRIS_ASSERT(idim == 2 && abs(bary[0] + bary[1] + bary[2] - 1) < 1.0e-12
               || idim == 3 && abs(bary[0] + bary[1] + bary[2] + bary[3] - 1) < 1.0e-12);

    if(nderiv == 0){
      rfld_intp(inode, 0) = sol.getSolBary(idim,ielem,bary);
    }else if(nderiv == 1){
      rfld_intp(inode, 0) = sol.getSolBary(idim,ielem,bary,{d1fld_intp[inode]});
    }else{
      rfld_intp(inode, 0) = sol.getSolBary(idim,ielem,bary,{d1fld_intp[inode],d2fld_intp[inode]});
    }

    if(nderiv == 0) continue;

    // Interpolation error derivative is 
    //   sum_beta \Psi_alpha(\xi^\beta) grad f(F_K(\xi^\beta))
    // Same degree field, we simply need to multiply the d1fun, d2fun coeffs
    // we just computed by \Psi(bary). 
    double basfun;
    // Compute \Psi_\alpha(\xi^\beta)
    if(msh.getBasis() == FEBasis::Lagrange){
      basfun = eval_lagrangefunc<ideg,idim>(ordent[ideg][idof], bary, -1, NULL);
    }else{
      basfun = eval_bezierfunc<ideg,idim>(ordent[ideg][idof], bary, -1, NULL);
    }

    for(int ii = 0; ii < idim; ii++) d1fld_intp(inode,ii) *= basfun;

    if(nderiv == 1) continue;

    for(int ii = 0; ii < nhess; ii++) d2fld_intp(inode,ii) *= basfun*basfun;
  }
  // These derivatives validated.


  //printf("Debug interpolate :\n");
  //rfld_intp.print();


  // -- Compute the interpolation error
  static dblAr2 rfld_err_lag(nnode_err,1);
  static dblAr2 d1fld_err_lag(nnode_err,idim);
  static dblAr2 d2fld_err_lag(nnode_err,nhess);
  // Interpolate evaluation
  constexpr auto eval_fintp = idim == 1 ? eval1<1,ideg_intp> 
                            : idim == 2 ? eval2<1,ideg_intp> : eval3<1,ideg_intp>;
  constexpr auto eval_d1intp = idim == 1 ? eval1<idim,ideg_intp> 
                             : idim == 2 ? eval2<idim,ideg_intp> : eval3<idim,ideg_intp>;
  constexpr auto eval_d2intp = idim == 1 ? eval1<nhess,ideg_intp> 
                             : idim == 2 ? eval2<nhess,ideg_intp> : eval3<nhess,ideg_intp>;
  static double d1fun[idim];
  static double d2fun[nhess];
  static double d1intp[idim];
  static double d2intp[nhess];
  for(int inode = 0; inode < nnode_err; inode++){
    for(int ii = 0; ii < idim + 1; ii++) 
      bary[ii] = ordent[ideg_err][inode][ii] / (double) ideg_err;

    // f(F_K(\xi))
    double feval; 
    if(nderiv == 0){
      feval = sol.getSolBary(idim,ielem,bary);
    }else if(nderiv == 1){
      feval = sol.getSolBary(idim,ielem,bary,{d1fun});
    }else{
      feval = sol.getSolBary(idim,ielem,bary,{d1fun,d2fun});
    }

    // Pi_K (f \circ F_K)
    double intrp; 
    eval_fintp(rfld_intp, &lfld[0], FEBasis::Lagrange, DifVar::None, DifVar::None,
              bary, &intrp, NULL, NULL);

    rfld_err_lag(inode, 0) = pnorm == 2 ? feval - intrp : abs(feval - intrp);
    //rfld_err_lag(inode, 0) = intrp;

    if(nderiv == 0) continue;

    double basfun;
    // Compute \Psi_\alpha(\xi^\beta)
    if(msh.getBasis() == FEBasis::Lagrange){
      basfun = eval_lagrangefunc<ideg,idim>(ordent[ideg][idof], bary, -1, NULL);
    }else{
      basfun = eval_bezierfunc<ideg,idim>(ordent[ideg][idof], bary, -1, NULL);
    }

    // Derivatives are straightforward here.
    eval_d1intp(d1fld_intp, &lfld[0], FEBasis::Lagrange, DifVar::None, DifVar::None,
                bary, d1intp, NULL, NULL);

    //for(int ii = 0; ii < idim; ii++)
    //  d1fld_err_lag(inode, ii) = d1intp[ii];
    if(feval - intrp > 0 || pnorm%2 == 0){
      for(int ii = 0; ii < idim; ii++)
        d1fld_err_lag(inode, ii) = basfun*d1fun[ii] - d1intp[ii];
    }else{
      for(int ii = 0; ii < idim; ii++)
        d1fld_err_lag(inode, ii) = d1intp[ii] - basfun*d1fun[ii];
    }

    if(nderiv == 1) continue;
    eval_d2intp(d2fld_intp, &lfld[0], FEBasis::Lagrange, DifVar::None, DifVar::None,
                bary, d2intp, NULL, NULL);

    if(feval - intrp > 0 || pnorm%2 == 0){
      for(int ii = 0; ii < nhess; ii++)
        d2fld_err_lag(inode, ii) = basfun*basfun*d2fun[ii] - d2intp[ii];
    }else{
      for(int ii = 0; ii < nhess; ii++)
        d2fld_err_lag(inode, ii) = d2intp[ii] - basfun*basfun*d2fun[ii];
    }
  }
  //printf("Debug error \n");
  //rfld_err_lag.print();




  // -- Convert to Bézier
  constexpr auto lag2bez_err = idim == 1 ? lag2bez1<ideg_err,1> 
                             : idim == 2 ? lag2bez2<ideg_err,1> : lag2bez3<ideg_err,1>;
  constexpr auto lag2bez_d1err = idim == 1 ? lag2bez1<ideg_err,idim> 
                               : idim == 2 ? lag2bez2<ideg_err,idim> : lag2bez3<ideg_err,idim>;
  constexpr auto lag2bez_d2err = idim == 1 ? lag2bez1<ideg_err,nhess> 
                               : idim == 2 ? lag2bez2<ideg_err,nhess> : lag2bez3<ideg_err,nhess>;
  static dblAr2 rfld_err_bez(nnode_err,1);
  static dblAr2 d1fld_err_bez(nnode_err,idim);
  static dblAr2 d2fld_err_bez(nnode_err,nhess);
  lag2bez_err(&lfld[0], rfld_err_lag, rfld_err_bez);
  // Linear transformation transforms derivatives linearly
  if(nderiv >= 1) lag2bez_d1err(&lfld[0], d1fld_err_lag, d1fld_err_bez);
  if(nderiv >= 2) lag2bez_d2err(&lfld[0], d2fld_err_lag, d2fld_err_bez);
  //printf("Debug error bez\n");
  //rfld_err_bez.print();



  // -- Square the interpolation error if needed
  static dblAr2 rfld_errp;
  static dblAr2 d1fld_errp;
  static dblAr2 d2fld_errp;
  if constexpr (pnorm == 2){
    rfld_errp.allocate(nnode_errp,1);
    rfld_errp.set_n(nnode_errp);

    if(nderiv == 0){
      square_bernstein<idim, idim, ideg_err>(rfld_err_bez, rfld_errp);
    }else if(nderiv == 1){
      d1fld_errp.allocate(nnode_errp,idim);
      d1fld_errp.set_n(nnode_errp);
      square_bernstein<idim, idim, ideg_err>(rfld_err_bez, rfld_errp, 
                  {&d1fld_err_bez}, {&d1fld_errp});
    }else{
      d1fld_errp.allocate(nnode_errp,idim);
      d2fld_errp.allocate(nnode_errp,nhess);
      d1fld_errp.set_n(nnode_errp);
      d2fld_errp.set_n(nnode_errp);
      square_bernstein<idim, idim, ideg_err>(rfld_err_bez, rfld_errp, 
                  {&d1fld_err_bez, &d2fld_err_bez}, {&d1fld_errp, &d2fld_errp});
    }

  }else{
    rfld_errp = rfld_err_bez;
    d1fld_errp= d1fld_err_bez;
    d2fld_errp= d2fld_err_bez;
  }


  //// Debug interrupt and check derivatives of this. 
  //double retdbg = 0;
  //for(int inode = 0; inode < nnode_errp; inode++){
  //  retdbg += rfld_errp(inode, 0);
  //  if(nderiv == 0) continue;

  //  for(int ii = 0; ii < idim; ii++) d1err[ii] += d1fld_errp(inode,ii); 
  //  if(nderiv == 1) continue;

  //  for(int ii = 0; ii < nhess; ii++) d2err[ii] += d2fld_errp(inode,ii); 
  //}
  //return retdbg;

  //printf("Debug error ^pnorm\n");
  //rfld_errp.print();

  // -- Compute the Jacobian determinant
  static dblAr2 ccoef(nnode_jdet,1);
  static dblAr2 d1ccoef(nnode_jdet,idim);
  static dblAr2 d2ccoef(nnode_jdet,nhess);
  // Note: second derivatives are always 0.
  if(nderiv == 0){
    getccoef<idim, idim, ideg>(msh, ielem, NULL, ccoef[0]);
  }else{
    getccoef_dpoint<idim, ideg>(msh, ielem, idof, ccoef[0], d1ccoef);
    if(nderiv == 2) d2ccoef.fill(0);
  }
  //printf("Debug ccoef\n");
  //ccoef.print();

  // -- Compute the product of Jacobian determinant with e_Kf^p
  static dblAr2 rfld_intgd(nnode_intgd, 1);
  static dblAr2 d1fld_intgd(nnode_intgd, idim);
  static dblAr2 d2fld_intgd(nnode_intgd, nhess);

  if(nderiv == 0){
    prod_bernstein<idim, idim, ideg_errp, ideg_jdet>(rfld_errp, ccoef, rfld_intgd);
  }else if(nderiv == 1){
    prod_bernstein<idim, idim, ideg_errp, ideg_jdet>(rfld_errp, ccoef, rfld_intgd,
                {&d1fld_errp}, 
                {&d1ccoef}, 
                {&d1fld_intgd});
  }else{
    prod_bernstein<idim, idim, ideg_errp, ideg_jdet>(rfld_errp, ccoef, rfld_intgd,
                {&d1fld_errp, &d2fld_errp}, 
                {&d1ccoef, &d2ccoef}, 
                {&d1fld_intgd, &d2fld_intgd});
  }

  //printf("Debug product\n");
  //rfld_intgd.print();

  // At this stage, we have the integrand as a Bézier polynomial.
  // Now we just need to compute the integral. 
  // All the Bernsteins have the same integral
  // Hence nnode x int(any Bernstein) = int_(^K) 1 = 1/d!. 
  double err_tot = 0;
  for(int inode = 0; inode < nnode_intgd; inode++){
    //printf("debug inode %d coef %e \n",inode,rfld_intgd(inode,0));
    err_tot += rfld_intgd(inode,0);

    if(nderiv == 0) continue;
    for(int ii = 0; ii < idim; ii++){
      //printf("debug inode %d ii %d add %f \n",inode,ii,d1err[ii]);
      d1err[ii] += d1fld_intgd(inode,ii);
    }

    if(nderiv == 1) continue;
    for(int ii = 0; ii < nhess; ii++) d2err[ii] += d2fld_intgd(inode,ii);
  }
  err_tot /= (nnode_intgd*ifact<idim>()); //()

  if(nderiv == 0) return err_tot;
  for(int ii = 0; ii < idim; ii++) d1err[ii] /= (nnode_intgd*ifact<idim>());

  if(nderiv == 1) return err_tot;
  for(int ii = 0; ii < nhess; ii++) d2err[ii] /= (nnode_intgd*ifact<idim>());

  return err_tot;
}
#define BOOST_PP_LOCAL_MACRO(n)\
template double interpErr<2, 1, n, 1,false>(const SolutionFieldAnalytical &sol, int ielem,\
                               int idof, std::initializer_list<double*> derr);\
template double interpErr<2, 1, n, 2,false>(const SolutionFieldAnalytical &sol, int ielem,\
                               int idof, std::initializer_list<double*> derr);\
template double interpErr<2, 2, n, 1,false>(const SolutionFieldAnalytical &sol, int ielem,\
                               int idof, std::initializer_list<double*> derr);\
template double interpErr<2, 2, n, 2,false>(const SolutionFieldAnalytical &sol, int ielem,\
                               int idof, std::initializer_list<double*> derr);\
template double interpErr<2, 1, n, 1,true >(const SolutionFieldAnalytical &sol, int ielem,\
                               int idof, std::initializer_list<double*> derr);\
template double interpErr<2, 1, n, 2,true >(const SolutionFieldAnalytical &sol, int ielem,\
                               int idof, std::initializer_list<double*> derr);\
template double interpErr<2, 2, n, 1,true >(const SolutionFieldAnalytical &sol, int ielem,\
                               int idof, std::initializer_list<double*> derr);\
template double interpErr<2, 2, n, 2,true >(const SolutionFieldAnalytical &sol, int ielem,\
                               int idof, std::initializer_list<double*> derr);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




} //namespace
