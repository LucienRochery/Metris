//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "low_metqua_d.hxx"
#include "quafun.hxx"

#include "../low_geo/normal.hxx"
#include "../low_geo/misc.hxx"
#include "../linalg/symidx.hxx"
#include "../linalg/det.hxx"
#include "../linalg/matprods.hxx"
#include "../Mesh/Mesh.hxx"

#include "../low_geo/measure.hxx"
#include "../utils/aux_misc.hxx"
#include "../utils/aux_pp_inc.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"

#include "aux_volumeMeasure.hxx"

namespace Metris{

// START DEBUG

//#define DEBUG_MACRO(r,SEQ)  toto(SEQ )
//BOOST_PP_SEQ_FOR_EACH_PRODUCT(DEBUG_MACRO,(MFT_SEQ)(ASDEG_SEQ)(FTYPE_SEQ))

template <class MFT, int gdim, int tdim, QuaFun iquaf, typename ftype>
ftype d_metqua(Mesh<MFT> &msh, AsDeg asdmsh, AsDeg asdmet,
               int ientt,
               int ivar, FEBasis dofbas, DifVar idifmet,
               ftype*__restrict__ dquael, ftype*__restrict__ hquael,
               double difto){
  static_assert(gdim==2 || gdim==3);
  const int pnorm = msh.param->opt_pnorm;
  METRIS_ASSERT(pnorm > 0);
  constexpr int nhess = (gdim*(gdim+1))/2;

  const intAr2 &ent2poi = msh.ent2poi(tdim);

  double bary[tdim+1];

  ftype qutet = 0;
  ftype W = 1;
  double nordev = 0;
  bool do_nordev = tdim == 2 && gdim == 3
    && msh.CAD()
    && abs(msh.param->qua_surf_wt_normal) > 1.0e-9*abs(msh.param->qua_surf_wt_quality);


  constexpr auto d_quafun_xi = get_d_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
  constexpr auto ordelt = ORDELT(tdim);

  const int ideg = msh.curdeg;
  const int ideg_eff = asdmsh == AsDeg::P1 ? 1 : ideg;
  const int nnode = getnnode(tdim, ideg_eff);

  #ifdef STEPDISTANCE

  if constexpr(iquaf == QuaFun::StepDistance){

  METRIS_ASSERT(ideg_eff == 1);
  METRIS_ASSERT(pnorm == 1);
  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);

  constexpr int nnmet = (gdim*(gdim+1))/2;

  if(ivar >= 0){
    for(int ii = 0; ii < gdim; ii++) dquael[ii] = 0;
    if(hquael != NULL){
      for(int ii = 0; ii < nhess; ii++) hquael[ii] = 0;
    }
  }

  qutet = 0.;

  // Denominator of the P1 target-volume normalized element value used by
  // CavityTargetAverage. With the target metric frozen, this denominator and
  // the resulting normalized quadrature weights are constant in the local
  // derivative model. For high order, the full point-dependent theta and its
  // quotient derivatives belong here.
  double target_average_denominator = 0.;
  const bool use_target_average =
      msh.param->step_distance_cavity_target_average;
  METRIS_ENFORCE_MSG(
      !(msh.param->step_distance_shape_volume && use_target_average),
      "Step Distance Shape Volume is a distinct integration variant");

  constexpr int nquad = tdim + 2; // vertices + barycenter
  for (int iquad = 0; iquad < nquad; iquad++){

    // ------------------------------------------------------------
    // Integration scheme: vertices: 0,...,tdim, + barycenter: tdim + 1
    // ------------------------------------------------------------

    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 0.0;

    if (iquad < tdim + 1){
      bary[iquad] = 1.;
    }else{
      for (int ii = 0; ii < tdim+1; ii++){
        bary[ii] = 1./(tdim+1);
      }
    }

    const ftype wquad = (ftype)1./nquad;

    // ------------------------------------------------------------
    // Geometry: coopr and canonical-reference Jacobian.
    //
    // jmat is tdim x gdim: d_{canonical ref} F.
    // ------------------------------------------------------------
    double coopr[gdim];
    double jmat[tdim*gdim];

    if constexpr(tdim == 2){
      eval2<gdim,1>(msh.coord, ent2poi[ientt], msh.getBasis(),
                    DifVar::Bary, DifVar::None,
                    bary, coopr, jmat, NULL);
    }else{
      eval3<gdim,1>(msh.coord, ent2poi[ientt], msh.getBasis(),
                    DifVar::Bary, DifVar::None,
                    bary, coopr, jmat, NULL);
    }

    // ------------------------------------------------------------
    // Jreg_T = Jreg^T = J0^{-T} Jcanonical^T, J0: jac from cannonical reference to ideal
    //
    // Constants::invtJ_0 is J0^{-T} in this transposed convention.
    // Shape: tdim x gdim.
    // ------------------------------------------------------------
    double Jreg_T[tdim*gdim];

    for(int i = 0; i < tdim; i++){
      for(int a = 0; a < gdim; a++){
        Jreg_T[i*gdim+a] = 0.0;
        for(int k = 0; k < tdim; k++){
          Jreg_T[i*gdim+a] +=
            Constants::invtJ_0[hana::type_c<double>][tdim][i*tdim+k]
            *jmat[k*gdim+a];
        }
      }
    }

    // ------------------------------------------------------------
    // Frozen metric at quadrature point.
    //
    // Important: DifVar::None. Metric is evaluated at this point,
    // but treated as constant w.r.t. node motion.
    // ------------------------------------------------------------
    double met[nnmet];

    auto get_frozen_metric_at_quad = [&](int iquad, const double* bary,
                                     const double* coopr,
                                     double* met){

      if(iquad < tdim + 1){
        // Vertex quadrature sample
        const int ipoin = ent2poi(ientt, iquad);
        for(int im = 0; im < nnmet; im++){
          met[im] = msh.met(ipoin, im);
        }
        return;
      }

      // Barycenter quadrature sample
      if constexpr(std::is_same<MFT, MetricFieldAnalytical>::value){
        msh.met.getMetPhys(DifVar::None, msh.met.getSpace(),
                          coopr, met, NULL);
      }else{
        msh.met.getMetBary(asdmet,
                          DifVar::None,
                          msh.met.getSpace(),
                          ent2poi[ientt],
                          tdim,
                          bary,
                          met,
                          NULL);
      }
    };
    get_frozen_metric_at_quad(iquad, bary, coopr, met);

    // ------------------------------------------------------------
    // Shape data for active local P1 vertex.
    // ------------------------------------------------------------
    double gradN[tdim];

    if(ivar >= 0){
      METRIS_ASSERT(ivar >= 0 && ivar < tdim + 1);

      for(int i = 0; i < tdim; i++) gradN[i] = 0.0;

      if(ivar == 0){
        for(int i = 0; i < tdim; i++){
          for(int k = 0; k < tdim; k++){
            gradN[i] -=
              Constants::invtJ_0[hana::type_c<double>][tdim][i*tdim+k];
          }
        }
      }else{
        const int column = ivar - 1;
        for(int i = 0; i < tdim; i++){
          gradN[i] =
            Constants::invtJ_0[hana::type_c<double>][tdim][i*tdim+column];
        }
      }
    }

    // ------------------------------------------------------------
    // Phi handled by d_quafun_xi.
    //
    // For frozen metric, phi derivatives should also ignore metric variation.
    // That means d_quafun_xi should interpret idifmet = DifVar::None for this
    // branch, or the step-distance quafun should not use metric derivatives.
    // ------------------------------------------------------------
    ftype phi;
    ftype dphi[gdim];
    ftype hphi[nhess];

    if(ivar < 0){
      phi = d_quafun_xi(msh, asdmsh, asdmet,
                        ent2poi[ientt], bary, met,
                        ivar, dofbas, DifVar::None,
                        NULL, NULL);
    }else if(hquael == NULL){
      phi = d_quafun_xi(msh, asdmsh, asdmet,
                        ent2poi[ientt], bary, met,
                        ivar, dofbas, DifVar::None,
                        dphi, NULL);
    }else{
      phi = d_quafun_xi(msh, asdmsh, asdmet,
                        ent2poi[ientt], bary, met,
                        ivar, dofbas, DifVar::None,
                        dphi, hphi);
    }

    if(msh.param->step_distance_shape_volume
       && phi >= ftype(0.5*step_distance_shape_volume_rejection_quality)){
      return ftype(step_distance_shape_volume_rejection_quality);
    }

    if(use_target_average){
      const double target_density =
          VolumeMeasureHelpers::eval_target_metric_volume_density<
              gdim,double>(met);
      const ftype target_weight = wquad*(ftype)target_density;

      qutet += target_weight*phi;
      target_average_denominator += wquad*target_density;

      if(ivar < 0) continue;

      for(int i = 0; i < gdim; i++){
        dquael[i] += target_weight*dphi[i];
      }

      if(hquael == NULL) continue;

      for(int i = 0; i < gdim; i++){
        for(int j = i; j < gdim; j++){
          hquael[sym2idx(i,j)] +=
              target_weight*hphi[sym2idx(i,j)];
        }
      }
      continue;
    }

    // ------------------------------------------------------------
    // Theta value and optional derivatives.
    // ------------------------------------------------------------
    double theta_d;
    double dtheta_d[gdim];
    double htheta[nhess];

    if(ivar < 0){
      VolumeMeasureHelpers::eval_theta_fixed_metric_grad<gdim,tdim,double>(
          Jreg_T, met, NULL,
          &theta_d, NULL);
    }else if(msh.param->step_distance_shape_volume){
      // Step Distance Shape Volume deliberately keeps theta frozen.  Its
      // pointwise SPD distance supplies the collapse coercivity.
      VolumeMeasureHelpers::eval_theta_fixed_metric_grad<gdim,tdim,double>(
          Jreg_T,met,NULL,&theta_d,NULL);
      for(int i = 0; i < gdim; i++) dtheta_d[i] = 0.0;
      if(hquael != NULL){
        for(int i = 0; i < nhess; i++) htheta[i] = 0.0;
      }
    }else{
      #ifdef STEPDISTANCE_INCLUDE_GEOM_THETA_DERIV

      VolumeMeasureHelpers::eval_theta_fixed_metric_grad<gdim,tdim,double>(
          Jreg_T, met, gradN,
          &theta_d, dtheta_d);

      if(hquael != NULL){
        VolumeMeasureHelpers::eval_theta_fixed_metric_hess_by_surreal<gdim,tdim>(
            Jreg_T, met, gradN,
            htheta);
      }

      #else

      VolumeMeasureHelpers::eval_theta_fixed_metric_grad<gdim,tdim,double>(
          Jreg_T, met, NULL,
          &theta_d, NULL);

      for(int i = 0; i < gdim; i++) dtheta_d[i] = 0.0;
      if(hquael != NULL){
        for(int i = 0; i < nhess; i++) htheta[i] = 0.0;
      }

      #endif
    }

    const ftype theta = (ftype)theta_d;

    // ------------------------------------------------------------
    // Unweighted metric-volume barrier.
    // rho = sqrt(det(Jreg^T M Jreg)); M is frozen, but J derivatives are
    // always retained because preventing geometric collapse is its purpose.
    // ------------------------------------------------------------
    double rho_d;
    double barrier_d;
    double dbarrier_d[gdim];
    double hbarrier_d[nhess];
    const double barrier_beta = msh.param->step_distance_shape_volume
                              ? 0.0
                              : msh.param->step_distance_barrier_beta;

    if(ivar < 0){
      VolumeMeasureHelpers::eval_metric_volume_barrier_fixed_metric_grad<
          gdim,tdim,double>(
              Jreg_T,met,NULL,
              msh.param->step_distance_barrier_rho0,
              barrier_beta,
              &rho_d,&barrier_d,NULL);
    }else{
      VolumeMeasureHelpers::eval_metric_volume_barrier_fixed_metric_grad<
          gdim,tdim,double>(
              Jreg_T,met,gradN,
              msh.param->step_distance_barrier_rho0,
              barrier_beta,
              &rho_d,&barrier_d,dbarrier_d);
      if(hquael != NULL){
        VolumeMeasureHelpers::
            eval_metric_volume_barrier_fixed_metric_hess_by_surreal<
                gdim,tdim>(
                    Jreg_T,met,gradN,
                    msh.param->step_distance_barrier_rho0,
                    barrier_beta,
                    hbarrier_d);
      }
    }

    // ------------------------------------------------------------
    // Value.
    // ------------------------------------------------------------
    qutet += wquad*(phi*theta + (ftype)barrier_d);

    if(ivar < 0) continue;

    // ------------------------------------------------------------
    // First derivative:
    //
    // d(phi theta + B) = theta dphi + phi dtheta + dB.
    // ------------------------------------------------------------
    for(int i = 0; i < gdim; i++){
      dquael[i] += wquad*(
          theta*dphi[i] + phi*(ftype)dtheta_d[i]
          + (ftype)dbarrier_d[i]);
    }

    if(hquael == NULL) continue;

    // ------------------------------------------------------------
    // Hessian:
    //
    // H(phi theta + B) =
    //   theta Hphi
    // + phi Htheta
    // + dtheta dphi^T
    // + dphi dtheta^T
    // + HB.
    // ------------------------------------------------------------
    for(int i = 0; i < gdim; i++){
      for(int j = i; j < gdim; j++){
        hquael[sym2idx(i,j)] += wquad*(
            theta*hphi[sym2idx(i,j)]
          + phi*(ftype)htheta[sym2idx(i,j)]
          + (ftype)dtheta_d[i]*dphi[j]
          + dphi[i]*(ftype)dtheta_d[j]
          + (ftype)hbarrier_d[sym2idx(i,j)]
        );
      }
    }
  }

  if(use_target_average){
    METRIS_ENFORCE_MSG(target_average_denominator > 0.0,
                       "Nonpositive StepDistance target-volume denominator");
    qutet /= target_average_denominator;

    if(ivar >= 0){
      for(int i = 0; i < gdim; i++){
        dquael[i] /= target_average_denominator;
      }
      if(hquael != NULL){
        for(int i = 0; i < nhess; i++){
          hquael[i] /= target_average_denominator;
        }
      }
    }
  }

  return qutet;
  }

#endif


  #ifdef TESTQUALITYALGO
  // Assumptions for quality algo:
  METRIS_ASSERT(ideg_eff == 1);
  METRIS_ASSERT(pnorm == 1);
  #endif

  // Accumulate normal error at the nodes (depending on asdmsh)
  if(do_nordev){
    double result[18];
    double norCAD[gdim], norelt[gdim];
    double *du = &result[3];
    double *dv = &result[6];
    const int iref = msh.fac2ref[ientt];
    const ego obj  = msh.CAD.cad2fac[iref];

    // Even if we use CAD normals at all vertices, we can compute this one just once.
    if(ideg == 1){
      getnorfacP1(ent2poi[ientt], msh.coord, norelt);
      if(normalize_vec<gdim>(norelt)){
        METRIS_THROW_MSG( "Normal (elt) vanishes");
      }
    }

    for(int inode = 0; inode < nnode; inode++){
      int ipoin = ent2poi(ientt, inode);
      int ibpoi = msh.poi2ebp(ipoin, tdim, ientt, iref);
      METRIS_ASSERT_MSG(ibpoi >= 0 && ibpoi < msh.nbpoi,
        "iface = {} iref = {} inode = {} ipoin = {} ibpoi = {}",
        ientt,iref,inode,ipoin,ibpoi);
      int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
      METRIS_ENFORCE_MSG(ierro == 0, "metqua0 EG_evaluate error {}", ierro);
      vecprod(du,dv,norCAD);
      if(normalize_vec<gdim>(norCAD)){
        // legitimate if e.g. cone tip.
        METRIS_ENFORCE_MSG(msh.getpoitdim(ipoin) == 0,
          "Normal (CAD) vanishes at non-corner point {}", ipoin)
        nordev += 0;
        continue;
      }

      if(ideg > 1){
        for(int ii = 0; ii < tdim + 1; ii++)
          bary[ii] = ordelt[ideg_eff][inode][ii]/((double) (ideg_eff));
        getnorfac(msh, ientt, bary, asdmsh, norelt);
        if(normalize_vec<gdim>(norelt)){
          GETVDEPTH(msh.param);
          MPRINTF("norelt vanished ientt = {} node {} point {} nodes {}\n",
                  ientt,inode,ipoin,intAr1(nnode, ent2poi[ientt]));
          for(int ii = 0; ii < gdim; ii++) MPRINTF("{}: {:23.15e}\n",ii,norelt[ii]);
          METRIS_THROW_MSG( "Normal (elt) vanishes");
        }
      }

      double dtprd = getprdl2<gdim>(norelt, norCAD);
      double tmp = 1 - abs(dtprd);
      METRIS_ASSERT(tmp >= 0);
      nordev += tmp*tmp;
    }
    nordev /= nnode;
    nordev = sqrt(nordev);
  }

  if(ideg_eff > 1){

    const auto ordelt = ORDELT(tdim);

    if(hquael != NULL && ivar >= 0)
      for(int ii = 0; ii < nhess; ii++) hquael[ii] = 0;
    if(ivar >= 0)
      for(int ii = 0; ii < gdim; ii++) dquael[ii] = 0;

    METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);

    ftype qua0, dqua0[gdim], hqua0[nhess];
    for(int iquad = 0; iquad < nnode; iquad++){

      for(int ii = 0; ii < tdim + 1; ii++)
        bary[ii] = ordelt[ideg][iquad][ii]/((double) (ideg));

      const int ipoin = ent2poi(ientt,iquad);
      if(hquael == NULL){
        qua0 = d_quafun_xi(msh,asdmsh,asdmet,
                           ent2poi[ientt],bary,msh.met[ipoin],
                           ivar,dofbas,idifmet,
                           dqua0,NULL);
      }else{
        qua0 = d_quafun_xi(msh,asdmsh,asdmet,
                           ent2poi[ientt],bary,msh.met[ipoin],
                           ivar,dofbas,idifmet,
                           dqua0,hqua0);
      }
      ftype powm1 = pow(abs(qua0 - difto),pnorm-1);
      qutet += abs(qua0 - difto)*abs(powm1)/nnode;

      if(ivar < 0) continue;

      // If p%2 = 1 and f = qua0 - difto < 0,
      // notice d|f|^p = d(-f)^p = p(-f')(-f)^(p-1) = p(-1)^p f' f^(p-1)
      // but since p%2 = 1, (-1)^p = -1, so d|f|^p = -pf'f^(p-1) = -d(f)^p
      // note f^(p-1) = |f|^(p-1) as p-1 = 0 [2]
      int sg = 1;
      if(qua0 - difto < 0) sg = -1;
      for(int ii = 0; ii < gdim; ii++){
        dquael[ii] += sg*pnorm*dqua0[ii]*powm1/nnode;
      }

      if(hquael == NULL) continue;

      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          hquael[sym2idx(ii,jj)] += sg*pnorm*hqua0[sym2idx(ii,jj)]*powm1;
        }
      }

      if(pnorm < 2) continue;

      ftype powm2 = pow(abs(qua0 - difto),pnorm-2);
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          // Only the gradient gets a sign difference
          hquael[sym2idx(ii,jj)] +=
            pnorm*(pnorm - 1)*dqua0[ii]*dqua0[jj]*powm2;
        }
      }

    }// for iquad
  }else{

  #if defined(ONEPOINTQUAL)

    #if 0
    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0/(tdim  + 1);
    qutet = d_quafun_xi(msh,asdmsh,asdmet,
                        ent2poi,bary,
                        ivar,dofbas,idifmet,
                        dquael,hquael);
    ftype powm1 = pow(abs(qutet - difto),pnorm-1);
    qutet = powm1*abs(qutet - difto);
    #else
    METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
    constexpr int nnmet = (gdim*(gdim+1))/2;
    double met[nnmet];
    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim  + 1);
    for(int jj = 0; jj < nnmet; jj++) met[jj] = 0;
    for(int ii = 0; ii < tdim + 1; ii++){
      int ipoin = ent2poi(ientt,ii);
      for(int jj = 0; jj < nnmet; jj++){
        met[jj] += msh.met(ipoin,jj) / (tdim + 1);
      }
    }
    qutet = d_quafun_xi(msh,asdmsh,asdmet,
                        ent2poi[ientt],bary,met,
                        ivar,dofbas,idifmet,
                        dquael,hquael);
    // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
    // time to run pow() here even if pnorm = 2 or 1.
    ftype quaerr = qutet - difto;
    int sg = quaerr < 0 ? -1 : 1;
    qutet = abs(quaerr);
    ftype powm1 = 1; // case pnorm 1
    #ifdef TESTQUALITYALGO
    METRIS_ASSERT(pnorm == 1);
    ftype meas;
    isvalideltP1<gdim,tdim>(msh,ientt,NULL,&meas);
    W = meas;
    #ifdef INTQUALINRIEMSPACE
    if (tdim == 2){
      W *= sqrt(met[0]*met[2] - met[1]*met[1]);
    }else {
      const double m11 = met[0];
      const double m12 = met[1];
      const double m22 = met[2];
      const double m13 = met[3];
      const double m23 = met[4];
      const double m33 = met[5];

      W *= sqrt(
        m11*(m22*m33 - m23*m23)
      - m12*(m12*m33 - m13*m23)
      + m13*(m12*m23 - m13*m22) );
    }
    #endif
    #else
    if(pnorm == 2){
      powm1 = qutet;
      qutet *= qutet;
    }else if(pnorm > 2){
      powm1 = pow(qutet,pnorm-1);
      qutet = qutet*powm1;
    }
    #endif
    #endif

    if(ivar < 0) return W * qutet;

    for(int ii = 0; ii < gdim; ii++){
      dquael[ii] = W * sg*pnorm*dquael[ii]*powm1;
    }

    if(hquael == NULL) return W * qutet;

    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        hquael[sym2idx(ii,jj)] = W * sg*pnorm*hquael[sym2idx(ii,jj)]*powm1;
      }
    }
    if(pnorm >= 2){
      ftype powm2 = pow(abs(qutet - difto),pnorm-2);
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          // Only the gradient gets a sign difference
          hquael[sym2idx(ii,jj)] += pnorm*(pnorm - 1)*dquael[ii]*dquael[jj]*powm2;
        }
      }
    }

  #elif defined(TDIM1POINTSQUAL)

    METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
    constexpr int nnmet = (gdim*(gdim+1))/2;

    // zero outputs
    if (ivar >= 0) {
      for (int ii=0; ii<gdim; ++ii) dquael[ii] = 0;
      if (hquael) for (int ii=0; ii<nhess; ++ii) hquael[ii] = 0;
    }

    auto sqrtDetM_of = [&](const double* metptr) -> double {
      if constexpr (tdim == 2){
        const double m11 = metptr[0];
        const double m12 = metptr[1];
        const double m22 = metptr[2];
        const double det = m11*m22 - m12*m12;
        return det > 0 ? std::sqrt(det) : 0.0;
      } else { // tdim==3
        const double m11 = metptr[0];
        const double m12 = metptr[1];
        const double m22 = metptr[2];
        const double m13 = metptr[3];
        const double m23 = metptr[4];
        const double m33 = metptr[5];
        const double det =
            m11*(m22*m33 - m23*m23)
          - m12*(m12*m33 - m13*m23)
          + m13*(m12*m23 - m13*m22);
        return det > 0 ? std::sqrt(det) : 0.0;
      }
    };

    double meas = 1.0;
    isvalideltP1<gdim,tdim>(msh, ientt, NULL, &meas);

    // common weight factor = |T|/5
    const ftype W0 = (ftype)(meas / 5.0);

    ftype val_sum = 0;

    // workspace for one sample
    ftype q0, dq0[gdim], hq0[nhess];

    auto accumulate_one_sample =
      [&](const double* bary_in, const double* met_in)
    {
      // copy bary for the call
      for(int ii=0; ii<tdim+1; ++ii) bary[ii] = bary_in[ii];

      // eval q and its derivatives
      if (hquael == NULL) {
        q0 = d_quafun_xi(msh, asdmsh, asdmet,
                        ent2poi[ientt], bary, met_in,
                        ivar, dofbas, idifmet,
                        dq0, NULL);
      } else {
        q0 = d_quafun_xi(msh, asdmsh, asdmet,
                        ent2poi[ientt], bary, met_in,
                        ivar, dofbas, idifmet,
                        dq0, hq0);
      }

      const ftype e  = q0 - (ftype)difto;
      const ftype ae = abs(e);
      const int sg   = (e < 0) ? -1 : 1;

      // value contribution: W0 * sqrt(detM) * |e|^p
      ftype epow = ae;                   // p=1
      if (pnorm == 2) epow = ae*ae;
      else if (pnorm > 2) epow = ae*pow(ae, pnorm-1);

      const ftype Wk = W0 * (ftype)sqrtDetM_of(met_in);

      val_sum += Wk * epow;

      if (ivar < 0) return;

      // powm1 = |e|^(p-1) (with p=1 => 1)
      ftype powm1 = 1;
      if (pnorm == 2) powm1 = ae;
      else if (pnorm > 2) powm1 = pow(ae, pnorm-1);

      for (int ii=0; ii<gdim; ++ii) {
        dquael[ii] += Wk * (ftype)(sg * pnorm) * dq0[ii] * powm1;
      }

      if (!hquael) return;

      // Hessian term from Hq
      for (int ii=0; ii<gdim; ++ii) {
        for (int jj=ii; jj<gdim; ++jj) {
          hquael[sym2idx(ii,jj)] += Wk * (ftype)(sg * pnorm) * hq0[sym2idx(ii,jj)] * powm1;
        }
      }

      if (pnorm < 2) return;

      // Outer-product term: p(p-1)|e|^(p-2) dq dq^T  (no sign)
      ftype powm2 = 1;
      if (pnorm == 2) powm2 = 1;
      else powm2 = pow(ae, pnorm-2);

      for (int ii=0; ii<gdim; ++ii) {
        for (int jj=ii; jj<gdim; ++jj) {
          hquael[sym2idx(ii,jj)] += Wk * (ftype)(pnorm*(pnorm-1)) * dq0[ii]*dq0[jj] * powm2;
        }
      }
    };

    // ---- 4 vertex samples ----
    for (int iver=0; iver<tdim+1; ++iver) {
      double bary_v[tdim+1] = {0};
      bary_v[iver] = 1.0;

      const int ipoin = ent2poi(ientt, iver);

      double met_v[nnmet];
      for (int jj=0; jj<nnmet; ++jj) met_v[jj] = msh.met(ipoin, jj);

      accumulate_one_sample(bary_v, met_v);
    }

    // ---- centroid sample ----
    double bary_c[tdim+1];
    for (int ii=0; ii<tdim+1; ++ii) bary_c[ii] = 1.0/(tdim+1);

    double met_c[nnmet];
    // compute physical centroid coordinates
    double coordC[3] = {0.0, 0.0, 0.0};
    for (int iver = 0; iver < tdim + 1; ++iver) {
      const int ip = ent2poi(ientt, iver);
      coordC[0] += bary_c[iver] * msh.coord(ip, 0);
      coordC[1] += bary_c[iver] * msh.coord(ip, 1);
      coordC[2] += bary_c[iver] * msh.coord(ip, 2);
    }

    // evaluate metric at centroid (analytic)
    if constexpr(std::is_same<MFT, MetricFieldAnalytical>::value) {
      msh.met.getMetPhys(DifVar::None, msh.met.getSpace(), coordC, met_c, NULL);
    } else {
      METRIS_THROW_MSG("TDIM1POINTSQUAL deriv: metric eval at centroid not implemented for MetricFieldFE");
    }

    accumulate_one_sample(bary_c, met_c);

    qutet = val_sum;

  #elif defined(KEAST4QUAL)

    METRIS_ASSERT_MSG(tdim == 3, "Keast degree-4 rule implemented for tetrahedra only (tdim=3).");
    METRIS_ASSERT(gdim==3);

    METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);

    constexpr int nnmet = (gdim*(gdim+1))/2;

    // zero outputs
    if (ivar >= 0) {
      for (int ii=0; ii<gdim; ii++) dquael[ii] = 0;
      if (hquael) for (int ii=0; ii<nhess; ii++) hquael[ii] = 0;
    }

    auto sqrtDetM_of = [&](const double* metptr) -> double {
      const double m11 = metptr[0];
      const double m12 = metptr[1];
      const double m22 = metptr[2];
      const double m13 = metptr[3];
      const double m23 = metptr[4];
      const double m33 = metptr[5];

      const double det =
          m11*(m22*m33 - m23*m23)
        - m12*(m12*m33 - m13*m23)
        + m13*(m12*m23 - m13*m22);

      return det > 0 ? std::sqrt(det) : 0.0;
    };

    double meas = 1.0;
    isvalideltP1<gdim,tdim>(msh, ientt, NULL, &meas);

    // Keast degree-4 (11-pt) rule on reference tetrahedron.
    // Renormalize weights by 6 so they sum to 1 (to match "meas * weighted average").
    constexpr double w0  =  0.013155555555555555;
    constexpr double w1  =  0.007622222222222222;
    constexpr double w2  =  0.024888888888888888;
    constexpr double wn0 = 6.0*w0;
    constexpr double wn1 = 6.0*w1;
    constexpr double wn2 = 6.0*w2;

    // [3,1] orbit parameter
    constexpr double a  = 0.100596423833200785;
    constexpr double b  = 1.0 - 3.0*a;

    // [2,2] orbit parameters
    constexpr double bb = 0.0714285714285714285;
    constexpr double cc = 0.5 - bb;

    // workspace
    ftype q0 = 0;
    ftype dq0[gdim];
    ftype hq0[nhess];

    ftype val_sum = 0;

    auto accumulate_one =
      [&](const double lam[4], double wn)
    {
      // bary for the call
      for (int ii=0; ii<4; ii++) bary[ii] = lam[ii];

      // physical coordinates of quadrature point
      double coordQ[3] = {0.0, 0.0, 0.0};
      for (int iv = 0; iv < 4; iv++) {
        const int ip = ent2poi(ientt, iv);
        coordQ[0] += lam[iv] * msh.coord(ip, 0);
        coordQ[1] += lam[iv] * msh.coord(ip, 1);
        coordQ[2] += lam[iv] * msh.coord(ip, 2);
      }

      // metric at quadrature point (analytic)
      double metq[nnmet];
      if constexpr(std::is_same<MFT, MetricFieldAnalytical>::value) {
        msh.met.getMetPhys(DifVar::None, msh.met.getSpace(), coordQ, metq, NULL);
      } else {
        METRIS_THROW_MSG("KEAST4QUAL deriv: metric eval at quadrature points not implemented for MetricFieldFE");
      }

      const ftype Wk = (ftype)meas * (ftype)wn * (ftype)sqrtDetM_of(metq);

      // evaluate q and derivatives
      if (hquael == NULL) {
        q0 = d_quafun_xi(msh, asdmsh, asdmet,
                        ent2poi[ientt], bary, metq,
                        ivar, dofbas, idifmet,
                        dq0, NULL);
      } else {
        q0 = d_quafun_xi(msh, asdmsh, asdmet,
                        ent2poi[ientt], bary, metq,
                        ivar, dofbas, idifmet,
                        dq0, hq0);
      }

      // pnorm==1: contribution is |q-difto|
      const ftype e  = q0 - (ftype)difto;
      const ftype ae = abs(e);
      const int sg   = (e < 0) ? -1 : 1;

      val_sum += Wk * ae;

      if (ivar < 0) return;

      // d|e| = sg * de
      for (int ii=0; ii<gdim; ii++) {
        dquael[ii] += Wk * (ftype)sg * dq0[ii];
      }

      if (!hquael) return;

      // d2|e| = sg * d2e
      for (int ii=0; ii<gdim; ii++) {
        for (int jj=ii; jj<gdim; jj++) {
          hquael[sym2idx(ii,jj)] += Wk * (ftype)sg * hq0[sym2idx(ii,jj)];
        }
      }
    };

    // centroid
    {
      const double lam[4] = {0.25, 0.25, 0.25, 0.25};
      accumulate_one(lam, wn0);
    }

    // point orbit: permutations of (b,a,a,a)
    for (int ibig=0; ibig<4; ++ibig) {
      double lam[4] = {a, a, a, a};
      lam[ibig] = b;
      accumulate_one(lam, wn1);
    }

    // 6-point orbit: permutations of (bb,bb,cc,cc)
    for (int i0=0; i0<4; ++i0) {
      for (int i1=i0+1; i1<4; ++i1) {
        double lam[4] = {cc, cc, cc, cc};
        lam[i0] = bb;
        lam[i1] = bb;
        accumulate_one(lam, wn2);
      }
    }

    qutet = val_sum;

  #else
    static_assert(tdim == 3, "No rule for quality integration defined!");
  #endif
  }

  if(do_nordev){
    METRIS_ASSERT(msh.param->qua_surf_wt_quality >= 0);
    METRIS_ASSERT(msh.param->qua_surf_wt_normal  >= 0);
    qutet = msh.param->qua_surf_wt_quality*qutet*W
          + msh.param->qua_surf_wt_normal*pow(nordev, pnorm); // for homogeneity
  }
  return W * qutet;
}

//// While cumbersome, this replaces a bunch of manual instantiations, about to
//// be made worse the day we add tdimn as a template argument.
//#define EXPAND_TEMPLATE(z,gdim,SEQ) \
//                  INSTANTIATE(gdim,BOOST_PP_SEQ_ELEM(0, SEQ),\
//                                   BOOST_PP_SEQ_ELEM(1, SEQ),\
//                                   BOOST_PP_SEQ_ELEM(2, SEQ))
//#define REPEAT_GDIM(r,SEQ) BOOST_PP_REPEAT(2,EXPAND_TEMPLATE,SEQ)
//#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical)
//#define QUAFUN_SEQ (QuaFun::Distortion)(QuaFun::Unit)
//#define INSTANTIATE(gdim,MFT_VAL,QUAFUN,FTYPE)\
//template FTYPE d_metqua< MFT_VAL , 2+gdim, QUAFUN,FTYPE>\
//                  (Mesh< MFT_VAL > &msh, AsDeg asdmsh, AsDeg asdmet,\
//                   int ientt, \
//                   int ivar, \
//                   FEBasis dofbas, \
//                   DifVar idifmet, \
//                   FTYPE*__restrict__ dquael, FTYPE*__restrict__ hquael, \
//                   double difto);
//BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_GDIM,\
//                              (MFT_SEQ)(QUAFUN_SEQ)(QUA_FTYPE_SEQ))
//#undef INSTANTIATE
//#undef EXPAND_TEMPLATE
//#undef REPEAT_GDIM



// While cumbersome, this replaces a bunch of manual instantiations, about to
// be made worse the day we add tdimn as a template argument.
#define EXPAND_TEMPLATE(r,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ),\
                              BOOST_PP_SEQ_ELEM(2, SEQ))
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical)
#define QUAFUN_SEQ (QuaFun::Distortion)(QuaFun::Unit)(QuaFun::SizeShape)(QuaFun::StepDistance)
#define INSTANTIATE(MFT_VAL,QUAFUN,FTYPE)\
template FTYPE d_metqua< MFT_VAL , 2, 2, QUAFUN,FTYPE>\
                  (Mesh< MFT_VAL > &msh, AsDeg asdmsh, AsDeg asdmet,\
                   int ientt, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, FTYPE*__restrict__ hquael, \
                   double difto);\
template FTYPE d_metqua< MFT_VAL , 3, 2, QUAFUN,FTYPE>\
                  (Mesh< MFT_VAL > &msh, AsDeg asdmsh, AsDeg asdmet,\
                   int ientt, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, FTYPE*__restrict__ hquael, \
                   double difto);\
template FTYPE d_metqua< MFT_VAL , 3, 3, QUAFUN,FTYPE>\
                  (Mesh< MFT_VAL > &msh, AsDeg asdmsh, AsDeg asdmet,\
                   int ientt, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, FTYPE*__restrict__ hquael, \
                   double difto);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,\
                              (MFT_SEQ)(QUAFUN_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE
#undef REPEAT_GDIM


/*
------------------------- Full derivatives (all control points)
*/


/*
D_quafun_distortion : computes derivatives (up to second) wrt all control points in the
element.
- dquael is ftype* size N = gdim*entnpps[ideg] (entnpps -> gdim)
- hquael is ftype* size sym NxN. Ordering:
  hquael(symidx(gdim*ip1 + i,gdim*ip2 + j)) = d_{ij}^{ip1,ip2}
Refer to d_quafun_distortion for additional comments within the routine.
*/
template <class MFT, int gdim, int ideg, AsDeg asdmet, typename ftype>
ftype D_quafun_distortion(Mesh<MFT> &msh,
                  const int* ent2poi,
                  const double*__restrict__ bary,
                  FEBasis dofbas,
                  DifVar idifmet,
                  ftype*__restrict__ dquael,
                  ftype*__restrict__ hquael){

  const int power = msh.param->opt_power;
  static_assert(gdim == 2 || gdim == 3);
  METRIS_ASSERT(gdim == msh.idim);
  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Log)
  // Differentiate or don't, but there is no barycentric derivative in this context
  METRIS_ASSERT(idifmet == DifVar::None || idifmet == DifVar::Phys);
  if(idifmet != DifVar::None) METRIS_THROW_MSG(
                           "Metric field derivative not implemented in quality")
  //METRIS_ASSERT( !(idiff != DifVar::None && dofbas == FEBasis::Undefined) );
  if(dofbas == FEBasis::Bezier && idifmet != DifVar::None) METRIS_THROW_MSG(
    "Ctrl pt dof not implemented -> do lag2bez derivatives of metric")
  //METRIS_ASSERT( !(idiff != DifVar::None && dquael == NULL) );

  constexpr int tdim  = gdim;
  constexpr int nnmet = (gdim*(gdim+1))/2;
  constexpr int nhess = nnmet;


  ftype quael;


  // Get Jacobian matrix at xi
  // Derivatives are not needed, we compute them ourselves, as they greatly simplify
  // see docs/quality/qualityiff.pdf
  double jmat[tdim*gdim],coopr[gdim];
  if constexpr(tdim == 2){
    eval2<gdim,ideg>(msh.coord,ent2poi,msh.getBasis(),DifVar::Bary,DifVar::None,
                                                          bary,coopr,jmat,NULL);
  }else{
    eval3<gdim,ideg>(msh.coord,ent2poi,msh.getBasis(),DifVar::Bary,DifVar::None,
                                                          bary,coopr,jmat,NULL);
  }

  // Compute J_0^{-T} J_K^T M J_K J_0^{-1}
  // Note that J_K is stored transposed w.r.t. above ! -> jmat[i,j] = d_i F_j
  // whereas (J_K)_{ij} = d_j F_i ..

  // Get J_0^{-T} J_K^T
  ftype invtJ0_tJ[tdim][gdim];
  matXmat<tdim,tdim,gdim>(Constants::invtJ_0[hana::type_c<ftype>][tdim],
                                                              jmat,invtJ0_tJ[0]);
  // Get metric interpolated at xi
  // The metric field class fetches geometric dimension from the mesh
  double met[nnmet],dmet[gdim][nnmet];
  msh.met.getMetFullinfo(asdmet,idifmet,MetSpace::Exp,
                         ent2poi,tdim,bary,coopr,met,dmet[0]);

  //printf("Debug metric = \n");
  //dblAr1(nnmet,met).print();

  // Compute the trace
  ftype tra;
  //ftype invtJ0_tJ_M_J_invJ0iag[tdim];
  tra = tra_matXsymXtmat<gdim, double, ftype, ftype>(met, invtJ0_tJ[0]);
  //tra = invtJ0_tJ_M_J_invJ0iag[0] + invtJ0_tJ_M_J_invJ0iag[1];
  //if constexpr (tdim == 3) tra += invtJ0_tJ_M_J_invJ0iag[2];
  // This is an actual exception that should never theoretically happen.
  if(tra < 1.0e-16) METRIS_THROW_MSG( "NEGATIVE J^TMJ trace = {:e}", tra);



  ftype det ;
  ftype detM, det_invtJ0_tJ;
  if constexpr(tdim == gdim){
    det_invtJ0_tJ = detmat<gdim,ftype>(invtJ0_tJ[0]); // Error of 10^{-19} = 10^-14 relative compared to Matlab on wonky case
    //ftype tmp = detsym<gdim>(met); // Error of 10^5 ... // Matlab yields 7.346639223296765e+09 we get 7.346911345315383e+09 // Even in relative this is terrible
    detM = detsym2<gdim>(met); // Also wrong, same error. Note Matlab gets 10^-9 final quality relative error ! Our determinant is terribly bad
    det = det_invtJ0_tJ*det_invtJ0_tJ*detM;
  }else{
    static_assert(tdim == 2);
    ftype J0tJtMJJ0[3];
    matXsymXtmat<2,3,double,ftype,ftype>(met,invtJ0_tJ[0],J0tJtMJJ0);
    det = detsym2<2>(J0tJtMJJ0);
  }

  if(abs(det) < 1.0e-16 && power > 0)
     METRIS_THROW_MSG( "Singular J^TMJ det = {:e}", det);

  //constexpr auto irpow = [&]<int n>(ftype x) -> ftype{
  //  if constexpr(std::is_same<ftype,double>::value) return idpow<n>(x);
  //  if constexpr(std::is_same<ftype,float4>::value) return id4pow<n>(x);
  //  if constexpr(std::is_same<ftype,float8>::value) return id8pow<n>(x);
  //};

  // This is used later on -> store it
  int dpowd = iipow<tdim>(tdim);
  ftype trapowdm2 = irpow<tdim-2,ftype>(tra);//irpow.template operator()<tdim-2>(tra);
  ftype trapowdm1 = trapowdm2*tra;
  ftype trapowd   = trapowdm1*tra;

  if(power > 0){
    quael = trapowd/(det*dpowd);
  }else{
    quael = (det*dpowd)/trapowd;
  }
  ftype quae1 = quael; // for derivatives
  quael = pow(quael, abs(power));



  //// From here, we compute derivatives.
  //if(idiff == DifVar::None) return quael;
  //// See docs/quality/qualityiff.pdf for details

  constexpr int nnode = tdim == 1 ? getnnod1(ideg) :
                        tdim == 2 ? getnnod2(ideg) : getnnod3(ideg) ;
  constexpr auto ordent = ORDELT(tdim);


  // Compute A^TM
  ftype invtJ0_tJ_M[tdim][gdim];
  matXsym<gdim>(invtJ0_tJ[0],met,invtJ0_tJ_M[0]);
  ftype J_invJ0[tdim][tdim];
  for(int ii = 0; ii < tdim; ii++){
    for(int jj = 0; jj < tdim; jj++){
      J_invJ0[jj][ii] = invtJ0_tJ[ii][jj];
    }
  }

  // Compute D_J_invJ0 which depends on each node,
  // then dtra and ddet
  ftype D_J_invJ0[nnode][tdim], ddet[nnode][gdim], ddetA[nnode][gdim], dtra[nnode][gdim];
  for(int inode = 0; inode < nnode; inode++){
    // Get the derivatives (d_k+1 - d_1) f
    double dfun[tdim];
    // multi-index of inode:
    int idx[tdim+1];
    for(int ii = 0; ii < tdim+1; ii++) idx[ii] = ordent[ideg][inode][ii];

    // This is what we called \psi_k in the pdf document
    if(dofbas == FEBasis::Bezier){
      eval_bezierfunc<ideg,tdim>(idx,bary,1,dfun);
    }else{
      eval_lagrangefunc<ideg,tdim>(idx,bary,1,dfun);
    }
    // The derivative is much simplified, see pdf doc
    // d_i(invtJ0_tJ)_{ij} = \sum_k \psi_k C_kj^T
    // In these notations, C_kj^T = Constants::invtJ_0[hana::type_c<ftype>][tdim][k][j]
    // BECAUSE THE JACOBIAN MATRICES ARE TRANSPOSED!
    for(int ii = 0; ii < tdim; ii++){
      D_J_invJ0[inode][ii] = 0;
      for(int kk = 0; kk < tdim; kk++){
        D_J_invJ0[inode][ii] += dfun[kk]*Constants::invtJ_0[hana::type_c<ftype>][tdim][tdim*ii+kk];
      }
    }

    for(int ii = 0; ii < gdim; ii++){
      dtra[inode][ii] = 0;
      for(int jj = 0; jj < gdim; jj++){
        dtra[inode][ii] += 2*invtJ0_tJ_M[jj][ii]
                            *D_J_invJ0[inode][jj];
      }
    }

    if constexpr (gdim == 2){
      ddetA[inode][0] = detvec2<ftype>(D_J_invJ0[inode],  J_invJ0[1]);
      ddetA[inode][1] = detvec2<ftype>(  J_invJ0[0]    ,D_J_invJ0[inode]);
    }else{
      ddetA[inode][0] = detvec3<ftype>(D_J_invJ0[inode],  J_invJ0[1]    ,  J_invJ0[2]    );
      ddetA[inode][1] = detvec3<ftype>(  J_invJ0[0]    ,D_J_invJ0[inode],  J_invJ0[2]    );
      ddetA[inode][2] = detvec3<ftype>(  J_invJ0[0]    ,  J_invJ0[1]    ,D_J_invJ0[inode]);
    }
    for(int ii = 0; ii < gdim;ii++){
      ddet[inode][ii] = 2*ddetA[inode][ii]*detM*det_invtJ0_tJ;
    }

    if(power > 0){
      for(int ii = 0; ii < gdim; ii++){
        dquael[inode*gdim + ii] = trapowdm1*( tdim*dtra[inode][ii]*det
                               - tra*ddet[inode][ii])
                    /(det*det*dpowd);
        //dquael[inode*gdim + ii] = trapowdm1*( tdim*dtra[inode][ii]*det_invtJ0_tJ
        //                                    - 2.0*ddetA[inode][ii]*tra )/(det_invtJ0_tJ*det*dpowd);
      }
    }else{
      for(int ii = 0; ii < gdim; ii++){
        dquael[inode*gdim+ii] = dpowd*(     2* tra   *ddetA[inode][ii]
                                      - tdim*dtra[inode][ii]*det_invtJ0_tJ)*detM*det_invtJ0_tJ/(trapowd*tra);
      }
    }
    if(abs(power) != 1){
      // Q^(abs(power)) , quael is Q^p, quae1 is just Q
      for(int ii = 0; ii < gdim; ii++){
        dquael[inode*gdim+ii] = abs(power)*dquael[inode*gdim+ii]*quael/quae1;
      }
    }

  }


  // Sized for two points local Hessian
  constexpr int nhess2 = ( (2*gdim)*(2*gdim + 1) )/ 2;
  // Compute Hessian
  //const int dpowdpowp = pow(dpowd,power);
  for(int inode = 0; inode < nnode; inode++){
    ftype sumDk2 = 0; // needed for second derivative of trace
    for(int ii = 0; ii < tdim; ii++){
      sumDk2 += D_J_invJ0[inode][ii]*D_J_invJ0[inode][ii];
    }

    ftype htra[nhess]; // These are temp
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        htra[sym2idx(ii,jj)] = 2.0*met[sym2idx(ii,jj)]*sumDk2;
      }
    }

    ftype hdet[nhess] = {0};//, hdet2[nhess2] = {0};
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        hdet[sym2idx(ii,jj)] = 2.0*detM*ddetA[inode][ii]*ddetA[inode][jj];
      }
    }

    if(power > 0){
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          hquael[sym2idx(gdim*inode+ii,gdim*inode+jj)] =
            ( tdim*trapowdm2*det*det*( tra*htra[sym2idx(ii,jj)]
                                     + (tdim - 1)*dtra[inode][ii]*dtra[inode][jj] )
            - trapowdm1*det*( tra*hdet[sym2idx(ii,jj)]
                            + tdim*dtra[inode][ii]*ddet[inode][jj]
                            + tdim*dtra[inode][jj]*ddet[inode][ii])
            + 2.0*trapowd*ddet[inode][ii]*ddet[inode][jj]
            )/(det*det*det);
          hquael[sym2idx(gdim*inode+ii,gdim*inode+jj)] /= dpowd;
        }
      }
    }else{
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          hquael[sym2idx(gdim*inode+ii,gdim*inode+jj)]  =
            (-(tdim+1)*dtra[inode][jj]*(ddet[inode][ii]*tra - tdim*det*dtra[inode][ii])
            + tra*( hdet[sym2idx(ii,jj)]*tra + ddet[inode][ii]*dtra[inode][jj]
                  - tdim*ddet[inode][jj]*dtra[inode][ii]
                  - tdim*det*htra[sym2idx(ii,jj)])
            )/(trapowd*tra*tra);
          hquael[sym2idx(gdim*inode+ii,gdim*inode+jj)]  *= dpowd;
        }
      }
    }// if power > 0
  }


  for(int inod1 = 0; inod1 < nnode; inod1++){
    for(int inod2 = inod1+1; inod2 < nnode; inod2++){
      ftype sumDkDk = 0;
      for(int ii = 0; ii < tdim; ii++){
        sumDkDk += D_J_invJ0[inod1][ii]*D_J_invJ0[inod2][ii];
      }

      ftype htra[nhess]; // These are temp
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          htra[sym2idx(ii,jj)] = 2.0*met[sym2idx(ii,jj)]*sumDkDk;
        }
      }

      // some extradiagonal terms are added due to d^2 detA != 0 now
      ftype hdet[nhess2] = {0};//, hdet2[nhess2] = {0};
      const ftype detA = det_invtJ0_tJ; // Just for clarity
      if constexpr (gdim == 2){
        // The crossed-point terms (also crossed coordinate)
        hdet[sym2idx(0*gdim+0,1*gdim+1)] = detvec2<ftype>(D_J_invJ0[inod1], D_J_invJ0[inod2])*detA;
        hdet[sym2idx(0*gdim+1,1*gdim+0)] = - hdet[sym2idx(0*gdim+0,1*gdim+1)];//detvec2<ftype>(D_J_invJ0[inod1], D_J_invJ0[inod2])*detA;
      }else{
        hdet[sym2idx(0*gdim+0,1*gdim+1)] =
          detvec3<ftype>(D_J_invJ0[inod1], D_J_invJ0[inod2],  J_invJ0[2]      )*detA;
        hdet[sym2idx(0*gdim+0,1*gdim+2)] =
          detvec3<ftype>(D_J_invJ0[inod1],   J_invJ0[1]    ,  D_J_invJ0[inod2])*detA;
        hdet[sym2idx(0*gdim+1,1*gdim+2)] =
          detvec3<ftype>(  J_invJ0[0]    , D_J_invJ0[inod1],  D_J_invJ0[inod2])*detA;

        hdet[sym2idx(1*gdim+0,0*gdim+1)] = - hdet[sym2idx(0*gdim+0,1*gdim+1)];
        hdet[sym2idx(1*gdim+0,0*gdim+2)] = - hdet[sym2idx(0*gdim+0,1*gdim+2)];
        hdet[sym2idx(1*gdim+1,0*gdim+2)] = - hdet[sym2idx(0*gdim+1,1*gdim+2)];
      }

      for(int ii = 0; ii < gdim; ii++){
        for(int jj = 0; jj < gdim; jj++){
          hdet[sym2idx(0*gdim+ii,1*gdim+jj)] += ddetA[inod1][ii]*ddetA[inod2][jj];
          hdet[sym2idx(0*gdim+ii,1*gdim+jj)] *= 2.0*detM;
        }
      }

      if(power > 0){
        for(int ii = 0; ii < gdim; ii++){
          for(int jj = 0; jj < gdim; jj++){
            hquael[sym2idx(gdim*inod1+ii,gdim*inod2+jj)] =
              ( tdim*trapowdm2*det*det*( tra*htra[sym2idx(ii,jj)]
                                       + (tdim - 1)*dtra[inod1][ii]*dtra[inod2][jj] )
              - trapowdm1*det*( tra*hdet[sym2idx(0*gdim+ii,1*gdim+jj)]
                              + tdim*dtra[inod1][ii]*ddet[inod2][jj]
                              + tdim*dtra[inod2][jj]*ddet[inod1][ii])
              + 2.0*trapowd*ddet[inod1][ii]*ddet[inod2][jj]
              )/(det*det*det);
            hquael[sym2idx(gdim*inod1+ii,gdim*inod2+jj)] /= dpowd;
          }
        }
      }else{
        for(int ii = 0; ii < gdim; ii++){
          for(int jj = 0; jj < gdim; jj++){
            hquael[sym2idx(gdim*inod1+ii,gdim*inod2+jj)] =
              (-(tdim+1)*dtra[inod2][jj]*(ddet[inod1][ii]*tra - tdim*det*dtra[inod1][ii])
              + tra*( hdet[sym2idx(0*gdim+ii,1*gdim+jj)]*tra
                    + ddet[inod1][ii]*dtra[inod2][jj]
                    - tdim*ddet[inod2][jj]*dtra[inod1][ii]
                    - tdim*det*htra[sym2idx(ii,jj)])
              )/(trapowd*tra*tra);
            hquael[sym2idx(gdim*inod1+ii,gdim*inod2+jj)] *= dpowd;
          }
        }
      }
      //if(power > 0){
      //  for(int ii = 0; ii < gdim; ii++){
      //    for(int jj = 0; jj < gdim; jj++){
      //      hquael[sym2idx(gdim*inod1+ii,gdim*inod2+jj)] =
      //        (- hquael[sym2idx(gdim*inod1+ii,gdim*inod2+jj)]*quael
      //        + 2.0*dquael[gdim*inod1+ii]*dquael[gdim*inod2+jj]) /(quael*quael*quael);
      //    }
      //  }
      //}

    }
  }



  return quael;
}


#define INSTANTIATE(gdim,ideg,MFT_VAL,ASDEG_VAL,FTYPE)\
template FTYPE D_quafun_distortion< MFT_VAL , 2+gdim, 1+ideg, ASDEG_VAL, FTYPE>\
                  (Mesh< MFT_VAL > &msh, \
                  const int* ent2poi,\
                  const double*__restrict__ bary, \
                  FEBasis dofbas, \
                  DifVar idifmet, \
                  FTYPE*__restrict__ dquael, \
                  FTYPE*__restrict__ hquael);

#define EXPAND_TEMPLATE(z,gdim,SEQ) \
                  INSTANTIATE(gdim,BOOST_PP_SEQ_ELEM(0, SEQ),\
                                   BOOST_PP_SEQ_ELEM(1, SEQ),\
                                   BOOST_PP_SEQ_ELEM(2, SEQ),\
                                   BOOST_PP_SEQ_ELEM(3, SEQ))

#define REPEAT_GDIM(z,n,SEQ) BOOST_PP_REPEAT(2,EXPAND_TEMPLATE,(n)SEQ)
#define REPEAT_IDEG(r,SEQ)   BOOST_PP_REPEAT(METRIS_MAX_DEG,REPEAT_GDIM,SEQ)

#define ASDEG_SEQ (AsDeg::P1)(AsDeg::Pk)
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_IDEG,\
                             (MFT_SEQ)(ASDEG_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE




template <class MFT, int gdim, int ideg, AsDeg asdmet, typename ftype>
ftype D_metqua(Mesh<MFT> &msh, int ielem,
               FEBasis dofbas, DifVar idifmet,
               ftype*__restrict__ dquael, ftype*__restrict__ hquael,
               double difto){
  constexpr int tdim = gdim;
  int* ent2poi = tdim == 2 ? msh.fac2poi[ielem] : msh.tet2poi[ielem];
  return D_metqua<MFT,gdim,ideg,asdmet,ftype>(msh,ent2poi,
                                                dofbas,idifmet,
                                                dquael,hquael,difto);
}


template <class MFT, int gdim, int ideg, AsDeg asdmet, typename ftype>
ftype D_metqua(Mesh<MFT> &msh, const int* ent2poi,
              FEBasis dofbas, DifVar idifmet,
              ftype*__restrict__ dquael,
              ftype*__restrict__ hquael,
              double difto){

  static_assert(gdim == 2 || gdim == 3);
  constexpr int tdim = gdim;
  const int pnorm = msh.param->opt_pnorm;

  double bary[tdim+1];

  ftype qutet = 0;

  for(int ii = 0; ii < gdim; ii++) dquael[ii] = 0;

  if constexpr(asdmet == AsDeg::Pk && ideg > 1){

    constexpr int idegj = SMOO_DEGJ(ideg);
    constexpr int nnodj = tdim == 2 ? getnnod2(idegj) : getnnod3(idegj);

    constexpr int nnode = tdim == 1 ? getnnod1(ideg) :
                          tdim == 2 ? getnnod2(ideg) : getnnod3(ideg) ;

    constexpr auto ordelt = ORDELT(tdim);

    if(dquael != NULL){
      for(int ii = 0; ii < gdim*nnode; ii++) dquael[ii] = 0;
    }
    if(hquael != NULL){
      for(int ii = 0; ii < gdim*nnode; ii++){
        for(int jj = ii; jj < gdim*nnode; jj++){
          hquael[sym2idx(ii,jj)] = 0;
        }
      }
    }
    constexpr int nhess = (nnode*gdim*(nnode*gdim + 1))/2;
    ftype dqua0[nnode*gdim], hqua0[nhess];
    for(int iquad = 0; iquad < nnodj; iquad++){
      for(int ii = 0; ii < tdim + 1; ii++){
        bary[ii] = ordelt[idegj][iquad][ii]/((double) (idegj));
      }
      ftype qua0 = D_quafun_distortion<MFT,gdim,ideg,asdmet,ftype>(msh,ent2poi,bary,
                                                       dofbas,idifmet,
                                                       dqua0,hqua0);
      ftype powm1 = pow(qua0 - difto,pnorm-1);
      qutet += (qua0 - difto)*powm1/nnodj;

      int sg = 1;
      if(qua0 - difto < 0 && pnorm % 2 == 1) sg = -1;


      if(hquael != NULL){

        for(int ii = 0; ii < gdim*nnode; ii++){
          for(int jj = ii; jj < gdim*nnode; jj++){
            hquael[sym2idx(ii,jj)] +=
             sg*pnorm*hqua0[sym2idx(ii,jj)]*powm1;
          }
        }

        if(pnorm >= 2){
          ftype powm2 = pow(qua0 - difto,pnorm-2);
          for(int ii = 0; ii < gdim*nnode; ii++){
            for(int jj = ii; jj < gdim*nnode; jj++){
              hquael[sym2idx(ii,jj)] +=
                sg*pnorm*(pnorm - 1)
               *dqua0[ii]*dqua0[jj]*powm2;
            }
          }
        } // if(pnorm >= 2)

      } // if(hquael != NULL)



      for(int ii = 0; ii < gdim*nnode; ii++){
        dquael[ii] += sg*pnorm*dqua0[ii]*powm1/nnodj;
      }


    }
  }else{
    constexpr int nnode = tdim + 1;

    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0/(tdim  + 1);
    ftype qua0 = D_quafun_distortion<MFT,gdim,1,asdmet,ftype>(msh,ent2poi,bary,
                                                      dofbas,idifmet,
                                                      dquael,hquael);

    int sg = 1;
    if(qua0 - difto < 0 && pnorm % 2 == 1) sg = -1;


    ftype powm1 = pow(qua0 - difto,pnorm-1);
    qutet = sg*powm1*(qua0 - difto);


    if(hquael != NULL){
      for(int ii = 0; ii < gdim*nnode; ii++){
        for(int jj = ii; jj < gdim*nnode; jj++){
          hquael[sym2idx(ii,jj)] =
            sg*pnorm*hquael[sym2idx(ii,jj)]*powm1;
        }
      }
      if(pnorm >= 2){
        ftype powm2 = pow(qua0 - difto,pnorm-2);
        for(int ii = 0; ii < gdim*nnode; ii++){
          for(int jj = ii; jj < gdim*nnode; jj++){
            hquael[sym2idx(ii,jj)] +=
                  sg*pnorm*(pnorm - 1)
                 *dquael[ii]*dquael[jj]*powm2;
          }
        }
      }// endif pnorm
    }

    for(int ii = 0; ii < gdim*nnode; ii++){
      dquael[ii] = sg*pnorm*dquael[ii]*powm1;
    }

  }
  return qutet;
}

#define EXPAND_TEMPLATE(z,gdim,SEQ) \
                  INSTANTIATE(gdim,BOOST_PP_SEQ_ELEM(0, SEQ),\
                                   BOOST_PP_SEQ_ELEM(1, SEQ),\
                                   BOOST_PP_SEQ_ELEM(2, SEQ),\
                                   BOOST_PP_SEQ_ELEM(3, SEQ))
#define REPEAT_GDIM(z,n,SEQ) BOOST_PP_REPEAT(2,EXPAND_TEMPLATE,(n)SEQ)
#define REPEAT_IDEG(r,SEQ)   BOOST_PP_REPEAT(METRIS_MAX_DEG,REPEAT_GDIM,SEQ)

#define INSTANTIATE(gdim,ideg,MFT_VAL,ASDEG_VAL,FTYPE)\
template FTYPE D_metqua< MFT_VAL , 2+gdim, 1+ideg, ASDEG_VAL, FTYPE>\
                  (Mesh< MFT_VAL > &msh, \
                   int ielem, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, FTYPE*__restrict__ hquael, \
                   double difto);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_IDEG,(MFT_SEQ)(ASDEG_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE

#define INSTANTIATE(gdim,ideg,MFT_VAL,ASDEG_VAL,FTYPE)\
template FTYPE D_metqua< MFT_VAL , 2+gdim, 1+ideg, ASDEG_VAL, FTYPE>\
                  (Mesh< MFT_VAL > &msh, \
                   const int* ent2poi, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, FTYPE*__restrict__ hquael, \
                   double difto);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_IDEG,(MFT_SEQ)(ASDEG_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE

#undef EXPAND_TEMPLATE
#undef REPEAT_GDIM
#undef REPEAT_IDEG
#undef MFT_SEQ // note these two could go into headers
#undef ASDEG_SEQ

} // End namespace
