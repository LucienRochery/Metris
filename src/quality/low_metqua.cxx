// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "low_metqua.hxx"

#include "../ho_constants.hxx"
#include "../libs/SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "../aux_topo.hxx"
#include "../low_geo/misc.hxx"
#include "../low_geo/measure.hxx"
#include "../low_geo/normal.hxx"
#include "../low_geo/lenedg.hxx"
#include "../io_libmeshb.hxx"
#include "../Mesh/Mesh.hxx"
#include "../linalg/det.hxx"

#include "../utils/aux_pp_inc.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"

#include "objective_quadrature_value.hxx"
#include "simplex_quadrature.hxx"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Metris
{

  template <class MFT, int gdim, QuaFun iquaf>
  double metqua1_length(Mesh<MFT> &msh, const int *edg2pol)
  {
    static_assert(gdim == 2 || gdim == 3);
    static_assert(iquaf == QuaFun::SizeShape
                  || iquaf == QuaFun::StepDistance);

    double endpoint_sizes[2];
    const double length = getlenedg_geosz<MFT,gdim,1>(
        msh,edg2pol,endpoint_sizes);
    if(!(length > 0.0) || !std::isfinite(length)){
      return std::numeric_limits<double>::infinity();
    }

    const double log_length = std::log(length);
    if constexpr(iquaf == QuaFun::SizeShape){
      // In one dimension A = L_M^2 and the natural SizeShape quality is
      // (A + A^{-1})/2 = cosh(2 log(L_M)).
      const double direct_quality = std::cosh(2.0*log_length);
      const double raw_quality = msh.param->opt_power > 0
                               ? direct_quality
                               : 1.0/direct_quality;
      double error = std::abs(raw_quality - 1.0);
      const int pnorm = msh.param->opt_pnorm;
      if(pnorm == 2) error *= error;
      else if(pnorm > 2) error = std::pow(error,pnorm);
      return error;
    }else{
      double distance2;
      if(msh.param->step_distance_shape_volume){
        // tdim = 1 has no shape coordinate.  Its ShapeVolume coordinate is
        // (A - A^{-1})/2, with A = L_M^2.
        const double volume_coordinate = 2.0*std::sinh(2.0*log_length);
        distance2 = 0.25*volume_coordinate*volume_coordinate;
      }else{
        const double log_metric_ratio = 2.0*log_length;
        distance2 = log_metric_ratio*log_metric_ratio;
      }
      const double regularization =
          msh.param->step_distance_regularization;
      return std::pow(distance2 + regularization*regularization,
                      msh.param->objective_p/2.0)
           - std::pow(regularization,msh.param->objective_p);
    }
  }

  template <class MFT, int gdim, int tdim>
  double step_distance_element_target_weight(Mesh<MFT> &msh,
                                             AsDeg asdmet,
                                             int ientt)
  {
    (void)msh;
    (void)asdmet;
    (void)ientt;
    return 1.0;
  }

  template <class MFT, int gdim, int tdim>
  StepDistanceObjectiveState step_distance_global_objective_state(
      Mesh<MFT> &msh,
      AsDeg asdmsh,
      AsDeg asdmet)
  {
    StepDistanceObjectiveState state;
    METRIS_ENFORCE(msh.param->step_distance_cavity_target_average);
    const intAr2 &ent2poi = msh.ent2poi(tdim);
    for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
      if(isdeadent(ientt,ent2poi)) continue;
      const double quality = metqua<MFT,gdim,tdim,QuaFun::StepDistance>(
          msh,asdmsh,asdmet,ientt,1.0);
      state.numerator += quality;
      state.element_count++;
    }
    METRIS_ENFORCE(state.element_count > 0);
    state.target_weight = state.element_count;
    return state;
  }

  template <class MFT, int gdim, int tdim, QuaFun iquaf, typename ftype>
  ftype metqua(Mesh<MFT> &msh, AsDeg asdmsh, AsDeg asdmet,
               int ientt, ftype difto)
  {
    static_assert(tdim == 2 || tdim == 3);

    const intAr2 &ent2poi = msh.ent2poi(tdim);

    METRIS_ASSERT(!isdeadent(ientt, ent2poi));

    static_assert(gdim == 2 || gdim == 3);

    const int ideg = msh.curdeg;
    const int ideg_eff = asdmsh == AsDeg::P1 ? 1 : ideg;

    constexpr bool objective_driven
        = iquaf == QuaFun::SizeShape || iquaf == QuaFun::StepDistance;
    if constexpr(objective_driven){
      if constexpr(iquaf == QuaFun::StepDistance){
        // StepDistance remains P1 until its separate high-order phase.
        METRIS_ASSERT(ideg_eff == 1);
        METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);

        const SimplexQuadratureView<tdim> quadrature
            = get_objective_quadrature<tdim>(
                  msh.param->objective_quadrature_order);
        return integrate_objective_quadrature_value<
            MFT,gdim,tdim,1,iquaf,ftype>(
                msh,asdmsh,asdmet,ent2poi[ientt],quadrature);
      }else{
        // SizeShape value evaluation follows the actual geometric degree.
        // The runtime AsDeg choice still owns the explicit compatibility
        // contract: AsDeg::P1 selects the corner map, while AsDeg::Pk selects
        // the complete high-order map. Metric sampling remains governed by
        // asdmet and is independent of this geometry dispatch.
        METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
        const SimplexQuadratureView<tdim> quadrature
            = get_objective_quadrature<tdim>(
                  msh.param->objective_quadrature_order);
        if(ideg_eff == 1){
          return integrate_objective_quadrature_value<
              MFT,gdim,tdim,1,iquaf,ftype>(
                  msh,asdmsh,asdmet,ent2poi[ientt],quadrature);
        }
        ftype objective_value = ftype(0);
        bool degree_was_dispatched = false;
        CT_FOR0_INC(2,METRIS_MAX_DEG,geometry_degree){
          if(geometry_degree == ideg_eff){
            objective_value = integrate_objective_quadrature_value<
                MFT,gdim,tdim,geometry_degree,iquaf,ftype>(
                    msh,asdmsh,asdmet,ent2poi[ientt],quadrature);
            degree_was_dispatched = true;
          }
        }CT_FOR1(geometry_degree);
        METRIS_ENFORCE_MSG(
            degree_was_dispatched,
            "Unsupported SizeShape geometry degree {}",ideg_eff);
        return objective_value;
      }
    }

    double bary[tdim + 1];
    ftype qutet = 0;
    if (tdim == 1 && gdim >= 2)
      METRIS_THROW_MSG("TODO: TODO: Edge quality with normal dev")

    // Performance impact should be zero
    constexpr auto quafun_xi = get_quafun_xi<MFT, gdim, tdim, iquaf, ftype>();
    constexpr int nnmet = (gdim * (gdim + 1)) / 2;

    //// Compute normal at the nodes. This is then used to interpolate a normal
    //// within the element. Fewer EG_evaluate calls needed and more robust as
    //// (u,v) interpolation followed by evaluation is not necessarily very stable.
    // double norfld[getnnode(tdim,METRIS_MAX_DEG)][gdim];
    // if(do_nordev){
    //   double result[18];
    //   double *du = &result[3];
    //   double *dv = &result[6];
    //   const int nnode = getnnode(tdim, asdmsh == AsDeg::P1 ? 1 : ideg);
    //   const int iref = msh.fac2ref[ientt];
    //   const ego obj  = msh.cad2fac[iref];

    //  for(int inode = 0; inode < nnode; inode++){
    //    int ipoin = ent2poi(ientt, inode);
    //    int ibpoi = msh.poi2ebp(ent2poi[ientt], tdim, ientt, iref);
    //    METRIS_ASSERT(ibpoi >= 0 && ibpoi < msh.nbpoi);
    //    int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
    //    METRIS_ENFORCE_MSG(ierro == 0, "metqua0 EG_evaluate error " << ierro);
    //    vecprod(du,dv,norfld[inode]);
    //    if(normalize_vec<gdim>(norfld[inode])) METRIS_THROW_MSG( "Normal vanishes");
    //  }
    //}

    const int pnorm = msh.param->opt_pnorm;
    double nordev = 0;
    const bool do_nordev = tdim == 2 && gdim == 3
                        && msh.CAD()
                        && abs(msh.param->qua_surf_wt_normal)
                           > 1.0e-9*abs(msh.param->qua_surf_wt_quality);
    constexpr auto ordelt = ORDELT(tdim);
    const int nnode = getnnode(tdim,ideg_eff);

    // Accumulate normal error at the nodes (depending on asdmsh)
    if (do_nordev)
    {
      if (asdmsh == AsDeg::P1)
      {
        nordev = getnordev<1>(msh, ientt);
      }
      else
      {
        CT_FOR0_INC(1, METRIS_MAX_DEG, ideg)
        {
          if (ideg == msh.curdeg)
          {
            nordev = getnordev<ideg>(msh, ientt);
          }
        }
        CT_FOR1(ideg);
      }
    }

#ifdef TESTQUALITYALGO
    // The historical Classical quality integration assumes opt_pnorm == 1.
    // P1 objective paths return above and do not use this scalar transform.
    METRIS_ASSERT(pnorm == 1);
#endif

    if (ideg_eff > 1)
    {
      qutet = 0.0;
      nordev = 0.0;

      // const int idegj = SMOO_DEGJ(ideg);
      // const int nnodj = tdim == 2 ? getnnod2(idegj) : getnnod3(idegj);

      METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);

      // double metl[nnmet];

      for (int iquad = 0; iquad < nnode; iquad++)
      {
        for (int ii = 0; ii < tdim + 1; ii++)
        {
          bary[ii] = ordelt[ideg][iquad][ii] / ((double)(ideg));
        }
        const int ipoin = ent2poi(ientt, iquad);
        ftype qua0 = quafun_xi(msh, asdmet, asdmsh, ent2poi[ientt], bary, msh.met[ipoin]);
        // About 6x slower if MetSpace::Log: leave an assert here.
        // if(msh.met.getSpace() == MetSpace::Exp){
        //  qua0 = quafun_xi(msh,asdmet,asdmsh,ent2poi[ientt],bary,msh.met[ipoin]);
        //}else{
        //  for(int ii = 0; ii < nnmet; ii++) metl[ii] = msh.met(ipoin,ii);
        //  getexpmet_inp<gdim>(metl);
        //  qua0 = quafun_xi(msh,asdmet,asdmsh,ent2poi[ientt],bary,metl);
        //}

        // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
        // time to run pow() here even if pnorm = 2 or 1.
        qua0 = abs(qua0 - difto);
        if (pnorm == 2)
        {
          qua0 *= qua0;
        }
        else if (pnorm > 2)
        {
          qua0 = pow(qua0, pnorm);
        }
        qutet += qua0 / nnode;

        // qutet += pow(abs(qua0 - difto),pnorm) / nnode;
      }
    }
    else
    {

      // Classical P1 sampling: evaluate once at the barycenter. FE metric
      // fields interpolate there log-Euclideanly; analytical fields evaluate
      // the metric function directly at the mapped physical point.

      METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
      double met[nnmet];
      for (int ii = 0; ii < tdim + 1; ii++)
        bary[ii] = 1.0 / (tdim + 1);
      msh.met.getMetBary(asdmet,DifVar::None,MetSpace::Exp,
                         ent2poi[ientt],tdim,bary,met,NULL);
      qutet = quafun_xi(msh, AsDeg::P1, asdmet, ent2poi[ientt], bary, met);
      // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
      // time to run pow() here even if pnorm = 2 or 1.
      qutet = abs(qutet - difto);
#ifdef TESTQUALITYALGO
      double meas;
      isvalideltP1<gdim, tdim>(msh, ientt, NULL, &meas);
      qutet *= meas;
#ifdef INTQUALINRIEMSPACE
      double sqrtDetM = 1;
      if (tdim == 2)
      {
        sqrtDetM = sqrt(met[0] * met[2] - met[1] * met[1]);
      }
      else
      {
        const double m11 = met[0];
        const double m12 = met[1];
        const double m22 = met[2];
        const double m13 = met[3];
        const double m23 = met[4];
        const double m33 = met[5];

        sqrtDetM = sqrt(
            m11 * (m22 * m33 - m23 * m23) - m12 * (m12 * m33 - m13 * m23) + m13 * (m12 * m23 - m13 * m22));
      }
      qutet *= sqrtDetM;
#endif
#endif
      if (pnorm == 2)
      {
        qutet *= qutet;
      }
      else if (pnorm > 2)
      {
        qutet = pow(qutet, pnorm);
      }

    }

    if (do_nordev)
    {
      METRIS_ASSERT(msh.param->qua_surf_wt_quality >= 0);
      METRIS_ASSERT(msh.param->qua_surf_wt_normal >= 0);
      qutet = msh.param->qua_surf_wt_quality * qutet + msh.param->qua_surf_wt_normal * pow(nordev, pnorm); // for homogeneity
    }

    return qutet;
  }

#define EXPAND_TEMPLATE(r, SEQ)          \
  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ), \
              BOOST_PP_SEQ_ELEM(1, SEQ), \
              BOOST_PP_SEQ_ELEM(2, SEQ))
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical)
#define QUAFUN_SEQ (QuaFun::Distortion)(QuaFun::Unit)(QuaFun::SizeShape)(QuaFun::StepDistance)

#define INSTANTIATE(MFT_VAL, QUAFUN, FTYPE)                                              \
  template FTYPE metqua<MFT_VAL, 2, 2, QUAFUN, FTYPE>(Mesh<MFT_VAL> & msh, AsDeg, AsDeg, \
                                                      int ielem, FTYPE difto);           \
  template FTYPE metqua<MFT_VAL, 3, 2, QUAFUN, FTYPE>(Mesh<MFT_VAL> & msh, AsDeg, AsDeg, \
                                                      int ielem, FTYPE difto);           \
  template FTYPE metqua<MFT_VAL, 3, 3, QUAFUN, FTYPE>(Mesh<MFT_VAL> & msh, AsDeg, AsDeg, \
                                                      int ielem, FTYPE difto);
  BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE, (MFT_SEQ)(QUAFUN_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE

#undef QUAFUN_SEQ
#undef MFT_SEQ // note these two could go into headers
#undef EXPAND_TEMPLATE

#define INSTANTIATE_STEPDISTANCE_TARGET_WEIGHT(MFT_VAL)                    \
  template double step_distance_element_target_weight<MFT_VAL,2,2>(       \
      Mesh<MFT_VAL>&,AsDeg,int);                                           \
  template double step_distance_element_target_weight<MFT_VAL,3,2>(       \
      Mesh<MFT_VAL>&,AsDeg,int);                                           \
  template double step_distance_element_target_weight<MFT_VAL,3,3>(       \
      Mesh<MFT_VAL>&,AsDeg,int);                                           \
  template StepDistanceObjectiveState                                     \
  step_distance_global_objective_state<MFT_VAL,2,2>(                      \
      Mesh<MFT_VAL>&,AsDeg,AsDeg);                                        \
  template StepDistanceObjectiveState                                     \
  step_distance_global_objective_state<MFT_VAL,3,2>(                      \
      Mesh<MFT_VAL>&,AsDeg,AsDeg);                                        \
  template StepDistanceObjectiveState                                     \
  step_distance_global_objective_state<MFT_VAL,3,3>(                      \
      Mesh<MFT_VAL>&,AsDeg,AsDeg);
INSTANTIATE_STEPDISTANCE_TARGET_WEIGHT(MetricFieldFE)
INSTANTIATE_STEPDISTANCE_TARGET_WEIGHT(MetricFieldAnalytical)
#undef INSTANTIATE_STEPDISTANCE_TARGET_WEIGHT

#define INSTANTIATE_METQUA1_LENGTH(MFT_VAL, GDIM_VAL, QUAFUN_VAL)          \
  template double metqua1_length<MFT_VAL,GDIM_VAL,QUAFUN_VAL>(            \
      Mesh<MFT_VAL>&,const int*);
INSTANTIATE_METQUA1_LENGTH(MetricFieldFE,2,QuaFun::SizeShape)
INSTANTIATE_METQUA1_LENGTH(MetricFieldFE,2,QuaFun::StepDistance)
INSTANTIATE_METQUA1_LENGTH(MetricFieldAnalytical,2,QuaFun::SizeShape)
INSTANTIATE_METQUA1_LENGTH(MetricFieldAnalytical,2,QuaFun::StepDistance)
#undef INSTANTIATE_METQUA1_LENGTH

} // End namespace
