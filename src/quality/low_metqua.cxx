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
#include "../io_libmeshb.hxx"
#include "../Mesh/Mesh.hxx"
#include "../linalg/det.hxx"

#include "../utils/aux_pp_inc.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"

namespace Metris
{

  template <class MFT, int gdim, int tdim, QuaFun iquaf, typename ftype>
  ftype metqua(Mesh<MFT> &msh, AsDeg asdmsh, AsDeg asdmet,
               int ientt, ftype difto)
  {
    static_assert(tdim == 2 || tdim == 3);

    const intAr2 &ent2poi = msh.ent2poi(tdim);

    METRIS_ASSERT(!isdeadent(ientt, ent2poi));

    static_assert(gdim == 2 || gdim == 3);

    double bary[tdim + 1];

    const int pnorm = msh.param->opt_pnorm;

    ftype qutet = 0;
    double nordev = 0;
    bool do_nordev = tdim == 2 && gdim == 3 && msh.CAD() && abs(msh.param->qua_surf_wt_normal) > 1.0e-9 * abs(msh.param->qua_surf_wt_quality);
    if (tdim == 1 && gdim >= 2)
      METRIS_THROW_MSG("TODO: TODO: Edge quality with normal dev")

    // Performance impact should be zero
    constexpr auto quafun_xi = get_quafun_xi<MFT, gdim, tdim, iquaf, ftype>();
    constexpr auto ordelt = ORDELT(tdim);

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
    const int ideg = msh.curdeg;
    const int ideg_eff = asdmsh == AsDeg::P1 ? 1 : ideg;
    const int nnode = getnnode(tdim, ideg_eff);
    constexpr int nnmet = (gdim * (gdim + 1)) / 2;

#ifdef TESTQUALITYALGO
    // Assumptions for quality algo:
    METRIS_ASSERT(pnorm == 1);
    METRIS_ASSERT(ideg_eff == 1);
#endif

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

#if 0
      // Value doesn't matter. (P1)
      for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim  + 1);
      METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
      for(int iquad = 0; iquad < tdim + 1; iquad++){
        int ipoin = ent2poi(ientt, iquad);
        ftype qua0 = quafun_xi(msh,asdmet,asdmsh,ent2poi[ientt],bary,msh.met[ipoin]);
        // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
        // time to run pow() here even if pnorm = 2 or 1.
        qua0 = abs(qua0 - difto);
        if(pnorm == 2){
          qua0 *= qua0;
        }else if(pnorm > 2){
          qua0 = pow(qua0, pnorm);
        }
        qutet += qua0 / (tdim + 1);
      }

#elif 0
      for (int ii = 0; ii < tdim + 1; ii++)
        bary[ii] = 1.0 / (tdim + 1);
#ifndef NDEBUG
      try
      {
#endif
        qutet = quafun_xi(msh, asdmet, asdmsh, ent2poi[ientt], bary, NULL);
#ifndef NDEBUG
      }
      catch (const MetrisExcept &e)
      {
        printf("## metqua ent2pol {}\n",
               intAr1(getnnode(tdim, ideg), ent2poi[ientt]));
        throw(e);
      }
#endif
      // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
      // time to run pow() here even if pnorm = 2 or 1.
      qutet = abs(qutet - difto);
      if (pnorm == 2)
      {
        qutet *= qutet;
      }
      else if (pnorm > 2)
      {
        qutet = pow(qutet, pnorm);
      }

#elif defined(ONEPOINTQUAL)

      METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
      double met[nnmet];
      for (int ii = 0; ii < tdim + 1; ii++)
        bary[ii] = 1.0 / (tdim + 1);
      for (int jj = 0; jj < nnmet; jj++)
        met[jj] = 0;
      for (int ii = 0; ii < tdim + 1; ii++)
      {
        int ipoin = ent2poi(ientt, ii);
        for (int jj = 0; jj < nnmet; jj++)
        {
          met[jj] += msh.met(ipoin, jj) / (tdim + 1);
        }
      }
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

#elif defined(TDIM1POINTSQUAL)

      METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);

      auto sqrtDetM_of = [&](const double *metptr) -> double
      {
        if constexpr (tdim == 2)
        {
          const double m11 = metptr[0];
          const double m12 = metptr[1];
          const double m22 = metptr[2];
          const double det = m11 * m22 - m12 * m12;
          return det > 0 ? std::sqrt(det) : 0.0;
        }
        else
        { // tdim == 3
          const double m11 = metptr[0];
          const double m12 = metptr[1];
          const double m22 = metptr[2];
          const double m13 = metptr[3];
          const double m23 = metptr[4];
          const double m33 = metptr[5];

          const double det =
              m11 * (m22 * m33 - m23 * m23) - m12 * (m12 * m33 - m13 * m23) + m13 * (m12 * m23 - m13 * m22);

          return det > 0 ? std::sqrt(det) : 0.0;
        }
      };

      double meas;
      isvalideltP1<gdim, tdim>(msh, ientt, NULL, &meas);

      // first sample integrand at element vertices
      ftype integrand = 0;
      for (int iver = 0; iver < tdim + 1; iver++)
      {

        for (int ii = 0; ii < tdim + 1; ii++)
          bary[ii] = 0.0;
        bary[iver] = 1.0;

        const int ipoin = ent2poi(ientt, iver);
        double met_v[nnmet];
        for (int jj = 0; jj < nnmet; ++jj)
          met_v[jj] = msh.met(ipoin, jj);

        ftype qua0 = quafun_xi(msh, AsDeg::P1, asdmet, ent2poi[ientt], bary, met_v);

        qua0 = abs(qua0 - difto);
        if (pnorm == 2)
        {
          qua0 *= qua0;
        }
        else if (pnorm > 2)
        {
          qua0 = pow(qua0, pnorm);
        }

        integrand += qua0 * (ftype)sqrtDetM_of(met_v);
      }

      // now integrand at centroid
      double met_c[nnmet];
      for (int jj = 0; jj < nnmet; jj++)
        met_c[jj] = 0.0;
      for (int ii = 0; ii < tdim + 1; ii++)
        bary[ii] = 1.0 / (tdim + 1);

      // actually evaluate metric at barycenter
      double coordBary[tdim];
      for (int icoord = 0; icoord < tdim; icoord++)
      {
        coordBary[icoord] = 0;
        for (int ibary = 0; ibary < tdim + 1; ibary++)
        {
          coordBary[icoord] += bary[ibary] * msh.coord(ent2poi(ientt, ibary), icoord);
        }
      }
      if constexpr (std::is_same<MFT, MetricFieldAnalytical>::value)
      {
        msh.met.getMetPhys(DifVar::None, msh.met.getSpace(), coordBary, met_c, NULL);
      }
      else
      {
        METRIS_THROW_MSG("metric evaluation in low_metqua not implemented for MetricFieldFE");
      }

      ftype qua_c = quafun_xi(msh, AsDeg::P1, asdmet, ent2poi[ientt], bary, met_c);

      qua_c = abs(qua_c - difto);
      if (pnorm == 2)
      {
        qua_c *= qua_c;
      }
      else if (pnorm > 2)
      {
        qua_c = pow(qua_c, pnorm);
      }

      integrand += qua_c * (ftype)sqrtDetM_of(met_c);

      qutet = meas / 5. * integrand;

#elif defined(KEAST4QUAL)

      METRIS_ASSERT_MSG(tdim == 3, "Keast degree-4 rule implemented for tetrahedra only (tdim=3).");
      METRIS_ASSERT(gdim == 3);

      METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);

      auto sqrtDetM_of = [&](const double *metptr) -> double
      {
        const double m11 = metptr[0];
        const double m12 = metptr[1];
        const double m22 = metptr[2];
        const double m13 = metptr[3];
        const double m23 = metptr[4];
        const double m33 = metptr[5];

        const double det =
            m11 * (m22 * m33 - m23 * m23) - m12 * (m12 * m33 - m13 * m23) + m13 * (m12 * m23 - m13 * m22);

        return det > 0 ? std::sqrt(det) : 0.0;
      };

      double meas = 1.0;
      isvalideltP1<gdim, tdim>(msh, ientt, NULL, &meas);

      // Keast degree-4 (11-pt) rule on reference tetrahedron.
      // We renormalize weights by 6 so they sum to 1 (to keep the "meas * weighted average" convention).

      // constexpr double w0 = -0.013155555555555555;
      constexpr double w0 = 0.013155555555555555;
      constexpr double w1 = 0.007622222222222222;
      constexpr double w2 = 0.024888888888888888;

      constexpr double wn0 = 6.0 * w0;
      constexpr double wn1 = 6.0 * w1;
      constexpr double wn2 = 6.0 * w2;

      constexpr double a = 0.100596423833200785;
      constexpr double b = 1.0 - 3.0 * a;
      constexpr double bb = 0.0714285714285714285;
      constexpr double cc = 0.5 - bb;

      // helper: evaluate analytic metric at barycentric point
      auto eval_met_at_lam = [&](const double lam[4], double metq[nnmet])
      {
        // compute physical coordinates of quadrature point
        double coordQ[3] = {0.0, 0.0, 0.0};
        for (int iv = 0; iv < 4; iv++)
        {
          const int ip = ent2poi(ientt, iv);
          coordQ[0] += lam[iv] * msh.coord(ip, 0);
          coordQ[1] += lam[iv] * msh.coord(ip, 1);
          coordQ[2] += lam[iv] * msh.coord(ip, 2);
        }

        if constexpr (std::is_same<MFT, MetricFieldAnalytical>::value)
        {
          msh.met.getMetPhys(DifVar::None, msh.met.getSpace(), coordQ, metq, NULL);
        }
        else
        {
          METRIS_THROW_MSG("Keast: metric eval at quadrature points not implemented for MetricFieldFE");
        }
      };

      ftype sum = 0;

      // centroid point
      {
        const double lam[4] = {0.25, 0.25, 0.25, 0.25};
        for (int ii = 0; ii < 4; ii++)
          bary[ii] = lam[ii];

        double metq[nnmet];
        eval_met_at_lam(lam, metq);

        ftype q = quafun_xi(msh, AsDeg::P1, asdmet, ent2poi[ientt], bary, metq);
        q = abs(q - difto); // pnorm==1

        sum += (ftype)wn0 * q * (ftype)sqrtDetM_of(metq);
      }

      // 4-point orbit: permutations of (b,a,a,a)
      for (int ibig = 0; ibig < 4; ibig++)
      {
        double lam[4] = {a, a, a, a};
        lam[ibig] = b;

        for (int ii = 0; ii < 4; ii++)
          bary[ii] = lam[ii];

        double metq[nnmet];
        eval_met_at_lam(lam, metq);

        ftype q = quafun_xi(msh, AsDeg::P1, asdmet, ent2poi[ientt], bary, metq);
        q = abs(q - difto);

        sum += (ftype)wn1 * q * (ftype)sqrtDetM_of(metq);
      }

      // 6-point orbit: permutations of (bb,bb,cc,cc)
      // Choose which 2 indices get bb; the rest get cc. There are C(4,2)=6 choices.
      for (int i0 = 0; i0 < 4; i0++)
      {
        for (int i1 = i0 + 1; i1 < 4; i1++)
        {
          double lam[4] = {cc, cc, cc, cc};
          lam[i0] = bb;
          lam[i1] = bb;

          for (int ii = 0; ii < 4; ++ii)
            bary[ii] = lam[ii];

          double metq[nnmet];
          eval_met_at_lam(lam, metq);

          ftype q = quafun_xi(msh, AsDeg::P1, asdmet, ent2poi[ientt], bary, metq);
          q = abs(q - difto);

          sum += (ftype)wn2 * q * (ftype)sqrtDetM_of(metq);
        }
      }

      qutet = (ftype)meas * sum;

#else
      static_assert(tdim == 3, "No rule for quality integration defined!");
#endif
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
#define QUAFUN_SEQ (QuaFun::Distortion)(QuaFun::Unit)(QuaFun::SizeShape)

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

} // End namespace
