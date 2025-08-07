//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "measure.hxx"
#include "misc.hxx"
#include "normal.hxx"

#include "../linalg/det.hxx"
#include "../utils/mprintf.hxx"

#include "../Mesh/MeshMetric.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

namespace Metris{




template<int gdim> 
double getmeasentP1(const int *ent2pol, const dblAr2& coord){
  static_assert(gdim >= 1 || gdim <= 3);
  if constexpr(gdim == 1){
    return coord[ent2pol[1]][0] - coord[ent2pol[0]][0];
  }else if(gdim == 2){
    return detvdif2(coord[ent2pol[1]],coord[ent2pol[0]]
                    ,coord[ent2pol[2]],coord[ent2pol[0]])/2;
  }else if(gdim == 3){
    return detvdif3(coord[ent2pol[1]],coord[ent2pol[0]]
                    ,coord[ent2pol[2]],coord[ent2pol[0]]
                    ,coord[ent2pol[3]],coord[ent2pol[0]])/6;
  }
}

template double getmeasentP1<1>(const int *ent2pol, const dblAr2 &coord);
template double getmeasentP1<2>(const int *ent2pol, const dblAr2 &coord);
template double getmeasentP1<3>(const int *ent2pol, const dblAr2 &coord);



template <int gdim, int tdim>
bool isvalideltP1(const MeshBase&__restrict__ msh, int ientt, const double*__restrict__ norref, 
                  double*__restrict__ meas, double nordev_tol){
  bool iflat;
  double meas1 = getmeasentP1<gdim,tdim>(msh, ientt, norref, &iflat, nordev_tol);
  if(meas != NULL) *meas = meas1;
  return !iflat;
}

template bool isvalideltP1<2,2>(const MeshBase&__restrict__ msh, int ientt, const double*__restrict__ norref, 
                                double*__restrict__ meas, double nordev_tol);
template bool isvalideltP1<3,2>(const MeshBase&__restrict__ msh, int ientt, const double*__restrict__ norref, 
                                double*__restrict__ meas, double nordev_tol);
template bool isvalideltP1<3,3>(const MeshBase&__restrict__ msh, int ientt, const double*__restrict__ norref, 
                                double*__restrict__ meas, double nordev_tol);


template <int gdim, int tdim>
bool isvalideltP1(const MeshBase&__restrict__ msh, int ientt, const int*__restrict__ nod2bpo, 
                  const double*__restrict__ norref, double*__restrict__ meas, double nordev_tol){
  bool iflat;
  double meas1 = getmeasentP1<gdim,tdim>(msh, msh.ent2poi(tdim)[ientt], nod2bpo, norref, &iflat, nordev_tol);
  if(meas != NULL) *meas = meas1;
  return !iflat;
}

template bool isvalideltP1<2,2>(const MeshBase&__restrict__ msh, int ientt, const int*__restrict__ nod2bpo, 
                                const double*__restrict__ norref, double*__restrict__ meas, double nordev_tol);
template bool isvalideltP1<3,2>(const MeshBase&__restrict__ msh, int ientt, const int*__restrict__ nod2bpo, 
                                const double*__restrict__ norref, double*__restrict__ meas, double nordev_tol);
template bool isvalideltP1<3,3>(const MeshBase&__restrict__ msh, int ientt, const int*__restrict__ nod2bpo, 
                                const double*__restrict__ norref, double*__restrict__ meas, double nordev_tol);

template <int gdim, int tdim>
bool isvalideltP1(const MeshBase&__restrict__ msh, const int*__restrict__ ent2pol, const int*__restrict__ nod2bpo, 
                  const double*__restrict__ norref, double*__restrict__ meas, double nordev_tol){
  bool iflat;
  double meas1 = getmeasentP1<gdim,tdim>(msh, ent2pol, nod2bpo, norref, &iflat, nordev_tol);
  if(meas != NULL) *meas = meas1;
  return !iflat;
}

template bool isvalideltP1<2,2>(const MeshBase&__restrict__ msh, const int*__restrict__ ent2pol, const int*__restrict__ nod2bpo, 
                                const double*__restrict__ norref, double*__restrict__ meas, double nordev_tol);
template bool isvalideltP1<3,2>(const MeshBase&__restrict__ msh, const int*__restrict__ ent2pol, const int*__restrict__ nod2bpo, 
                                const double*__restrict__ norref, double*__restrict__ meas, double nordev_tol);
template bool isvalideltP1<3,3>(const MeshBase&__restrict__ msh, const int*__restrict__ ent2pol, const int*__restrict__ nod2bpo, 
                                const double*__restrict__ norref, double*__restrict__ meas, double nordev_tol);


template<int gdim, int tdim>
double getmeasentP1(const MeshBase &msh, int ientt, 
                    const double* norref, bool*__restrict__ iflat, double nordev_tol){
  int nod2bpo[3];
  if constexpr(tdim == 2 && gdim == 3){
    nod2bpo[0] = msh.poi2ebp(msh.fac2poi(ientt,0), 2, ientt, msh.fac2ref[ientt]);
    nod2bpo[1] = msh.poi2ebp(msh.fac2poi(ientt,1), 2, ientt, msh.fac2ref[ientt]);
    nod2bpo[2] = msh.poi2ebp(msh.fac2poi(ientt,2), 2, ientt, msh.fac2ref[ientt]);
  }
  return getmeasentP1<gdim,tdim>(msh, msh.ent2poi(tdim)[ientt], nod2bpo, norref, iflat, nordev_tol);
}
template double getmeasentP1<2,2>(const MeshBase&__restrict__ msh, int ientt, 
                                  const double*__restrict__ norref, bool*__restrict__ iflat, double nordev_tol);
template double getmeasentP1<3,2>(const MeshBase&__restrict__ msh, int ientt, 
                                  const double*__restrict__ norref, bool*__restrict__ iflat, double nordev_tol);
template double getmeasentP1<3,3>(const MeshBase&__restrict__ msh, int ientt, 
                                  const double*__restrict__ norref, bool*__restrict__ iflat, double nordev_tol);

// This variant returns whether above or below specified tolerance
// nrmal only required if tdim == 2 and gdim == 3 (surface), can be NULL otherwise
// norCAD can be computed discretely, it is just a reference normal pointing inwards
template<int gdim, int tdim>
double getmeasentP1(const MeshBase&__restrict__ msh, const int*__restrict__ ent2pol, 
                    const int*__restrict__ nod2bpo,
                    const double*__restrict__ norref, bool*__restrict__ iflat, double nordev_tol){
  GETVDEPTH(msh.param);
  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim);

  //const intAr2 &ent2poi = msh.ent2poi(tdim);

  double fac, det;
  if constexpr(tdim == 2){ 

    double nrm1 = geterrl2<gdim>(msh.coord[ent2pol[0]],msh.coord[ent2pol[1]]);
    double nrm2 = geterrl2<gdim>(msh.coord[ent2pol[0]],msh.coord[ent2pol[2]]);
    double nrm3 = geterrl2<gdim>(msh.coord[ent2pol[1]],msh.coord[ent2pol[2]]);

    fac = 2*std::cbrt(nrm1*nrm2*nrm3); // cubic root, homo to h^2

    if constexpr(gdim == 2){
      det = detvdif2(msh.coord[ent2pol[1]],msh.coord[ent2pol[0]],
                     msh.coord[ent2pol[2]],msh.coord[ent2pol[0]]);
    }else{
      // Measure of the face projected in the plane norCAD ^ orth. Could be zero
      // Notice there exists rotation R st edges l1, l2 verify 
      // l1^(flat) = Rl1 = (0 l1^(2D)), 
      // l2^(flat) = Rl2 = (0 l2^(2D))
      // Furthermore, notice that 
      // l1^(flat) x l2^(flat) = R(l1 x l2)
      // where x is the vector product. 
      // Now the vector product of the "flattened" edges is simply 
      // l1^(flat) x l2^(flat) = (det(l1^(2D) l2^(2D)) 0 0)
      // Thus we simply replace the 2D determinant with the norm of the normal 
      double norfac[3];

      getnorfacP1(ent2pol,msh.coord,norfac);
      double nrm = getnrml2<3>(norfac);
      if(nrm < Constants::vecNrmTol){
        *iflat = true;
        return 0;
      }

      
      double norCAD[3];
      if(norref == NULL){
        norref = norCAD;
        if(msh.CAD()){
          int ierro = getnorfacCAD(msh,ent2pol,nod2bpo,norCAD);
          METRIS_ASSERT_MSG(ierro == 0, "getnorfacCAD failed ierro = " << ierro);
        }else{
          getnorfacP1(ent2pol,msh.coord,norCAD);
        }
      }
      
      double nrm2 = getnrml2<3>(norref);
      if(nrm2 < Constants::vecNrmTol){
        *iflat = true;
        return 0;
      }

      // norfac is l1 x l2 is already homo h^2 despite norref O(1)
      det = getprdl2<3>(norfac,norref)/sqrt(nrm2);

      // Additionally, check normal deviation
      if(nordev_tol >= 0){
        double nordev = getnordev<1>(msh,ent2pol,nod2bpo,norfac);
        CPRINTF1(" - getmeasentP1 got nordev = {}, tol = {}\n",nordev,nordev_tol);
        if(nordev > nordev_tol){
          //printf("## DEBUG nordev = {} > tol = {}\n",nordev,nordev_tol);
          //printf("WAIT HERE nordev \n");
          //wait();
          *iflat = true;
          return det/2;
        }
      }

    }
    det /= 2;

  }else if(tdim == 3){

    double nrm1 = geterrl2<gdim>(msh.coord[ent2pol[0]],msh.coord[ent2pol[1]]);
    double nrm2 = geterrl2<gdim>(msh.coord[ent2pol[0]],msh.coord[ent2pol[2]]);
    double nrm3 = geterrl2<gdim>(msh.coord[ent2pol[0]],msh.coord[ent2pol[3]]);
    double nrm4 = geterrl2<gdim>(msh.coord[ent2pol[1]],msh.coord[ent2pol[2]]);
    double nrm5 = geterrl2<gdim>(msh.coord[ent2pol[1]],msh.coord[ent2pol[3]]);
    double nrm6 = geterrl2<gdim>(msh.coord[ent2pol[2]],msh.coord[ent2pol[3]]);
    // full prod is homo h^12; det only h^3
    fac = 6*sqrt(sqrt(nrm1*nrm2*nrm3*nrm4*nrm5*nrm6));

    det = detvdif3(msh.coord[ent2pol[1]],msh.coord[ent2pol[0]],
                   msh.coord[ent2pol[2]],msh.coord[ent2pol[0]],
                   msh.coord[ent2pol[3]],msh.coord[ent2pol[0]]);
    det /= 6; 

  } 
  *iflat = (det < msh.param->vtol * fac) || fac < 1.0e-16;
  //CPRINTF1("-- END getmeasentP1 got meas = {}, iflat = {}\n",det,*iflat);
  return det;
}

template double getmeasentP1<2,2>(const MeshBase&__restrict__ msh, const int*__restrict__ ent2pol,
                                  const int*__restrict__ nod2bpo, 
                                  const double*__restrict__ norref, bool*__restrict__ iflat, double nordev_tol);
template double getmeasentP1<3,2>(const MeshBase&__restrict__ msh, const int*__restrict__ ent2pol,
                                  const int*__restrict__ nod2bpo, 
                                  const double*__restrict__ norref, bool*__restrict__ iflat, double nordev_tol);
template double getmeasentP1<3,3>(const MeshBase&__restrict__ msh, const int*__restrict__ ent2pol,
                                  const int*__restrict__ nod2bpo, 
                                  const double*__restrict__ norref, bool*__restrict__ iflat, double nordev_tol);



template <class MFT, int gdim, int ideg>
double getmeasent_aniso(const MeshMetric<MFT> &msh, int ientt){

  constexpr int tdim = gdim;
  constexpr int nnode = getnnode(tdim, ideg);
  constexpr auto eval = tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;
  constexpr auto ordent = ORDELT(tdim);
  const intAr2 &ent2poi = tdim == 3 ? msh.tet2poi : msh.fac2poi;

  double bary[tdim+1];
  double dum[gdim],jmat[tdim*gdim];
  const double* metl;

  const double dq = 1.0 / nnode;
  double volM = 0;
  if constexpr (ideg == 1){
    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1 / (tdim + 1.0);
    eval(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::Bary,DifVar::None,
         bary, dum, jmat, NULL);
  }
  for(int inode = 0; inode < nnode; inode++){
    for(int ii = 0; ii < tdim + 1; ii++)
      bary[ii] = ordent[ideg][inode][ii] / (double) ideg;
      
    int ipoin = ent2poi(ientt, inode);

    metl = msh.met[ipoin];

    // and Jacobian 
    if constexpr(ideg > 1){
      eval(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::Bary,DifVar::None,
           bary, dum, jmat, NULL);
    }

    double detJ = detmat<gdim>(jmat) / tdim / (tdim - 1); // Works in 2D and 3D
    // Determinant of M^{1/2}
    double det12;
    if(msh.met.getSpace() == MetSpace::Exp){
      det12 = sqrt(detsym2<gdim>(metl));
    }else{ // In log format, simply the trace
      det12 = metl[0] + metl[2];
      if constexpr(gdim == 3) det12 += metl[5];
      det12 = exp(det12/2);
    } 
    volM += det12 * detJ * dq ;
  }

  return volM;
}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double getmeasent_aniso<MetricFieldAnalytical,2,n>(const MeshMetric<MetricFieldAnalytical> &msh, int ientt);\
template double getmeasent_aniso<MetricFieldFE        ,2,n>(const MeshMetric<MetricFieldFE        > &msh, int ientt);\
template double getmeasent_aniso<MetricFieldAnalytical,3,n>(const MeshMetric<MetricFieldAnalytical> &msh, int ientt);\
template double getmeasent_aniso<MetricFieldFE        ,3,n>(const MeshMetric<MetricFieldFE        > &msh, int ientt);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




template <int gdim>
void getmeasentP1grad(const int *ent2pol, const dblAr2& coord, int idof, double *grad){
  if constexpr(gdim == 3){
    METRIS_THROW_MSG(TODOExcept(),"Implement getmeasentP1grad with idim == 3 (use subdetvec)");
  }else{

    int inxt1 = (idof + 1) % 3;
    int inxt2 = (inxt1 +1) % 3;

    const double* coop1 = coord[ent2pol[inxt1]];
    const double* coop2 = coord[ent2pol[inxt2]];
    grad[0] =   coop1[1] - coop2[1];
    grad[1] = -(coop1[0] - coop2[0]);
  }
}
template void getmeasentP1grad<2>(const int *ent2pol, const dblAr2& coord, int idof, double *grad);
template void getmeasentP1grad<3>(const int *ent2pol, const dblAr2& coord, int idof, double *grad);





}// namespace Metris