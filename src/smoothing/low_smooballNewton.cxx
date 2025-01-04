//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../smoothing/low_smooballNewton.hxx"


namespace Metris{



template<class MetricFieldType, int idim, int ideg>
int smooballNewton(Mesh<MetricFieldType>& msh, int ipoin, int nball, const int*__restrict__ lball, double*__restrict__ qball, 
                   double*__restrict__ qavg0, double*__restrict__ qmax0, double*__restrict__ qavg1, double*__restrict__ qmax1, 
                   double qpower, int qpnorm, double difto, double maxwt, int inorm, int iverb, int ithread){

  constexpr int nnmet = (idim*(idim+1))/2;
  constexpr auto getverent = idim == 2 ? getverfac<1> : getvertet<1>;
  const double vtol = Defaults::vtol;
  const int ipower = MAX(inorm, 1);
  intAr2& ent2poi = idim == 2 ? msh.fac2poi : msh.tet2poi;


  double coor0[idim];
  double met0[nnmet];
  for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord[ipoin][ii];
  for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met[ipoin][ii];
 
  *qavg0 = 0;
  *qmax0 = -1.0e30;

  double coonw[idim] = {0};

  for(int iball = 0; iball < nball; iball++){
    int ientt = lball[iball];
    double quael;
    metqua<MetricFieldType,idim,idim,ideg,AsDeg::Pk,double>(msh,ientt,qpower,&quael,qpnorm,difto);
    qball[iball] = quael;

    *qavg0 += pow(quael,ipower);
    *qmax0  = MAX(quael,*qmax0);
  }

    //iflgp = iflag;
    //optim_newton_drivertype(idim    ,
    //                        &xcur[1] ,&fcur  ,rhs   ,hess ,
    //                        xtol ,stpmin, 1 , 
    //                        wlfc1,wlfc2 ,ratnew ,
    //                        &niter,maxit ,iprt > 1 ? 2 : 0 ,
    //                        &iflag,&ihess , 
    //                        nrwrk,rwrkN ,
    //                        niwrk,iwrkN ,
    //                        xopt ,&fopt ,&ierro);



}





} // end namespace
