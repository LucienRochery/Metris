//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#define BOOST_TEST_MODULE MyTest 


#include <egads.h>

#include <random>
#include <string>
#include <regex>

#include "common_setup.hxx"
#include "../src/utils/aux_pp_inc.hxx"
#include "../src/Localization/low_localization.hxx"
#include "../src/low_ccoef.hxx"
#include "../src/linalg/invmat.hxx"

#include <boost/hana.hpp> 
#include <nlopt.hpp>

namespace hana = boost::hana;
using namespace hana::literals;

namespace utf = boost::unit_test;

namespace Metris {

typedef MetricFieldAnalytical MFT;




#if 0
template<int gdim, int ideg>
void invevaltest(const MeshBase &msh,
                 const int* ent2pol,
                 const double* coor0, 
                 double* __restrict__ coopr, 
                 double* __restrict__ bary,
                 double tol){

  const double qstep = 0.95;
  
  constexpr int nhess = (gdim*(gdim+1))/2;
  static_assert(gdim == 1 || gdim == 2 || gdim == 3);
  constexpr int tdim = gdim;
  constexpr auto evalf = tdim == 1 ? eval1<gdim,ideg> : 
                         tdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;
  const int iverb = msh.param->iverb;

  inventP1<gdim>(ent2pol, msh.coord, coor0, bary);

  double gcur[gdim], hess[nhess];

  //double wlfc1 = 0.3;
  double wlfc1 = 1e-4;
  double wlfc2 = 0.9;

  for(int niter = 0; niter < 20; niter++){

    double fun = invevalfun<gdim,ideg>(msh, ent2pol, msh.coord, coor0, bary, 1, coopr, gcur, hess);

    double desc[gdim+1];
    symXvec<gdim>(hess,gcur,&desc[1]);
    desc[0] = 0;
    for(int ii = 1; ii < gdim+1; ii++){
      desc[ii] *= -1;
      desc[0] -= desc[ii];
    }

    double eigval[gdim], eigvec[gdim][gdim];
    geteigsym<gdim>(hess, eigval, eigvec[0]);
    //printf("Hessian %23.15e %23.15e %23.15e \n",hess[0],hess[1],hess[2]);
    printf("Hess eigval ");
    dblAr1(gdim,eigval).print();
    printf("-- eval f %f bary  ",fun);
    dblAr1(gdim+1,bary).print();

    double bar0[gdim+1];
    for(int ii = 0; ii < gdim+1; ii++) bar0[ii] = bary[ii];

    double dtpgd0 = getprdl2<gdim>(&desc[1],gcur);

    double step = 1.0;
    printf("-- LS start at %f dtpgd0 = %f \n",step,dtpgd0);

    double stpopt = -1;
    while(abs(step) >= 1.0e-6){
      double dum[gdim];
      for(int ii = 0; ii < gdim+1; ii++){
        bary[ii] = bar0[ii] + step * desc[ii];
      }
      double funLS = invevalfun<gdim,ideg>(msh, ent2pol, msh.coord, coor0, bary, 0, coopr, dum, NULL);
      double dtpgd1 = getprdl2<gdim>(&desc[1],dum);
      double nrm1 = getnrml2<gdim>(dum);
      bool cdt1 = (funLS <= fun + wlfc1*step*dtpgd0);
      bool cdt2 = (dtpgd1  >= wlfc2*dtpgd0);
      bool cdt3 = (abs(dtpgd1) <= wlfc2*abs(dtpgd0));
      // dtpgd0 is negative 
      double c1lim = (funLS - fun) / (step * dtpgd0);
      double c2lim = dtpgd1 / dtpgd0;
      printf(" - %12.7e : ratio %9.4e cdt %d %d %d c1 mx %9.4e c2 %9.4e nxt g nrm %e\n"
        ,step,funLS/fun,cdt1,cdt2,cdt3,c1lim,c2lim,nrm1);
      step *= qstep; 
      if(cdt1 && cdt3 && stpopt < 0){
        printf(" -> picked this one\n");
        stpopt = step;
        // but we carry on
      }

      //if(cdt1 && cdt3){
      //  stpopt = step;
      //  break;
      //}
    }
    wait();

    if(stpopt < 0){
      printf("Failed LS\n");
      break;
    }

    for(int ii = 0; ii < gdim+1; ii++){
      bary[ii] = bar0[ii] + stpopt * desc[ii];
    }


  }

  return;
}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void invevaltest<2,n>(const MeshBase &msh,\
                                 const int* ent2pol,\
                                 const double* coor0, \
                                 double* __restrict__ coopr, \
                                 double* __restrict__ bary,\
                                 double tol);\
template void invevaltest<3,n>(const MeshBase &msh,\
                 const int* ent2pol,\
                 const double* coor0, \
                 double* __restrict__ coopr, \
                 double* __restrict__ bary,\
                 double tol);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()
#endif



template<int gdim,int ideg>
double invevalfun_nlointf1(const std::vector<double> &x, std::vector<double> &grad, void* f_data){
  invevalfun_data *mydata = (invevalfun_data*)(f_data);
  double buf[gdim];
  double bary[gdim+1];
  for(int ii = 0; ii < gdim; ii++) bary[1+ii] = x[ii];
  bary[0] = 1;
  for(int ii = 1; ii < gdim + 1; ii++) bary[0] -= bary[ii];
  double val = invevalfun<gdim,ideg>(mydata->msh,mydata->ent2pol,mydata->coord,mydata->coor0,bary, 0, mydata->coopr, buf, NULL);
  //printf(" x = %f %f f = %f g %f %f req %d \n",x[0],x[1],val,buf[0],buf[1],grad.size()>0);
  if(grad.size() > 0){
    for(int ii = 0; ii < gdim; ii++) grad[ii] = buf[ii];
  }
  return val;
}
#define BOOST_PP_LOCAL_MACRO(n)\
template double invevalfun_nlointf1<2,n>(const std::vector<double> &x, std::vector<double> &grad, void* f_data);\
template double invevalfun_nlointf1<3,n>(const std::vector<double> &x, std::vector<double> &grad, void* f_data);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


// Using high level nlopt function
template<int gdim, int ideg> 
int inveval_nloptD(MeshBase &msh, int ientt, const double*__restrict__ coor0, 
                   double*__restrict__ coopr, double*__restrict__ bary,
                   double tol){

  intAr2& ent2poi = msh.ent2poi(gdim);
  constexpr auto evalf = gdim == 1 ? eval1<gdim,ideg> : 
                         gdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;

  inventP1<gdim>(ent2poi[ientt], msh.coord, coor0, bary);

  // NLOPT_LD_TNEWTON_PRECOND_RESTART.
  // Second, simplified versions NLOPT_LD_TNEWTON_PRECOND (same without restarting), 
  // NLOPT_LD_TNEWTON_RESTART (same without preconditioning), 
  // and NLOPT_LD_TNEWTON (same without restarting or preconditioning).
  //for(auto algo : {nlopt::LD_TNEWTON_PRECOND_RESTART,
  //                 nlopt::LD_TNEWTON_RESTART,
  //                 nlopt::LD_TNEWTON}){
  auto algo = nlopt::LD_TNEWTON_PRECOND_RESTART;
    try{
      nlopt::opt optimizer(algo, gdim);


      invevalfun_data mydata(msh,ent2poi[ientt],msh.coord,coor0,coopr);

      optimizer.set_min_objective(invevalfun_nlointf1<gdim,ideg>, &mydata);

      std::vector<double> x(gdim);
      for(int ii = 0; ii < gdim; ii++) x[ii] = bary[ii+1];
      double opt_f;
      nlopt::result res = optimizer.optimize(x, opt_f);

      bary[0] = 1.0;
      for(int ii = 0; ii < gdim; ii++){
        bary[0] -= x[ii];
        bary[ii+1] = x[ii];
      }

      evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
            DifVar::None,bary,coopr,NULL,NULL);

    }catch(const std::exception& e){
      return 2;
    }
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
template int inveval_nloptD<2,n>(MeshBase &msh, int ientt, const double*__restrict__ coor0,\
                    double*__restrict__ coopr, double*__restrict__ bary,\
                    double tol);\
template int inveval_nloptD<3,n>(MeshBase &msh, int ientt, const double*__restrict__ coor0,\
                    double*__restrict__ coopr, double*__restrict__ bary,\
                    double tol);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




BOOST_AUTO_TEST_CASE(test_inveval) 
{
 
// Case invevalP2 To loc:    6.492172113352532e-01   -1.913581807459048e-02
  std::vector<std::string> meshes = {
    //"../cases/square_adap.p2.meshb",
    "../cases/2D/square.circmet.50.curved.meshb",
    "../cases/2D/square.circmet.5k.curved.meshb",
    "../cases/invevalP2.mesh",
    //"../cases/invevalP2_2",
  };

  const int nsamp = 100;
  dblAr2 bar1(nsamp,2),bar2(nsamp,3),bar3(nsamp,4);
  dblAr2 bar1_out(nsamp,2),bar2_out(nsamp,3),bar3_out(nsamp,4);
  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(0);
  for(int ii = 0; ii < nsamp; ii++){
    double sum1 = 0, sum2 = 0, sum3 = 0;
    do{
      for(int jj = 0; jj < 4; jj++){
        if(jj < 2) bar1(ii,jj) = unif(rng);
        if(jj < 3) bar2(ii,jj) = unif(rng);
                   bar3(ii,jj) = unif(rng);
        if(jj < 2) sum1 += bar1(ii,jj);
        if(jj < 3) sum2 += bar2(ii,jj);
                   sum3 += bar3(ii,jj);
      }
    }while(abs(sum1) < 1.0e-16 || abs(sum2) < 1.0e-16 || abs(sum3) < 1.0e-16);

    for(int jj = 0; jj < 4; jj++){
      if(jj < 2) bar1(ii,jj) /= sum1;
      if(jj < 3) bar2(ii,jj) /= sum2;
                 bar3(ii,jj) /= sum3;
    }
  }
  std::uniform_real_distribution<double> unif_out(-100.0,100.0);
  std::default_random_engine rng_out(0);
  for(int ii = 0; ii < nsamp; ii++){
    double sum1 = 0, sum2 = 0, sum3 = 0;
    do{
      for(int jj = 0; jj < 4; jj++){
        if(jj < 2) bar1_out(ii,jj) = unif_out(rng_out);
        if(jj < 3) bar2_out(ii,jj) = unif_out(rng_out);
                   bar3_out(ii,jj) = unif_out(rng_out);
        if(jj < 2) sum1 += bar1_out(ii,jj);
        if(jj < 3) sum2 += bar2_out(ii,jj);
                   sum3 += bar3_out(ii,jj);
      }
    }while(abs(sum1) < 1.0e-16 || abs(sum2) < 1.0e-16 || abs(sum3) < 1.0e-16);

    for(int jj = 0; jj < 4; jj++){
      if(jj < 2) bar1_out(ii,jj) /= sum1;
      if(jj < 3) bar2_out(ii,jj) /= sum2;
                 bar3_out(ii,jj) /= sum3;
    }
  }


  
  for(auto s : meshes)
  {
    std::cout << "Mesh " << s << std::endl;
    cargHandler arg("-in " + s + "  -anamet 1 -verb 1");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);

    msh.cleanup();

    const double xtol = 1.0e-8;

    double maxErrPrecond = -1.0e30;
    double maxErr = -1.0e30;

    const double tol   = 1.0e-12;
    const int mdx = 256;
    const double dx0   = 1.0e-1;
    const double dx1   = 1.0e-8;
    const double qdx   = 10.0;
    const double minsl = 0.5;
    int ndx = 0;
    for(double dx = dx0; dx > dx1; dx /= qdx){
      ndx++;
    }
    if(ndx > mdx) METRIS_THROW_MSG(SMemExcept(),"Increase mdx")
    double errgdx[mdx], errhdx[mdx], logdx[mdx];



    int naliv = 0;
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      CT_FOR0_INC(2,3,idim){if(idim == msh.get_tdim()){
        constexpr int gdim = idim;
        constexpr int tdim = idim;
        constexpr int nhess = (gdim*(gdim + 1))/2;
        double coopr[idim];


        intAr2& ent2poi = msh.ent2poi(idim);
        int nentt = msh.nentt(idim);
        double coor0[idim];
        constexpr auto evalf = idim == 1 ? eval1<idim,ideg> : 
                               idim == 2 ? eval2<idim,ideg> : 
                                           eval3<idim,ideg>;


        dblAr2& bary = idim == 1 ? bar1 :
                       idim == 2 ? bar2 : bar3;
        dblAr2& bary_out = idim == 1 ? bar1_out :
                           idim == 2 ? bar2_out : bar3_out;

        printf("-- TEST 0: test derivatives of cost function\n");
        for(int ientt = 0; ientt < nentt; ientt++){
          if(isdeadent(ientt,ent2poi)) continue;

          bool iinva;
          double ccoef[getnnod3(idim*(ideg-1))];
          getsclccoef<idim,idim,ideg>(msh,ientt,NULL,ccoef,&iinva);
          BOOST_REQUIRE(!iinva);

          naliv++;
          double eps = getepsent<idim>(msh,idim,ientt);
          for(int isamp = 0; isamp < nsamp; isamp++){
            evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                  DifVar::None,bary[isamp],coor0,NULL,NULL);

            double evals[1+idim], grads[1+idim][idim];
            double grad_disc[idim], hess_disc[nhess], hess[nhess];
            evals[0] = invevalfun<gdim,ideg>(msh, ent2poi[ientt], 
                                             msh.coord, 
                                             coor0, bary[isamp], 
                                             1  , coopr, grads[0], hess);

            int idx = 0;
            for(double dx = dx0; dx > dx1; dx /= qdx){
              for(int ii = 0; ii < idim; ii++){
                double bar1[idim+1];
                for(int jj = 0; jj < idim + 1; jj++) bar1[jj] = bary[isamp][jj];
                  bar1[ii + 1] += dx;
                bar1[0] -= dx;
                evals[ii+1] = invevalfun<gdim,ideg>(msh, ent2poi[ientt],msh.coord,
                                                    coor0, bar1, 
                                                    0  , coopr, grads[ii+1], NULL);
              }
              for(int idif1 = 0; idif1 < idim; idif1++){
                grad_disc[idif1] = (evals[idif1+1] - evals[0])/dx;
                for(int idif2 = 0; idif2 < idim; idif2++){
                  hess_disc[sym2idx(idif1,idif2)] 
                     = (grads[1+idif1][idif2] - grads[0][idif2])/dx;
                }
              }
              double errg = sqrt(geterrl2<idim >(grad_disc,grads[0]));
              double errh = sqrt(geterrl2<nhess>(hess_disc,hess));
              errgdx[idx] = log(errg);
              errhdx[idx] = log(errh);
              logdx[idx] = log(dx);
              idx++;
            }// for dx
            double slopg = linearRegression(ndx,logdx,errgdx);
            double sloph = linearRegression(ndx,logdx,errhdx);
            BOOST_CHECK_MESSAGE(minsl < slopg || exp(errgdx[0]) < tol,
              " Grad slopg "<<slopg<<" under minimum "<<minsl<<"\n");
            BOOST_CHECK_MESSAGE(minsl < sloph || exp(errhdx[0]) < tol,
              " Grad sloph "<<sloph<<" under minimum "<<minsl<<"\n");
            //printf("Debug grad slope %f hess %f \n",slopg,sloph);
          }// for isamp
        }// for ientt


        int itest = 1;
        for(const dblAr2* bary_p : {&bary,&bary_out}){
          const dblAr2 &bary_ = *bary_p;
          printf(" - TEST %d: sample in and locate ",itest);
          if(itest == 1) printf("in      element\n");
          else           printf("outside element\n");
          itest++;

          int ifun = 0;
          using InvevalFunc = int (*)(MeshBase&, int, const double*, double*, double*, double);
          std::vector<InvevalFunc> funcs = {&inveval_badNewton<idim, ideg>,
                                            &inveval<idim, ideg>,
                                            &inveval_nloptD<idim, ideg>};
          for(auto inveval_fun : funcs){
            //printf("-- Testing function %d ",ifun);
            if(ifun == 0) printf(" - inveval_badNewton:");
            if(ifun == 1) printf(" - inveval          :");
            if(ifun == 2) printf(" - inveval_nloptD   :");

            double t0 = get_wall_time();
            int nerro = 0;
            int nsucc = 0;
            for(int ientt = 0; ientt < nentt; ientt++){
              if(isdeadent(ientt,ent2poi)) continue;

              bool iinva;
              double ccoef[getnnod3(idim*(ideg-1))];
              getsclccoef<idim,idim,ideg>(msh,ientt,NULL,ccoef,&iinva);
              BOOST_REQUIRE(!iinva);

              naliv++;
              double eps = getepsent<idim>(msh,idim,ientt);
              for(int isamp = 0; isamp < nsamp; isamp++){
                evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                      DifVar::None,bary_[isamp],coor0,NULL,NULL);

                double bar1[idim+1] = {1.0 / (idim + 1)};
                int ncall;
                int ierro = inveval_fun(msh,ientt,coor0,coopr,bar1,xtol);
                double err = sqrt(geterrl2<idim>(coor0,coopr));
                //#ifndef NDEBUG
                //BOOST_TEST(err / eps < xtol);
                //#endif
                if(err < xtol*eps){
                  nsucc++;
                }else{
                  nerro++;
                }
                maxErr = MAX(err/eps,maxErr);
                //#if 0
                if(err >= xtol*eps && ifun == 1 && msh.param->dbgfull){
                  msh.param->iverb = 10;
                  double meas = getmeasentP1<gdim>(ent2poi[ientt],msh.coord);
                  meas = sqrt(meas);
                  printf("Failed with err = %e xtol %e eps %e meas %e\n",err,xtol,eps,meas);
                  double eps2 = 0; 
                  for(int ii = 0; ii < tdim; ii++){
                    for(int jj = ii+1;jj < tdim; jj++){
                      eps2 += geterrl2<gdim>(msh.coord[ent2poi(ientt,ii)],msh.coord[ent2poi(ientt,jj)]);
                    }
                  }
                  eps2 = sqrt(eps2);

                  double dum[gdim], jmat[gdim][gdim];
                  evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::Bary,
                        DifVar::None,bary_[isamp],dum,jmat[0],NULL);
                  double jmat2[gdim][gdim];
                  matXtmat<gdim>(jmat[0],jmat[0],jmat2[0]);
                  double eps3 = 0;
                  for(int ii = 0; ii < gdim; ii++) eps3 += jmat2[ii][ii];
                  eps3 = sqrt(eps3);

                  // (Frob J^-{1})^{-1}
                  invmat<gdim>(jmat2[0]);
                  double eps4 = 0;
                  for(int ii = 0; ii < gdim; ii++) eps4 += jmat2[ii][ii];
                  eps4 = 1.0 / sqrt(eps4);

                  double eigval[gdim], eigvec[gdim*gdim];
                  double symm[(gdim*(gdim+1))/2];
                  for(int ii = 0; ii < gdim; ii++){
                    for(int jj = ii ; jj < gdim; jj++){
                      symm[sym2idx(ii,jj)] = jmat2[ii][jj];
                    }
                  }
                  geteigsym<gdim>(symm,eigval,eigvec);
                  printf("min eigval %e max %e one over %e %e \n",
                    eigval[0], eigval[gdim-1],
                    1.0/eigval[0], 1.0/eigval[gdim-1]);


                  printf(" eps2 = %e eps3 = %e eps4 %e crit2 %d crit3 %d\n",eps2,eps3,eps4,
                    err < xtol*eps2,err < xtol*eps3);
                  printf(" max eps %e\n",err/xtol);

                  printf("\n\nWith prints:\n\n\n");
                  ierro = inveval_fun(msh,ientt,coor0,coopr,bar1,xtol);
                  printf("## EXIT after %d success \n",nsucc);
                  printf("Returns bary ");
                  dblAr1(idim+1,bar1).print();
                  printf("Was seeking bary = "); 
                  dblAr1(idim+1,bary_[isamp]).print();
                  printf("Error = ");
                  for(int ii = 0; ii < idim + 1; ii++)printf(" %e ",abs(bary_[isamp][ii] - bar1[ii]) );
                  printf("\n");
                  err = sqrt(geterrl2<idim>(coor0,coopr))/eps;
                  printf("Got error = %15.7e\n",err);

                  //printf("## TRY WITH NLOPT\n");
                  //inveval_nloptD<gdim,ideg>(msh,ent2poi[ientt],coor0, coopr);
                  wait();

                  //printf("## TRY with adhoc newton\n");
                 // invevaltest<gdim,ideg>(msh,ent2poi[ientt],coor0,coopr,bar1,xtol);

                } // if error 
                //#endif

              } // for isamp 
            } // for ientt
            double t1 = get_wall_time();
            double pct_err = nerro / (double) (nsucc + nerro) * 100.0;
            printf(" nerro %7d nsucc %7d pct err %4.1f maxErr = %.2e total time %f = %dk op/s\n",
                   nerro,nsucc,pct_err,maxErr, t1-t0,(int)(nsamp*naliv/(1000*(t1-t0))));

            ifun++;
          } // for function 
        } // for bary
      }}CT_FOR1(idim);
    }}CT_FOR1(ideg);

    

  }
}

} //namespace Metris 
