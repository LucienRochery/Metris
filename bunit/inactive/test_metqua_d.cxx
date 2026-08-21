//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include <bunit/common_setup.hxx>

#include <boost/timer/progress_display.hpp>

//#include "utils/aux_misc.hxx"
#include "utils/CT_loop.hxx"
#include "quality/low_metqua.hxx"
#include "quality/low_metqua_d.hxx"

#include "SANS/Surreal/SurrealS.h"
#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
//#include <src/msh_metric.hxx>

namespace utf = boost::unit_test;


namespace Metris{


typedef MetricFieldAnalytical MFT;


BOOST_AUTO_TEST_CASE(test_metqua_d) 
{//METRIS_MAX_DEG


  using ftyp1 = float8;
  //using ftyp2 = double;
  using ftyp2 = float8;

  const int nstrmsh = 3;
  std::vector<std::string> meshes = 
  {"../cases/debug/4tri.mesh",
   "../cases/2D/square.p1.10.meshb",
   //"../cases/2D/square.circmet.50.curved.meshb",
   "../cases/2D/square.circmet.5k.curved.meshb",
   //"../cases/1200_p1.meshb", 
   //"../cases/1200_p2.meshb",
   //"../cases/curved_p2.bez.meshb",
   //"../cases/curved_p2.meshb",
   //#if METRIS_MAX_DEG >= 3
   //  "../cases/1200_p3.meshb",
   //  "../cases/curved_p3.bez.meshb",
   //  "../cases/curved_p3.meshb"
   //#endif
  };
  const ftyp1 dx0   = 1.0e-2;
  const ftyp1 dx1   = 1.0e-8;
  const ftyp1 qdx   = 5.0;
  const double minsl = 0.1;
  const double maxsl = 2.5;
  const double tol   = 1.0e-5;
  const double tol0  = 1.0e-5; // close to 0


  double maxerr_eval = -1.0;
  double avgerr_eval = 0;
  unsigned long long int nerr_eval = 0; 

  double maxerr_diff_rel = -1.0;
  double maxerr_diff_abs = -1.0;
  double maxerr_diff_abs_hess = -1.0;
  double maxerr_diff_abs_Hess = -1.0;
  double avgerr_diff_rel = 0;
  double avgerr_diff_abs = 0;
  double avgerr_diff_abs_hess = 0;
  double avgerr_diff_abs_Hess = 0;
  unsigned long long int nerr_diff = 0; 
  unsigned long long int nerr_diff2 = 0; 

  int ndx = 0;
  for(ftyp1 dx = dx0; dx > dx1; dx /= qdx){
    ndx++;
  }
  double errdx_rel[512],errdx_abs[512],logdx[512];
  double errdx_abs_hess[512];

  double t0,t1,t2;
  ftyp1 dqua_disc[3], hqua_disc[6];
  ftyp1 dqua_surr[3], hqua_surr[6];
  ftyp1 Hqua_disc[6*(6+1)/2]; // size for 2*gdim (max 3) Hess
  ftyp1 Hqua_surr[6*(6+1)/2]; // size for 2*gdim (max 3) Hess

  const int power = 1;

  char **argv = (char**) malloc(256*sizeof(char*));
  int argc;
  
  std::vector<std::string> isomsh = 
  {"cases/isotri.meshb",
   "cases/isotet.meshb",
   "cases/isotet_p2.bez.meshb"};

  int imsh = 0;
  for(auto s : meshes)
  {

    imsh++;
    std::cout<<"Mesh "<<s<<std::endl;
    cargHandler arg("-i "+s+" -verb 0 -anamet 1 -sclmet 1");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);

    msh.cleanup();

    //msh.met.setSpace(MetSpace::Log);

    constexpr AsDeg asdmet = AsDeg::Pk;
    constexpr FEBasis dofbas = FEBasis::Lagrange; 
    constexpr DifVar idifmet = DifVar::None;

    msh.setBasis(dofbas);


    CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){


      ftyp1 qutets[4][2];
      ftyp1 dqutets[4][2][gdim];
      ftyp1 dqutets11[4][2][gdim],
            dqutets12[4][2][gdim],
            dqutets21[4][2][gdim],
            dqutets22[4][2][gdim];


      constexpr int tdim = gdim;
      double bary[tdim+1];
      constexpr int nhess = (gdim*(gdim+1))/2;
      constexpr int nHess2 = (gdim*2 * (gdim*2 + 1))/2; // Hess size for 2 pt
      ftyp2 qutet,dqutet[gdim],hqutet[nhess];
      ftyp1 hdum[nhess];
      
      intAr2 &ent2poi = msh.ent2poi(tdim);
      int nentt = msh.nentt(tdim);

      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
        constexpr int nrfld = tdim == 2 ? facnpps[ideg] : tetnpps[ideg]; 
        constexpr int nHess  = (gdim*nrfld*(gdim*nrfld+1))/2;
        #if 0
        for(int ientt = 0; ientt < nentt; ientt++){
          if(isdeadent(ientt,ent2poi)) continue;

          for(int ibary = 0; ibary < gdim+1; ibary++){
            for(int ii = 0; ii < gdim + 1; ii++){
              bary[ii] = 0;
            }
            bary[ibary] = 1;
            qutet = d_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp2>(msh,ent2poi[ientt],bary,power,-1,dofbas,idifmet,NULL);

            ftyp1 qutet0;
            quafun_distortion<MFT,gdim,gdim,ideg,asdmet,ftyp1>(msh,ent2poi[ientt],power,bary,&qutet0);
            
            double err = (double)abs((ftyp1)qutet - qutet0);
            double nrm = (double)abs(qutet0);
            maxerr_eval  = maxerr_eval > err/nrm ? maxerr_eval : err/nrm;
            avgerr_eval += err/nrm;
            nerr_eval++;


            BOOST_CHECK_MESSAGE(err < (double)(tol*(double)qutet0),
              "Large quality value error at xi:"<<abs(qutet - qutet0)/(double)qutet0<<"\n");

            if(err >= (double)(tol*(double)qutet0)){
              std::cout<<"Large quality value error at xi:"<<abs(qutet - qutet0)/(double)qutet0<<"\n";
              std::cout<<"ielem = "<<ientt<<"\n";
              std::cout<<"metqua   gives = "<<qutet0<<"\n";
              std::cout<<"metqua_d gives = "<<qutet<<"\n";
              std::cout<<"bary was ";
              dblAr1(gdim+1,bary).print();
              for(int ii = 0; ii < nrfld; ii++){
                int ipoin = ent2poi(ientt,ii);
                printf(" node %d = %d coords = ",ii,ipoin);
                dblAr1(gdim,msh.coord[ipoin]).print();
              }
            }

            if(err > 1.0*qutet0){
              printf("Very large error %15.7e qutet0 = %15.7e \n",(double)err,(double)qutet0);
              wait();
            }
          }   
        }

        printf("\n\n Tests passed for _xi \n\n");


        for(int ientt = 0; ientt < nentt; ientt++){
          if(isdeadent(ientt,ent2poi)) continue;

          qutet = d_metqua<MFT,gdim,ideg,asdmet,ftyp2>(msh,ientt,power,-1,dofbas,idifmet,NULL);

          ftyp1 qutet0;
          metqua<MFT,gdim,gdim,ideg,asdmet,ftyp1>(msh,ientt,power,&qutet0);
          
          double err = (double)abs((ftyp1)qutet - qutet0);
          double nrm = (double)abs(qutet0);
          maxerr_eval  = maxerr_eval > err/nrm ? maxerr_eval : err/nrm;
          avgerr_eval += err/nrm;
          nerr_eval++;


          BOOST_CHECK_MESSAGE(err < (double)(tol*(double)qutet0),
            "Large quality value error:"<<abs(qutet - qutet0)/(double)qutet0<<"\n");

          if(err >= (double)(tol*(double)qutet0)){
            std::cout<<"Large quality value error:"<<abs(qutet - qutet0)/(double)qutet0<<"\n";
            std::cout<<"ielem = "<<ientt<<"\n";
            std::cout<<"metqua   gives = "<<qutet0<<"\n";
            std::cout<<"metqua_d gives = "<<qutet<<"\n";

          }

          if(err > 1.0*qutet0){
            printf("Very large error %15.7e qutet0 = %15.7e \n",(double)err,(double)qutet0);
            wait();
          }
        }   
        #endif
        
        printf("-- Value tests passed.\n");
        printf("-- Start diff tests.\n");
        // Test derivatives
        // Linear function: exactly given by diffs
        boost::timer::progress_display progress(nentt);
        for(int ii = 0; ii < tdim+1; ii++){
          bary[ii] = 1.0 / (tdim + 1);
        }
        for(int ientt = 0; ientt < nentt; ientt++){
          if(isdeadent(ientt,ent2poi)) continue;

          for(int ivar = 0; ivar < nrfld; ivar++){
            int ipoin = ent2poi(ientt,ivar);

            int idx = 0;
            //qutet = d_metqua<MFT,gdim,ideg,asdmet,ftyp2>(msh,ientt,power,ivar,dofbas,idifmet,dqutet);
            qutet = d_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp2>(msh,ent2poi[ientt],bary,
                                              power,ivar,dofbas,idifmet,dqutet,hqutet);
            //printf("Debug ivar = %d hess = \n",ivar);
            //MeshArray1D<ftyp2>(3,hqutet).print();
            //printf(" grad = \n");
            //MeshArray1D<ftyp2>(gdim,dqutet).print();

            qutets[0][0] = qutet;

            double minerr_abs = 1.0e30;
            double minerr_abs_hess = 1.0e30;
            double minerr_rel = 1.0e30;
            bool ismall_grad = false;
            bool ismall_hess = false;
            ftyp1 cc = (ftyp1) getepsent<gdim>(msh,gdim,ientt);
            for(ftyp1 dx = dx0*cc; dx > dx1*cc; dx /= qdx){
              double coor0[gdim];
              for(int i = 0; i < gdim; i++) coor0[i] = msh.coord(ipoin,i);
              
  
              for(int idim = 0 ; idim < gdim; idim++){
                int jsig = 0;
                for(int isig = -1; isig <= 1; isig += 2){
                  for(int ii = 0; ii < gdim; ii++) msh.coord(ipoin,ii) = coor0[ii];
                  msh.coord(ipoin,idim) = coor0[idim] + (double)(isig*dx);
                  //metqua<MFT,gdim,gdim,ideg,asdmet,ftyp1>(msh,ientt,power,&qutets[idim+1][jsig]);
                  //quafun_distortion<MFT,gdim,gdim,ideg,asdmet,ftyp1>(msh,ent2poi[ientt],power,
                  //                                  bary,&qutets[idim+1][jsig]);
                  qutets[idim+1][jsig]= d_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp1>(
                                                  msh,ent2poi[ientt],bary,power,
                                                  ivar,dofbas,idifmet,dqutets[idim+1][jsig],hdum);
                  jsig++;
                }
              }

              for(int ii = 0; ii < gdim; ii++) msh.coord(ipoin,ii) = coor0[ii];
              for(int ii = 0; ii < gdim; ii++){
                dqua_disc[ii] = (qutets[ii+1][1] - qutets[ii+1][0])/(2*dx);
                //dqua_disc[i] = (qutets[i+1][1] - qutets[0][0])/(dx);
                dqua_surr[ii] = (ftyp1)dqutet[ii];//(ftyp1)(qutet.deriv(i));
                for(int jj = 0; jj < gdim; jj++){
                  hqua_disc[sym3idx(ii,jj)] = (dqutets[ii+1][1][jj] - dqutets[ii+1][0][jj])/(2*dx);
                  hqua_surr[sym3idx(ii,jj)] = hqutet[sym3idx(ii,jj)];
                }
              }

              double nrm = (double)getnrml2<gdim,ftyp1>(dqua_surr);
              if(nrm < tol0) nrm = tol0;
              double err_rel = (double)sqrt(geterrl2<gdim,ftyp1>(dqua_disc,dqua_surr)/nrm);
              double err_abs = (double)sqrt(geterrl2<gdim,ftyp1>(dqua_disc,dqua_surr));

              double err_abs_hess = (double)sqrt(geterrl2<nhess,ftyp1>(hqua_disc,hqua_surr));

              logdx[idx] = (double)(dx);
              errdx_rel[idx] = (double)(err_rel);
              errdx_abs[idx] = (double)(err_abs);
              errdx_abs_hess[idx] = (double)(err_abs_hess);
              if(err_abs < tol0) ismall_grad = true;
              if(err_abs_hess < tol0) ismall_hess = true;
              minerr_rel = METRIS_MIN(minerr_rel,errdx_rel[idx]);
              minerr_abs = METRIS_MIN(minerr_abs,errdx_abs[idx]);
              minerr_abs_hess = METRIS_MIN(minerr_abs_hess,errdx_abs_hess[idx]);
              idx++;
              //printf("Dbg dx = %15.7e error = %15.7e\n",dx,errdx[idx]);
              //printf("Debug dqua_disc = %f %f %f \n",dqua_disc[0],dqua_disc[1],dqua_disc[2]);
              //printf("Debug dqua_surr = %f %f %f \n",dqua_surr[0],dqua_surr[1],dqua_surr[2]);
              //printf("Debug diff      = %f %f %f \n",dqua_surr[0]-dqua_disc[0],dqua_surr[1]-dqua_disc[1],dqua_surr[2]-dqua_disc[2]);
            }

            maxerr_diff_rel = METRIS_MAX(maxerr_diff_rel,minerr_rel);
            maxerr_diff_abs = METRIS_MAX(maxerr_diff_abs,minerr_abs);
            maxerr_diff_abs_hess = METRIS_MAX(maxerr_diff_abs_hess,minerr_abs_hess);
            avgerr_diff_rel += minerr_rel;
            avgerr_diff_abs += minerr_abs;
            avgerr_diff_abs_hess += minerr_abs_hess;
            nerr_diff++;

            bool iok = true; 
            double slope_grad, slope_hess;
            if(!ismall_grad){
              //bool ismall = false;
              //for(int ii = 0; ii < ndx; ii++){
              //  logdx[ii] = log(logdx[ii]);
              //  if(errdx_rel[ii] < tol0){
              //    // Basically 0
              //    ismall = true;
              //    break;
              //  }
              //  errdx_rel[ii] = log(errdx_rel[ii]);
              //}
              //double slope;
              //if(!ismall) slope = linearRegression(ndx,logdx,errdx_abs);
              //iok = iok && (minsl < slope || ismall); 
              slope_grad = linearRegression(ndx,logdx,errdx_abs);
              iok = iok && minsl < slope_grad;
            }

            if(!ismall_hess){
              slope_hess = linearRegression(ndx,logdx,errdx_abs_hess);
              iok = iok && minsl < slope_hess;
            }
            //#ifndef NDEBUG
            if(!iok){
              printf("First test failure \n");
              printf("Failed test %d/%d.\n",ientt*nrfld+ivar,nentt*nrfld);

              printf("-- Gradient errors:\n");
              printf("min err log = %15.7e exp = %15.7e\n",log(minerr_rel),minerr_rel);
              MeshArray1D<double>(ndx,errdx_rel).print();
              printf("absolute errors: ");
              MeshArray1D<double>(ndx,errdx_abs).print();

              printf("-- Hessian errors:\n");
              printf("absolute errors: ");
              MeshArray1D<double>(ndx,errdx_abs_hess).print();

              printf("Slope grad = %f hess = %f \n",slope_grad,slope_hess);
              printf("logdx = \n");
              dblAr1(ndx,logdx).print();
              
              printf("Analytical gradient:\n");
              MeshArray1D<ftyp1>(gdim,dqua_surr).print();
              printf("Diff gradient:\n");
              MeshArray1D<ftyp1>(gdim,dqua_disc).print();

              printf("Analytical Hessian:\n");
              MeshArray1D<ftyp1>(nhess,hqua_surr).print();
              printf("Diff Hessian:\n");
              MeshArray1D<ftyp1>(nhess,hqua_disc).print();


              std::cout<<"ielem = "<<ientt<<"\n";
              std::cout<<"bary was ";
              dblAr1(gdim+1,bary).print();
              for(int ii = 0; ii < nrfld; ii++){
                int ipoin = ent2poi(ientt,ii);
                printf(" node %d = %d coords = ",ii,ipoin);
                dblAr1(gdim,msh.coord[ipoin]).print();
              }

              qutet = d_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp2>(msh,ent2poi[ientt],bary,
                                                power,ivar,dofbas,idifmet,dqutet,hqutet);

              wait();
            }
            //#endif
            BOOST_CHECK_MESSAGE(iok,
              " Slope grad "<<slope_grad<<" >? "<<minsl<<"\n"<<
              " Slope Hess "<<slope_hess<<" >? "<<minsl<<"\n"<<             
              " ivar = "<<ivar<<" last error "<<
              errdx_abs[ndx-1]<<" firsst "<<errdx_abs[0]<<
              " minimum "<<minerr_abs<<"\n");
          } // for int ivar 

          //printf("Passed grad and single point hessian tests (1elt)\n");

          // 2 nodes 
          ftyp2 Dqutet[nrfld][gdim];
          ftyp2 Hqutet[nHess];

          qutet = D_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp2>(msh,ent2poi[ientt],bary,
                                            power,dofbas,idifmet,Dqutet[0],Hqutet);
          for(int inod1 = 0; inod1 < nrfld; inod1++){
            int ipoi1 = ent2poi(ientt,inod1);

            for(int inod2 = inod1; inod2 < nrfld; inod2++){
              int ipoi2 = ent2poi(ientt,inod2);
              
              double coor01[gdim], coor02[gdim];
              for(int i = 0; i < gdim; i++) coor01[i] = msh.coord(ipoi1,i);
              for(int i = 0; i < gdim; i++) coor02[i] = msh.coord(ipoi2,i);

              // Move ipoi2 and use gradient with respect to ipoi1 
              ftyp1 cc = (ftyp1) getepsent<gdim>(msh,gdim,ientt);
              int idx = 0;
              bool ismall = false;
              double minerr_abs = 1.0e30;
              for(ftyp1 dx = dx0*cc; dx > dx1*cc; dx /= qdx){

                for(int idim = 0 ; idim < gdim; idim++){
                  int jsig = 0;
                  for(int isig = -1; isig <= 1; isig += 2){
                    for(int ii = 0; ii < gdim; ii++) msh.coord(ipoi2,ii) = coor02[ii];
                    msh.coord(ipoi2,idim) += (double)(isig*dx);
                    d_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp1>(
                               msh,ent2poi[ientt],bary,power,
                               inod1,dofbas,idifmet,dqutets12[idim+1][jsig],hdum);
                    d_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp1>(
                               msh,ent2poi[ientt],bary,power,
                               inod2,dofbas,idifmet,dqutets22[idim+1][jsig],hdum);
                    jsig++;
                  } // for isig 
                } // for idim 

                for(int idim = 0 ; idim < gdim; idim++){
                  int jsig = 0;
                  for(int isig = -1; isig <= 1; isig += 2){
                    for(int ii = 0; ii < gdim; ii++) msh.coord(ipoi1,ii) = coor01[ii];
                    msh.coord(ipoi1,idim) += (double)(isig*dx);
                    d_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp1>(
                               msh,ent2poi[ientt],bary,power,
                               inod1,dofbas,idifmet,dqutets11[idim+1][jsig],hdum);
                    d_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp1>(
                               msh,ent2poi[ientt],bary,power,
                               inod2,dofbas,idifmet,dqutets21[idim+1][jsig],hdum);
                    jsig++;
                  } // for isig 
                } // for idim 

                for(int ii = 0; ii < gdim; ii++) msh.coord(ipoi2,ii) = coor01[ii];
                for(int ii = 0; ii < gdim; ii++) msh.coord(ipoi2,ii) = coor02[ii];

                for(int ii = 0; ii < gdim; ii++){
                  for(int jj = 0; jj < gdim; jj++){
                    // dxx
                    //std::cout<<" using values "<<dqutets11[ii+1][1][jj]<<" and "
                    //<<dqutets11[ii+1][0][jj]<<" \n";
                    Hqua_disc[sym3idx(gdim*0+ii,gdim*0+jj)] 
                         = (dqutets11[ii+1][1][jj] - dqutets11[ii+1][0][jj])/(2*dx);
                    // dxy 
                    Hqua_disc[sym3idx(gdim*0+jj,gdim*1+ii)] 
                         = (dqutets12[ii+1][1][jj] - dqutets12[ii+1][0][jj])/(2*dx);
                         //- (dqutets21[ii+1][1][jj] - dqutets21[ii+1][0][jj])/(2*dx); // confirm finidiff ok
                    // dyy 
                    Hqua_disc[sym3idx(gdim*1+ii,gdim*1+jj)] 
                         = (dqutets22[ii+1][1][jj] - dqutets22[ii+1][0][jj])/(2*dx);

                    Hqua_surr[sym3idx(gdim*0+ii,gdim*0+jj)] = Hqutet[sym3idx(gdim*inod1+ii,gdim*inod1+jj)];
                    Hqua_surr[sym3idx(gdim*1+ii,gdim*0+jj)] = Hqutet[sym3idx(gdim*inod2+ii,gdim*inod1+jj)];
                    Hqua_surr[sym3idx(gdim*1+ii,gdim*1+jj)] = Hqutet[sym3idx(gdim*inod2+ii,gdim*inod2+jj)];
                  }
                }

                double err_abs = (double)sqrt(geterrl2<nHess2,ftyp1>(Hqua_disc,Hqua_surr));

                logdx[idx]     = (double)(dx);
                errdx_abs[idx] = (double)(err_abs);
                if(err_abs < tol0) ismall = true;
                minerr_abs = METRIS_MIN(minerr_abs,errdx_abs[idx]);
                idx++;
              } // for dx 

              maxerr_diff_abs_Hess = METRIS_MAX(maxerr_diff_abs_Hess,minerr_abs);
              avgerr_diff_abs_Hess += minerr_abs;
              nerr_diff2++;

              bool iok = true; 
              double slope;
              if(!ismall){
                slope = linearRegression(ndx,logdx,errdx_abs);
                iok = iok && minsl < slope;
              }
              if(!iok){
                printf("Failed test 2 ientt = %d inod1 = %d inod2 = %d.\n",ientt,inod1,inod2);
                printf("Hessian ordering (points X, Y):\n");
                for(int ii = 0; ii < gdim; ii++){
                  for(int jj = 0; jj < gdim; jj++){
                    printf("x%d x%d = %d\n",ii,jj,sym3idx(gdim*0+ii,gdim*0+jj));
                    printf("x%d y%d = %d\n",ii,jj,sym3idx(gdim*0+ii,gdim*1+jj));
                    printf("y%d x%d = %d\n",jj,ii,sym3idx(gdim*1+jj,gdim*0+ii));
                    printf("y%d y%d = %d\n",ii,jj,sym3idx(gdim*1+ii,gdim*1+jj));
                  }
                }

                printf("-- Hessian errors:\n");
                printf("absolute errors: ");
                MeshArray1D<double>(ndx,errdx_abs).print();

                printf("Slope hess = %f ismall = %d\n",slope,ismall);
                printf("logdx = \n");
                dblAr1(ndx,logdx).print();
                
                printf("Analytical Hessian:\n");
                MeshArray1D<ftyp1>(nHess2,Hqua_surr).print();
                printf("Diff Hessian:\n");
                MeshArray1D<ftyp1>(nHess2,Hqua_disc).print();

                printf("Analytical gradients: \n");
                MeshArray1D<ftyp2>(gdim,Dqutet[inod1]).print();
                MeshArray1D<ftyp2>(gdim,Dqutet[inod2]).print();

                std::cout<<"ielem = "<<ientt<<"\n";
                std::cout<<"bary was ";
                dblAr1(gdim+1,bary).print();
                for(int ii = 0; ii < nrfld; ii++){
                  int ipoin = ent2poi(ientt,ii);
                  printf(" node %d = %d coords = ",ii,ipoin);
                  dblAr1(gdim,msh.coord[ipoin]).print();
                }

                qutet = d_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp2>(msh,ent2poi[ientt],bary,
                                                  power,inod1,dofbas,idifmet,dqutet,hqutet);

                wait();
              }
              BOOST_CHECK_MESSAGE(iok,
                " Slope Hess "<<slope<<" >? "<<minsl<<"\n"<<             
                " inod1 = "<<inod1<<" last error "<<
                errdx_abs[ndx-1]<<" firsst "<<errdx_abs[0]<<
                " minimum "<<minerr_abs<<"\n");
            } // for inod2 


            //#ifndef NDEBUG
            //if(!iok){
            //  printf("Failed test %d/%d.\n",ientt*nrfld+inod1,nentt*nrfld);

            //  printf("-- Gradient errors:\n");
            //  printf("min err log = %15.7e exp = %15.7e\n",log(minerr_rel),minerr_rel);
            //  MeshArray1D<double>(ndx,errdx_rel).print();
            //  printf("absolute errors: ");
            //  MeshArray1D<double>(ndx,errdx_abs).print();

            //  printf("-- Hessian errors:\n");
            //  printf("absolute errors: ");
            //  MeshArray1D<double>(ndx,errdx_abs_hess).print();

            //  printf("Slope grad = %f hess = %f \n",slope_grad,slope_hess);
            //  printf("logdx = \n");
            //  dblAr1(ndx,logdx).print();
            //  
            //  printf("Analytical gradient:\n");
            //  MeshArray1D<ftyp1>(gdim,dqua_surr).print();
            //  printf("Diff gradient:\n");
            //  MeshArray1D<ftyp1>(gdim,dqua_disc).print();

            //  printf("Analytical Hessian:\n");
            //  MeshArray1D<ftyp1>(nhess,hqua_surr).print();
            //  printf("Diff Hessian:\n");
            //  MeshArray1D<ftyp1>(nhess,hqua_disc).print();


            //  std::cout<<"ielem = "<<ientt<<"\n";
            //  std::cout<<"bary was ";
            //  dblAr1(gdim+1,bary).print();
            //  for(int ii = 0; ii < nrfld; ii++){
            //    int ipoin = ent2poi(ientt,ii);
            //    printf(" node %d = %d coords = ",ii,ipoin);
            //    dblAr1(gdim,msh.coord[ipoin]).print();
            //  }

            //  qutet = d_quafun_distortion<MFT,gdim,ideg,asdmet,ftyp2>(msh,ent2poi[ientt],bary,
            //                                    power,inod1,dofbas,idifmet,dqutet,hqutet);

            //  wait();
            //}
            //#endif
          }// end inod1
          ++progress;
        } // end for ielem

        #ifdef NDEBUG
        double Dqutet[nrfld][gdim];
        double Hqutet[nHess];
        double dqutetD[gdim],hqutetD[nhess];
        ftyp1 hdum[nhess];
        double t0 = get_cpu_time();
        for(int ii = 0; ii < tdim+1; ii++){
          bary[ii] = 1.0 / (tdim + 1);
        }
        double dum = 0;
        for(int ientt = 0; ientt < nentt; ientt++){
          if(isdeadent(ientt,ent2poi)) continue;
          for(int inode = 0; inode < nrfld; inode++){
            dum += d_quafun_distortion<MFT,gdim,ideg,asdmet,double>(msh,ent2poi[ientt],bary,
                                              power,inode,dofbas,idifmet,dqutetD,hqutetD);
          }
        }
        double t1 = get_cpu_time();
        int ps = (int) ((double)(nentt*nrfld) / (t1 - t0) / 1000);
        printf("d_quafun_distortion benchmark t = %f op/s = %dk\n",t1-t0,ps);
        for(int ientt = 0; ientt < nentt; ientt++){
          if(isdeadent(ientt,ent2poi)) continue;
          dum += D_quafun_distortion<MFT,gdim,ideg,asdmet,double>(msh,ent2poi[ientt],bary,
                          power,dofbas,idifmet,Dqutet[0],Hqutet);
        }
        t0 = get_cpu_time();
        ps = (int) ((double)(nentt) / (t0 - t1) / 1000);
        printf("D_quafun_distortion benchmark t = %f op/s = %dk\n",t1-t0,ps);
        printf("dum = %f\n",dum);
        printf("\n\n");

        #endif


        // END body
      }}CT_FOR1(ideg);
      std::cout<<"-- All tests passed for mesh "<<s<<std::endl;
      //wait();
    }}CT_FOR1(gdim);
  } // end for auto meshes 

  avgerr_eval /= nerr_eval;
  avgerr_diff_abs /= nerr_diff;
  avgerr_diff_rel /= nerr_diff;
  avgerr_diff_abs_hess /= nerr_diff;
  printf("Eval err max = %15.7e avg = %15.7e\n",maxerr_eval,avgerr_eval);
  printf("Grad err max = %15.7e avg = %15.7e (rel)\n",maxerr_diff_rel,avgerr_diff_rel);
  printf("Grad err max = %15.7e avg = %15.7e (abs)\n",maxerr_diff_abs,avgerr_diff_abs);
  printf("Hess err max = %15.7e avg = %15.7e (abs)\n",maxerr_diff_abs_hess,avgerr_diff_abs_hess);

}


} // end namesapce
