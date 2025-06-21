//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include "common_setup.hxx"

#include "../src/quality/quafun_tradet.hxx"

namespace Metris{

typedef MetricFieldFE MFT;

template <class MFT, int gdim, int tdim, typename ftype>
void quafun_tradet_nodet(Mesh<MFT> &msh,AsDeg asdmsh, AsDeg asdmet,
                         const int*__restrict__ ent2pol,  
                         const double*__restrict__ bary,
                         const double*__restrict__ met_,
                         ftype*__restrict__ tra,
                         ftype*__restrict__ det);

BOOST_AUTO_TEST_CASE(test_eigen) 
{

  std::vector<std::string> meshes = 
  {METRIS_CASES_DIR "/2D/square.p1.10.meshb -sclmet 0.05 -adapt 20",
   METRIS_CASES_DIR "/1200_p1.meshb",
   METRIS_CASES_DIR "/2D/square.circmet.5k.curved.meshb  -sclmet 0.5 -adapt 20 -prefix tmp/ -out out",
  };


  for(auto testcase : meshes){
    std::string mesh_name = testcase;
    std::cout<<"-- Test case mesh = "<<mesh_name<<std::endl;
    cargHandler arg("-in "+mesh_name+" -verb 0 -t 2");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);

    run.adaptMesh();

    msh.cleanup();
    msh.met.setSpace(MetSpace::Log);

    for(int ideg = 0; ideg < 2; ideg++){
      if(ideg == 1) run.degElevate();

      // i = gdim, j = tdim
      double opps_full[4][4],opps_nodet[4][4];
      double dum = 0;
      CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
        CT_FOR0_INC(2,gdim,tdim){if(tdim <= msh.get_tdim()){
          double bary[tdim + 1];
          for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1/(tdim + 1.0);
          auto quafun = quafun_tradet<MFT,gdim,tdim,double>;
          auto quafun_nodet = quafun_tradet_nodet<MFT,gdim,tdim,double>;

          int nentt = msh.nentt(tdim);
          const intAr2& ent2poi = msh.ent2poi(tdim);

          double t0_full = get_wall_time();
          for(int ientt = 0; ientt < nentt; ientt++){
            double tra, det;
            quafun(msh, AsDeg::Pk, AsDeg::Pk, ent2poi[ientt], bary, NULL, &tra, &det);
            dum += tra*det;
          }// for ientt
          double t1_full = get_wall_time();
          opps_full[gdim][tdim] = nentt / (t1_full - t0_full);

          double t0_nodet = get_wall_time();
          for(int ientt = 0; ientt < nentt; ientt++){
            double tra, det;
            quafun_nodet(msh, AsDeg::Pk, AsDeg::Pk, ent2poi[ientt], bary, NULL, &tra, &det);
            dum += tra*det;
          }// for ientt
          double t1_nodet = get_wall_time();
          opps_nodet[gdim][tdim] = nentt / (t1_nodet - t0_nodet);

          printf("-- gdim %d tdim %d ideg %d full %dk/s nodet %dk/s\n",gdim,tdim,msh.curdeg,
            (int) (opps_full[gdim][tdim]/1000),(int)(opps_nodet[gdim][tdim]/1000));

        }}CT_FOR1(tdim);
      }}CT_FOR1(gdim);



    }

  }// for testcase


}//BOOST_AUTO_TEST_CASE




// For some special barys (nodes), met is already known -> pass it in
template <class MFT, int gdim, int tdim, typename ftype>
void quafun_tradet_nodet(Mesh<MFT> &msh,AsDeg asdmsh, AsDeg asdmet,
                         const int*__restrict__ ent2pol,  
                         const double*__restrict__ bary,
                         const double*__restrict__ met_,
                         ftype*__restrict__ tra,
                         ftype*__restrict__ det){

  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim); 

  METRIS_ASSERT(gdim == msh.idim);



  // Only if metric interpolation is needed
  if(msh.met.getSpace() != MetSpace::Log && met_ == NULL) METRIS_THROW_MSG(WArgExcept(),
      "## SET MESH METRIC TO LOG BEFORE CALLING metqua2_xi");

  constexpr int nnmet = (gdim*(gdim+1))/2;

  double jmat[tdim*gdim],coopr[gdim];
  double met[nnmet];


  // Get Jacobian matrix at xi
  if(asdmsh == AsDeg::P1 || msh.curdeg == 1){
    if constexpr(tdim == 2){  
      eval2<gdim,1>(msh.coord,ent2pol,msh.getBasis(),
                       DifVar::Bary,DifVar::None,
                       bary,coopr,jmat,NULL);
    }else{
      eval3<gdim,1>(msh.coord,ent2pol,msh.getBasis(),
                       DifVar::Bary,DifVar::None,
                       bary,coopr,jmat,NULL);
    }
  }else{
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){

      if constexpr(tdim == 2){  
        eval2<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),
                         DifVar::Bary,DifVar::None,
                         bary,coopr,jmat,NULL);
      }else{
        eval3<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),
                         DifVar::Bary,DifVar::None,
                         bary,coopr,jmat,NULL);
      }

    }}CT_FOR1(ideg);
  }

  if(met_ == NULL){
    // Get metric interpolated at xi (or eval if analytical)
    msh.met.getMetFullinfo(asdmet,DifVar::None,MetSpace::Exp,
                           ent2pol,tdim,bary,coopr,met,NULL);
  }else{
    for(int ii = 0; ii < nnmet; ii++) met[ii] = met_[ii];
  }


  // Compute J_0^{-T} J_K^T M J_K J_0^{-1}
  // Starting with J_K J_0^{-1}
  // Note that J_K is stored transposed w.r.t. above ! -> jmat[i,j] = d_i F_j 
  // whereas (J_K)_{ij} = d_j F_i .. 


  // Get J_0^{-T} J_K^T
  ftype invtJ0_tJK[tdim*gdim];

  matXmat<tdim,tdim,gdim>(Constants::invtJ_0[hana::type_c<ftype>][tdim],
                          jmat,invtJ0_tJK);


  ftype J0tJtMJJ0_diag[tdim];
  matXsymXtmat_diag<gdim, tdim, double, ftype, ftype>(met, invtJ0_tJK, J0tJtMJJ0_diag);
  *tra = J0tJtMJJ0_diag[0] + J0tJtMJJ0_diag[1];
  if constexpr (tdim == 3) *tra += J0tJtMJJ0_diag[2];
  
  // This is an actual exception that should never theoretically happen. 
  if(*tra < 1.0e-16) METRIS_THROW_MSG(GeomExcept(),"Zero trace of spd matrix? "<<*tra);


  if constexpr(tdim == gdim){
    *det = 0;
  }else{
    static_assert(tdim == 2 && gdim == 3);
    ftype tJ0_tJK_M_JK_J0[3];
    matXsymXtmat<2,3,double,ftype,ftype>(met,invtJ0_tJK,tJ0_tJK_M_JK_J0);
    *det = tJ0_tJK_M_JK_J0[0];
  }

   return;
}



}//namespace