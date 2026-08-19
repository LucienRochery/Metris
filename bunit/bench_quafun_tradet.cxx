//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 
#include <boost/test/included/unit_test.hpp> 

#ifndef USE_MULTIPRECISION

#warning "Unit test bench_quafun_tradet not compiled if USE_MULTIPRECISION=OFF"

// main
namespace Metris{
BOOST_AUTO_TEST_CASE(bench_quafun_tradet) 
{
}
}
#else



#include "common_setup.hxx"
#include "gen_bary.hxx"

#include "quality/quafun_tradet.hxx"
#include "low_eval.hxx"
#include "linalg/det.hxx"


namespace Metris{

typedef MetricFieldFE MFT;


template<int gdim>
void generate_frame(float8* eigvec,
                    std::uniform_real_distribution<double>& unif, 
                    std::default_random_engine& rng);


template <class MFT, int gdim, int tdim, typename ftype>
void quafun_tradet_nodet(Mesh<MFT> &msh,AsDeg asdmsh, AsDeg asdmet,
                         const int*__restrict__ ent2pol,  
                         const double*__restrict__ bary,
                         const double*__restrict__ met_,
                         ftype*__restrict__ tra,
                         ftype*__restrict__ det);

BOOST_AUTO_TEST_CASE(bench_quafun_tradet) 
{

  std::vector<std::string> meshes = 
  {METRIS_CASES_DIR "/unit/2D/square/iso.p1.100k",
   METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k",
   METRIS_CASES_DIR "/unit/2D/square/circmet.p2.500",
  };

  const double aniso_max = 2e7;
  const double aniso_mul = 10;
  const int ntrans = 2; // anisotropic linear transformations
  const int nsamp = 10; // bary samples


  dblAr2 barys[4];
  genBary(nsamp, 2, barys[2]);
  genBary(nsamp, 3, barys[3]);

  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(0);

  MeshArray2D<float8> transfo[4];
  for(int gdim = 2; gdim <= 3; gdim++){
    transfo[gdim].allocate(ntrans,gdim*gdim);
    transfo[gdim].set_n(ntrans);
  }
  for(int itrans = 0; itrans < ntrans; itrans++){
    generate_frame<2>(transfo[3][itrans], unif, rng);
    generate_frame<3>(transfo[3][itrans], unif, rng);
  }
  for(auto testcase : meshes){
    std::string mesh_name = testcase;
    std::cout<<"-- Test case mesh = "<<mesh_name<<std::endl;
    cargHandler arg("-in "+mesh_name+" -verb 0 -t 2");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);

    run.adaptMesh();

    msh.cleanup();
    msh.met.setSpace(MetSpace::Log);


    CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
      CT_FOR0_INC(2,gdim,tdim){if(tdim <= msh.get_tdim()){
        constexpr int nnmet = (gdim*(gdim+1))/2;
        auto quafun = quafun_tradet<MFT,gdim,tdim,double>;
        auto quafun_nodet = quafun_tradet_nodet<MFT,gdim,tdim,double>;

        int nentt = msh.nentt(tdim);
        const intAr2& ent2poi = msh.ent2poi(tdim);

        // ---------------------------------------------------------------

        for(int idegelev = 0; idegelev <= 1; idegelev++){
          if(idegelev == 1) run.degElevate();

          CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){

            const int nnode = getnnode(tdim,ideg);

            // --- Get function error, only the determinant term
            MeshArray2D<float8> coordf8(nnode,gdim);
            MeshArray2D<double> eltcrd(nnode,gdim);

            for(double aniso = 2; aniso <= aniso_max + 1; aniso *= aniso_mul){

              MinMaxAvg errsmat;

              for(int ientt = 0; ientt < nentt; ientt++){

                for(int inode = 0; inode < nnode; inode++){
                  int ipoin = ent2poi(ientt, inode);
                  for(int ii = 0; ii < gdim; ii++){
                    coordf8(inode, ii) = (float8) msh.coord(ipoin,ii);
                    eltcrd(inode, ii) = msh.coord(ipoin,ii);
                  }
                }
                for(int itrans = 0; itrans < ntrans; itrans++){

                  float8 f8_trans0[gdim*gdim], f8_trans1[gdim*gdim];

                  // Transfo is an orthogonal matrix R
                  // Make it anisotropic by multiplying by D = diag(1/aniso, 1, ...)
                  // Map coord by L = RD, metric by L^-T M L-1 = R D-1MD-1 R^T
                  for(int ii = 0; ii < gdim; ii++){
                    for(int jj = 0; jj < gdim; jj++){
                      f8_trans0[gdim*ii+jj] = (float8) transfo[gdim](itrans,gdim*ii+jj);
                      f8_trans1[gdim*ii+jj] = (float8) transfo[gdim](itrans,gdim*ii+jj);
                    }
                  }

                  // Multiply by diag(1/aniso, 1, ...)
                  for(int ii = 0; ii < gdim; ii++){
                    f8_trans1[ii] /= aniso;
                  }

                  for(int inode = 0; inode < nnode; inode++){
                    int ipoin = ent2poi(ientt, inode);
                    float8 coop[gdim];
                    // Transform the coordinates by transfo[itrans]
                    matXvec<gdim>(f8_trans1, coordf8[inode], coop);
                    for(int ii = 0; ii < gdim; ii++) msh.coord(ipoin,ii) = (double) coop[ii];
                  }

                  for(int isamp = 0; isamp < nsamp; isamp++){

                    double *bary = barys[tdim][isamp];


                    // 
                    double f2_jmat[tdim*gdim], f2_coopr[gdim];
                    float8 f8_jmat[tdim*gdim], f8_coopr[gdim];
                    if constexpr (tdim == 2){
                      eval2<gdim,ideg>(msh.coord,ent2poi[ientt],msh.getBasis(),
                          DifVar::Bary,DifVar::None,
                          bary,f2_coopr,f2_jmat,NULL);
                    }else{
                      eval3<gdim,ideg>(msh.coord,ent2poi[ientt],msh.getBasis(),
                          DifVar::Bary,DifVar::None,
                          bary,f2_coopr,f2_jmat,NULL);
                    }
                    for(int ii = 0; ii < gdim; ii++) f8_coopr[ii] = (float8) f2_coopr[ii];
                    for(int ii = 0; ii < tdim*gdim; ii++) f8_jmat[ii] = (float8) f2_jmat[ii];


                    double f2_invtJ0_tJK[tdim*gdim];
                    float8 f8_invtJ0_tJK[tdim*gdim];
                    matXmat<tdim,tdim,gdim>(Constants::invtJ_0[hana::type_c<double>][tdim],
                                            f2_jmat,f2_invtJ0_tJK);
                    matXmat<tdim,tdim,gdim>(Constants::invtJ_0[hana::type_c<float8>][tdim],
                                            f8_jmat,f8_invtJ0_tJK);

                    double f2_met[nnmet], f2_met0[nnmet];
                    float8 f8_met[nnmet], f8_met0[nnmet];
                    msh.met.getMetBary(AsDeg::Pk, DifVar::None, MetSpace::Exp, 
                                       ent2poi[ientt], tdim, bary, f2_met0, NULL);
                    // Transform the metric in quadruple precision. 
                    for(int ii = 0; ii < nnmet; ii++) f8_met0[ii] = (float8) f2_met0[ii];
                    // Get D^-1 M D-1
                    f8_met0[sym2idx(0,0)] *= (float8) (aniso*aniso);
                    for(int ii = 1; ii < gdim; ii++)
                      f8_met0[sym2idx(ii,0)] *= (float8) aniso;

                    // Get RD-1MD-1R^T
                    matXsymXtmat<gdim,gdim,float8,float8,float8>(f8_met0,f8_trans0,f8_met);
                    for(int ii = 0; ii < nnmet; ii++) f2_met[ii] = (double) f8_met[ii];


                    constexpr int nmat = (tdim*(tdim+1))/2;
                    double f2_tJ0_tJK_M_JK_J0[nmat];
                    float8 f8_tJ0_tJK_M_JK_J0[nmat];
                    matXsymXtmat<tdim,gdim,double,double,double>(f2_met,f2_invtJ0_tJK,f2_tJ0_tJK_M_JK_J0);
                    matXsymXtmat<tdim,gdim,float8,float8,float8>(f8_met,f8_invtJ0_tJK,f8_tJ0_tJK_M_JK_J0);
                    float8 errmat = 0;
                    for(int ii = 0; ii < nmat; ii++)
                      errmat += abs(f8_tJ0_tJK_M_JK_J0[ii] - (float8)f2_tJ0_tJK_M_JK_J0[ii]);

                    errsmat += (double) errmat;

                  }// for isamp
                }// for itrans
              }// for ientt

              printf("  -- DONE with aniso ratio %5.1e gdim %d tdim %d\n",aniso,gdim,tdim);
              printf("   - error min = %e avg = %e max = %e (Naive)\n",errsmat.min(),errsmat.avg(),errsmat.max());


            }// for aniso

            // --- Benchmark function without determinant calls
            // i = gdim, j = tdim
            double opps_full[4][4],opps_nodet[4][4];
            double dum = 0;
            double bary[tdim + 1];
            for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1/(tdim + 1.0);

            double t0_full = get_cpu_time();
            for(int ientt = 0; ientt < nentt; ientt++){
              double tra, det;
              quafun(msh, AsDeg::Pk, AsDeg::Pk, ent2poi[ientt], bary, NULL,
                      &tra, &det, QualitySingularityPolicy::Reject);
              dum += tra*det;
            }// for ientt
            double t1_full = get_cpu_time();
            opps_full[gdim][tdim] = nentt / (t1_full - t0_full);

            double t0_nodet = get_cpu_time();
            for(int ientt = 0; ientt < nentt; ientt++){
              double tra, det;
              quafun_nodet(msh, AsDeg::Pk, AsDeg::Pk, ent2poi[ientt], bary, NULL, &tra, &det);
              dum += tra*det;
            }// for ientt
            double t1_nodet = get_cpu_time();
            opps_nodet[gdim][tdim] = nentt / (t1_nodet - t0_nodet);

            printf("-- gdim %d tdim %d ideg %d full %dk/s nodet %dk/s\n",gdim,tdim,msh.curdeg,
              (int) (opps_full[gdim][tdim]/1000),(int)(opps_nodet[gdim][tdim]/1000));

          }}CT_FOR1(ideg);
        }// for idegelev
      }}CT_FOR1(tdim);
    }}CT_FOR1(gdim);
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
  if(msh.met.getSpace() != MetSpace::Log && met_ == NULL) METRIS_THROW_MSG(
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
  if(*tra < 1.0e-16) METRIS_THROW_MSG("Zero trace of spd matrix? "<<*tra);


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

template<int ndim>
void generate_frame(float8 *eigvec,
                    std::uniform_real_distribution<double>& unif, 
                    std::default_random_engine& rng){
  // Generate frame
  while(true){
    for(int jj = 0; jj < ndim; jj++){
     eigvec[0*ndim+jj] = unif(rng);
    }
    float8 nrm = sqrt(getnrml2<ndim>(&eigvec[0*ndim]));

    if(nrm < 1.0e-12) continue;

    for(int jj = 0; jj < ndim; jj++){
     eigvec[0*ndim+jj] = eigvec[0*ndim+jj] / nrm;
    }

    bool ifnd = false;
    for(int jj = 0; jj < ndim; jj++){
      if( abs(eigvec[0*ndim+jj]) > 1.0e-6 ){
        for(int kk = 0; kk < ndim; kk++){
          eigvec[1*ndim+kk] = 0;
        }
        eigvec[1*ndim+(jj+1)%ndim] = -eigvec[0*ndim+jj];
        eigvec[1*ndim+jj]          = eigvec[0*ndim+(jj+1)%ndim];
        nrm = sqrt(getnrml2<ndim>(&eigvec[1*ndim]));
        for(int jj = 0; jj < ndim; jj++){
         eigvec[1*ndim+jj] = eigvec[1*ndim+jj] / nrm;
        }

        ifnd = true;
        break;
      }
    }

    if(!ifnd) continue;

    METRIS_ENFORCE( abs(getprdl2<ndim>(&eigvec[0*ndim],&eigvec[1*ndim])) < 1.0e-12);

    if constexpr(ndim == 3){
      vecprod(&eigvec[0*ndim],&eigvec[1*ndim],&eigvec[2*ndim]);
      float8 nrm = sqrt(getnrml2<ndim>(&eigvec[2*ndim]));
      for(int jj = 0; jj < ndim; jj++){
       eigvec[2*ndim+jj] = eigvec[2*ndim+jj] / nrm;
      }
      METRIS_ENFORCE_MSG( abs(getnrml2<ndim>(&eigvec[2*ndim]) - 1) < 1.0e-12,
        "Not orthogonormal, nrm dif to 1 = "<<abs(getnrml2<ndim>(&eigvec[2*ndim]) - 1));
      METRIS_ENFORCE_MSG( abs(getprdl2<ndim>(&eigvec[0*ndim],&eigvec[1*ndim])) < 1.0e-12,
        "Not orthogonormal, prd = "<<abs(getprdl2<ndim>(&eigvec[0*ndim],&eigvec[1*ndim])));
      METRIS_ENFORCE_MSG( abs(getprdl2<ndim>(&eigvec[2*ndim],&eigvec[1*ndim])) < 1.0e-12,
        "Not orthogonormal, prd = "<<abs(getprdl2<ndim>(&eigvec[2*ndim],&eigvec[1*ndim])));
      METRIS_ENFORCE_MSG( abs(getprdl2<ndim>(&eigvec[2*ndim],&eigvec[0*ndim])) < 1.0e-12,
        "Not orthogonormal, prd = "<<abs(getprdl2<ndim>(&eigvec[2*ndim],&eigvec[0*ndim])));
    }

    break;
  }
}

template 
void generate_frame<2>(float8* eigvec,
                    std::uniform_real_distribution<double>& unif, 
                    std::default_random_engine& rng);
template 
void generate_frame<3>(float8* eigvec,
                    std::uniform_real_distribution<double>& unif, 
                    std::default_random_engine& rng);



}//namespace


#endif
