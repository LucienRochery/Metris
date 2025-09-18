//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_surf

#include "common_setup.hxx"

#include <random>

using namespace Metris;

typedef MetricFieldAnalytical MFT;

void test_match(Mesh<MFT> &msh2D, Mesh<MFT> &msh3D);
void test_normals(Mesh<MFT> &msh3D);
void test_getnordev(Mesh<MFT> &msh3D);
void test_measure(Mesh<MFT> &msh2D, Mesh<MFT> &msh3D);
void test_increase_cavity(Mesh<MFT> &msh2D, Mesh<MFT> &msh3D);

// Take a 2D mesh, elevate it to 3D by adding 0 along z. 
// Then test various functions give same outputs in 2D and surf.
BOOST_AUTO_TEST_CASE(test_surf) 
{

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char** argv = boost::unit_test::framework::master_test_suite().argv;
  MetrisOptions cmdopt;
  cmdopt.parse(argc,argv);
  MetrisParameters cmdparam(cmdopt);
  GETVDEPTH((&cmdparam));
  
  const int nsamp = 100;
  dblAr2 bary(nsamp,3);
  dblAr2 bary_out(nsamp,3);
  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(0);
  for(int ii = 0; ii < nsamp; ii++){
    double sum = 0;
    do{
      for(int jj = 0; jj < 3; jj++){
        bary(ii,jj) = unif(rng);
        sum += bary(ii,jj);
      }
    }while(abs(sum) < 1.0e-16);

    for(int jj = 0; jj < 3; jj++){
      bary(ii,jj) /= sum;
    }
  }
  std::uniform_real_distribution<double> unif_out(-100.0,100.0);
  std::default_random_engine rng_out(0);
  for(int ii = 0; ii < nsamp; ii++){
    double sum = 0;
    do{
      for(int jj = 0; jj < 3; jj++){
       bary_out(ii,jj) = unif_out(rng_out);
       sum += bary_out(ii,jj);
      }
    }while(abs(sum) < 1.0e-16);

    for(int jj = 0; jj < 3; jj++){
      bary_out(ii,jj) /= sum;
    }
  }

  std::string mesh = METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k";
  std::string CADf = METRIS_CASES_DIR "/unit/2D/square/square.egads";

  // Using boundary-layer mesh along x centered at 0.5. This is 3 in 2D and 10 in 3D.
  cargHandler arg2D("-in " + mesh + " -cad " + CADf + " -anamet 3 -sclmet 1 -vdepth 0 -verb 0 -opt-niter 0 -adp-opt-niter 0 -adapt 0");
  MetrisRunner run2D(arg2D.c, arg2D.v);
  //run2D.param->ivdepth = MAX(run2D.param->ivdepth, cmdparam.ivdepth);
  //run2D.param->iverb   = MAX(run2D.param->iverb  , cmdparam.iverb  );
  Mesh<MFT> &msh2D = *((Mesh<MFT>*) run2D.msh_g);

  msh2D.cleanup();


  std::cout<<"\n\n------------------------------------------------\n";
  std::cout<<"Mesh "<<mesh<<"\n";
  std::cout<<"------------------------------------------------\n";


  // We can add a dimension to coord, it shouldn't change anything as 
  // long as all we don't flatten (IO only)
  msh2D.coord.set_stride(3);
  for(int ipoin = 0; ipoin < msh2D.npoin; ipoin++){
    msh2D.coord(ipoin,2) = 0;
  }

  msh2D.idim = 3;
  msh2D.set_nelem(0);
  writeMesh(mesh+".3D",msh2D);


  cargHandler arg3D("-in " + mesh + ".3D -cad " + CADf + " -anamet 10 -sclmet 1 -vdepth 0 -verb 0 -opt-niter 0 -adp-opt-niter 0 -adapt 1");
  MetrisRunner run3D(arg3D.c, arg3D.v);
  //run3D.param->ivdepth = MAX(run3D.param->ivdepth, cmdparam.ivdepth);
  //run3D.param->iverb   = MAX(run3D.param->iverb  , cmdparam.iverb  );
  Mesh<MFT> &msh3D = *((Mesh<MFT>*) run3D.msh_g);

  msh2D.coord.set_stride(2);
  msh2D.idim = 2;

  check_topo(msh2D, 0);
  check_topo(msh3D, 0);

  test_match(msh2D,msh3D);

  test_normals(msh3D);
  test_getnordev(msh3D);
  test_measure(msh2D,msh3D);
  test_increase_cavity(msh2D,msh3D);

}// end boost test case


void test_match(Mesh<MFT> &msh2D, Mesh<MFT> &msh3D){
  const int ideg = msh2D.curdeg;
  BOOST_REQUIRE(ideg == msh3D.curdeg);
  BOOST_REQUIRE(msh2D.idim == 2);
  BOOST_REQUIRE(msh3D.idim == 3);
  const int nnode = msh2D.nnode(2);
  BOOST_REQUIRE(msh2D.nface == msh3D.nface);

  for(int iface = 0; iface < msh2D.nface; iface++){
    if(isdeadent(iface, msh2D.fac2poi)) continue;
    for(int ii = 0; ii < nnode; ii++){
      int ipoin = msh2D.fac2poi(iface,ii);
      BOOST_REQUIRE(ipoin == msh3D.fac2poi(iface,ii));
      BOOST_REQUIRE(msh2D.coord(ipoin,0) == msh3D.coord(ipoin,0));
      BOOST_REQUIRE(msh2D.coord(ipoin,1) == msh3D.coord(ipoin,1));
      BOOST_REQUIRE(msh3D.coord(ipoin,2) == 0.0);
    }
  }

  // Check the 3D mesh is not tagged as periodic or non-manifold
  fmt::print("Debug msh3D.isperiodic_face = {} nperiodic_face = {}\n", msh3D.isperiodic_face,msh3D.nperiodic_face);
  BOOST_REQUIRE(msh3D.isperiodic_face.get_n() == 1);
  BOOST_REQUIRE(!msh3D.isperiodic_face[0]);
  BOOST_REQUIRE(msh3D.nperiodic_face == 0);
  
}


void test_measure(Mesh<MFT> &msh2D, Mesh<MFT> &msh3D){
  
  // Check function getmeasentP1 gives same results.
  // With very low nordev tol, the surface version should still return valid
  // as the surface is flat.
  for(int iface = 0; iface < msh2D.nface; iface++){
    if(isdeadent(iface, msh2D.fac2poi)) continue;
    bool iflat2D, iflat3D;
    double meas2D = getmeasentP1<2,2>(msh2D,iface,NULL,&iflat2D,-1);
    double meas3D = getmeasentP1<3,2>(msh3D,iface,NULL,&iflat3D,1.0e-12);
    BOOST_CHECK_CLOSE(meas2D, meas3D, Defaults::vtol);
    BOOST_CHECK(iflat2D == iflat3D);
  }
  fmt::print("-- SUCCESS getmeasentP1, tested {} faces\n",msh2D.nface);

  // Check isvalideltP1, same tests.
  for(int iface = 0; iface < msh2D.nface; iface++){
    if(isdeadent(iface, msh2D.fac2poi)) continue;
    double meas2D, meas3D;
    bool iflat2D = !isvalideltP1<2,2>(msh2D,iface,NULL,&meas2D,-1);
    bool iflat3D = !isvalideltP1<3,2>(msh3D,iface,NULL,&meas3D,1.0e-12);
    BOOST_CHECK_CLOSE(meas2D, meas3D, Defaults::vtol);
    BOOST_CHECK(iflat2D == iflat3D);
  }
  fmt::print("-- SUCCESS isvalideltP1, tested {} faces\n",msh2D.nface);

}

// Check discrete and CAD normals coincide (up to norm).
void test_normals(Mesh<MFT> &msh3D){
  double norfac[3], norCAD[3];
  for(int iface = 0; iface < msh3D.nface; iface++){
    if(isdeadent(iface, msh3D.fac2poi)) continue; 
    getnorfacP1(msh3D.fac2poi[iface], msh3D.coord, norfac);
    getnorfacCAD(msh3D, iface, norCAD);
    BOOST_REQUIRE(!normalize_vec<3>(norfac));
    BOOST_REQUIRE(!normalize_vec<3>(norCAD));
    BOOST_CHECK_CLOSE(norfac[0], norCAD[0], Defaults::vtol);
    BOOST_CHECK_CLOSE(norfac[1], norCAD[1], Defaults::vtol);
    BOOST_CHECK_CLOSE(norfac[2], norCAD[2], Defaults::vtol);
  }
  fmt::print("-- SUCCESS normals, tested {} faces\n", msh3D.nface);
}


// Check getnordev returns 0 (flat surface).
void test_getnordev(Mesh<MFT> &msh3D){
  double max_nordev = -1;
  for(int iface = 0; iface < msh3D.nface; iface++){
    if(isdeadent(iface, msh3D.fac2poi)) continue;
    double nordev;
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh3D.curdeg){
      nordev = getnordev<ideg>(msh3D, iface, NULL);
    }}CT_FOR1(ideg);
    BOOST_CHECK_CLOSE(nordev, 0.0, Defaults::vtol);
    max_nordev = MAX(nordev, max_nordev);
  }
  fmt::print("-- SUCCESS getnordev, tested {} faces, max nordev = {:e}\n", msh3D.nface, max_nordev);
}

// Test increase_cavity yields same results in 2D and flat surf
void test_increase_cavity(Mesh<MFT> &msh2D, Mesh<MFT> &msh3D){

  GETVDEPTH(msh2D.param);

  std::mt19937 rng(0); // seed the generator
  std::uniform_int_distribution<> distr(0, msh3D.nface-1);

  const int ntest = 100;

  MshCavity cav2D(0,100,0), cav2D_validity(0,100,0), cav3D(0,100,0);
  double result[18];
  int ntotc = 0, ierro;
  for(int itest = 0; itest < ntest; itest++){
    int iface = -1;
    for(; iface < 0 || isdeadent(iface,msh3D.fac2poi); iface = distr(rng)){}
    int ifac2 = -1;
    for(; ifac2 == iface || ifac2 < 0; ifac2 = distr(rng)){}
    // Generate a point in iface, then seed cavity with ifac2.
    int ip3D = msh3D.newpoitopo(PointType::Vertex,2, iface);
    int ib3D = msh3D.newbpotopo(Vertex{ip3D}, 2, iface);

    int ib1 = msh3D.poi2ebp(msh3D.fac2poi(iface,0), 2, iface, msh3D.fac2ref[iface]);
    int ib2 = msh3D.poi2ebp(msh3D.fac2poi(iface,1), 2, iface, msh3D.fac2ref[iface]);
    int ib3 = msh3D.poi2ebp(msh3D.fac2poi(iface,2), 2, iface, msh3D.fac2ref[iface]);
    BOOST_REQUIRE(ib1 >= 0 && ib2 >= 0 && ib3 >= 0);

    for(int ii = 0; ii < nrbi; ii++) 
      msh3D.bpo2rbi(ib3D, ii) = ( msh3D.bpo2rbi(ib1,ii)
                                + msh3D.bpo2rbi(ib2,ii)
                                + msh3D.bpo2rbi(ib3,ii))/3;

    ego obj = msh3D.CAD.cad2fac[msh3D.fac2ref[iface]];
    ierro = EG_evaluate(obj, msh3D.bpo2rbi[ib3D], result);
    BOOST_REQUIRE(ierro == EGADS_SUCCESS);


    int ip2D = msh2D.newpoitopo(PointType::Vertex,2, iface);
    msh2D.coord(ip2D,0) = msh3D.coord(ip3D,0) = result[0];
    msh2D.coord(ip2D,1) = msh3D.coord(ip3D,1) = result[1];
                          msh3D.coord(ip3D,2) = result[2];
    BOOST_REQUIRE(msh3D.coord(ip3D,2) == 0.0);


    msh2D.interpMetBack(ip2D);
    msh3D.interpMetBack(ip3D);
    BOOST_CHECK_CLOSE(msh2D.met(ip2D,0), msh3D.met(ip3D,0), 1.0e-15);
    BOOST_CHECK_CLOSE(msh2D.met(ip2D,1), msh3D.met(ip3D,1), 1.0e-15);
    BOOST_CHECK_CLOSE(msh2D.met(ip2D,2), msh3D.met(ip3D,2), 1.0e-15);

    // Now seed a cavity with ifac2.
    cav2D.reset();
    cav2D_validity.reset();
    cav3D.reset();

    cav2D.ipins = ip2D;
    cav2D_validity.ipins = ip2D;
    cav3D.ipins = ip3D;

    cav2D.lcfac.stack(ifac2);
    cav2D_validity.lcfac.stack(ifac2);
    cav3D.lcfac.stack(ifac2);

    
    increase_cavity_validity(msh2D, cav2D_validity, 0);
    increase_cavity_validity(msh3D, cav3D, 0);

    BOOST_REQUIRE(cav2D_validity.lcfac.get_n() == cav3D.lcfac.get_n());
    BOOST_REQUIRE(cav2D_validity.lctet.get_n() == cav3D.lctet.get_n());
    BOOST_REQUIRE(cav2D_validity.lcedg.get_n() == cav3D.lcedg.get_n());

    for(int ii = 0; ii < cav2D_validity.lcfac.get_n(); ii++){
      BOOST_REQUIRE(cav2D_validity.lcfac[ii] == cav3D.lcfac[ii]);
    }
    ntotc += cav2D_validity.lcfac.get_n();


    // Reset cavity to seed, then test other increase_cav routine
    cav2D.lcfac.set_n(1);
    cav3D.lcfac.set_n(1);

    increase_cavity(msh2D, cav2D, false, 0, 1);
    increase_cavity(msh3D, cav3D, false, 0, 1);

    BOOST_REQUIRE(cav2D.lcfac.get_n() == cav3D.lcfac.get_n());
    BOOST_REQUIRE(cav2D.lctet.get_n() == cav3D.lctet.get_n());
    BOOST_REQUIRE(cav2D.lcedg.get_n() == cav3D.lcedg.get_n());

    for(int ii = 0; ii < cav2D.lcfac.get_n(); ii++){
      BOOST_REQUIRE(cav2D.lcfac[ii] == cav3D.lcfac[ii]);
    }

    ntotc += cav2D.lcfac.get_n();

    // Check these two coincide, given that idelaunay was false.
    BOOST_REQUIRE(cav2D.lcfac.get_n() == cav2D_validity.lcfac.get_n());
    for(int ii = 0; ii < cav2D.lcfac.get_n(); ii++){
      BOOST_REQUIRE(cav2D.lcfac[ii] == cav2D_validity.lcfac[ii]);
    }

    
    // Test Delaunay 2D <==> 3D
    cav2D.lcfac.set_n(1);
    cav3D.lcfac.set_n(1);


    CPRINTF1("\n\n-- Call 2D Delaunay\n");
    {
      INCVDEPTH(msh2D.param);
      increase_cavity_Delaunay(msh2D, cav2D, 2, 10, 1);
    }
    CPRINTF1("\n\n-- Call 3D Delaunay\n");
    {
      INCVDEPTH(msh3D.param);
      increase_cavity_Delaunay(msh3D, cav3D, 2, 10, 1);
    }
    
    if(cav2D.lcfac.get_n() != cav3D.lcfac.get_n()){
      PRINTF("Debug cav2D:\n");
      cav2D.print(msh2D, 10);
      PRINTF("Debug cav3D:\n");
      cav3D.print(msh3D, 10);
    }
    BOOST_REQUIRE(cav2D.lcfac.get_n() == cav3D.lcfac.get_n());
    BOOST_REQUIRE(cav2D.lctet.get_n() == cav3D.lctet.get_n());
    BOOST_REQUIRE(cav2D.lcedg.get_n() == cav3D.lcedg.get_n());

    for(int ii = 0; ii < cav2D.lcfac.get_n(); ii++){
      BOOST_REQUIRE(cav2D.lcfac[ii] == cav3D.lcfac[ii]);
    }


    // Test Delaunay + validity 2D <==> 3D
    cav2D.lcfac.set_n(1);
    cav3D.lcfac.set_n(1);

    increase_cavity(msh2D, cav2D, true, 0, 1);
    increase_cavity(msh3D, cav3D, true, 0, 1);
    
    BOOST_REQUIRE(cav2D.lcfac.get_n() == cav3D.lcfac.get_n());
    BOOST_REQUIRE(cav2D.lctet.get_n() == cav3D.lctet.get_n());
    BOOST_REQUIRE(cav2D.lcedg.get_n() == cav3D.lcedg.get_n());

    for(int ii = 0; ii < cav2D.lcfac.get_n(); ii++){
      BOOST_REQUIRE(cav2D.lcfac[ii] == cav3D.lcfac[ii]);
    }



    //fmt::print("Debug cav2D:\n");
    //cav2D.print(msh2D, 10);
    //fmt::print("Debug cav3D:\n");
    //cav3D.print(msh3D, 10);
    //wait();
    

    msh3D.killpoint(ip3D);
    msh3D.set_npoin(ip3D);
    msh3D.set_nbpoi(ib3D);
  }
  fmt::print("-- SUCCESS increase_cavity(idelaunay = true/false), increase_cavity_validity, increase_cavity_Delaunay, tested {} faces, total cav ent {}\n",ntest,ntotc);
}