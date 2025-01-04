//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../quality/low_surfqua.hxx"

#include "../linalg/matprods.hxx"
#include "../linalg/explogmet.hxx"

namespace Metris{
// Using tr(J^TMJ)^3 / det(J^TMJ)
//       tr(J^TMJ)   / det(J^TMJ)^{1/3} is another common choice
// as well as tr(J^TMJ)^{3/2} / det(J^TMJ)^{1/2}
template <int ideg, int ilag, typename ftype>
void metquaS_xi(const Mesh &msh, int ielem, int power, 
	              const double* bary, ftype* qufac){
	assert(power != 0);
	if(msh.ilogmet != 1)METRIS_THROW_MSG(WArgExcept(),
		"## SET MESH METRIC TO LOG BEFORE CALLING metquaS_xi");
	double jmat[9],met[6],dum[3],lmet[6];
	// Get Jacobian matrix at xi
	eval3<3,ideg,ilag,1,0>(msh.coord,msh.tet2poi[ielem],bary,dum,jmat,NULL);

	// Get metric interpolated at xi
	if(msh.ilagmet == 1){
		eval3<6,ideg,1,0,0>(msh.met,msh.tet2poi[ielem],bary,lmet,NULL,NULL);
	}else{
		eval3<6,ideg,0,0,0>(msh.met,msh.tet2poi[ielem],bary,lmet,NULL,NULL);
	}

	getexpmet_cpy<3>(lmet,met);

	// Compute J_0^{-T} J_K^T M J_K J_0^{-1}
	// Starting with J_K J_0^{-1}
	// Note that J_K is stored transposed w.r.t. above 
	// The below J_0^{-1} is also transposed. 
  ftype invtJ_0[3][3] = {
      -0.57735026918962562,-0.57735026918962562 , 0                 ,
      1                   ,-1                   , 0                 ,
      -0.40824829046386302,-0.4082482904638630  , 1.2247448713915889};

  // Get J_0^{-T} J_K^T
  ftype invtJ_0tJ_K[9];
  matXmat<3>(invtJ_0[0],jmat,invtJ_0tJ_K);

  // 2. Get whole product
	//double J0tJtMJJ0[6];
	//matXsymXtmat<3>(met, invtJ_0tJ_K, J0tJtMJJ0);
	//double tra = J0tJtMJJ0[0] + J0tJtMJJ0[2] + J0tJtMJJ0[5];

	//ftype J0tJtMJJ0_diag[3];
	ftype tra = tra_matXsymXtmat<3, double, ftype, ftype>(met, invtJ_0tJ_K);
	//ftype tra = J0tJtMJJ0_diag[0] + J0tJtMJJ0_diag[1] + J0tJtMJJ0_diag[2];
	// Unfortunately, the only numerically stable possibility is this. 
	// As assessed using octuple precision. Following methods were tested (example values):
	// detsym<3>(J0tJtMJJ0) = -2.4724652e+01
	// detmat<3>(J0tJtMJJ0) = -2.4724652e+01 (having copied to a 3x3 matrix)
	// -> det1*det1*detsym<3>(met) = 8.0969526e-01
	// detmat<3>(A) =  7.2640309e-01 
	//   with A being J0tJtMJJ0 computed using matXmat<3> then matXtmat<3> instead of matXsymXtmat
	// -> detsym<3>(J0tJtMJJ0) = 8.0969526e-01 (using octuple precision)
//	double det = detsym<3>(J0tJtMJJ0);
  ftype det1 = detmat<3>(invtJ_0tJ_K);
	ftype det = det1*det1*detsym<3>(met);



//
//
//	double tra =    jmat[3*0+0]*jmat[3*0+0]*met[0] 
//	            +   jmat[3*0+1]*jmat[3*0+1]*met[2]
//	            +   jmat[3*0+2]*jmat[3*0+2]*met[5]
//	            + 2*jmat[3*0+0]*jmat[3*0+1]*met[1]
//	            + 2*jmat[3*0+0]*jmat[3*0+2]*met[3]
//	            + 2*jmat[3*0+1]*jmat[3*0+2]*met[4]
//       
//              +   jmat[3*1+0]*jmat[3*1+0]*met[0] 
//	            +   jmat[3*1+1]*jmat[3*1+1]*met[2]
//	            +   jmat[3*1+2]*jmat[3*1+2]*met[5]
//	            + 2*jmat[3*1+0]*jmat[3*1+1]*met[1]
//	            + 2*jmat[3*1+0]*jmat[3*1+2]*met[3]
//	            + 2*jmat[3*1+1]*jmat[3*1+2]*met[4]
//       
//              +   jmat[3*2+0]*jmat[3*2+0]*met[0] 
//	            +   jmat[3*2+1]*jmat[3*2+1]*met[2]
//	            +   jmat[3*2+2]*jmat[3*2+2]*met[5]
//	            + 2*jmat[3*2+0]*jmat[3*2+1]*met[1]
//	            + 2*jmat[3*2+0]*jmat[3*2+2]*met[3]
//	            + 2*jmat[3*2+1]*jmat[3*2+2]*met[4];
	if(tra < 1.0e-16) METRIS_THROW_MSG(GeomExcept(),
		"NEGATIVE J^TMJ trace "<<tra);
	if(det < 1.0e-16){

		ftype J0tJtMJJ0[6];
		matXsymXtmat<3>(met, invtJ_0tJ_K, J0tJtMJJ0);

    ftype tmpMat[3][3];
    tmpMat[0][0] = J0tJtMJJ0[0];

    tmpMat[0][1] = J0tJtMJJ0[1];
    tmpMat[1][0] = J0tJtMJJ0[1];

    tmpMat[1][1] = J0tJtMJJ0[2];

    tmpMat[0][2] = J0tJtMJJ0[3];
    tmpMat[2][0] = J0tJtMJJ0[3];

    tmpMat[1][2] = J0tJtMJJ0[4];
    tmpMat[2][1] = J0tJtMJJ0[4];

    tmpMat[2][2] = J0tJtMJJ0[5];


    // Recompute using matXmat<3> instead of matXsymXtmat
    ftype met33[3][3];
    met33[0][0] = met[0];

    met33[0][1] = met[1];
    met33[1][0] = met[1];

    met33[1][1] = met[2];

    met33[0][2] = met[3];
    met33[2][0] = met[3];

    met33[1][2] = met[4];
    met33[2][1] = met[4];

    met33[2][2] = met[5];


    ftype J0tJtMJJ02[9],tmpMat2[3][3];
    matXmat<3>(invtJ_0tJ_K,met33[0],tmpMat2[0]);
    matXtmat<3>(tmpMat2[0],invtJ_0tJ_K,J0tJtMJJ02);

    ftype det1 = detmat<3,ftype>(invtJ_0tJ_K);
    printf("Got det = %15.7e as 3x3 = %15.7e",(double)det,(double)detmat<3,ftype>(tmpMat[0]));
    printf(" as product of dets = %15.7e as 3x3 using matXmat<3> %15.7e\n",
    	(double)(det1*det1*detsym<3>(met)),(double)detmat<3,ftype>(J0tJtMJJ02));

    printf("J0tJtMJJ0 original:\n");
    MeshArray2D<ftype>(3,3,tmpMat[0]).print();
    printf("J0tJtMJJ0 matmat\n");
    MeshArray2D<ftype>(3,3,J0tJtMJJ02).print();
    printf("Debug diff between original and recomputed \n");
    //for(int i = 0; i < 3; i++){
    //	for(int j = 0 ; j < 3; j++){
    //		printf("%15.7e \n",tmpMat[i][j] - J0tJtMJJ02[3*i+j]);
    //	}
    //}

    float8 jmat_8[9], invtJ_0_8[9], met_8[6];
    for(int i = 0; i < 3 ; i++){
    	for(int j = 0; j < 3 ; j++){
    		jmat_8[3*i+j] = jmat[3*i+j];
    		invtJ_0_8[3*i+j] = invtJ_0[i][j];
    		if(3*i+j < 6) met_8[3*i+j] = met[3*i+j];
    	}
    }
  	float8 invtJ_0tJ_K_8[9];
    matXmat<3,float8,float8,float8>(invtJ_0_8,jmat_8,invtJ_0tJ_K_8);
		float8 J0tJtMJJ0_8[6];
		matXsymXtmat<3,float8,float8,float8>(met_8, invtJ_0tJ_K_8, J0tJtMJJ0_8);
	
		float8 tra_8 = J0tJtMJJ0_8[0] + J0tJtMJJ0_8[2] + J0tJtMJJ0_8[5];
		float8 det_8 = detsym<3,float8>(J0tJtMJJ0_8);

		std::cout << "Using oct precision: "<<
		std::setprecision(std::numeric_limits<float8>::max_digits10)
       << "tra = " << tra_8 << " det = " << det_8 << std::endl;


		METRIS_THROW_MSG(GeomExcept(),
    "NEGATIVE met det  "<<det<<" met det "<<detsym<3>(met)
    <<" prod det "<<det1<<" met = "<<
    met[0]<<" "<<met[1]<<" "<<met[2]<<" "<<
    met[3]<<" "<<met[4]<<" "<<met[5]<<" ielem = "<<ielem<<"\n");
	} 
	//matXsymXtmat<3>(met, jmat, tJMJ);

	// Now we compute the trace and determinant 
	//double tra = tJMJ[0] + tJMJ[2] + tJMJ[5];
	//double det = detsym<3>(tJMJ);

	// The trace is the arithmetic average of eigvals ; tr =  3 at the min
	//       det        geometric                      det =  1 
	// Thus tr^3 / det > 27
	// And the magic factor 2  !

	if(power > 0){
		*qufac = pow(tra*tra*tra/det/27,power);
	}else{
		*qufac = pow(27*det/(tra*tra*tra),-power);
	}
}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void metquaS_xi< n ,0, double>(const Mesh&,int,int,const double*,double*);\
template void metquaS_xi< n ,1, double>(const Mesh&,int,int,const double*,double*);\
template void metquaS_xi< n ,0, float8>(const Mesh&,int,int,const double*,float8*);\
template void metquaS_xi< n ,1, float8>(const Mesh&,int,int,const double*,float8*);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


} // End namespace
