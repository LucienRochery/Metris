//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __MSH_LAG2BEZ__
#define __MSH_LAG2BEZ__


#include "low_eval.hxx"


namespace Metris{

class MeshBase;


// Work arrays needed of size at least 3npoin
// Also trying to use Mesh's own rwks but these are too small in back mesh
// which is allocated to fit exactly 
// i.e. for front mesh don't bother supplying, only for back
//template<int ideg>
//void setMeshBezier(MeshBase &msh, int nwork = 0, double *rwork = NULL);
//template<int ideg>
//void setMeshLagrange(MeshBase &msh, int nwork = 0, double *rwork = NULL);
//template<int ideg>
//void setMetricBezier(Mesh &msh, int nwork = 0, double *rwork = NULL);
//template<int ideg>
//void setMetricLagrange(Mesh &msh, int nwork = 0, double *rwork = NULL);


// Avoid calling these directly.
template<int ideg, int szfld>
void setFieldBezier(MeshBase &msh, dblAr2 &rfld);
template<int ideg, int szfld>
void setFieldLagrange(MeshBase &msh, dblAr2 &rfld);


// Low-level routines.
template<int ideg, int szfld>
void bez2lag1(const int* __restrict__ lfld,
 							const dblAr2& __restrict__ rfld0,
 							dblAr2& __restrict__ rfld1){
	int nn = edgnpps[ideg];
	double bary[2];
	for(int i = 0; i < nn; i++){
		bary[0] = ordedg.s[ideg][i][0]/((double)ideg);
		bary[1] = ordedg.s[ideg][i][1]/((double)ideg);
		eval1_bezier<szfld,ideg>(rfld0,lfld,DifVar::None,bary,rfld1[lfld[i]],NULL);
	}
}
template<int ideg, int szfld>
void bez2lag2(const int* __restrict__ lfld,
 								const dblAr2& __restrict__ rfld0,
 								dblAr2& __restrict__ rfld1){
	int nn = facnpps[ideg];
	double bary[3];
	for(int i = 0; i < nn; i++){
		bary[0] = ordfac.s[ideg][i][0]/((double)ideg);
		bary[1] = ordfac.s[ideg][i][1]/((double)ideg);
		bary[2] = ordfac.s[ideg][i][2]/((double)ideg);
		eval2_bezier<szfld,ideg>(rfld0,lfld,DifVar::None,bary,rfld1[lfld[i]],NULL);
	}
}
template<int ideg, int szfld>
void bez2lag3(const int* __restrict__ lfld,
 								const dblAr2& __restrict__ rfld0,
 								dblAr2& __restrict__ rfld1){
	int nn = tetnpps[ideg];
	double bary[4];
	for(int i = 0; i < nn; i++){
		bary[0] = ordtet.s[ideg][i][0]/((double)ideg);
		bary[1] = ordtet.s[ideg][i][1]/((double)ideg);
		bary[2] = ordtet.s[ideg][i][2]/((double)ideg);
		bary[3] = ordtet.s[ideg][i][3]/((double)ideg);
		eval3_bezier<szfld,ideg>(rfld0,lfld,DifVar::None,DifVar::None,bary,rfld1[lfld[i]],NULL,NULL);
	}
}


} // End namespace

#endif
