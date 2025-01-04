//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



//#include "low_subdivide.hxx"


///*
//	Store existing subdivision in spl2bar, barys. 
//	Split one of these subdivided entities: isplx. 
//*/
//template <int ideg>
//int subdivref(int *nsplx, intAr2 &spl2bar, 
//	            int *nbary, dblAr2 &barys, 
//	            int isplx){
//	int melem = spl2bar.size()/spl2bar.get_stride();
//	int nnew = subref.ntet[ideg];
//
//	if(*nsplx + nnew > melem) return 1;
//
//	int mpoin = barys.size()/barys.get_stride();
//	if(*nbary + tetnpps[ideg] > mpoin) return 2;
//
//	double *bary1 = barys[spl2bar[isplx][0]];
//	double *bary2 = barys[spl2bar[isplx][1]];
//	double *bary3 = barys[spl2bar[isplx][2]];
//	double *bary4 = barys[spl2bar[isplx][3]];
//	int npoi0 = (*npoin);
//
//	for(int ipnew = 0; ipnew < tetnpps[ideg]; ipnew++){
//		for(int ii = 0; ii < 3; ii++){
//			barys[*nbary][ii] = ordtet.s[ideg][ipnew][0]*bary1[ii]/ideg
//												+ ordtet.s[ideg][ipnew][1]*bary2[ii]/ideg
//												+ ordtet.s[ideg][ipnew][2]*bary3[ii]/ideg
//												+ ordtet.s[ideg][ipnew][3]*bary4[ii]/ideg;
//		}
//		(*nbary)++;
//	}
//
//	for(int itnew = 0; itnew < nnew; itnew++){
//		spl2bar[*nsplx][0] = npoi0 + subdiv.tet[ideg][itnew][0];
//		spl2bar[*nsplx][1] = npoi0 + subdiv.tet[ideg][itnew][1];
//		spl2bar[*nsplx][2] = npoi0 + subdiv.tet[ideg][itnew][2];
//		spl2bar[*nsplx][3] = npoi0 + subdiv.tet[ideg][itnew][3];
//		(*nsplx)++;
//	}
//
//	return 0;
//}
