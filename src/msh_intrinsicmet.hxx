//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_INTRINSICMET__
#define __METRIS_MSH_INTRINSICMET__


#include "aux_utils.hxx"
#include "Mesh/MeshMetric.hxx"

#include "types.hxx"
//#include "ho_constants.hxx"
#include "aux_topo.hxx"
#include "low_lenedg.hxx"


namespace Metris{

// -----------------------------------------------------------------------------
template<class MetricFieldType,int ideg>
void getMetMesh(const MetrisParameters &param, MeshMetric<MetricFieldType> &msh);


template<class MetricFieldType, int gdim, int tdim, int ideg>
void getMetMesh0_lplib(int ient0, int ient1,int ithread, MeshMetric<MetricFieldType> *msh_, int poitag);

//// This reoutine re-interpolates at edge extremities from the shell
//template<int ideg,int ilag>
//void get_met_shell(Mesh &msh, int ip1, int ip2, int ielem);



// These are old routines that probably don't work as intended (pre metric refactor, gdim)
// -----------------------------------------------------------------------------
template <class MetricFieldType,int ideg>
void getmshedglen(MeshMetric<MetricFieldType> &msh, HshTabDbl2 &rlened, int nquad){
	rlened.reserve((int)(6.5*msh.nelem/5.5));

	int edg2pol[ideg+1];

	for(int ielem = 0; ielem < msh.nelem; ielem++){
		if(isdeadent(ielem,msh.tet2poi)) continue;
		for(int ied = 0; ied < 3; ied++){
			int ip1 = msh.tet2poi(ielem,lnoed3[ied][0]);
			int ip2 = msh.tet2poi(ielem,lnoed3[ied][1]);

			auto key = stup2(ip1,ip2);
			auto t = rlened.find(key);
			if(t != rlened.end()) continue;
			edg2pol[0] = ip1;
			edg2pol[1] = ip2;
			int idx0 = 4 + ied*(ideg-1);
			for(int i = 0; i < ideg-1; i++){
				edg2pol[2+i] = msh.tet2poi[ielem][idx0+i];
			}

			double len = getlenedg_quad<MetricFieldType,3,ideg>(edg2pol,msh.coord,msh.met,nquad);

			rlened.insert({key,len});
		}
	}
}



// -----------------------------------------------------------------------------
template <class MetricFieldType,int ideg>
void getmshedglen_shell(MeshMetric<MetricFieldType> &msh, HshTabDbl2 &rlened, int nquad){
	rlened.reserve((int)(6.5*msh.nelem/5.5));
	int edg2pol[ideg+1];

	for(int ielem = 0; ielem < msh.nelem; ielem++){
		if(isdeadent(ielem,msh.tet2poi)) continue;
		for(int ied = 0; ied < 3; ied++){
			int ip1 = msh.tet2poi(ielem,lnoed3[ied][0]);
			int ip2 = msh.tet2poi(ielem,lnoed3[ied][1]);

			auto key = stup2(ip1,ip2);
			auto t = rlened.find(key);
			if(t != rlened.end()) continue;

			//ierro = get_met_shell<ideg>(msh,nei,ip1,ip2,ielem);
			//if(ierro > 0){
			//	printf("## get_met_shell FAILED %d \n",ierro);
			//	return 1;
			//}

			// Debug

			edg2pol[0] = ip1;
			edg2pol[1] = ip2;
			int idx0 = 4 + ied*(ideg-1);
			for(int i = 0; i < ideg-1; i++){
				edg2pol[2+i] = msh.tet2poi[ielem][idx0+i];
			}
			for(int j = 0;j < 6; j++){
				msh.met[ip1][j] = msh.met[msh.tet2poi(ielem,idx0)][j];
				msh.met[ip2][j] = msh.met[msh.tet2poi(ielem,idx0)][j];
			}

			double len = getlenedg_quad<MetricFieldType,3,ideg>(edg2pol,msh.coord,msh.met,nquad);

			rlened.insert({key,len});
		}
	}
}

//// Compute anisotropic quality tr/det of J_K^T M J_K. 
//template <int ideg>
//void getmshqua(MeshBase &msh, std::vector<double> &quael){
//	quael.resize(msh.nelem);
//	for(int ielem = 0; ielem < msh.nelem; ielem++){
//
//	}
//}
} // End namespace


#endif

