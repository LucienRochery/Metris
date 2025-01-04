//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "common_includes.hxx"

#include "../low_geosurf.hxx"

namespace Metris{

/*
Metric-based surface optimization routines. 
Explicit interpolation of derivatives is carried out by optimizing with constraints. 
---
EDGES
--- 
Correct CAD-supplied tangents to be coplanar (project in bisection plane passing by straight edge)
P2: intersect tangents; if 1 pt, set there. If +inf, minimize Jacobian conformity i.e. length. 
P3: control points along tangents, minimize Jacobian conformity
P4+: vertex-adjacent control points along tangents. Inner control points free. 

---
FACES
--- 
Vertex-adjacent control points in vertex tangent plane. 
P2: control point either along line or whole plane if surface locally planar
P3+: vertex-adjacent control points along planes, inner completely free

--- 
Jacobian conformity 
--- 
We consider surface entities like traces of volume entities. 
A unit tetrahedron has J_K^TJ_K = M^{-1} everywhere. 
This is relaxed to minimizing tr(J_K^TMJ_K) / det(J_K^T M J_K)
which, if minimized pointwise, does yield a perfectly unit element. 

A tetrahedron face therefore verifies


*/


/*
	Constrained surface smoothing: 
	 - Each control point is constrained by neighbouring vertex tangent planes
	 - Jacobian conformity to the metric is maximized 

	For a surface element, the Jacobian conformity is computed using the Jacobian components in the 
	local tangent plane. 
*/
template<int ideg>
void smoo_surfqua_G1(Mesh &msh){
	int ilag = 0;
	if(msh.ilag == 1){
		ilag = 1;
    setMeshBezier<ideg>(msh);
	}

	int *lbpoi = &msh.poi2iwk[0];
	if(6*msh.nbpoi > msh.poi2rwk.size())
		METRIS_THROW_MSG(DMemExcept(),"Increase size of poi2rwk");
	double *du = &msh.poi2rwk[0];
	double *dv = &msh.poi2rwk[3*msh.nbpoi];






	// Compute normals at vertices. Also list vertices. 
	msh.tag[0]++;
	int nbpof = 0; 
	for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
		if(msh.poi2tag(0,ibpoi) == msh.tag[0]) continue;
		msh.poi2tag(0,ibpoi) = msh.tag[0];

		// Unroll linked list and tag along the way
		int ibpo2 = ibpoi;
		int ibpo1 = -1;
		int nn = 0; 
		do{
			if(nn++ > METRIS_MAX_WHILE) METRIS_THROW_MSG(TopoExcept(), 
				"ill-formed linked list of bpois");

			msh.poi2tag(0,ibpo2) = msh.tag[0];
			if(msh.bpo2ibi[ibpo2][1] == 2) ibpo1 = ibpo2;
		
			ibpo2 = msh.bpo2ibi[ibpo2][3];
		}while(ibpo2 != ibpoi && ibpo2 != -1);

		if(ibpo1 < 0)	METRIS_THROW_MSG(TopoExcept(),"LINE MESH?");
		// At this stage, ibpo1 is a triangle-attached boundary point. 

		// Store for quicker reference (and shorter code, really)
		lbpoi[nbpof] = ibpo1;

		// Store du, dv. 
		double nrmal[3],du[3],dv[3];
		bpo2CADnormal(msh,ibpoi,&du[3*nbpof],&dv[3*nbpof],nrmal);
		nbpof ++;
		// Point allowed to move along planes (du,dv). 

	}

	if(ilag == 1){
    setMeshLagrange<ideg>(msh);
	}
}


} // End namespace
