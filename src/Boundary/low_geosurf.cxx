//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "low_geosurf.hxx"

#include "../msh_structs.hxx"
#include "../Mesh/Mesh.hxx"

#include "../aux_utils.hxx"
#include "../aux_exceptions.hxx"
#include "../low_geo.hxx"

#include <egadsTypes.h>
#include <egads.h>

#include <assert.h>
#include <stddef.h>


namespace Metris{

// If tol > 0 and the distance to the projected is larger than tol, return error.
void projptsurf(MeshBase &msh, int ibpoi, double *coop, double tol){
	assert(ibpoi >=0  && ibpoi < msh.nbpoi);

	int ipoin = msh.bpo2ibi(ibpoi,0);
	assert("bpo2ibi valid" && ipoin >= 0 && ipoin < msh.npoin);

	int ityp = msh.bpo2ibi(ibpoi,1);
	if(ityp == 0){
		for(int i=0; i<3; i++) coop[i] = msh.coord(ipoin,i);
		return;
	}
	int ientt = msh.bpo2ibi(ibpoi,2);
	assert(ientt >= 0);
	assert(ientt < msh.nface && ityp == 2 
     || ientt < msh.nedge && ityp == 1  );

	int iref;
	ego obj;
	if(ityp == 1){
		iref = msh.edg2ref[ientt]; 
		if(iref < 0){
			projptsurf_disc(msh, ibpoi, coop);
			return;
		} 
		obj = msh.CAD.cad2edg[iref];
		assert(obj!=NULL);
	}else{
		iref = msh.fac2ref[ientt]; 
		if(iref < 0){
			projptsurf_disc(msh, ibpoi, coop);
			return;
		} 
		obj = msh.CAD.cad2fac[iref];
		assert(obj!=NULL);
	}

	double result[18];
	int ierro = EG_invEvaluate(obj, msh.coord[ipoin], msh.bpo2rbi[ibpoi], result);
	if(ierro != 0){
   print_EGADS_error("EG_getBodyTopos (EDGE)",ierro);
		METRIS_THROW(TopoExcept());
	}
	double dist = geterrl2<3>(msh.coord[ipoin],result);
	if(tol > 0.0 && dist > tol*tol) METRIS_THROW(RealExcept());

	coop[0] = result[0];
	coop[1] = result[1];
	coop[2] = result[2];
}

void bpo2CADnormal(MeshBase &msh, int ibpoi, double *du, double *dv, double *nrmal){
  if(ibpoi < 0 || ibpoi >= msh.nbpoi) METRIS_THROW(WArgExcept());
  if(!msh.CAD()) METRIS_THROW_MSG(TopoExcept(), "CAD not initialized");
  if(msh.bpo2ibi(ibpoi,1) != 2) METRIS_THROW_MSG(WArgExcept(),"Point not attached to CAD face")

  int iface = msh.bpo2ibi(ibpoi,2);
  if(iface < 0 || iface >= msh.nface) METRIS_THROW(TopoExcept());
  int iref = msh.fac2ref[iface];
  if(iref < 0 || iref >= msh.CAD.ncadfa) METRIS_THROW(TopoExcept());

  ego obj = msh.CAD.cad2fac[iref];

  double result[18];
  int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
  if(ierro != 0){
   print_EGADS_error("EG_getBodyTopos (EDGE)",ierro);
   METRIS_THROW(TopoExcept());
  }

  du[0] = result[3];
  du[1] = result[4];
  du[2] = result[5];
  dv[0] = result[6];
  dv[1] = result[7];
  dv[2] = result[8];

  nrmal[0] = du[1]*dv[2];
  nrmal[1] = du[2]*dv[0];
  nrmal[2] = du[0]*dv[1];
}

int projveredg(MeshBase &msh, int iedge, int iver, double *coop);
int projverfac(MeshBase &msh, int iface, int iver, double *coop);

} // End namespace
