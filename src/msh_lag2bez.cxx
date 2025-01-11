//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_lag2bez.hxx"
#include "aux_topo.hxx"
#include "codegen_lag2bez.hxx"
#include "Mesh/MeshBase.hxx"


namespace Metris{



//template<> void setFieldBezier<1,2>(MeshBase &msh, dblAr2 &rfld){}
//template<> void setFieldBezier<1,3>(MeshBase &msh, dblAr2 &rfld){}
//template<> void setFieldBezier<1,6>(MeshBase &msh, dblAr2 &rfld){}
//template<> void setFieldLagrange<1,2>(MeshBase &msh, dblAr2 &rfld){}
//template<> void setFieldLagrange<1,3>(MeshBase &msh, dblAr2 &rfld){}
//template<> void setFieldLagrange<1,6>(MeshBase &msh, dblAr2 &rfld){}



/*
// Array rwork of size nwork
// Needed size 3*npoin, realloc if needed
// Field is indexed by tet2poi etc. Basically coord or met.
template<int ideg>
void setMeshBezier(Mesh &msh, int nwork, double *rwork){
	METRIS_THROW_MSG(TODOExcept(),"Replace this call with msh.setBezier");
	if(msh.ilag == 0) return;
  printf("-- Converting mesh to Bézier ideg = %d idim = %d .\n",ideg,msh.idim);

  if(msh.idim == 2){
		setFieldBezier<ideg,2>(msh,msh.coord,nwork,rwork);
  }else{
		setFieldBezier<ideg,3>(msh,msh.coord,nwork,rwork);
  }
	msh.ilag = 0;
}

template<int ideg>
void setMeshLagrange(Mesh &msh, int nwork, double *rwork){
	METRIS_THROW_MSG(TODOExcept(),"Replace this call with msh.setBezier");
	if(msh.ilag == 1) return;
  printf("-- Converting mesh to Lagrange ideg = %d.\n",ideg);
  if(msh.idim == 2){
		setFieldLagrange<ideg,2>(msh,msh.coord,nwork,rwork);
  }else{
		setFieldLagrange<ideg,3>(msh,msh.coord,nwork,rwork);
  }
	msh.ilag = 1;
}
#define BOOST_PP_LOCAL_MACRO(n)\
template void setMeshBezier<n>(Mesh &msh, int nwork, double *rwork);\
template void setMeshLagrange<n>(Mesh &msh, int nwork, double *rwork);\
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()
*/

/*
template<int ideg>
void setMetricBezier(Mesh &msh, int nwork, double *rwork){
	if(msh.ilagmet == 0) return;
  if(msh.idim == 2){
		setFieldBezier<ideg,3>(msh,msh.met,nwork,rwork);
  }else{
		setFieldBezier<ideg,6>(msh,msh.met,nwork,rwork);
  }
	msh.ilagmet = 0;
}

template<int ideg>
void setMetricLagrange(Mesh &msh, int nwork, double *rwork){
	if(msh.ilagmet  == 1) return;
  if(msh.idim == 2){
		setFieldLagrange<ideg,3>(msh,msh.met,nwork,rwork);
  }else{
		setFieldLagrange<ideg,6>(msh,msh.met,nwork,rwork);
  }
	msh.ilagmet = 1;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void setMeshBezier<n>(Mesh &msh, int nwork, double *rwork);\
template void setMeshLagrange<n>(Mesh &msh, int nwork, double *rwork);\
template void setMetricBezier<n>(Mesh &msh, int nwork, double *rwork);\
template void setMetricLagrange<n>(Mesh &msh, int nwork, double *rwork);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()
*/







// Avoid calling this routine directly. Call the specialized ones to ensure the 
// tags are properly handled. 
// Array rwork of size nwork
// Needed size szfld*npoin, realloc if needed
// Field is indexed by tet2poi etc. Basically coord or met.
template<int ideg, int szfld>
void setFieldBezier(MeshBase &msh, dblAr2 &rfld){

  if constexpr(ideg == 1) return;

  dblAr2 rwrk(msh.npoin,szfld);
  rwrk.set_n(msh.npoin);


  int nentt = msh.nelem > 0 ? msh.nelem :
              msh.nface > 0 ? msh.nface :
                              msh.nedge ;
  auto lag2bez = msh.nelem > 0 ? lag2bez3<ideg,szfld> :
                 msh.nface > 0 ? lag2bez2<ideg,szfld> :
                                 lag2bez1<ideg,szfld> ;
  intAr2 &ent2poi =  msh.nelem > 0 ? msh.tet2poi :
                     msh.nface > 0 ? msh.fac2poi :
                                     msh.edg2poi ;

	// Low dimension entities are always attached to high dimension ones if those exist
	// Therefore, we can simply work from the highest dim entities directly, and that will
	// cover all points in the mesh. 
  for(int ientt = 0; ientt < nentt; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;
    lag2bez(ent2poi[ientt],rfld,rwrk);
  }

	// Simply replace now
	for(int ipoin=0; ipoin<msh.npoin; ipoin++){
		if(msh.poi2ent(ipoin,0) < 0) continue;
		for(int j = 0; j < szfld; j++){
			rfld[ipoin][j] = rwrk[ipoin][j];
		}
	}

}



template<int ideg, int szfld>
void setFieldLagrange(MeshBase &msh, dblAr2 &rfld){

  if constexpr(ideg == 1) return;

  dblAr2 rwrk(msh.npoin,szfld);
  rwrk.set_n(msh.npoin);

  int nentt = msh.nelem > 0 ? msh.nelem :
              msh.nface > 0 ? msh.nface :
                              msh.nedge ;
  auto bez2lag = msh.nelem > 0 ? bez2lag3<ideg,szfld> :
                 msh.nface > 0 ? bez2lag2<ideg,szfld> :
                                 bez2lag1<ideg,szfld> ;
  intAr2 &ent2poi =  msh.nelem > 0 ? msh.tet2poi :
                     msh.nface > 0 ? msh.fac2poi :
                                     msh.edg2poi ;
                                     
  // Low dimension entities are always attached to high dimension ones if those exist
  // Therefore, we can simply work from the highest dim entities directly, and that will
  // cover all points in the mesh. 
  for(int ientt = 0; ientt < nentt; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;
    bez2lag(ent2poi[ientt],rfld,rwrk);
  }

  // Simply replace now
  for(int ipoin=0; ipoin<msh.npoin; ipoin++){
    // Detached point (neither corner nor attached to anything)
    if(msh.poi2ent(ipoin,0) < 0) continue;
    // Corners can be loose despite >= poi2ent and should also be left alone 
    // (they are always vertices anyways)
    int ibpoi = msh.poi2bpo[ipoin];
    if(ibpoi >= 0){
      if(msh.bpo2ibi(ibpoi,1) == 0) continue;
    }

    for(int j = 0; j < szfld; j++) rfld[ipoin][j] = rwrk[ipoin][j];
  }

}

// Metric fields
#define BOOST_PP_LOCAL_MACRO(n)\
template void setFieldBezier<n,2>(MeshBase &msh, dblAr2 &rfld);\
template void setFieldBezier<n,3>(MeshBase &msh, dblAr2 &rfld);\
template void setFieldBezier<n,6>(MeshBase &msh, dblAr2 &rfld);\
template void setFieldLagrange<n,2>(MeshBase &msh, dblAr2 &rfld);\
template void setFieldLagrange<n,3>(MeshBase &msh, dblAr2 &rfld);\
template void setFieldLagrange<n,6>(MeshBase &msh, dblAr2 &rfld);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

} // End namespace

