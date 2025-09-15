//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "io_libmeshb.hxx"

#include "aux_exceptions.hxx"                          // for METRIS_THROW_MSG
#include "ho_constants.hxx"                            // for edgnpps, facnpps
#include "low_eval.hxx"                                // for hana
#include "utils/mprintf.hxx"                                // for hana
#include "utils/CT_loop.hxx"
#include "metris_constants.hxx"                        // for METRIS_MAXTAGS
#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "cavity/msh_cavity.hxx"
#include "Boundary/msh_inisurf.hxx"
#include "msh_lag2bez.hxx"

#include <cstdlib>

namespace Metris{



std::string correctExtension_meshb(const std::string &s){
  std::string ret = s;
  if(ret.find(".mesh") == std::string::npos){
    //std::cout<<"## No filename extension provided, defaulting to .meshb\n";
    ret += ".meshb";
  }
  return ret;
}

std::string correctExtension_solb(const std::string &s){
  std::string ret = s;
  if(ret.find(".sol") == std::string::npos){
    //std::cout<<"## No filename extension provided, defaulting to .solb\n";
    ret += ".solb";
  }
  return ret;
}

// If no extension, default to .egads. Leaves .legads 
std::string correctExtension_egads(const std::string &s){
  std::string ret = s;
  if(ret.find(".legads") != std::string::npos) return ret; 
  if(ret.find(".egads") == std::string::npos){
    //std::cout<<"## No filename extension provided, defaulting to .egads\n";
    ret += ".egads";
  }
  return ret;
}

/*
Refer to docs/libMeshb7.pdf 
Block reads and writes are very convenient as they are the only routines implemented 
with arrays of arguments, as well as likely faster. 
For reference, file versions:
  | int size  | real size  | max file size 
1 |     32    |     32     |     2Go
2 |     32    |     64     |     2Go
3 |     32    |     64     |     8Exa Octet
4 |     64    |     64     |     8Exa Octet

If we are going the HPC route, we need to support 64bit integers. 
For now, let's throw an error for any version not 2 or 3.
*/


template<> int64_t MetrisOpenMeshFile<GmfWrite>(std::string meshName, int meshDim){
  int libVer = 3; // Set for when rwtyp == GmfWrite

  int64_t libIdx = GmfOpenMesh(meshName.c_str(), GmfWrite,  libVer, meshDim);

  METRIS_ENFORCE_MSG(libIdx && (libVer == 2 || libVer == 3),
    "FILE COULDNT BE OPENED OR WRONG VERSION: {}, name = {}", libVer, meshName);
  
  return libIdx;
}

template<> int64_t MetrisOpenMeshFile<GmfRead>(std::string meshName, int *meshDim){
  int libVer; // Set for when rwtyp == GmfWrite

  int64_t libIdx = GmfOpenMesh(meshName.c_str(), GmfRead, &libVer, meshDim);

  METRIS_ENFORCE_MSG(libIdx && (libVer == 2 || libVer == 3),
    "FILE COULDNT BE OPENED OR WRONG VERSION: {}, name = {}", libVer, meshName);

  METRIS_ENFORCE_MSG(*meshDim == 3 || *meshDim == 2, "Unsupported dimension {}", *meshDim);

  return libIdx;
}
/*
  Entities are 1-indexed in file, but 0-indexed in code. 
  Element refs are kept as-is, they relate to CAD entities. 
  Point refs relate to CAD nodes. We don't store those explicitly. 
  Thus, they become the first "ibpois". 
  We will need to find the correct ref by interrogating the CAD. 
*/




void writeMeshCavity(std::string meshName_, MeshBase &msh, const MshCavity& cav){
  GETVDEPTH(msh.param);

  std::string meshName = msh.param->outmPrefix + meshName_;

  MPRINTF("-- Write cavity filename = {}\n",meshName.c_str());

  int ierro = 0;
  if(msh.curdeg < 1){
    MPRINTF("## Invalid maximum degree {} in readMesh !",msh.curdeg);
    ierro = 1;
  }
  if(msh.curdeg > __MAX_LIBMESHB_DEG__){
    MPRINTF("## Maximum degree supported by libmeshb = {} passed {} \n",
             __MAX_LIBMESHB_DEG__,msh.curdeg);
    ierro = 1;
  }
  if(meshName.find(".mesh") == std::string::npos){
    //std::cout<<"## No filename extension provided, defaulting to .meshb\n";
    meshName += ".meshb";
  }
  // this ensures all error messages are printed
  METRIS_ENFORCE(ierro <= 0);



  FEBasis ibas0 = msh.getBasis();
  msh.setBasis(msh.param->outbasis);
  

  int64_t libIdx;
  libIdx = MetrisOpenMeshFile<GmfWrite>(meshName.c_str(), msh.idim);


  if(msh.getBasis() == FEBasis::Bezier) GmfSetKwd( libIdx, GmfBezierBasis, 1);

  const int ncedg = cav.lcedg.get_n();
  const int ncfac = cav.lcfac.get_n();
  const int nctet = cav.lctet.get_n();

  intAr2 ent2po2;
  intAr1 ent2re2;

  if(ncedg > 0){
  //   Note on edges: for some reason, libmeshb uses 1 index for edges
  // but the usual 3 for triangles, 4 for tets, etc.     
    constexpr int mppe = getnnod1(METRIS_MAX_DEG);

    int fKwd = libmeshb::edgeKwds[msh.curdeg];
    int npp  = getnnod1(msh.curdeg);
    if(msh.curdeg > 1){

      int myOrd[mppe]; // see remark above
      for(int i = 0; i < npp; i++){
        for(int j = 0; j < 1; j++){
          myOrd[i+j] = ordedg.s[msh.curdeg][i][j];
        }
      }

      GmfSetKwd(libIdx, libmeshb::edgeOrdKwds[msh.curdeg], npp); 
      GmfSetBlock(libIdx, libmeshb::edgeOrdKwds[msh.curdeg], 1, npp, 0, NULL, NULL,
        GmfIntVec, npp, &myOrd[0], &myOrd[(npp-1)]);
    }

    ent2po2.set_n(0); // No copy 
    ent2po2.allocate(ncedg,npp);
    ent2po2.set_n(ncedg);

    ent2re2.set_n(0); // No copy 
    ent2re2.allocate(ncedg);
    ent2re2.set_n(ncedg);

    //intAr2 edg2po2(ncedg,npp,&buffp[0]);
    //intAr1 edg2re2(ncedg,&buffr[0]);

    for(int ii = 0; ii < ncedg; ii++){
      int iedge = cav.lcedg[ii];
      for(int jj = 0; jj < npp; jj++){
        int ip = msh.edg2poi(iedge,jj);
        ent2po2[ii][jj] = ip + 1;
      }
      ent2re2[ii] = msh.edg2ref[iedge] + 1; //iedge + 1;
    }

    GmfSetKwd( libIdx, fKwd, ncedg);
    GmfSetBlock(libIdx, fKwd, 1, ncedg, 0, NULL, NULL,
      GmfIntVec, npp, &ent2po2[0][0], &ent2po2[ncedg-1][0],
      GmfInt   ,      &ent2re2[0   ], &ent2re2[ncedg-1  ]);

  }



  if(ncfac > 0){
    constexpr int mppf = getnnod2(METRIS_MAX_DEG);

    int fKwd = libmeshb::faceKwds[msh.curdeg];
    int npp = getnnod2(msh.curdeg);


    if(msh.curdeg > 1){
      int myOrd[3*mppf];
      for(int i = 0; i < npp; i++){
        for(int j = 0; j < 3; j++){
          myOrd[i*3+j] = ordfac.s[msh.curdeg][i][j];
        }
      }
      GmfSetKwd(libIdx, libmeshb::faceOrdKwds[msh.curdeg], npp); 
      GmfSetBlock(libIdx, libmeshb::faceOrdKwds[msh.curdeg], 1, npp, 0, NULL, NULL,
        GmfIntVec, npp, &myOrd[0], &myOrd[3*(npp-1)]);

    }

    ent2po2.set_n(0); // No copy 
    ent2po2.allocate(ncfac,npp);
    ent2po2.set_n(ncfac);

    ent2re2.set_n(0); // No copy 
    ent2re2.allocate(ncfac);
    ent2re2.set_n(ncfac);

    //intAr2 fac2po2(ncfac,npp,&buffp[0]);
    //intAr1 fac2re2(ncfac,&buffr[0]);

    //printf("Debug cav.lcfac = ");
    //cav.lcfac.print();
    //for(int iface : cav.lcfac){
    //  printf(" {} : ",iface);
    //  intAr1(npp, msh.fac2poi[iface]).print();
    //}

    for(int ii = 0; ii < ncfac; ii++){
      int iface = cav.lcfac[ii];
      for(int jj = 0; jj < npp; jj++){
        int ip = msh.fac2poi(iface,jj);
        ent2po2[ii][jj] = ip + 1;
      }
      ent2re2[ii] = msh.fac2ref[iface] + 1; //iface + 1;
    }

    GmfSetKwd( libIdx, fKwd, ncfac);
    GmfSetBlock(libIdx, fKwd, 1,  ncfac, 0, NULL, NULL,
      GmfIntVec, npp, &ent2po2[0][0], &ent2po2[ncfac-1][0],
      GmfInt   ,      &ent2re2[0   ], &ent2re2[ncfac-1  ]);

  }

  // Redundant but safer
  if(nctet > 0 && msh.idim >= 3){
    constexpr int mppt = getnnod3(METRIS_MAX_DEG);

    int eKwd = libmeshb::elemKwds[msh.curdeg];
    int npp = getnnod3(msh.curdeg);

    if(msh.curdeg > 1){
      int myOrd[4*mppt];
      for(int i = 0; i < npp; i++){
        for(int j = 0; j < 4; j++){
          myOrd[i*4+j] = ordtet.s[msh.curdeg][i][j];
        }
      }

      GmfSetKwd(libIdx, libmeshb::elemOrdKwds[msh.curdeg], npp); 
      GmfSetBlock(libIdx, libmeshb::elemOrdKwds[msh.curdeg], 1, npp, 0, NULL, NULL,
        GmfIntVec, npp, &myOrd[0], &myOrd[4*(npp-1)]);
    }

    //intAr2 tet2po2(nctet,npp,&buffp[0]);
    //intAr1 tet2re2(nctet,&buffr[0]);

    ent2po2.set_n(0); // No copy 
    ent2po2.allocate(nctet,npp);
    ent2po2.set_n(nctet);

    ent2re2.set_n(0); // No copy 
    ent2re2.allocate(nctet);
    ent2re2.set_n(nctet);

    for(int ii = 0; ii < nctet; ii++){
      int itete = cav.lctet[ii];
      for(int jj = 0; jj < npp; jj++){
        int ip = msh.tet2poi(itete,jj);
        ent2po2[ii][jj] = ip + 1;
      }
      ent2re2[ii] = msh.tet2ref[itete] + 1; //itete + 1;
    }


    GmfSetKwd( libIdx, eKwd, nctet);
    GmfSetBlock(libIdx, eKwd, 1, nctet, 0, NULL, NULL,
      GmfIntVec, npp, &ent2po2[0][0], &ent2po2[nctet-1][0],
             GmfInt , &ent2re2[0   ], &ent2re2[nctet-1  ]); //idx0 + nface[iDeg]-1
  }

  // Only vertices can be corners
  int mcorn = 2 * ncedg * msh.isboundary_edges();
  if(msh.isboundary_faces()) mcorn = MAX(mcorn, 3*ncfac);
  // But there can be no more than points, obviously. 

  mcorn++; // ipins

  intAr1 lcorn(mcorn);

  bool iipns = false;

  // Quadratic search in the corners so we don't have to use tags.
  for(int ii = 0; ii < ncedg; ii++){
    int iedge = cav.lcedg[ii];
    if(isdeadent(iedge,msh.edg2poi)) continue;
    for(int jj = 0; jj < 2; jj++){
      int ip = msh.edg2poi(iedge,jj);
      int ib = msh.poi2bpo[ip];
      METRIS_ASSERT(ib >= 0);

      if(msh.bpo2ibi(ib,1) != 0) continue;

      bool iskip = false;
      for(int icorn : lcorn){
        if(icorn != ip+1) continue;
        iskip = true;
        break;
      }

      if(iskip) continue;

      lcorn.stack(ip + 1);
      if(ip == cav.ipins) iipns = true;
    }
  }

  if(msh.isboundary_faces()){
    for(int ii = 0; ii < ncfac; ii++){
      int iface = cav.lcfac[ii];
      if(isdeadent(iface,msh.fac2poi)) continue;
      for(int jj = 0; jj < 3; jj++){
        int ip = msh.fac2poi(iface,jj);
        int ib = msh.poi2bpo[ip];
        METRIS_ASSERT(ib >= 0);

        if(msh.bpo2ibi(ib,1) != 0) continue;

        bool iskip = false;
        for(int icorn : lcorn){
          if(icorn != ip+1) continue;
          iskip = true;
          break;
        }

        if(iskip) continue;

        lcorn.stack(ip + 1);
        if(ip == cav.ipins) iipns = true;
      }
    }
  }

  if(!iipns){
    lcorn.stack(cav.ipins + 1);
  }


  int ncorn = lcorn.get_n();
  METRIS_ASSERT(ncorn >= 1);

  GmfSetKwd(libIdx, GmfCorners, ncorn);
  GmfSetBlock(libIdx, GmfCorners, 1, ncorn, 0, NULL, NULL, 
              GmfInt, &lcorn[0], &lcorn[ncorn-1]);

  GmfSetKwd(libIdx, GmfVertices, msh.npoin);
  if(msh.idim == 3){
    for(int ii = 0; ii < msh.npoin; ii++){
      GmfSetLin(libIdx, GmfVertices, msh.coord(ii,0), msh.coord(ii,1), 
                                     msh.coord(ii,2),0);
    }
  }else if(msh.idim == 2){
    for(int ii = 0; ii < msh.npoin; ii++){
      GmfSetLin(libIdx, GmfVertices, msh.coord(ii,0), msh.coord(ii,1), 0);
    }
  }

  GmfCloseMesh( libIdx );

  msh.setBasis(ibas0);
}


void debugInveval(std::string meshName_, MeshBase &msh, int tdim, int* ent2pol, double *coop){
  METRIS_ASSERT(tdim > 1);


  std::string meshName = msh.param->outmPrefix + meshName_;

  FEBasis ibas0 = msh.getBasis();
  msh.setBasis(FEBasis::Lagrange); 

  if(meshName.find(".mesh") == std::string::npos){
    //std::cout<<"## No filename extension provided, defaulting to .meshb\n";
    meshName += ".meshb";
  }

  if(msh.param->iverb >= 1) std::cout<<"-- Write file "<<meshName<<std::endl;



  int nnode = tdim == 1 ? getnnod1(msh.curdeg) :
              tdim == 2 ? getnnod2(msh.curdeg) : getnnod3(msh.curdeg);

  for(int ii = 0; ii < nnode; ii++) ent2pol[ii] += 1;

  int ipnew = msh.newpoitopo(PointType::Vertex,0,-1);
  msh.newbpotopo(Vertex{ipnew},0,-1);

  for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipnew,ii) = coop[ii];

  int64_t libIdx;
  libIdx = MetrisOpenMeshFile<GmfWrite>(meshName.c_str(), msh.idim);

  int eltKwd = tdim == 2 ? libmeshb::faceKwds[msh.curdeg] 
                         : libmeshb::elemKwds[msh.curdeg];


  if(msh.curdeg > 1){
    constexpr int mnode = getnnod3(METRIS_MAX_DEG); // largest possible for static alloc
    int myOrd[4*mnode];
    for(int ii = 0; ii < nnode; ii++){
      for(int jj = 0; jj < tdim+1; jj++){
        if(tdim == 2) myOrd[ii*(tdim+1)+jj] = ordfac.s[msh.curdeg][ii][jj];
        else          myOrd[ii*(tdim+1)+jj] = ordtet.s[msh.curdeg][ii][jj];
      }
    }
    int ordKwd = tdim == 2 ? libmeshb::faceOrdKwds[msh.curdeg] : libmeshb::elemOrdKwds[msh.curdeg];
    GmfSetKwd(libIdx, ordKwd, nnode); 
    GmfSetBlock(libIdx, ordKwd, 1, nnode, 0, NULL, NULL,
               GmfIntVec, nnode, &myOrd[0], &myOrd[(tdim+1)*(nnode-1)]);
  }

  GmfSetKwd( libIdx, eltKwd, 1);
  int iref = 1;
  GmfSetBlock(libIdx, eltKwd, 1, 1, 0, NULL, NULL,
    GmfIntVec, nnode, &ent2pol[0]      , &ent2pol[0],
    GmfInt   ,        &iref            , &iref     );




  GmfSetKwd(libIdx, GmfCorners, 1);
  int icorn = ipnew + 1;
  GmfSetBlock(libIdx, GmfCorners, 1, 1, 0, NULL, NULL, 
              GmfInt, &icorn, &icorn);



  GmfSetKwd(libIdx, GmfVertices, msh.npoin);
  if(msh.idim == 3){
    for(int ii = 0; ii < msh.npoin; ii++){
      GmfSetLin(libIdx, GmfVertices, msh.coord(ii,0), msh.coord(ii,1), msh.coord(ii,2),0);
    }
  }else if(msh.idim == 2){
    for(int ii = 0; ii < msh.npoin; ii++){
      GmfSetLin(libIdx, GmfVertices, msh.coord(ii,0), msh.coord(ii,1), 0);
    }
  }


  for(int ii = 0; ii < nnode; ii++) ent2pol[ii] -= 1;
  msh.setBasis(ibas0);
  msh.killpoint(ipnew);

  GmfCloseMesh( libIdx );

}


void writeMesh(std::string meshName, const MeshBase &msh, bool iprefix,
               int iedg0, int ifac0, int iele0,
               int iedg1, int ifac1, int iele1){
  GETVDEPTH(msh.param);

  if(iedg1 < 0) iedg1 = msh.nedge;
  if(ifac1 < 0) ifac1 = msh.nface;
  if(iele1 < 0) iele1 = msh.nelem;

  METRIS_ASSERT(iedg0 >= 0 && ifac0 >= 0 && iele0 >= 0);
  METRIS_ASSERT(iedg1 >= iedg0 && ifac1 >= ifac0 && iele1 >= iele0);

  std::string eff_meshName = iprefix ? msh.param->outmPrefix + meshName
                                     : meshName;

  MPRINTF("-- Write file {} npoin {} nedge {} nface {} tetra {} \n",
          eff_meshName.c_str(),msh.npoin,iedg1,ifac1,iele1);

  dblAr2 coord;
  if(msh.getBasis() == msh.param->outbasis || msh.curdeg == 1){
    coord = msh.coord;
  }else{
    // We need to convert the coordinates to the output basis.
    // This is a bit of a hack, but it works for now.
    coord.allocate(msh.npoin, msh.idim);
    coord.set_n(msh.npoin);
    for(int ii = 0; ii < msh.npoin; ii++){
      for(int jj = 0; jj < msh.idim; jj++){
        coord(ii,jj) = msh.coord(ii,jj);
      }
    }
    int tdim = msh.get_tdim();
    const intAr2 &ent2poi = msh.ent2poi(tdim);
    MSH_DIM_DEG0(msh){ // defines ideg and idim as constexpr
      if(msh.param->outbasis == FEBasis::Bezier){
        auto lag2bez = tdim == 1 ? lag2bez1<ideg,gdim> :
                       tdim == 2 ? lag2bez2<ideg,gdim> :
                                   lag2bez3<ideg,gdim> ;
        for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
          if(isdeadent(ientt,ent2poi)) continue;
          lag2bez(ent2poi[ientt], msh.coord, coord);
        }
      }else{
        auto bez2lag = tdim == 1 ? bez2lag1<ideg,gdim> :
                       tdim == 2 ? bez2lag2<ideg,gdim> :
                                   bez2lag3<ideg,gdim> ;
        for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
          if(isdeadent(ientt,ent2poi)) continue;
          bez2lag(ent2poi[ientt], msh.coord, coord);
        }
      }
    }MSH_DIM_DEG1();
  }


  int ierro = 0;
  if(msh.curdeg < 1){
    MPRINTF("## Invalid maximum degree in readMesh ! \n");
    ierro = 1;
  }
  if(msh.curdeg > __MAX_LIBMESHB_DEG__){
    MPRINTF("## Maximum degree supported by libmeshb = {} passed {}\n",
             __MAX_LIBMESHB_DEG__,msh.curdeg);
    ierro = 1;
  }
  if(eff_meshName.find(".mesh") == std::string::npos){
    //std::cout<<"## No filename extension provided, defaulting to .meshb\n";
    eff_meshName += ".meshb";
  }
  // this ensures all error messages are printed
  METRIS_ENFORCE(ierro <= 0);

  intAr1 edg2ref(iedg1-iedg0);
  intAr1 fac2ref(ifac1-ifac0);
  intAr1 tet2ref(iele1-iele0);
  intAr2 edg2poi(iedg1-iedg0, msh.edg2poi.get_stride());
  intAr2 fac2poi(ifac1-ifac0, msh.fac2poi.get_stride());
  intAr2 tet2poi(iele1-iele0, msh.tet2poi.get_stride());


  if(msh.idim >= 3){
    for(int ii = iele0; ii < iele1;ii++){
      tet2ref[ii-iele0] = msh.tet2ref[ii] + 1;
      for(int jj = 0; jj < getnnod3(msh.curdeg); jj++)
        tet2poi(ii-iele0,jj) = msh.tet2poi(ii,jj) + 1;
    }
  }
  for(int ii = ifac0; ii < ifac1;ii++){
    fac2ref[ii-ifac0] = msh.fac2ref[ii] + 1;
    for(int jj = 0; jj < getnnod2(msh.curdeg); jj++)
      fac2poi(ii-ifac0,jj) = msh.fac2poi(ii,jj) + 1;
  }
  for(int ii = iedg0; ii < iedg1; ii++){
    edg2ref[ii-iedg0] = msh.edg2ref[ii] + 1;
    for(int jj = 0; jj < getnnod1(msh.curdeg); jj++)
      edg2poi(ii-iedg0,jj) = msh.edg2poi(ii,jj) + 1;
  }

  int64_t libIdx;
  libIdx = MetrisOpenMeshFile<GmfWrite>(eff_meshName.c_str(), msh.idim);


  if(msh.getBasis() == FEBasis::Bezier) GmfSetKwd( libIdx, GmfBezierBasis, 1);

  if(iedg1 - iedg0 > 0){
  //   Note on edges: for some reason, libmeshb uses 1 index for edges
  // but the usual 3 for triangles, 4 for tets, etc.     
    //CPRINTF2(" - START writing edges: {} -> {} \n",iedg0, iedg1);
    constexpr int mppe = getnnod1(METRIS_MAX_DEG);

    int fKwd = libmeshb::edgeKwds[msh.curdeg];
    int npp  = getnnod1(msh.curdeg);
    if(msh.curdeg > 1){

      int myOrd[mppe]; // see remark above
      for(int i = 0; i < npp; i++){
        for(int j = 0; j < 1; j++){
          myOrd[i+j] = ordedg.s[msh.curdeg][i][j];
        }
      }

      GmfSetKwd(libIdx, libmeshb::edgeOrdKwds[msh.curdeg], npp); 
      GmfSetBlock(libIdx, libmeshb::edgeOrdKwds[msh.curdeg], 1, npp, 0, NULL, NULL,
        GmfIntVec, npp, &myOrd[0], &myOrd[(npp-1)]);
    }
    GmfSetKwd( libIdx, fKwd, iedg1- iedg0);
    GmfSetBlock(libIdx, fKwd, 1, iedg1- iedg0, 0, NULL, NULL,
      GmfIntVec, npp,  edg2poi[0],  edg2poi[iedg1-iedg0-1],
      GmfInt   ,      &edg2ref[0], &edg2ref[iedg1-iedg0-1]);

    //CPRINTF2(" - DONE writing edges\n");
  }



  if(ifac1 - ifac0 > 0){
    //CPRINTF2(" - START writing triangles: {} \n",ifac1);
    constexpr int mppf = getnnod2(METRIS_MAX_DEG);

    int fKwd = libmeshb::faceKwds[msh.curdeg];
    int nppf = getnnod2(msh.curdeg);


    if(msh.curdeg > 1){
      int myOrd[3*mppf];
      for(int i = 0; i < nppf; i++){
        for(int j = 0; j < 3; j++){
          myOrd[i*3+j] = ordfac.s[msh.curdeg][i][j];
        }
      }
      GmfSetKwd(libIdx, libmeshb::faceOrdKwds[msh.curdeg], nppf); 
      GmfSetBlock(libIdx, libmeshb::faceOrdKwds[msh.curdeg], 1, nppf, 0, NULL, NULL,
        GmfIntVec, nppf, &myOrd[0], &myOrd[3*(nppf-1)]);

    }

    GmfSetKwd( libIdx, fKwd, ifac1 - ifac0);
    GmfSetBlock(libIdx, fKwd, 1,  ifac1 - ifac0, 0, NULL, NULL,
      GmfIntVec, nppf,  fac2poi[0], fac2poi[ifac1-ifac0-1],
      GmfInt   ,       &fac2ref[0],&fac2ref[ifac1-ifac0-1]);

    //CPRINTF2(" - DONE writing triangles\n");
  }

  // Redundant but safer
  if(iele1 - iele0 > 0 && msh.idim >= 3){
    //CPRINTF2(" - START writing tetrahedra: {} \n",iele1);
    constexpr int mppt = getnnod3(METRIS_MAX_DEG);

    int eKwd = libmeshb::elemKwds[msh.curdeg];
    int nppt = getnnod3(msh.curdeg);

    if(msh.curdeg > 1){
      int myOrd[4*mppt];
      for(int i = 0; i < nppt; i++){
        for(int j = 0; j < 4; j++){
          myOrd[i*4+j] = ordtet.s[msh.curdeg][i][j];
        }
      }

      GmfSetKwd(libIdx, libmeshb::elemOrdKwds[msh.curdeg], nppt); 
      GmfSetBlock(libIdx, libmeshb::elemOrdKwds[msh.curdeg], 1, nppt, 0, NULL, NULL,
        GmfIntVec, nppt, &myOrd[0], &myOrd[4*(nppt-1)]);
    }
    GmfSetKwd( libIdx, eKwd, iele1 - iele0);
    GmfSetBlock(libIdx, eKwd, 1, iele1 - iele0, 0, NULL, NULL,
                GmfIntVec, nppt, tet2poi[0], tet2poi[iele1-iele0-1],
                GmfInt   ,      &tet2ref[0],&tet2ref[iele1-iele0-1]); //idx0 + nface[iDeg]-1
    //CPRINTF2(" - DONE writing tetrahedra\n");
  }





  //CPRINTF2(" - START writing point: {} \n",msh.npoin);
  //  printf("Debug now bpo2ibi:\n");
  //  msh.bpo2ibi.print(msh.nbpoi);
  /*
  Build tables for VerticesOnGeometricEntities
  */

  //if(msh.CAD()){
    try{

      // stream CAD into file 

      //size_t nbyte;
      //char* stream;
      //ierro = EG_exportModel(msh.CAD.EGADS_model, &nbyte, &stream);
      //if(ierro != 0){
      //  print_EGADS_error("EG_exportModel",ierro);
      //  METRIS_THROW_MSG("Failed to export model to stream.");
      //}
      //if(iverb > 0) printf(" - Stream of size {}b \n",nbyte);
      //GmfWriteByteFlow(libIdx, stream, (int) nbyte);

      //CADInfo CAD2;
      //CAD2.setModel(nbyte, stream);
      //printf("here\n");
      //wait();



      //int* dbgptr = (int*)stream;
      //for(int ii = 0; ii < MIN(10,nbyte / sizeof(int)); ii++){
      //  printf(" DEBUG STREAM AS INT: {} = {} \n",ii,*dbgptr);
      //  dbgptr++;
      //}

      //ego geom;
      //int oclass,mtype,nbody,*dum;
      //ego *bodies;
      //ierro = EG_getTopology(msh.CAD.EGADS_model,&geom,&oclass,&mtype,NULL,&nbody,&bodies,&dum);
      //printf("## DEBUG WRITING CHECK NBODY = {} \n",nbody);

      //int *nbodyptr = (int*) (&stream[0] + 3*sizeof(int) + 6*sizeof(double));
      //printf("## DEBUG NBODY FROM STREAM {} \n",*nbodyptr);


      int mcorn = getNumCorners(msh);
      int mgpoe = msh.nbpoi;
      int mgpof = msh.nbpoi;
      intAr1 lcorn(mcorn);
      intAr1 lpoic(msh.npoin);
      intAr2 lgpoe(mgpoe,2);
      intAr2 lgpof(mgpof,2);
      dblAr2 rgpoe(mgpoe,2);
      dblAr2 rgpof(mgpof,3);
      genOnGeometricEntLists(msh, lcorn, lpoic,
                                  lgpoe, rgpoe,
                                  lgpof, rgpof, 1);
      int ncorn = lcorn.get_n();
      int ngpoe = lgpoe.get_n();
      int ngpof = lgpof.get_n();


      GmfSetKwd(libIdx, GmfVertices, msh.npoin);
      GmfSetBlock(libIdx, GmfVertices, 1, msh.npoin, 0, NULL, NULL,
        GmfDoubleVec, msh.idim, &coord(0,0), &coord[msh.npoin-1][0],
        GmfInt                , &lpoic[0]  , &lpoic[msh.npoin-1]);

      if(ngpof > 0){
        GmfSetKwd(libIdx, GmfVerticesOnGeometricTriangles, ngpof);
        GmfSetBlock(libIdx, GmfVerticesOnGeometricTriangles, 1, ngpof, 0, NULL, NULL, 
                    GmfIntVec   , 2, &lgpof[0][0], &lgpof[ngpof-1][0],
                    GmfDoubleVec, 3, &rgpof[0][0], &rgpof[ngpof-1][0]);
      }
//      for(int ii = 0; ii < ngpof; ii++){
//        GmfSetLin(libIdx, GmfVerticesOnGeometricTriangles, 
//          lgpof[ii][0],lgpof[ii][1],
//          rgpof[ii][0],rgpof[ii][1],rgpof[ii][2]);
//      }


      if(ngpoe > 0){
        GmfSetKwd(libIdx, GmfVerticesOnGeometricEdges, ngpoe);
        GmfSetBlock(libIdx, GmfVerticesOnGeometricEdges, 1, ngpoe, 0, NULL, NULL, 
                    GmfIntVec   , 2, &lgpoe[0][0], &lgpoe[ngpoe-1][0],
                    GmfDoubleVec, 2, &rgpoe[0][0], &rgpoe[ngpoe-1][0]);
      }

      //GmfSetKwd(libIdx, GmfVerticesOnGeometricVertices, ncorn);
      //GmfSetBlock(libIdx, GmfVerticesOnGeometricVertices, 1, ncorn, 0, NULL, NULL 
      //           ,GmfIntVec   , 2, &lcorn[0],  &lcorn[2*ncorn-1]);
      if(ncorn > 0){
        GmfSetKwd(libIdx, GmfCorners, ncorn);
        GmfSetBlock(libIdx, GmfCorners, 1, ncorn, 0, NULL, NULL, 
                    GmfInt, &lcorn[0], &lcorn[ncorn-1]);
      }
    }catch(const MetrisExcept &e){
      std::cout<<"## Type: "<<e.what()<<std::endl;

      #ifndef NO_BOOST_EXCEPT
        if(std::string const * ms=boost::get_error_info<excMessage>(e) )
          std::cout<<"## Message: "<<*ms; 
        if(boost::stacktrace::stacktrace const * tr=boost::get_error_info<excStackTrace>(e) )
          std::cerr << "## Call stack: \n" << *tr;
      #endif
      METRIS_THROW_MSG("TODO: Manage corners (lpoic) in this context");
    }

  #if 0
  }else{
    intAr1 lpoic(msh.npoin);
    intAr1 lcorn(10);
    lcorn.set_n(0);
    for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
      int ipoin = msh.bpo2ibi(ibpoi,0);
      if(ipoin < 0) continue;
      if(msh.bpo2ibi(ibpoi,1) != 0) continue;
      lcorn.stack(ipoin+1);
      lpoic[ipoin] = lcorn.get_n();
    }
    int ncorn = lcorn.get_n();
    CPRINTF1(" - counted {} corners in mesh\n",ncorn);

    GmfSetKwd(libIdx, GmfVertices, msh.npoin);
    GmfSetBlock(libIdx, GmfVertices, 1, msh.npoin, 0, NULL, NULL,
      GmfDoubleVec, msh.idim, &coord(0,0), &coord[msh.npoin-1][0],
      GmfInt                , &lpoic[0]      , &lpoic[msh.npoin-1]);

    if(ncorn > 0){
      GmfSetKwd(libIdx, GmfCorners, ncorn);
      GmfSetBlock(libIdx, GmfCorners, 1, ncorn, 0, NULL, NULL, 
                  GmfInt, &lcorn[0], &lcorn[ncorn-1]);
    }
    //GmfSetKwd(libIdx, GmfVertices, msh.npoin);
    //if(msh.idim == 3){
    //  for(int ii = 0; ii < msh.npoin; ii++){
    //    GmfSetLin(libIdx, GmfVertices, msh.coord(ii,0), msh.coord(ii,1), msh.coord(ii,2),0);
    //  }
    //}else if(msh.idim == 2){
    //  for(int ii = 0; ii < msh.npoin; ii++){
    //    GmfSetLin(libIdx, GmfVertices, msh.coord(ii,0), msh.coord(ii,1), 0);
    //  }
    //}
  }
  #endif
 
 


  GmfCloseMesh( libIdx );

}




void writeField(std::string outname, const MeshBase &msh, SolTyp stype, dblAr1 &rfld, int ndim){
  GETVDEPTH(msh.param);

  if(stype != SolTyp::P0Elt && stype != SolTyp::CG) 
    METRIS_THROW_MSG("TODO: Implement other SolTyps in writeField");

  int tdimn = msh.get_tdim();
  METRIS_ENFORCE(tdimn == 1 || tdimn == 2 || tdimn == 3);
 


  std::string metName = correctExtension_solb(outname);
  metName = msh.param->outmPrefix + metName;

  int szfld;
  if(ndim == 1){
    szfld = GmfSca;
  }else{
    szfld = GmfVec;
  }

  int ndof = -1;
  if(stype == SolTyp::CG)         ndof = msh.npoin;
  else if(stype == SolTyp::P0Elt) ndof = msh.nentt(msh.get_tdim());
  else METRIS_THROW_MSG("TODO: Sol typ not implemented ndof")
  if(rfld.get_n() / ndim != ndof){
    CPRINTF1("## WARNING rfld.get_n / ndim = {} ndof {} \n",
             rfld.get_n() / ndim, ndof);
    return;
  }



  int64_t libIdx;

  libIdx = MetrisOpenMeshFile<GmfWrite>(metName, tdimn);


  MPRINTF("-- Write file {}\n",metName.c_str());

  int solKwd = 0;
  bool iHO = false;
  if(stype == SolTyp::P0Elt){
    int solPosKwd;
    if(msh.curdeg == 1){
      solKwd = tdimn == 1 ? GmfSolAtEdges :
               tdimn == 2 ? GmfSolAtTriangles : GmfSolAtTetrahedra;
    }else if(msh.curdeg == 2){
      iHO = true;
      solKwd = tdimn == 1 ? GmfHOSolAtEdgesP2 :
               tdimn == 2 ? GmfHOSolAtTrianglesP2 : GmfHOSolAtTetrahedraP2;
      solPosKwd = tdimn == 1 ? GmfHOSolAtEdgesP2NodesPositions :
                  tdimn == 2 ? GmfHOSolAtTrianglesP2NodesPositions : 
                               GmfHOSolAtTetrahedraP2NodesPositions;
      GmfSetKwd(libIdx, solPosKwd, 1);
      if(tdimn == 1){
        GmfSetLin(libIdx, solPosKwd, 0.5, 0.5);
      }else if(tdimn == 2){
        GmfSetLin(libIdx, solPosKwd, 0.5, 0.5, 0.5);
      }else{
        GmfSetLin(libIdx, solPosKwd, 0.5, 0.5, 0.5, 0.5);
      }
    }else{
      METRIS_THROW_MSG("TODO: SoltAtElt kwd deg > 2 not implemented");
    }
  }else if(stype == SolTyp::CG){
    solKwd = GmfSolAtVertices;
  }

  if(!iHO) GmfSetKwd(libIdx, solKwd, ndof, 1, &szfld);
  else     GmfSetKwd(libIdx, solKwd, ndof, 1, &szfld, 0, 1);

  GmfSetBlock(libIdx, solKwd, 1, ndof, 0, NULL, NULL,
    GmfDoubleVec, ndim, &rfld[0], &rfld[ndim*(ndof-1)]);

  //for(int ipoin = 0; ipoin < nfld; ipoin++){
  //  if(ndim == 1){
  //    GmfSetLin(libIdx, solKwd, rfld[ipoin]);
  //  }else if(ndim == 2){
  //    GmfSetLin(libIdx, solKwd, &rfld[2*ipoin]); // , rfld[2*ipoin+1]
  //  }else{
  //    GmfSetLin(libIdx, solKwd, rfld[3*ipoin+0], rfld[3*ipoin+1], rfld[3*ipoin+2]);
  //  }
  //}
  //GmfSetBlock(libIdx, solKwd, 1, nfld, 0, NULL, NULL,
  //            GmfDoubleVec, ndim, &rfld[0], &rfld[ndim*(nfld-1)]);

  CPRINTF2("-- Done  writing field\n");


  GmfCloseMesh( libIdx );

}

static int warning_print_writeBackLinks = 0;

template<class MFT>
void writeBackLinks(std::string solName, Mesh<MFT>& msh){
  GETVDEPTH(msh.param);

  if(msh.met.metricClass() != MetricClass::MetricFieldFE){
    if(warning_print_writeBackLinks++ < 2){
      CPRINTF1(" # writeBackLinks disabled in analytical metric mode\n");
      if(warning_print_writeBackLinks >= 2) CPRINTF1(" # suppressing print\n");
    }
    return;
  }

  dblAr2 field(msh.npoin,msh.idim);

  double bary[4], coom[3];

  for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
    if(msh.isdeadpoint(ipoin)){
      for(int ii = 0; ii < msh.idim; ii++) field(ipoin, ii) = 0.0;
      continue;
    } 

    // Get centroid of back element 
    int pdim = msh.getpoitdim(ipoin);
    if(pdim == 0) continue;

    //for(int tdim = 1; tdim <= msh.get_tdim(); tdim++){
    int ientt = msh.poi2bak[ipoin];
    if(ientt < 0) continue;

    const intAr2& ent2poi = msh.bak->ent2poi(pdim);

    //int nentb = msh.bak->nentt(tdim);
    //if(ientt >= nentb){
    //  printf("## IENTT >= NENTB tdim = {} \n",tdim);
    //  for(int ii = 1; ii <= 3; ii++){
    //    printf("{} : {} {} \n", ii, msh.poi2bak(ientt,))
    //  }
    //  wait();
    //}

    for(int ii = 0; ii < pdim + 1; ii++) bary[ii] = 1.0 / (pdim + 1);

    if(msh.idim == 2){
      if(pdim == 1){
        eval1<2,1>(msh.bak->coord, ent2poi[ientt], msh.bak->getBasis(), 
          DifVar::None, DifVar::None, bary, coom, NULL, NULL);
      }else if(pdim == 2){
        eval2<2,1>(msh.bak->coord, ent2poi[ientt], msh.bak->getBasis(), 
          DifVar::None, DifVar::None, bary, coom, NULL, NULL);
      }else{
        METRIS_THROW_MSG("pdim 3 but gdim 2");
      }
    }else{
      if(pdim == 1){
        eval1<3,1>(msh.bak->coord, ent2poi[ientt], msh.bak->getBasis(), 
          DifVar::None, DifVar::None, bary, coom, NULL, NULL);
      }else if(pdim == 2){
        eval2<3,1>(msh.bak->coord, ent2poi[ientt], msh.bak->getBasis(), 
          DifVar::None, DifVar::None, bary, coom, NULL, NULL);
      }else{
        eval3<3,1>(msh.bak->coord, ent2poi[ientt], msh.bak->getBasis(), 
          DifVar::None, DifVar::None, bary, coom, NULL, NULL);
      }
    }

    for(int ii = 0; ii < msh.idim; ii++) 
      field(ipoin,ii) = coom[ii] - msh.coord(ipoin,ii);

    break;  
    //}

  }

  dblAr1 dum(field); // flatten to 1d array
  writeField(solName, msh, SolTyp::CG, dum, msh.idim);
}
template void writeBackLinks(std::string solNmae, Mesh<MetricFieldFE>& msh);
template void writeBackLinks(std::string solNmae, Mesh<MetricFieldAnalytical>& msh);


void writeEdgesLengths(const MeshBase &msh, std::string outnroot, 
                       intAr2 &edg2poi, const dblAr1 &rlned){
  std::string mshName = correctExtension_meshb(outnroot);
  std::string solName = correctExtension_solb(outnroot);

  const int idim = msh.idim;


  int nedge = edg2poi.get_n();
  METRIS_ENFORCE(nedge == rlned.get_n());

  for(int ii = 0; ii < nedge; ii++){
    for(int jj = 0; jj < 2; jj++){
      edg2poi(ii,jj) += 1;
    }
  }


  intAr1 edg2ref(nedge);
  edg2ref.fill(1);

  std::cout<<"-- Write file "<<mshName<<std::endl;

  int edgKwd = libmeshb::edgeKwds[1]; 
  int64_t libIdx = MetrisOpenMeshFile<GmfWrite>(mshName, idim);
  GmfSetKwd( libIdx, edgKwd, nedge);
  GmfSetBlock(libIdx, edgKwd, 1, nedge, 0, NULL, NULL,
    GmfIntVec, 2, &edg2poi(0,0), &edg2poi[nedge-1][0],
    GmfInt   ,    &edg2ref[0   ], &edg2ref[nedge-1]);

  GmfSetKwd(libIdx, GmfVertices, msh.npoin);
  if(msh.idim == 3){
    for(int ii = 0; ii < msh.npoin; ii++){
      GmfSetLin(libIdx, GmfVertices, msh.coord(ii,0), msh.coord(ii,1), 
                                     msh.coord(ii,2),0);
    }
  }else if(msh.idim == 2){
    for(int ii = 0; ii < msh.npoin; ii++){
      GmfSetLin(libIdx, GmfVertices, msh.coord(ii,0), msh.coord(ii,1), 0);
    }
  }

  GmfCloseMesh( libIdx );


  std::cout<<"-- Write file "<<solName<<std::endl;

  libIdx = MetrisOpenMeshFile<GmfWrite>(solName, idim);
  int szfld = GmfSca;
  GmfSetKwd(libIdx, GmfSolAtEdges, nedge, 1, &szfld);
  GmfSetBlock(libIdx, GmfSolAtEdges, 1, nedge, 0, NULL, NULL,
              GmfDouble, &rlned[0], &rlned[nedge-1]);
  GmfCloseMesh( libIdx );

  for(int ii = 0; ii < nedge; ii++){
    for(int jj = 0; jj < 2; jj++){
      edg2poi(ii,jj) -= 1;
    }
  }

}




#if 0
void writeMeshVecs(std::string meshName, MeshBase &msh, const dblAr2 &poi2vec){
  int iverb = msh.param->iverb;

  std::string eff_meshName = msh.param->outmPrefix + meshName;

  std::cout<<"-- Write file "<<eff_meshName<<std::endl;

  int ierro = 0;
  if(msh.curdeg < 1){
    std::cout<<"Invalid maximum degree in readMesh !"<<std::endl;
    ierro = 1;
  }
  if(msh.curdeg > __MAX_LIBMESHB_DEG__){
    std::cout<<"Maximum degree supported by libmeshb = "
    <<__MAX_LIBMESHB_DEG__<<" passed "<<msh.curdeg<<std::endl;
    ierro = 1;
  }
  if(eff_meshName.find(".mesh") == std::string::npos){
    //std::cout<<"## No filename extension provided, defaulting to .meshb\n";
    eff_meshName += ".meshb";
  }
  // this ensures all error messages are printed
  METRIS_ENFORCE(ierro <= 0);

  if(msh.idim >= 3){
    for(int i =0; i < msh.nelem;i++){
      msh.tet2ref[i] += 1;
      for(int j =0;j<getnnod3(msh.curdeg);j++)
        msh.tet2poi(i,j) += 1;
    }
  }
  for(int i = 0; i < msh.nface;i++){
    msh.fac2ref[i] += 1;
    for(int j = 0; j < getnnod2(msh.curdeg); j++)
      msh.fac2poi(i,j) += 1;
  }
  for(int i = 0; i < msh.nedge; i++){
    msh.edg2ref[i] += 1;
    for(int j = 0; j < getnnod1(msh.curdeg); j++)
      msh.edg2poi(i,j) += 1;
  }

  int64_t libIdx;
  libIdx = MetrisOpenMeshFile<GmfWrite>(eff_meshName.c_str(), msh.idim);


  if(msh.getBasis() == FEBasis::Bezier) GmfSetKwd( libIdx, GmfBezierBasis, 1);

  if(msh.nedge > 0){
  //   Note on edges: for some reason, libmeshb uses 1 index for edges
  // but the usual 3 for triangles, 4 for tets, etc.     
    if(iverb>0)std::cout<<"-- Start writing edges"<<std::endl;
    constexpr int mppe = getnnod1(METRIS_MAX_DEG);

    int fKwd = libmeshb::edgeKwds[msh.curdeg];
    int npp  = getnnod1(msh.curdeg);
    if(msh.curdeg > 1){

      int myOrd[mppe]; // see remark above
      for(int i = 0; i < npp; i++){
        for(int j = 0; j < 1; j++){
          myOrd[i+j] = ordedg.s[msh.curdeg][i][j];
        }
      }

      GmfSetKwd(libIdx, libmeshb::edgeOrdKwds[msh.curdeg], npp); 
      GmfSetBlock(libIdx, libmeshb::edgeOrdKwds[msh.curdeg], 1, npp, 0, NULL, NULL,
        GmfIntVec, npp, &myOrd[0], &myOrd[(npp-1)]);
    }

    GmfSetKwd( libIdx, fKwd, msh.nedge);
    GmfSetBlock(libIdx, fKwd, 1, msh.nedge, 0, NULL, NULL,
      GmfIntVec, npp, &msh.edg2poi(0,0), &msh.edg2poi[msh.nedge-1][0],
      GmfInt   ,      &msh.edg2ref[0   ], &msh.edg2ref[msh.nedge-1  ]);

    if(iverb>0)std::cout<<"-- Done  writing edges; nTot = "<<msh.nedge<<std::endl;
  }



  if(msh.nface > 0){
    if(iverb>0)std::cout<<"-- Start writing triangles"<<std::endl;
    constexpr int mppf = getnnod2(METRIS_MAX_DEG);

    int fKwd = libmeshb::faceKwds[msh.curdeg];
    int nppf = getnnod2(msh.curdeg);


    if(msh.curdeg > 1){
      int myOrd[3*mppf];
      for(int i = 0; i < nppf; i++){
        for(int j = 0; j < 3; j++){
          myOrd[i*3+j] = ordfac.s[msh.curdeg][i][j];
        }
      }
      GmfSetKwd(libIdx, libmeshb::faceOrdKwds[msh.curdeg], nppf); 
      GmfSetBlock(libIdx, libmeshb::faceOrdKwds[msh.curdeg], 1, nppf, 0, NULL, NULL,
        GmfIntVec, nppf, &myOrd[0], &myOrd[3*(nppf-1)]);

    }

    GmfSetKwd( libIdx, fKwd, msh.nface);
    GmfSetBlock(libIdx, fKwd, 1,  msh.nface, 0, NULL, NULL,
      GmfIntVec, nppf, &msh.fac2poi(0,0), &msh.fac2poi[msh.nface-1][0],
      GmfInt   ,       &msh.fac2ref[0   ], &msh.fac2ref[msh.nface-1  ]);

    if(iverb>0)std::cout<<"-- Done  writing triangles; nTot = "<<msh.nface<<std::endl;
  }

  // Redundant but safer
  if(msh.nelem > 0 && msh.idim >= 3){
    if(iverb>0)std::cout<<"-- Start writing tetrahedra"<<std::endl;
    constexpr int mppt = getnnod3(METRIS_MAX_DEG);

    int eKwd = libmeshb::elemKwds[msh.curdeg];
    int nppt = getnnod3(msh.curdeg);

    if(msh.curdeg > 1){
      int myOrd[4*mppt];
      for(int i = 0; i < nppt; i++){
        for(int j = 0; j < 4; j++){
          myOrd[i*4+j] = ordtet.s[msh.curdeg][i][j];
        }
      }

      GmfSetKwd(libIdx, libmeshb::elemOrdKwds[msh.curdeg], nppt); 
      GmfSetBlock(libIdx, libmeshb::elemOrdKwds[msh.curdeg], 1, nppt, 0, NULL, NULL,
        GmfIntVec, nppt, &myOrd[0], &myOrd[4*(nppt-1)]);
    }
    GmfSetKwd( libIdx, eKwd, msh.nelem);
    GmfSetBlock(libIdx, eKwd, 1, msh.nelem, 0, NULL, NULL,
      GmfIntVec, nppt, &msh.tet2poi(0,0), &msh.tet2poi[msh.nelem-1][0],
                GmfInt         , &msh.tet2ref[0   ], &msh.tet2ref[msh.nelem-1  ]); //idx0 + nface[iDeg]-1
    if(iverb>0)std::cout<<"-- Done  writing tetrahedra; nTot = "<<msh.nelem<<std::endl;
  }





  if(iverb>0)std::cout<<"-- Start writing points "<<msh.npoin<<std::endl;
  //  printf("Debug now bpo2ibi:\n");
  //  msh.bpo2ibi.print(msh.nbpoi);
  /*
  Build tables for VerticesOnGeometricEntities
  */

  if(msh.CAD()){
    try{

      // stream CAD into file 

      //size_t nbyte;
      //char* stream;
      //ierro = EG_exportModel(msh.CAD.EGADS_model, &nbyte, &stream);
      //if(ierro != 0){
      //  print_EGADS_error("EG_exportModel",ierro);
      //  METRIS_THROW_MSG("Failed to export model to stream.");
      //}
      //if(iverb > 0) printf(" - Stream of size {}b \n",nbyte);
      //GmfWriteByteFlow(libIdx, stream, (int) nbyte);

      //CADInfo CAD2;
      //CAD2.setModel(nbyte, stream);
      //printf("here\n");
      //wait();



      //int* dbgptr = (int*)stream;
      //for(int ii = 0; ii < MIN(10,nbyte / sizeof(int)); ii++){
      //  printf(" DEBUG STREAM AS INT: {} = {} \n",ii,*dbgptr);
      //  dbgptr++;
      //}

      //ego geom;
      //int oclass,mtype,nbody,*dum;
      //ego *bodies;
      //ierro = EG_getTopology(msh.CAD.EGADS_model,&geom,&oclass,&mtype,NULL,&nbody,&bodies,&dum);
      //printf("## DEBUG WRITING CHECK NBODY = {} \n",nbody);

      //int *nbodyptr = (int*) (&stream[0] + 3*sizeof(int) + 6*sizeof(double));
      //printf("## DEBUG NBODY FROM STREAM {} \n",*nbodyptr);


      int mcorn = getNumCorners(msh);
      int mgpoe = msh.nbpoi;
      int mgpof = msh.nbpoi;
      intAr1 lcorn(mcorn);
      intAr1 lpoic(msh.npoin);
      intAr2 lgpoe(mgpoe,2);
      intAr2 lgpof(mgpof,2);
      dblAr2 rgpoe(mgpoe,2);
      dblAr2 rgpof(mgpof,3);
      genOnGeometricEntLists(msh, lcorn, lpoic,
                                  lgpoe, rgpoe,
                                  lgpof, rgpof, 1);
      int ncorn = lcorn.get_n();
      int ngpoe = lgpoe.get_n();
      int ngpof = lgpof.get_n();


      GmfSetKwd(libIdx, GmfVertices, msh.npoin);
      GmfSetBlock(libIdx, GmfVertices, 1, msh.npoin, 0, NULL, NULL,
        GmfDoubleVec, msh.idim, &msh.coord(0,0), &msh.coord[msh.npoin-1][0],
        GmfInt                , &lpoic[0]       , &lpoic[msh.npoin-1]);

      if(ngpof > 0){
        GmfSetKwd(libIdx, GmfVerticesOnGeometricTriangles, ngpof);
        GmfSetBlock(libIdx, GmfVerticesOnGeometricTriangles, 1, ngpof, 0, NULL, NULL, 
                    GmfIntVec   , 2, &lgpof[0][0], &lgpof[ngpof-1][0],
                    GmfDoubleVec, 3, &rgpof[0][0], &rgpof[ngpof-1][0]);
      }
//      for(int ii = 0; ii < ngpof; ii++){
//        GmfSetLin(libIdx, GmfVerticesOnGeometricTriangles, 
//          lgpof[ii][0],lgpof[ii][1],
//          rgpof[ii][0],rgpof[ii][1],rgpof[ii][2]);
//      }


      if(ngpoe > 0){
        GmfSetKwd(libIdx, GmfVerticesOnGeometricEdges, ngpoe);
        GmfSetBlock(libIdx, GmfVerticesOnGeometricEdges, 1, ngpoe, 0, NULL, NULL, 
                    GmfIntVec   , 2, &lgpoe[0][0], &lgpoe[ngpoe-1][0],
                    GmfDoubleVec, 2, &rgpoe[0][0], &rgpoe[ngpoe-1][0]);
      }

      //GmfSetKwd(libIdx, GmfVerticesOnGeometricVertices, ncorn);
      //GmfSetBlock(libIdx, GmfVerticesOnGeometricVertices, 1, ncorn, 0, NULL, NULL 
      //           ,GmfIntVec   , 2, &lcorn[0],  &lcorn[2*ncorn-1]);
      if(ncorn > 0){
        GmfSetKwd(libIdx, GmfCorners, ncorn);
        GmfSetBlock(libIdx, GmfCorners, 1, ncorn, 0, NULL, NULL, 
                    GmfInt, &lcorn[0], &lcorn[ncorn-1]);
      }
    }catch(const MetrisExcept &e){
      std::cout<<"## Type: "<<e.what()<<std::endl;

      #ifndef NO_BOOST_EXCEPT
        if(std::string const * ms=boost::get_error_info<excMessage>(e) )
          std::cout<<"## Message: "<<*ms; 
        if(boost::stacktrace::stacktrace const * tr=boost::get_error_info<excStackTrace>(e) )
          std::cerr << "## Call stack: \n" << *tr;
      #endif
      METRIS_THROW_MSG("TODO: Manage corners (lpoic) in this context");
    }
  }else{
    GmfSetKwd(libIdx, GmfVertices, msh.npoin);
    if(msh.idim == 3){
      for(int ii = 0; ii < msh.npoin; ii++){
        GmfSetLin(libIdx, GmfVertices, msh.coord(ii,0), msh.coord(ii,1), msh.coord(ii,2),0);
      }
    }else if(msh.idim == 2){
      for(int ii = 0; ii < msh.npoin; ii++){
        GmfSetLin(libIdx, GmfVertices, msh.coord(ii,0), msh.coord(ii,1), 0);
      }
    }
  }



  GmfSetKwd(libIdx, GmfNormals, msh.npoin);
  GmfSetBlock(libIdx, GmfNormals, 1, msh.npoin, 0, NULL, NULL,
    GmfDoubleVec, msh.idim, &poi2vec[0][0], &poi2vec[msh.npoin-1][0]);

  int nentt = msh.nentt(msh.idim);
  intAr2 &ent2poi = msh.ent2poi(msh.idim);
  GmfSetKwd(libIdx, GmfNormalAtTriangleVertices, nentt);
  GmfSetBlock(libIdx, GmfNormalAtTriangleVertices, 1, nentt, 0, NULL, NULL,
    GmfIntVec, msh.idim + 1, &ent2poi(0,0), &ent2poi[nentt-1][0]);




  GmfCloseMesh( libIdx );

  for(int i =0; i < msh.nelem;i++){
    msh.tet2ref[i] -= 1;
    for(int j =0;j<getnnod3(msh.curdeg);j++)
      msh.tet2poi(i,j) -= 1;
  }
  for(int i =0; i < msh.nface;i++){
    msh.fac2ref[i] -= 1;
    for(int j =0;j<getnnod2(msh.curdeg);j++)
      msh.fac2poi(i,j) -= 1;
  }
  for(int i =0; i < msh.nedge;i++){
    msh.edg2ref[i] -= 1;
    for(int j =0;j<getnnod1(msh.curdeg);j++)
      msh.edg2poi(i,j) -= 1;
  }

  //std::cout<<" - Done  writing "<<meshName<<std::endl;
}
#endif


#if 0
void writeMesh(std::string meshName, int ideg, int ilag,
               int npoin, dblAr2 &coord,
               int nelem, intAr2 &tet2poi, intAr1 &tet2ref,
               int nface, intAr2 &fac2poi, intAr1 &fac2ref,
               int nedge, intAr2 &edg2poi, intAr1 &edg2ref){
  int iverb = 0;
  std::cout<<"-- Write file "<<meshName<<std::endl;

  int ierro = 0;
  if(ideg < 1){
    std::cout<<"Invalid maximum degree in readMesh !"<<std::endl;
    ierro = 1;
  }
  if(ideg > __MAX_LIBMESHB_DEG__){
    std::cout<<"Maximum degree supported by libmeshb = "
    <<__MAX_LIBMESHB_DEG__<<" passed "<<ideg<<std::endl;
    ierro = 1;
  }
  if(meshName.find(".mesh") == std::string::npos){
    //std::cout<<"## No filename extension provided, defaulting to .meshb\n";
    meshName += ".meshb";
  }
  // this ensures all error messages are printed
  METRIS_ENFORCE(ierro <= 0);

  int idim = coord.get_stride();

  if(idim != 2 && idim != 3) METRIS_THROW_MSG("Unsupported dimension " << idim);


  if(idim >= 3){
    for(int i =0; i < nelem;i++){
      tet2ref[i] += 1;
      for(int j =0;j<getnnod3(ideg);j++)
        tet2poi(i,j) += 1;
    }
  }
  for(int i =0; i < nface;i++){
    fac2ref[i] += 1;
    for(int j =0;j<getnnod2(ideg);j++)
      fac2poi(i,j) += 1;
  }
  for(int i =0; i < nedge;i++){
    edg2ref[i] += 1;
    for(int j =0;j<getnnod1(ideg);j++)
      edg2poi(i,j) += 1;
  }

  int64_t libIdx;
  libIdx = MetrisOpenMeshFile<GmfWrite>(meshName.c_str(), idim);


  if(!ilag) GmfSetKwd( libIdx, GmfBezierBasis, 1);

  if(nedge > 0){
  //   Note on edges: for some reason, libmeshb uses 1 index for edges
  // but the usual 3 for triangles, 4 for tets, etc.     
    if(iverb>0)std::cout<<"-- Start writing edges"<<std::endl;
    constexpr int mppe = getnnod1(METRIS_MAX_DEG);

    int fKwd = libmeshb::edgeKwds[ideg];
    int npp  = getnnod1(ideg);
    if(ideg > 1){

      int myOrd[mppe]; // see remark above
      for(int i = 0; i < npp; i++){
        for(int j = 0; j < 1; j++){
          myOrd[i+j] = ordedg.s[ideg][i][j];
        }
      }

      GmfSetKwd(libIdx, libmeshb::edgeOrdKwds[ideg], npp); 
      GmfSetBlock(libIdx, libmeshb::edgeOrdKwds[ideg], 1, npp, 0, NULL, NULL,
        GmfIntVec, npp, &myOrd[0], &myOrd[(npp-1)]);
    }

    GmfSetKwd( libIdx, fKwd, nedge);
    GmfSetBlock(libIdx, fKwd, 1, nedge, 0, NULL, NULL,
      GmfIntVec, npp, &edg2poi(0,0), &edg2poi[nedge-1][0],
      GmfInt   ,      &edg2ref[0   ], &edg2ref[nedge-1  ]);

    if(iverb>0)std::cout<<"-- Done  writing edges; nTot = "<<nedge<<std::endl;
  }



  if(nface > 0){
    if(iverb>0)std::cout<<"-- Start writing triangles"<<std::endl;
    constexpr int mppf = getnnod2(METRIS_MAX_DEG);

    int fKwd = libmeshb::faceKwds[ideg];
    int nppf = getnnod2(ideg);


    if(ideg > 1){
      int myOrd[3*mppf];
      for(int i = 0; i < nppf; i++){
        for(int j = 0; j < 3; j++){
          myOrd[i*3+j] = ordfac.s[ideg][i][j];
        }
      }
      GmfSetKwd(libIdx, libmeshb::faceOrdKwds[ideg], nppf); 
      GmfSetBlock(libIdx, libmeshb::faceOrdKwds[ideg], 1, nppf, 0, NULL, NULL,
        GmfIntVec, nppf, &myOrd[0], &myOrd[3*(nppf-1)]);

    }

    GmfSetKwd( libIdx, fKwd, nface);
    GmfSetBlock(libIdx, fKwd, 1,  nface, 0, NULL, NULL,
      GmfIntVec, nppf, &fac2poi(0,0), &fac2poi[nface-1][0],
      GmfInt   ,       &fac2ref[0   ], &fac2ref[nface-1  ]);

    if(iverb>0)std::cout<<"-- Done  writing triangles; nTot = "<<nface<<std::endl;
  }


  // Redundant but safer while we test 
  if(nelem > 0 && idim >= 3){
    if(iverb>0)std::cout<<"-- Start writing tetrahedra"<<std::endl;
    constexpr int mppt = getnnod3(METRIS_MAX_DEG);

    int eKwd = libmeshb::elemKwds[ideg];
    int nppt = getnnod3(ideg);

    if(ideg > 1){
      int myOrd[4*mppt];
      for(int i = 0; i < nppt; i++){
        for(int j = 0; j < 4; j++){
          myOrd[i*4+j] = ordtet.s[ideg][i][j];
        }
      }

      GmfSetKwd(libIdx, libmeshb::elemOrdKwds[ideg], nppt); 
      GmfSetBlock(libIdx, libmeshb::elemOrdKwds[ideg], 1, nppt, 0, NULL, NULL,
        GmfIntVec, nppt, &myOrd[0], &myOrd[4*(nppt-1)]);
    }
    GmfSetKwd( libIdx, eKwd, nelem);
    GmfSetBlock(libIdx, eKwd, 1, nelem, 0, NULL, NULL,
      GmfIntVec, nppt, &tet2poi(0,0), &tet2poi[nelem-1][0],
                GmfInt         , &tet2ref[0   ], &tet2ref[nelem-1  ]); //idx0 + nface[iDeg]-1
    if(iverb>0)std::cout<<"-- Done  writing tetrahedra; nTot = "<<nelem<<std::endl;
  }





  if(iverb>0)std::cout<<"-- Start writing points "<<npoin<<std::endl;
  /*
  Build tables for VerticesOnGeometricEntities
  */

  GmfSetKwd(libIdx, GmfVertices, npoin);
  if(idim == 2){
    for(int ii = 0; ii < npoin; ii++){
      GmfSetLin(libIdx, GmfVertices, coord(ii,0), coord(ii,1), 0);
    }
  }else if(idim == 3){
    for(int ii = 0; ii < npoin; ii++){
      GmfSetLin(libIdx, GmfVertices, coord(ii,0), coord(ii,1), coord(ii,2), 0);
    }
  }




  GmfCloseMesh( libIdx );

  if(idim >= 3){
    for(int i = 0; i < nelem; i++){
      tet2ref[i] -= 1;
      for(int j = 0; j < getnnod3(ideg); j++)
        tet2poi(i,j) -= 1;
    }
  }
  for(int i = 0; i < nface;i++){
    fac2ref[i] -= 1;
    for(int j = 0; j < getnnod2(ideg); j++)
      fac2poi(i,j) -= 1;
  }
  for(int i = 0; i < nedge; i++){
    edg2ref[i] -= 1;
    for(int j = 0; j < getnnod1(ideg); j++)
      edg2poi(i,j) -= 1;
  }

  //std::cout<<" - Done  writing "<<meshName<<std::endl;
}
#endif



} // End namespace

