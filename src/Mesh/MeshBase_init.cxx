//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "MeshBase.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_exceptions.hxx"
#include "../metris_constants.hxx"
#include "../ho_constants.hxx"
#include "../libs/libmeshb.hxx"
#include "../aux_utils.hxx"
#include "../CT_loop.hxx"
#include "../aux_topo.hxx"
#include "../aux_timer.hxx"
#include "../mprintf.hxx"
#include "../low_geo.hxx"
#include "../msh_inineigh.hxx"
#include "../ho_constants.hxx"
#include "../Boundary/msh_inisurf.hxx"
#include "../API/MetrisAPI.hxx"

namespace Metris{

void MeshBase::iniFromData(MetrisAPI &data, int usrTarDeg){
  readConstants(data, usrTarDeg);
  zeroArrays();
  readMeshData(data);
}

void MeshBase::iniFromFile(std::string fname, int usrTarDeg){
  int64_t libIdx = MetrisOpenMeshFile<GmfRead>(fname.c_str(), &idim);
  readConstants(libIdx, usrTarDeg); 
  zeroArrays();
  readMeshFile(libIdx);
  GmfCloseMesh(libIdx);
}

// data is freed 
void MeshBase::initialize(MetrisAPI *data, 
                          MetrisParameters &param){

  this->param = &param;
  const int iverb = param.iverb;

  int usrTarDeg = -1;
  if(this->meshClass() == MeshClass::Mesh){
    usrTarDeg = param.usrTarDeg;
  }else if(this->meshClass() == MeshClass::MeshBack){
    usrTarDeg = -1;
  }else{
    METRIS_THROW_MSG(WArgExcept(), "Unknown code path (1) MeshBase init with type" 
      << (int) this->meshClass());
  }
  //int usrTarDeg = imsh == whichMesh::Front ? param.usrTarDeg : -1;

  if(data == NULL){
    // If this is front, get front. 
    // Otherwise get back, but only if dedicated back, otherwise front. 
    std::string fname;
    if(this->meshClass() == MeshClass::Mesh){
      fname = param.meshFileName;
    }else if(this->meshClass() == MeshClass::MeshBack){
      fname = param.inpBack ? param.backFileName : param.meshFileName; 
    }else{
      METRIS_THROW_MSG(WArgExcept(), "Unknown code path (2) MeshBase init with type" 
        << (int) this->meshClass());
    }
    //std::string fname = imsh == whichMesh::Front ? param.meshFileName 
    //                  : param.inpBack == true ? param.backFileName : param.meshFileName; 
    iniFromFile(fname,usrTarDeg);
  }else{
    iniFromData(*data,usrTarDeg);
  }


  // Compute bounding box
  for(int i = 0; i < idim; i++){
    bb[i][0] = coord(0,i);
    bb[i][1] = coord(0,i);
  }
  for(int ipoin = 0; ipoin < npoin; ipoin++){
    for(int i = 0; i < idim; i++){
      double x = coord(ipoin,i);
      bb[i][0] = bb[i][0] < x ? bb[i][0] : x;
      bb[i][1] = bb[i][1] > x ? bb[i][1] : x;
    }
  }

  // Initialize hash tables for edges and triangles
  if(idim >= 3){
    facHshTab.reserve(nface);
    for(int iface = 0; iface < nface; iface++){
      if(isdeadent(iface,fac2poi)) continue;
      int i1 = fac2poi(iface,0);
      int i2 = fac2poi(iface,1);
      int i3 = fac2poi(iface,2);
      facHshTab.insert({stup3(i1,i2,i3),iface});
      fac2tet(iface,0) = -1;
      fac2tet(iface,1) = -1;
    }
  }
  edgHshTab.reserve(nedge);
  for(int iedge = 0;iedge < nedge;iedge++){
    if(isdeadent(iedge,edg2poi)) continue;
    int i1 = edg2poi(iedge,0);
    int i2 = edg2poi(iedge,1);
    edgHshTab.insert({stup2(i1,i2),iedge});
  }

  int nbpo0 = nbpoi; // Those before creation by neighbours, reconst

  iniNeighbours();
  iniBdryPoints();
  
  
  iniCADLink(nbpo0);


  // Orient edges so interior is to their left (manifold only, of course)
  // This is the same as saying the 1 -> 2 tangent's orthogonal in the 
  // clockwise direction is outgoing  
  if(idim == 2){

    for(int iedge = 0; iedge < nedge; iedge++){
      if(isdeadent(iedge,edg2poi)) continue;
      int ipoi1 = edg2poi(iedge,0);
      int ipoi2 = edg2poi(iedge,1);

      double tang[2] = {coord(ipoi2,0) - coord(ipoi1,0),
                        coord(ipoi2,1) - coord(ipoi1,1)};

      // ingoing normal iff positively oriented
      double normal[2] = {tang[1], -tang[0]};

      // get triangle adjacent 
      int iface = edg2fac[iedge];
      METRIS_ASSERT(iface >= 0);
      METRIS_ASSERT(!isdeadent(iface,fac2poi));

      int iver = -1;
      for(int ii = 0; ii < 3; ii++){
        int jpoi1 = fac2poi(iface,lnoed2[ii][0]);
        int jpoi2 = fac2poi(iface,lnoed2[ii][1]);
        if(ipoi1 == jpoi1 && ipoi2 == jpoi2 || 
           ipoi1 == jpoi2 && ipoi2 == jpoi1 ){
          iver = ii;
          break;
        }
      }
      METRIS_ASSERT(iver >= 0);
      int ipoi3 = fac2poi(iface,iver);

      // Any non-tangent vector does the trick (either 3 - 1 or 3 - 2)
      double du[2] = {coord(ipoi3,0) - coord(ipoi1,0),
                      coord(ipoi3,1) - coord(ipoi1,1)};
 
      double dtprd = getprdl2<2>(du,normal);
      if(dtprd <= 0) continue;

      // Positive dotprod: the normal is ingoing, not outgoing. 
      if(iverb >= 5) 
        printf(" - reorient edge %d was %d %d dtprd = %f\n",
               iedge,ipoi1,ipoi2,dtprd);
      //if(iverb >= 4){
      //  printf("tang = ");
      //  dblAr1(idim, tang).print();
      //  printf("normal = ");
      //  dblAr1(idim, normal).print();
      //  printf("du = ");
      //  dblAr1(idim, du).print();
      //}

      edg2poi(iedge,0) = ipoi2;
      edg2poi(iedge,1) = ipoi1;
      // Reverse ctrl pt order for degree 3+ 
      for(int ii = 1; ii <= (curdeg - 1) / 2; ii++){
        int jj = curdeg - ii;

        int iver1 = mul2nod(ii,jj);
        int iver2 = mul2nod(jj,ii);
        int tmp = edg2poi(iedge,iver1);
        edg2poi(iedge,iver1) = edg2poi(iedge,iver2);
        edg2poi(iedge,iver2) = tmp;
      }

      // Reverse neighbours 
      int tmp = edg2edg(iedge,0);
      edg2edg(iedge,0) = edg2edg(iedge,1);
      edg2edg(iedge,1) = tmp;

    }

  }else if(CAD()){ // dim == 3

    // Orient faces in 3D
    // The "natural" normal getnorfacP1 should be pointing outwards. 
    // Compute triangle normals, CAD normals, and (if exists) tet "normals"

    double nor_disc[3], norCAD[3], nor_tet[3];
    double dum[3];

    for(int iface = 0; iface < nface; iface++){
      if(isdeadent(iface,fac2poi)) continue;
      getnorfacP1(fac2poi[iface], coord, nor_disc);

      if(normalize_vec<3>(nor_disc)) METRIS_THROW_MSG(TODOExcept(), 
        "Error handling in face orientation disc normals.")

      // To compute the CAD normal, the safest is to average the vertex normals.
      // This is because taking the average of the (u,v)'s can send us just about
      // anywhere.
      int ierro = getnorfacCAD(*this,iface,norCAD);

      METRIS_ASSERT_MSG(ierro == 0, "Manage CAD normal errors. Stack elements"
          " with failures and deal with them in a second time.");

      if(normalize_vec<3>(norCAD)) METRIS_THROW_MSG(TODOExcept(), 
        "Error handling in face orientation CAD normals.")

      double dtprd = getprdl2<3>(norCAD, nor_disc);

      METRIS_ASSERT_MSG(abs(dtprd) >= Constants::dtprdMisAlign,
        "Check meaning of apparently very badly aligned CAD and face normal. "
        "dtprd = "<<dtprd)

      int iref = fac2ref[iface];
      if(iverb >= 4) printf("Debug iface %d iref %d dtprd = %f \n",
                             iface,iref,dtprd);

      if(dtprd < 0){
        // Misaligned, switch d00 and 0d0
        for(int i3 = 0; i3 < curdeg; i3++){
          for(int i2 = 0; i2 <= (curdeg - i3)/2; i2++){
            int i1   = curdeg - i3 - i2;
            int irnk1 = mul2nod(i1,i2,i3);
            int irnk2 = mul2nod(i2,i1,i3);
            int tmp = fac2poi(iface,irnk1);
            fac2poi(iface,irnk1) = fac2poi(iface,irnk2);
            fac2poi(iface,irnk2) = tmp;
          }
        }

        int tmp = fac2fac(iface,0);
        fac2fac(iface,0) = fac2fac(iface,1);
        fac2fac(iface,1) = tmp;
      }// endif dtprd

      #ifndef NDEBUG
      bool iflat;
      double meas =  getmeasentP1<3,2>(*this,fac2poi[iface],norCAD,&iflat);
      if(iflat || meas < 0){
        printf("## DEBUG meas = %15.7e iflat %d \n",meas,iflat);
        writeMesh("debugsurf",*this);
        METRIS_THROW(GeomExcept());
      }
      #endif

    }


    //METRIS_THROW_MSG(TODOExcept(), "Implement edge and triangle orientation in 3D")
  }




#if 0
  // For edges, the tangent deviation is computed at the vertices. 
  // check edges only belong to one loop
  geodev[0] = -1;
  int imax = -1;
  if(this->CAD()){

    for(int iedge = 0; iedge < nedge; iedge++){
      if(isdeadent(iedge,edg2poi)) continue;
      int iref = edg2ref[iedge];
      ego obj = CAD.cad2edg[iref];
      for(int ii = 0; ii < 2; ii++){
        int ipoin = edg2poi(iedge,ii);

        // Get CAD tangent 
        int ibpoi = poi2bpo[ipoin];
        METRIS_ASSERT(ibpoi >= 0);
        ibpoi = getref2bpo(*this,ibpoi,iref,1);
        double result[18];
        int ierro = EG_evaluate(obj, bpo2rbi[ibpoi], result);
        METRIS_ENFORCE(ierro == 0);
        double *tanCAD = &result[3];

        // Get elt tangent
        double tanedg[3], dum[3];
        double bary[2] = {1.0 - ii, (double)ii};
        CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == curdeg){
          if(idim == 2){
            eval1<2,ideg>(coord, edg2poi[iedge],
                          getBasis(), DifVar::Bary, DifVar::None,
                          bary, dum, tanedg, NULL);
          }else{
            eval1<3,ideg>(coord, edg2poi[iedge],
                          getBasis(), DifVar::Bary, DifVar::None,
                          bary, dum, tanedg, NULL);
          }
        }}CT_FOR1(ideg); 

        // Normalize both 
        double normCAD, normedg; 
        if(idim == 2){
          normCAD = getnrml2<2>(tanCAD);
          normedg = getnrml2<2>(tanedg);
        }else{
          normCAD = getnrml2<3>(tanCAD);
          normedg = getnrml2<3>(tanedg);
        }
        METRIS_ENFORCE(normCAD >= 1.0e-32);
        METRIS_ENFORCE(normedg >= 1.0e-32);

        for(int ii = 0; ii < idim; ii++) tanCAD[ii] /= sqrt(normCAD);
        for(int ii = 0; ii < idim; ii++) tanedg[ii] /= sqrt(normedg);

        double dtprd = idim == 2 ? getprdl2<2>(tanCAD,tanedg)
                                 : getprdl2<3>(tanCAD,tanedg);

        // Force dtprd to be positive. We can do this because edges have been 
        // oriented, so there's nothing to check here.
        dtprd = abs(dtprd);

        if(iverb >= 3){
          printf("  - iedge %d ipoin %d dtprd %f tanCAD = ",ipoin,iedge,dtprd);
          dblAr1(idim,tanCAD).print();
          printf("  - tanedg = ");
          dblAr1(idim,tanedg).print();
        }

        METRIS_ASSERT_MSG(dtprd >= 1.0e-16,"zero dtprd = "<<dtprd);

        double dev = 1 - abs(dtprd);

        if(dev >= geodev[0]){
          imax = iedge;
          geodev[0] = dev;
        }


      }
    }

  }// if msh.CAD()
  else{
    geodev[0] = 1;
  }
  
  if(iverb >= 2) printf("-- Computed max tandev edges, got %15.7e at %d\n",
                               geodev[0], imax);


  // Compute maximum normal deviation for localization tolerance 
  if(isboundary_faces()) METRIS_THROW_MSG(TODOExcept(), 
                                         "Implement geodev for triangles in 3D")
  geodev[1] = -1;
#endif

}
 


void MeshBase::iniNeighbours(){
  GETVDEPTH((*this));
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == this->curdeg){
    double t1 = get_wall_time(); 
    iniMeshNeighbours<ideg>(*this);
    double t2 = get_wall_time(); 
    CPRINTF2(" - Done neighbours, time = %7.3fs \n",t2-t1); 
  }}CT_FOR1(ideg);
}


void MeshBase::iniBdryPoints(){
  hana::while_(hana::less_equal.than(hana::int_c<METRIS_MAX_DEG>), 1_c, [&](auto ideg_c){
    constexpr int ideg = ideg_c;
    if(ideg == this->curdeg){

      if(param->iverb >= 1) printf("-- Update bdry point link to entities\n");
      int ncrea = iniMeshBdryPoints<ideg>(*this); 
      if(param->iverb >= 1) printf("   %d boundary points created\n",ncrea);
    }
  return ideg_c+1_c;});
}

void MeshBase::iniCADLink(int nbpo0){
  CAD.iniCADLink(*param,*this,nbpo0);
  cfa2tag.allocate(METRIS_MAXTAGS, CAD.ncadfa);
  cfa2tag.set_n(METRIS_MAXTAGS); 
  ced2tag.allocate(METRIS_MAXTAGS, CAD.ncaded);
  ced2tag.set_n(METRIS_MAXTAGS); 
  cno2tag.allocate(METRIS_MAXTAGS, CAD.ncadno);
  cno2tag.set_n(METRIS_MAXTAGS); 

  cfa2tag.fill(METRIS_MAXTAGS, CAD.ncadfa, 0);
  ced2tag.fill(METRIS_MAXTAGS, CAD.ncaded, 0);
  cno2tag.fill(METRIS_MAXTAGS, CAD.ncadno, 0);
}



void MeshBase::readConstants(int64_t libIdx, int usrMinDeg){

  set_npoin(GmfStatKwd( libIdx, GmfVertices ));
  if(npoin == 0) METRIS_THROW_MSG(TopoExcept(),"EMPTY MESH (NO VERTICES)");

  // We don't know yet. 
  nbpoi_ = 0;

  int isuppr    = 0;
  nedge_ = nface_ = nelem_ = 0;
  for(int iDeg = 1; iDeg <= __MAX_LIBMESHB_DEG__; iDeg++){
    int i1 = GmfStatKwd(libIdx, libmeshb::edgeKwds[iDeg]);
    int i2 = GmfStatKwd(libIdx, libmeshb::faceKwds[iDeg]);
    int i3 = GmfStatKwd(libIdx, libmeshb::elemKwds[iDeg]);
    if(i1 > 0 || i2 > 0 || i3 > 0){
      METRIS_ENFORCE(iDeg <= METRIS_MAX_DEG);
      if(param->iverb >= 1) std::cout<<"-- Mesh of degree "<<iDeg<<std::endl;

      // Degree used for storage. User wants usrMinDeg and up to METRIS_MAX_DEG if the mesh file expects more. 
      // Therefore this is min(max(usrMinDeg,iDeg),METRIS_MAX_DEG)
      int strMaxDeg = usrMinDeg > iDeg ? usrMinDeg : iDeg;
      strMaxDeg = strMaxDeg < METRIS_MAX_DEG ? strMaxDeg : METRIS_MAX_DEG;

      strdeg = strMaxDeg; 
      curdeg = iDeg;
    } 
    if((i1 > 0 && nedge > 0) || 
       (i2 > 0 && nface > 0) || 
       (i3 > 0 && nelem > 0)) METRIS_THROW_MSG(TopoExcept(),"SEVERAL DEGREES IN THE MESH");
    if(i1 > 0) set_nedge(i1);
    if(i2 > 0) set_nface(i2);
    if(i3 > 0) set_nelem(i3);
  }
  if(isuppr > 0)METRIS_THROW_MSG(TopoExcept(),"MESH HAS ELTS OF DEG > METRIS_MAX_DEG");

  METRIS_ENFORCE_MSG( idim >= 3 || nelem == 0, "TETRAHEDRA IN DIMENSION 2 FILE");
}


void MeshBase::readConstants(const MetrisAPI &data, int usrMinDeg){
  idim   = data.idim;
  curdeg = data.ideg;
  strdeg = MAX(usrMinDeg, curdeg);

  nbpoi_ = 0;
  // Only skip allocs if curdeg != strdeg. 
  set_nedge(data.nedge, true); //curdeg == strdeg
  set_nface(data.nface, true); //curdeg == strdeg
  set_nelem(data.nelem, true); //curdeg == strdeg
  set_npoin(data.npoin, true);
  METRIS_ENFORCE_MSG( idim >= 3 || nelem == 0, "TETRAHEDRA IN DIMENSION 2 DATA");
}


void MeshBase::copyConstants(const MeshBase &msh, int MAX_DEG){
  idim   = msh.idim;
  curdeg = msh.curdeg;
  strdeg = msh.curdeg; // We don't copy allocation sizes, we copy effective sizes.

  set_nbpoi(msh.nbpoi);
  set_npoin(npoin);
  set_nedge(msh.nedge);
  set_nface(msh.nface);
  set_nelem(msh.nelem);
}



void MeshBase::zeroArrays(){
  poi2ent.fill(npoin,2,-1);
  edg2fac.fill(nedge,-1);

  if(idim >= 3) fac2tet.fill(nface,2,-1);

  poi2tag.fill(METRIS_MAXTAGS,npoin,0);
  edg2tag.fill(METRIS_MAXTAGS,nedge,0);
  fac2tag.fill(METRIS_MAXTAGS,nface,0);
  if(idim >= 3) tet2tag.fill(METRIS_MAXTAGS,nelem,0);
  if(idim >= 3) tet2ftg.fill(nelem,false);

  for(int itag = 0; itag < METRIS_MAXTAGS; itag++) tag[itag] = 0;

  //if(nelem > 0){
  //  for(int iface = 0; iface < nface; iface++){
  //    for(int ii = 0; ii < 2; ii++){
  //      int ielem = fac2tet(iface,ii);
  //      if(ielem < 0) continue;
  //      tet2ftg[ielem] = true;
  //    }
  //  }
  //}
}



void MeshBase::readMeshFile(int64_t libIdx, int ithread){

  GETVDEPTH((*this));

  int ilag = 1 - GmfStatKwd(libIdx, GmfBezierBasis);
  if(ilag == 1){
    ibasis = FEBasis::Lagrange;
    CPRINTF2(" - Mesh in Lagrange format. depth\n");
  }else if(ilag == 0){
    ibasis = FEBasis::Bezier;
    CPRINTF2(" - Mesh in Bézier format.\n");
  }else{
    ibasis = FEBasis::Undefined;
    METRIS_THROW_MSG(WArgExcept(),"Invalid basis in mesh file");
  }


  if(idim != 2 && idim !=3) METRIS_THROW_MSG(WArgExcept(), "Dimension unsupported "<<idim)

  /* -------------------------------------------------------------------------------- */
  /* ----------------------------------- Points ------------------------------------- */
  /* -------------------------------------------------------------------------------- */

  CPRINTF2("-- Start reading %10d points\n",npoin);
  GmfGotoKwd( libIdx, GmfVertices );
  GmfGetBlock(libIdx, GmfVertices, 1, npoin, 0, NULL, NULL,
    GmfDoubleVec, idim, &coord(0,0), &coord[npoin-1][0],
    GmfInt            , &poi2bpo[0] , &poi2bpo[npoin-1]);
  CPRINTF2("-- Done reading %10d points\n",npoin);

  /* --------------------------------- Corners
  Point refs in file relate to corners. 
  The following is an upper bound that can never be exceeded:
  */
  //int nbpo_guess = edgnpps[strdeg]*nedge + facnpps[strdeg]*nface;
  for(int i = 0; i < nbpoi; i++){
    bpo2ibi(i,0) = -1;
    bpo2ibi(i,3) = -1;
    bpo2rbi(i,0) = 0;
    bpo2rbi(i,1) = 0;
  }

  // Initialize corners -> using corners instead 
  nbpoi_ = 0;
  for(int ipoin = 0; ipoin < npoin; ipoin++){
    poi2bpo[ipoin] = -1;
    ////poi2ent[ipoin] = -1;
    //poi2bpo[ipoin] -= 1;
    //if(poi2bpo[ipoin] < 0) continue;
    //int icorn = poi2bpo[ipoin];
    //poi2bpo[ipoin] = -1;
    //newbpotopo<0>(ipoin,icorn);
  }



  int ncorn = GmfStatKwd(libIdx, GmfCorners);
  if(ncorn > 0){
    CPRINTF2("-- Start reading %10d corners\n",ncorn);
    iwork.allocate(ncorn);
    iwork.set_n(ncorn);
    GmfGetBlock(libIdx, GmfCorners, 1, ncorn, 0, NULL, NULL,
      GmfInt, &iwork[0] , &iwork[ncorn-1]);
    // feflo.a has output some files with a bunch of 0 corners
    int ncor1 = 0; 
    for(int icorn = 0;icorn < ncorn; icorn++){
      INCVDEPTH((*this));
      int ipoin = iwork[icorn] - 1;
      if(ipoin < 0) continue;
      if(ipoin > npoin){
        CPRINTF1("## INVALID CORNERS TABLE IN FILE! %d > %d (max)\n",ipoin+1,npoin);
      }
      if(poi2bpo[ipoin] >= 0){
        CPRINTF1("## Warning: point %d already supplied as corner %d. Would have become %d \n",ipoin,
                 bpo2ibi[poi2bpo[ipoin]][2],icorn);
        continue;
      }
      newbpotopo<0>(ipoin,icorn);
      ncor1++;
    }
    //CPRINTF1(" - Added %d corners\n",ncor1);
  }

  ncorn = GmfStatKwd(libIdx, GmfVerticesOnGeometricVertices);
  if(ncorn > 0){
    CPRINTF2("-- Start reading %10d VerticesOnGeometricVertices (corners)\n",ncorn);
    iwork.allocate(2*ncorn);
    iwork.set_n(2*ncorn);
    //GmfGetBlock(libIdx, GmfVerticesOnGeometricVertices, 1, ncorn, 0, NULL, NULL,
    //            GmfIntVec, &iwork[0] , &iwork[2*ncorn-1]);
    GmfGotoKwd(libIdx, GmfVerticesOnGeometricVertices);
    for(int ii = 0; ii < ncorn; ii++){
      GmfGetLin(libIdx, GmfVerticesOnGeometricVertices, 
                &iwork[2*ii+0],&iwork[2*ii+1]);
    }

    int ncor1 = 0; 
    for(int icor0 = 0; icor0 < ncorn; icor0++){
      INCVDEPTH((*this));
      int ipoin = iwork[2*icor0    ] - 1;
      int icorn = iwork[2*icor0 + 1] - 1;
      if(ipoin < 0) continue;
      if(ipoin > npoin){
        CPRINTF1("## INVALID CORNERS TABLE IN FILE! %d > %d (max)\n",ipoin+1,npoin);
      }
      if(poi2bpo[ipoin] >= 0){
        CPRINTF1("## Warning: point %d already supplied as corner %d. Would have become %d \n",ipoin,
                 bpo2ibi[poi2bpo[ipoin]][2],icorn);
        continue;
      }
      newbpotopo<0>(ipoin,icorn);
      ncor1++;
    }
    //CPRINTF1(" - Added %d corners\n",ncor1);
  }



  bool redorefs[3] = {false, false, false};

  /* -------------------------------------------------------------------------------- */
  /* --------------------------------- Elements  ----------------------------------- */
  /* -------------------------------------------------------------------------------- */
  
  if(nelem > 0){

    CPRINTF2("-- Start reading %10d tetrahedra\n",nelem);
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == curdeg){
    //          if(nelem[iDeg] == 0) continue;
      CPRINTF2("-- Start reading %10d P%d tetrahedra\n",nelem,ideg);
  
      int eKwd = libmeshb::elemKwds[ideg];
  
  
      constexpr int nppe = tetnpps[ideg];
  
      GmfGotoKwd( libIdx, eKwd );
  
      if(ideg > 1){
        int FileOrdering[4*nppe];
        int myOrd[4*nppe];
        for(int i = 0; i < nppe; i++){
          for(int j = 0; j < 4; j++){
            myOrd[i*4+j] = ordtet.s[ideg][i][j];
          }
        }
  
        if(!GmfStatKwd(libIdx, libmeshb::elemOrdKwds[ideg])){
          gen_ordering_Vizir<ideg,3>(FileOrdering);
          std::cout<<"## No ordering given in file ! Use e.g. GmfTetrahedraP2Ordering."<<
          std::endl<<"Defaulting to Vizir4 (\"P.-L.\") ordering"<<std::endl;
        }else{
          GmfGetBlock(libIdx, libmeshb::elemOrdKwds[ideg], 1, nppe, 0, NULL,
            NULL, GmfIntTab, nppe, &FileOrdering[0], &FileOrdering[4*(nppe-1)]);
        }
        GmfSetHONodesOrdering(  libIdx, eKwd,
          &myOrd[0], FileOrdering );
      }
  
      GmfGetBlock(libIdx, eKwd, 1, nelem, 0, NULL, NULL,
        GmfIntVec, nppe, &tet2poi(0,0), &tet2poi[nelem-1][0],
        GmfInt         , &tet2ref[0   ], &tet2ref[nelem-1  ]);
    }}CT_FOR1(ideg);
  
    for(int i = 0; i< nelem;i++){
      tet2ref[i] -= 1;
      for(int j = 0; j < tetnpps[curdeg]; j++){
        tet2poi(i,j) -= 1;
      }
      if(isdeadent(i,tet2poi)) continue;
      for(int j = 0; j < tetnpps[curdeg]; j++){
        poi2ent[tet2poi(i,j)][0] = i ;
        poi2ent[tet2poi(i,j)][1] = 3 ;
      }
    }
    
    if(param->refineConventions){
      int nseen = 0;
      bool ineg = false;
      for(int ielem = 0; ielem < nelem; ielem++){
        if(isdeadent(ielem,tet2poi)) continue;
        nseen++;
        if(nseen >= 100 && !ineg) break;
        double meas;
        if(!ineg){
          meas = getmeasentP1<3>(tet2poi[ielem], coord);
        }
        if(ineg || meas < param->vtol){
          METRIS_ENFORCE_MSG(ineg || nseen == 1, "## FIRST NEGATIVE ELEMENT IS RANK "<<nseen
            << " meas "<<meas);
          ineg = true;
          int tmp = tet2poi(ielem,0);
          tet2poi(ielem,0) = tet2poi(ielem,1);
          tet2poi(ielem,1) = tmp;
        }
      }
      if(ineg) CPRINTF1("## FLIPPED ELEMENT SIGNS !\n");
    }

    CPRINTF2("-- Done reading %10d tetrahedra\n",nelem);
  }



  if(nface > 0){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == curdeg){
      CPRINTF2("-- Start reading %10d P%d triangles\n",nface,ideg);
  
      int fKwd = libmeshb::faceKwds[ideg];
  
  
      constexpr int nppf = facnpps[ideg];
  
      GmfGotoKwd( libIdx, fKwd );
  
      if(ideg > 1){
        int FileOrdering[3*nppf];
        int myOrd[3*nppf];
        for(int i = 0; i < nppf; i++){
          for(int j = 0; j < 3; j++){
            myOrd[i*3+j] = ordfac.s[ideg][i][j];
          }
        }
  
        if(!GmfStatKwd(libIdx, libmeshb::faceOrdKwds[ideg])){
          gen_ordering_Vizir<ideg,2>(FileOrdering);
          std::cout<<"!! No ordering given in file ! Use e.g. GmfTrianglesP2Ordering."<<
          std::endl<<"Defaulting to Vizir4 (\"P.-L.\") ordering"<<std::endl;
        }else{
          GmfGetBlock(libIdx, libmeshb::faceOrdKwds[ideg], 1, nppf, 0, NULL,
            NULL, GmfIntVec, nppf, &FileOrdering[0], &FileOrdering[3*(nppf-1)]);
        }
        GmfSetHONodesOrdering(  libIdx, fKwd,
          &myOrd[0], FileOrdering );
      }
  
      GmfGetBlock(libIdx, fKwd, 1,nface, 0, NULL, NULL,
        GmfIntVec, nppf, &fac2poi(0,0), &fac2poi[nface-1][0],
        GmfInt   ,       &fac2ref[0   ], &fac2ref[nface-1  ]);
    }}CT_FOR1(ideg);
    
    for(int i = 0; i < nface;i++){
      fac2ref[i] -= 1;
      for(int j = 0; j < facnpps[curdeg]; j++){
        fac2poi(i,j) -= 1;
      }
      if(fac2ref[i] < 0) redorefs[2] = true; 

      if(isdeadent(i,fac2poi)) continue;
      for(int jj = 0; jj < facnpps[curdeg]; jj++){
        poi2ent[fac2poi(i,jj)][0] = i;
        poi2ent[fac2poi(i,jj)][1] = 2;
      }
    }

    if(param->refineConventions){
      int nseen = 0;
      bool ineg = false;
      for(int iface = 0; iface < nface; iface++){
        if(isdeadent(iface,fac2poi)) continue;
        nseen++;
        if(nseen >= 100 && !ineg) break;
        double meas;
        if(!ineg){
          meas = getmeasentP1<2>(fac2poi[iface], coord);
        }
        if(ineg || meas < param->vtol){
          METRIS_ENFORCE_MSG(ineg || nseen == 1, "## FIRST NEGATIVE ELEMENT IS RANK "<<nseen
            <<" meas "<<meas);
          ineg = true;
          int tmp = fac2poi(iface,0);
          fac2poi(iface,0) = fac2poi(iface,1);
          fac2poi(iface,1) = tmp;
        }
      }
      if(ineg) CPRINTF1("## FLIPPED ELEMENT SIGNS !\n");
    }

    CPRINTF2("-- Done reading %10d triangles\n",nelem);
  }


  if(nedge > 0){


    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == curdeg){
      CPRINTF2("-- Start reading %10d P%d edges\n",nedge,ideg);
      int fKwd = libmeshb::edgeKwds[ideg];
  
  
      constexpr int npp = edgnpps[ideg];
  
      GmfGotoKwd( libIdx, fKwd );
  
      if(ideg > 1){
        int FileOrdering[3*npp];
        int myOrd[npp]; // see remark above
        for(int i = 0; i < npp; i++){
          for(int j = 0; j < 1; j++){
            myOrd[i+j] = ordedg.s[ideg][i][j];
          }
        }
        if(!GmfStatKwd(libIdx, libmeshb::edgeOrdKwds[ideg])){
          gen_ordering_Vizir<ideg,1>(FileOrdering);
          std::cout<<"!! No ordering given in file ! Use e.g. GmfEdgesP2Ordering."<<
          std::endl<<"Defaulting to Vizir4 (\"P.-L.\") ordering"<<std::endl;
        }else{
          GmfGetBlock(libIdx, libmeshb::edgeOrdKwds[ideg], 1, npp, 0, NULL,
            NULL, GmfIntVec, npp, &FileOrdering[0], &FileOrdering[npp-1]);
        }
        GmfSetHONodesOrdering(  libIdx, fKwd,
          &myOrd[0], FileOrdering );
      }
  
      GmfGetBlock(libIdx, fKwd, 1,nedge, 0, NULL, NULL,
        GmfIntVec, npp, &edg2poi(0,0),&edg2poi[nedge-1][0],
        GmfInt   ,      &edg2ref[0   ],&edg2ref[nedge-1  ]);
  
    }}CT_FOR1(ideg);

    for(int i = 0; i < nedge; i++){
      edg2ref[i] -= 1;
      for(int j = 0; j < edgnpps[curdeg]; j++){
        edg2poi(i,j) -= 1;
      }
      if(edg2ref[i] < 0) redorefs[1] = true; 
      if(isdeadent(i,edg2poi)) continue;
      for(int j = 0; j < edgnpps[curdeg]; j++){
        poi2ent[edg2poi(i,j)][0] = i;
        poi2ent[edg2poi(i,j)][1] = 1;
      }
    }
    CPRINTF2("-- Done reading %10d edges\n",nedge);
  }


  /* Correct references */

  intAr1 entstack(100); 
  for(int tdims = 1; tdims <= 2; tdims++){
    if(!redorefs[tdims]) continue;
    CPRINTF1(" - Correct dim %d references -> invalid provided ! \n",tdims);

    intAr2 &ent2tag_ = ent2tag(tdims);
    intAr2 &ent2ent_ = ent2ent(tdims);
    intAr1 &ent2ref_ = ent2ref(tdims);
    intAr2 &ent2poi_ = ent2poi(tdims);
    int nentt_ = nentt(tdims);

    tag[ithread]++;
    int ient0 = 0;
    entstack.set_n(0);
    int nref = 0;
    bool dosomething; 
    do{
      dosomething = false;
      for(int ientt = ient0; ientt < nentt_; ientt++){
        if(isdeadent(ientt,ent2poi_))continue;
        if(ent2tag_[ithread][ientt] >= tag[ithread]) continue;
        ient0 = ientt+1;
        entstack.stack(ientt);
        ent2tag_[ithread][ientt] = tag[ithread];
        nref++;
        dosomething = true;
        while(entstack.get_n() > 0){
          int ient1 = entstack.pop();
          ent2ref_[ient1] = nref;
          for(int ii = 0; ii < 2; ii++){
            int ientv = ent2ent_[ient1][ii];
            if(ientv < 0) continue;
            if(ent2tag_[ithread][ientv] >= tag[ithread]) continue;
            ent2tag_[ithread][ientv] = tag[ithread];
            entstack.stack(ientv);
          }
        }
      }
    }while(dosomething);
  }
  

  /* -------------------------------------------------------------------------------- */
  /* ---------------------------------- Bdry info ----------------------------------- */
  /* -------------------------------------------------------------------------------- */

  //// CAD 
  //int nbyte;
  //char *stream = GmfReadByteFlow(libIdx, &nbyte);
  //if(stream != NULL && nbyte > 0){
  //  if(iverb >= 1) printf(" - Read %db CAD stream\n",nbyte);
  //  CAD.setModel((size_t)nbyte, stream);
  //}

  int nwarn = 0;
  int mwarn = 5;
  int ngpoe = GmfStatKwd(libIdx, GmfVerticesOnGeometricEdges);
  if(ngpoe > 0){
    CPRINTF2(" - File has %d bdry pts -> edge links\n",ngpoe);
    if(param->refineConventions) CPRINTF1(" - Using refine convention\n");
    intAr2 lgpoe(ngpoe,2);
    dblAr2 rgpoe(ngpoe,2); 
    lgpoe.set_n(ngpoe);
    rgpoe.set_n(ngpoe);


    GmfGotoKwd(libIdx, GmfVerticesOnGeometricEdges);
    for(int ii = 0; ii < ngpoe; ii++){
      GmfGetLin(libIdx, GmfVerticesOnGeometricEdges, 
                &lgpoe(ii,0),&lgpoe(ii,1),
                &rgpoe(ii,0),&rgpoe(ii,1));
    }

    tag[0]++;
    int maxtag = tag[0];
    for(int igpoe = 0; igpoe < ngpoe; igpoe++){
      int ipoin = lgpoe(igpoe,0) - 1;
      // If in refineConvention, this will be a ref 1-n, otherwise an edge
      int iedge = param->refineConventions ? -lgpoe(igpoe,1) 
                                           : (lgpoe(igpoe,1) - 1);
      if(ipoin < 0 || iedge < 0 && !param->refineConventions){
        printf("## WARNING invalid entry %d/%d in GmfVerticesOnGeometricEdges: %d %d \n",
          igpoe,ngpoe,ipoin,iedge);
        continue;
      }
      METRIS_ASSERT_MSG(iedge >= 0 && iedge < nedge || param->refineConventions,
        "iedge = "<<iedge<<" refineConventions = "<<param->refineConventions
        << " lgpoe = "<<lgpoe(igpoe,1));
      //METRIS_ASSERT_MSG(iedge >= 0 && iedge < CAD.ncaded || !param->refineConventions,
      //                  "Invalid edge reference in refine convention iedge "
      //                  <<iedge<<" CAD.ncaded "<<CAD.ncaded<< " refine conv "
      //                  <<param->refineConventions)

      if(!param->refineConventions && isdeadent(iedge,edg2poi)){
        if(nwarn++ < mwarn) 
          printf("## FILE CONTAINS IBPOS POINTING TO DEAD ENTITIES\n");
        continue;
      }


      if(!param->refineConventions){
        if(edg2tag(0,iedge) < tag[0]) 
          edg2tag(0,iedge) = tag[0];
        else{
          edg2tag(0, iedge)++;
          maxtag = MAX(maxtag, edg2tag(0,iedge));
        }
        if(edg2tag(0, iedge) > 10){
          printf("## EDGE %d = %d %d REFERENCED > 10 VerticesOnGeometricEdges\n",
                 iedge,edg2poi(iedge,0),edg2poi(iedge,1));
          printf("Is this a refine mesh ? Kill and restart with -refine-conventions.\n");
          wait();
        }
      }

      // The link created in edg2bpo is temporary: if the point turns out 
      // not to be a corner then there is no need to keep the link. 
      // It will be deleted in iniMeshBdryPoints. 
      int ibpoi = newbpotopo<1>(ipoin,iedge);
      if(ibpoi < 0) continue;
      bpo2rbi(ibpoi,0) = rgpoe(igpoe,0);
      bpo2rbi(ibpoi,1) = 0.0;

    }
    tag[0] = maxtag;

    if(param->refineConventions){
      for(int iedge = 0; iedge < nedge; iedge++){
        INCVDEPTH((*this));
        if(isdeadent(iedge,edg2poi)) continue;
        for(int ii = 0; ii < 2; ii++){
          int ipoin = edg2poi(iedge,ii);
          int ibpoi = poi2bpo[ipoin];
          METRIS_ASSERT(ibpoi >= 0 && ibpoi < nbpoi);
          // Find first ibpoi entry that is associated to an edge of same ref
          // and does not have an entity yet. 
          for(;ibpoi >= 0; ibpoi = bpo2ibi(ibpoi,3)){
            int itype = bpo2ibi(ibpoi,1);
            if(itype != 1) continue;
            int ientt = bpo2ibi(ibpoi,2);
            if(ientt >= 0) continue;
            // In refine convention, the onGeometricEdges entry stores the ref
            // we put here - the entry. 
            int iref = - ientt - 1;
            //if(ipoin == 7){
            //  printf("ipoin 7 ibpoi %d type %d ientt %d iref %d edg2ref[iedge] %d\n",
            //    ibpoi,itype,ientt,iref,edg2ref[iedge]);
            //}
            if(iref != edg2ref[iedge]) continue;
            bpo2ibi(ibpoi,2) = iedge;
            CPRINTF1(" - create link ipoin %d ibpoi %d -> edge %d\n"
                                         ,ipoin, ibpoi, iedge);
            //break;
          }
        }
      }
    }

  }


  int ngpof = GmfStatKwd(libIdx, GmfVerticesOnGeometricTriangles);
  if(ngpof > 0){
    CPRINTF1(" - File contains %d boundary points -> face links\n",ngpof);
    intAr2 lgpof(ngpof,2);
    dblAr2 rgpof(ngpof,3);
    lgpof.set_n(ngpof);
    rgpof.set_n(ngpof);

    GmfGotoKwd(libIdx, GmfVerticesOnGeometricTriangles);  
    for(int ii = 0; ii < ngpof; ii++){
      GmfGetLin(libIdx, GmfVerticesOnGeometricTriangles, 
                &lgpof(ii,0),&lgpof(ii,1),
                &rgpof(ii,0),&rgpof(ii,1),&rgpof(ii,2));
    }

    for(int igpof = 0; igpof < ngpof; igpof++){
      int ipoin = lgpof(igpof,0) - 1;
      int iface = param->refineConventions ? -lgpof(igpof,1) 
                                           : (lgpof(igpof,1) - 1);
      if(!param->refineConventions && isdeadent(iface,fac2poi)){
        if(nwarn++ < mwarn) printf("## FILE CONTAINS IBPOS POINTING TO DEAD ENTITIES\n");
        continue;
      }


      int ibpoi = newbpotopo<2>(ipoin,iface);
      if(ibpoi < 0) continue;
      // Third value is unused
      bpo2rbi(ibpoi,0) = rgpof(igpof,0);
      bpo2rbi(ibpoi,1) = rgpof(igpof,1);
    }

    if(param->refineConventions && idim >= 3){
      for(int iface = 0; iface < nface; iface++){
        if(isdeadent(iface,fac2poi)) continue;
        for(int ii = 0; ii < 3; ii++){
          int ipoin = fac2poi(iface,ii);
          int ibpoi = poi2bpo[ipoin];
          METRIS_ASSERT(ibpoi >= 0 && ibpoi < nbpoi);
          // Find first ibpoi entry that is associated to an face of same ref
          // and does not have an entity yet. 
          for(;ibpoi >= 0; ibpoi = bpo2ibi(ibpoi,3)){
            int itype = bpo2ibi(ibpoi,1);
            if(itype != 2) continue;
            int ientt = bpo2ibi(ibpoi,2);
            if(ientt >= 0) continue;
            // In refine convention, the onGeometricfaces entry stores the ref
            // we put here - the entry. 
            int iref = - ientt - 1;
            if(iref != fac2ref[iface]) continue;
            bpo2ibi(ibpoi,2) = iface;
            break;
          }
        }
      }
    }

  }

}


void MeshBase::readMeshData(MetrisAPI &data){
  ibasis = data.mshbasis;
  METRIS_ASSERT(ibasis == FEBasis::Lagrange || ibasis == FEBasis::Bezier);

  if(idim != 2 && idim !=3) METRIS_THROW_MSG(WArgExcept(), "Dimension unsupported "<<idim)

  coord = std::move(data.coord);
  poi2bpo.fill(npoin,-1);
  //for(int ipoin = 0; ipoin < npoin; ipoin++){
  //  for(int ii = 0; ii < idim; ii++){
  //    coord(ipoin,ii) = data.coord(ipoin,ii);
  //  }
  //  poi2bpo[ipoin] = -1;
  //}

  CAD = std::move(data.CAD); 

  //nbpoi_ = data.ncorn;
  for(int icorn = 0; icorn < data.ncorn; icorn++){
    int ipoin = data.lcorn[icorn];
    // even metrisAPI deg change can't remove a corner
    METRIS_ASSERT(ipoin >= 0 && ipoin < npoin);
    newbpotopo<0>(ipoin,icorn);
  }
  data.lcorn.free();



  if(curdeg == strdeg){
    tet2poi = std::move(data.tet2poi);
    fac2poi = std::move(data.fac2poi);
    edg2poi = std::move(data.edg2poi);
  }else{
    data.tet2poi.copyTo(tet2poi); 
    data.tet2poi.free();
    data.fac2poi.copyTo(fac2poi); 
    data.fac2poi.free();
    data.edg2poi.copyTo(edg2poi); 
    data.edg2poi.free();
  }
  tet2ref = std::move(data.tet2ref);
  fac2ref = std::move(data.fac2ref);
  edg2ref = std::move(data.edg2ref);


  int nwarn = 0;
  int mwarn = 5;
  for(int igpoe = 0; igpoe < data.ngpoe; igpoe++){
    int ipoin = data.lgpoe[igpoe][0];
    // ipoin < 0 if dead and >= npoin if degree has been reduced in API
    if(ipoin < 0 || ipoin >= npoin) continue;
    int iedge = data.lgpoe[igpoe][1];
    if(isdeadent(iedge,edg2poi)){
      if(nwarn++ < mwarn) printf("## FILE CONTAINS IBPOS POINTING TO DEAD EDGES");
      continue;
    }
    int ibpon = newbpotopo<1>(ipoin,iedge);
    bpo2rbi(ibpon,0) = data.rgpoe[igpoe][0];
    bpo2rbi(ibpon,1) = 0.0;
  }
  data.lgpoe.free();
  data.rgpoe.free();

  nwarn = 0;
  for(int igpof = 0; igpof < data.ngpof; igpof++){
    int ipoin = data.lgpof[igpof][0];
    // ipoin < 0 if dead and >= npoin if degree has been reduced in API
    if(ipoin < 0 || ipoin >= npoin) continue;
    int iface = data.lgpof[igpof][1];
    if(isdeadent(iface,fac2poi)){
      if(nwarn++ < mwarn) printf("## FILE CONTAINS IBPOS POINTING TO DEAD TRIANGLES");
      continue;
    }
    int ibpon = newbpotopo<2>(ipoin,iface);
    bpo2rbi(ibpon,0) = data.rgpof[igpof][0];
    bpo2rbi(ibpon,1) = data.rgpof[igpof][1];
  }
  data.lgpof.free();
  data.rgpof.free();
  



  int nnode;
  nnode = tetnpps[curdeg];
  for(int ielem = 0; ielem < nelem; ielem++){
    //for(int ii = 0; ii < nnode; ii++){
    //  tet2poi(ielem,ii) = data.tet2poi(ielem,ii);
    //}
    //tet2ref[ielem] = data.tet2ref[ielem];
    if(isdeadent(ielem,tet2poi)) continue;
    for(int ii = 0; ii < nnode; ii++){
      poi2ent[tet2poi(ielem,ii)][0] = ielem;
      poi2ent[tet2poi(ielem,ii)][1] = 3;
    }
  }

  nnode = facnpps[curdeg];
  for(int iface = 0; iface < nface; iface++){
    //for(int ii = 0; ii < nnode; ii++){
    //  fac2poi(iface,ii) = data.fac2poi(iface,ii);
    //}
    //fac2ref[iface] = data.fac2ref[iface];
    if(isdeadent(iface,fac2poi)) continue;
    for(int ii = 0; ii < nnode; ii++){
      poi2ent[fac2poi(iface,ii)][0] = iface;
      poi2ent[fac2poi(iface,ii)][1] = 2;
    }
  }

  nnode = edgnpps[curdeg];
  for(int iedge = 0; iedge < nedge; iedge++){
    //for(int ii = 0; ii < nnode; ii++){
    //  edg2poi(iedge,ii) = data.edg2poi(iedge,ii);
    //}
    //edg2ref[iedge] = data.edg2ref[iedge];
    if(isdeadent(iedge,edg2poi)) continue;
    for(int ii = 0; ii < nnode; ii++){
      poi2ent[edg2poi(iedge,ii)][0] = iedge; 
      poi2ent[edg2poi(iedge,ii)][1] = 1; 
    }
  }


}



} // End namespace
