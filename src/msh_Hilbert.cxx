//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php




#if 0
template<int ideg>
void reoderHilbert(MeshBase &msh){

  double t0,t1;

//  int nthread = MAX(GetNumberOfCores(),4);
  int nthread = GetNumberOfCores();
  if(msh.nproc > 0) nthread = MIN(nthread,msh.nproc);
  printf("Running Hilbert with %d threads\n",nthread);
  
  uint64_t LPlibIdx = InitParallel(nthread);
  double bbLPlib[6];
  for(int ii = 0; ii < 3 ; ii++){
    bbLPlib[  ii] = msh.bb[ii][0];
    bbLPlib[3+ii] = msh.bb[ii][1];
  }
  uint64_t *lorder = NULL;
  bool ialloc = false;
  int fac = sizeof(uint64_t)/sizeof(int);
  int isear = 0;
  if(msh.poi2iwk.size() >= 2 * msh.npoin * fac){
    lorder = (uint64_t *)&msh.poi2iwk[0];
    isear = 1;
  }else if(msh.edg2iwk.size() >= 2 * msh.npoin * fac){
    lorder = (uint64_t *)&msh.edg2iwk[0];
    isear = 2;
  }else if(msh.fac2iwk.size() >= 2 * msh.npoin * fac){
    lorder = (uint64_t *)&msh.fac2iwk[0];
    isear = 3;
  }else if(msh.tet2iwk.size() >= 2 * msh.npoin * fac){
    lorder = (uint64_t *)&msh.tet2iwk[0];
    isear = 4;
  }else{
    lorder = new uint64_t[msh.npoin];
    if(lorder == NULL) METRIS_THROW(DMemExcept());
    ialloc = true;
  }

  double *crd = msh.coord[0];

  t0 = get_wall_time();
  HilbertRenumbering(LPlibIdx, msh.npoin, bbLPlib, &crd[-3], &lorder[-2]);
  t1 = get_wall_time();
  printf("HilbertRenumbering() call time %f \n",t1-t0);

  int *invord = NULL;
  bool iallo2 = false;
  if(msh.poi2iwk.size() >= msh.npoin && isear != 1){
    invord = &msh.poi2iwk[0];
  }else if(msh.edg2iwk.size() >= msh.npoin && isear != 2){
    invord = &msh.edg2iwk[0];
  }else if(msh.fac2iwk.size() >= msh.npoin && isear != 3){
    invord = &msh.fac2iwk[0];
  }else if(msh.tet2iwk.size() >= msh.npoin && isear != 4){
    invord = &msh.tet2iwk[0];
  }else{
    invord = new int[msh.npoin];
    if(invord == NULL) METRIS_THROW(DMemExcept());
    iallo2 = true;
  }

  // To clear the confusion on the order. 
  // Originally io0 = 0, io1 = 1
  constexpr int io0 = 1;
  constexpr int io1 = 0;

  // First link everyone to whoever is pointing to us
  for(int ii = 0; ii < msh.npoin; ii++){
    int idx1 = (int)lorder[2*ii + io1] - 1;
    invord[idx1] = ii;
  }


  msh.tag[0]++;
  for(int ii = 0; ii < msh.npoin; ii++){
    if(msh.poi2tag(0,ii) == msh.tag[0]) continue;
    msh.poi2tag(0,ii) = msh.tag[0];

    int idx0 = (int)lorder[2*ii + io0] - 1;
    int idx1 = (int)lorder[2*ii + io1] - 1; 

    if(idx0 == idx1) continue;

    double tmp[3]; 
    for(int kk = 0; kk < 3; kk++) tmp[kk] = msh.coord[idx1][kk]; 
    int itmp1 = msh.poi2bpo[idx1];
    int itmp2 = msh.poi2ent[idx1];

    for(int kk = 0; kk < 3; kk++) msh.coord[idx1][kk] = msh.coord[idx0][kk]; 
    msh.poi2bpo[idx1] = msh.poi2bpo[idx0];
    msh.poi2ent[idx1] = msh.poi2ent[idx0];
    //if(idx1 == 2) printf(" (beg) Point 2 replaced by %d \n",idx0);

    //printf("Debug chain start at %d idx0 id1 %d %d \n",ii,idx0,idx1);

    // We now want the ii that yields idx0 as its idx1. 
    int jp, jj = invord[idx0]; 

    //printf("Debug init ii idx0 idx1 %d %d %d first jj = %d \n",ii,idx0,idx1,jj);fflush(stdout);
    //int idx0_n, idx1_n;
    int nchain = 0;
    do{
      // Go to next
      if(msh.poi2tag(0,jj) >= msh.tag[0]) METRIS_THROW_MSG(AlgoExcept(),"ERROR PERMUTING POINTS");
      msh.poi2tag(0,jj) = msh.tag[0];

      //int idx1_n = (int)lorder[2*jj + 1] - 1;
      //for(int kk = 0; kk < 3; kk++) msh.coord[idx1][kk] = tmp[idx1_n]; 
//
      //idx1 = idx1_n;
      idx0 = (int)lorder[2*jj + io0] - 1;
      idx1 = (int)lorder[2*jj + io1] - 1;
      for(int kk = 0; kk < 3; kk++) msh.coord[idx1][kk] = msh.coord[idx0][kk]; 
      msh.poi2bpo[idx1] = msh.poi2bpo[idx0];
      msh.poi2ent[idx1] = msh.poi2ent[idx0];
      //if(idx1 == 2) printf(" (mid) Point 2 replaced by %d \n",idx0);

      //printf("Debug loop jj idx0 idx1 %d %d %d \n",jj,idx0,idx1);fflush(stdout);

      jp = jj;
      jj = invord[idx0]; 

      nchain++;
      if(nchain > msh.npoin) METRIS_THROW_MSG(AlgoExcept(),"INFINITE CHAIN OF POINTS?")
    }while(jj != ii);

    if( jp == ii) METRIS_THROW_MSG(AlgoExcept(),"Empty chain")


    idx1 = (int)lorder[2*jp + io1] - 1;
    for(int kk = 0; kk < 3; kk++) msh.coord[idx1][kk] = tmp[kk];
    msh.poi2bpo[idx1] = itmp1;
    msh.poi2ent[idx1] = itmp2;
  }


  // Now all entities: tetrahedra, etc. 

  // First, setup invord to point from idx0 rather than idx1. 
  // Indeed, it is idx0 that is replaced by idx1. 
  for(int ii = 0; ii < msh.npoin; ii++){
    int idx0 = (int)lorder[2*ii + io0] - 1;
    invord[idx0] = ii;
  }



  void (*tet_loop)(int,int,int,Mesh*,int *, uint64_t *)
  =[] (int ilo0, int ilo1, int , Mesh *msh, int *invord, uint64_t *lorder){
    int nnt = tetnpps[msh->curdeg];
    for(int ielem = ilo0-1; ielem < ilo1; ielem++){
      for(int ii = 0; ii < nnt; ii++){
        int ipoin = msh->tet2poi(ielem,ii);
        // This ipoin is someone's idx0:
        int iidx = invord[ipoin]; 
        int idx1 = (int)lorder[2*iidx + io1] - 1;
        msh->tet2poi(ielem,ii) = idx1;
//        msh->poi2ent[idx1] = ielem;
      }
    }
  };

  void (*fac_loop)(int,int,int,Mesh*,int *, uint64_t *)
  =[] (int ilo0, int ilo1, int , Mesh *msh, int *invord, uint64_t *lorder){
    int nnf = facnpps[msh->curdeg];
    for(int iface = ilo0-1; iface < ilo1; iface++){
      for(int ii = 0; ii < nnf; ii++){
        int ipoin = msh->fac2poi(iface,ii);
        // This ipoin is someone's idx0:
        int iidx = invord[ipoin]; 
        int idx1 = (int)lorder[2*iidx + io1] - 1;
        msh->fac2poi(iface,ii) = idx1;
//        msh->poi2ent[idx1] = iface;
      }
    }
  };

  void (*edg_loop)(int,int,int,Mesh*,int *, uint64_t *)
  =[] (int ilo0, int ilo1, int , Mesh *msh, int *invord, uint64_t *lorder){
    int nne = edgnpps[msh->curdeg];
    for(int iedge = ilo0-1; iedge < ilo1; iedge++){
      for(int ii = 0; ii < nne; ii++){
        int ipoin = msh->edg2poi(iedge,ii);
        // This ipoin is someone's idx0:
        int iidx = invord[ipoin]; 
        int idx1 = (int)lorder[2*iidx + io1] - 1;
        msh->edg2poi(iedge,ii) = idx1;
//        msh->poi2ent[idx1] = iedge;
      }
    }
  };

  void (*bpo_loop)(int,int,int,Mesh*)
  =[] (int ilo0, int ilo1, int , Mesh *msh){
    for(int ipoin = ilo0-1; ipoin < ilo1; ipoin++){
      if(msh->poi2ent[ipoin] < 0) continue;
      int ibpoi = msh->poi2bpo[ipoin];
      if(ibpoi < 0) continue;
      msh->bpo2ibi[ibpoi][0] = ipoin;
    }
    //for(int ibpoi = ilo0-1; ibpoi < ilo1; ibpoi++){
    //  int ipoin = msh->bpo2ibi[ibpoi][0];
    //  if(ipoin < 0 || ipoin >= msh->npoin) continue;
    //  int iidx = invord[ipoin]; 
    //  int idx1 = (int)lorder[2*iidx + io1] - 1;
    //  msh->bpo2ibi[ibpoi][0] = idx1;
    //}
  };

  // Now replace
  int LP_tet = NewType(LPlibIdx, msh.nelem);
  int LP_fac = NewType(LPlibIdx, msh.nface);
  int LP_edg = NewType(LPlibIdx, msh.nedge);
  int LP_poi = NewType(LPlibIdx, msh.npoin);
  float acc; 
  acc = LaunchParallelMultiArg(LPlibIdx, LP_tet, 0, (void*) tet_loop,
                               3, &msh, invord, lorder);
  printf("Acceleration factor (tet) %f \n",acc);
  //tet_loop(1,msh.nelem,0,&msh,invord,lorder);
  if(msh.nface > 100000){
    acc = LaunchParallelMultiArg(LPlibIdx, LP_fac, 0, (void*) fac_loop,
                                 3, &msh, invord, lorder);
    printf("Acceleration factor (fac) %f \n",acc);
  }else{
    fac_loop(1,msh.nface,0,&msh,invord,lorder);
  }
  if(msh.nedge > 100000){
    acc = LaunchParallelMultiArg(LPlibIdx, LP_edg, 0, (void*) edg_loop,
                                 3, &msh, invord, lorder);
    printf("Acceleration factor (edg) %f \n",acc);
  }else{
    edg_loop(1,msh.nedge,0,&msh,invord,lorder);
  }
  
  acc = LaunchParallelMultiArg(LPlibIdx, LP_poi, 0, (void*) bpo_loop,
                               1, &msh);
  printf("Acceleration factor (bpo) %f \n",acc);



  //int sztet = tetnpps[ideg]*sizeof(int);
  //int szfac = facnpps[ideg]*sizeof(int);
  //int szedg = edgnpps[ideg]*sizeof(int);




  //int (*cmp)(const void*, const void*) = [](const void* x, const void* y){
  //  int *ix = (int*) x;
  //  int *iy = (int*) y;
  //  if(*ix < *iy) return -1;
  //  else if (*ix == *iy) return 0;
  //  else return 1;
  //};

  auto cmpe = [&](int i, int j) -> bool {
    assert(i>= 0 && i < msh.nedge);
    assert(j>= 0 && j < msh.nedge);
    int ip1 = MIN(msh.edg2poi(i,0),msh.edg2poi(i,1));
    int ip2 = MIN(msh.edg2poi(j,0),msh.edg2poi(j,1));
    if(ip1 > ip2) return false;
    else return true;
  };
  auto cmpf = [&](int i, int j) -> bool {
    assert(i>= 0 && i < msh.nface);
    assert(j>= 0 && j < msh.nface);
    int ip1 = MIN(MIN(msh.fac2poi(i,0),msh.fac2poi(i,1)),msh.fac2poi(i,2));
    int ip2 = MIN(MIN(msh.fac2poi(j,0),msh.fac2poi(j,1)),msh.fac2poi(j,2));
    if(ip1 > ip2) return false;
    else return true;
  };
  auto cmpt = [&](int i, int j) -> bool {
    assert(i>= 0 && i < msh.nelem);
    assert(j>= 0 && j < msh.nelem);
    int ip1 = MIN(MIN(MIN(msh.tet2poi(i,0),msh.tet2poi(i,1)),msh.tet2poi(i,2)),msh.tet2poi(i,3));
    int ip2 = MIN(MIN(MIN(msh.tet2poi(j,0),msh.tet2poi(j,1)),msh.tet2poi(j,2)),msh.tet2poi(j,3));
    if(ip1 > ip2) return false;
    else return true;
  };


  //printf("Tet 1-10 pre\n");
  t0 = get_wall_time();
  //ParallelQsort(LPlibIdx, (void *)msh.tet2poi[0], msh.nelem, sztet, cmp);
  //ParallelQsort(LPlibIdx, (void *)msh.fac2poi[0], msh.nface, szfac, cmp);
  //ParallelQsort(LPlibIdx, (void *)msh.edg2poi[0], msh.nedge, szedg, cmp);
  int nm = MAX(MAX(MAX(msh.nelem,msh.nface),msh.nedge),msh.npoin);
  int *vbuf = (int*)lorder;
  std::vector<int, PreAllocator<int> > vect(nm, 0, PreAllocator<int>(vbuf, nm));

  //// Use lorder after vector buffer as invord array
  //// There is enough space: we originally sized for 2*npoin*2 (as int64)
  //invord = &(vbuf[msh.nelem]);
//
  //for(int ii = 0; ii < msh.nelem; ii++){
  //  invord[vect[ii]] = ii;
  //}

 
  auto apply_sort = [&]<int itype>(std::integral_constant<int,itype>){
    int nentt = (itype == 1 ? msh.nedge : itype == 2 ? msh.nface : msh.nelem);
    intAr1 &ent2ref = (itype == 1 ? msh.edg2ref : itype == 2 ? msh.fac2ref : msh.tet2ref); 
    intAr2 &ent2tag = (itype == 1 ? msh.edg2tag : itype == 2 ? msh.fac2tag : msh.tet2tag);
    intAr2 &ent2poi = (itype == 1 ? msh.edg2poi : itype == 2 ? msh.fac2poi : msh.tet2poi);
    constexpr int nnode = (itype == 1 ? edgnpps[ideg] :
                           itype == 2 ? facnpps[ideg] :
                                        tetnpps[ideg] ); 

    std::iota(vect.begin(), vect.begin() + nentt, 0);
    if(itype == 1){
      std::sort(vect.begin(), vect.begin() + nentt, cmpe);
    }else if(itype == 2){
      std::sort(vect.begin(), vect.begin() + nentt, cmpf);
    }else{
      std::sort(vect.begin(), vect.begin() + nentt, cmpt);
    }

//    printf("Debug first 10 entries in vect \n");
//    for(int i = 0; i < 10; i++) printf("%d %d (%d %d ...)\n",i,vect[i],ent2poi[vect[i]][0]
//      ,ent2poi[vect[i]][1]);

    msh.tag[0]++;
    int lnode[nnode]; // Sized for maximum 

    for(int ient0 = 0; ient0 < nentt; ient0++){
      if(ent2tag(0,ient0) >= msh.tag[0]) continue;
      ent2tag(0,ient0) = msh.tag[0];
  
  
      int ient1 = vect[ient0];
      if(ient1 == ient0) continue;
  
      int iref0 = ent2ref[ient0];
      for(int inode = 0; inode < nnode; inode++) lnode[inode] = ent2poi(ient0,inode);
  
      ent2ref[ient0] = ent2ref[ient1];
      for(int inode = 0; inode < nnode; inode++) 
        ent2poi(ient0,inode) = ent2poi(ient1,inode);
//      printf("Debug start chain from iele0 %d iele1 = %d (%d %d ...) \n",ient0,ient1,
//        ent2poi(ient1,0),ent2poi(ient1,1));
//      if(ient0 == 1) printf("Found 1 in start\n");
      int nchain = 0;
      int ientp = ient0;
      do{
        assert(ent2tag(0,ient1) != msh.tag[0]);
        ent2tag(0,ient1) = msh.tag[0];
        int ient2 = vect[ient1];
//        printf("   next chain from iele1 %d iele2 = %d \n",ient1,ient2);
//        if(ient1 == 1) printf("Found 1 in mid\n");
  
        ent2ref[ient1] = ent2ref[ient2];
        for(int inode = 0; inode < nnode; inode++) 
          ent2poi(ient1,inode) = ent2poi(ient2,inode);
    
        ientp = ient1;
        ient1 = ient2;
        if(nchain > msh.nelem) METRIS_THROW_MSG(AlgoExcept(), "Infinite chain (Hilbert elt sort)")
      }while(ient1 != ient0);
    
//      if(ientp == 1) printf("Found 1 in end\n");
//      printf("   end chain at %d \n",ientp);
      ent2ref[ientp] = iref0;
      for(int inode = 0; inode < nnode; inode++) ent2poi(ientp,inode) = lnode[inode];
    }

//    printf("Debug first 10 entities:\n");
//    ent2poi.print(10);
  };


  apply_sort(std::integral_constant<int,1>{});
  apply_sort(std::integral_constant<int,2>{});
  apply_sort(std::integral_constant<int,3>{});


//  int nnode = tetnpps[msh.curdeg];
//  int lnode[tetnpps[msh.curdeg]];
//  msh.tag[0]++;
//  for(int iele0 = 0; iele0 < msh.nelem; iele0++){
//    if(msh.tet2tag(0,iele0) >= msh.tag[0]) continue;
//    msh.tet2tag(0,iele0) = msh.tag[0];
//
//
//    int iele1 = vect[iele0];
//    if(iele1 == iele0) continue;
//
//    int iref0 = msh.tet2ref[iele0];
//    for(int inode = 0; inode < nnode; inode++) lnode[inode] = msh.tet2poi(iele0,inode);
//
//    msh.tet2ref[iele0] = msh.tet2ref[iele1];
//    for(int inode = 0; inode < nnode; inode++) 
//      msh.tet2poi(iele0,inode) = msh.tet2poi(iele1,inode);
//
//    int nchain = 0;
//    int ielep = iele0;
//    do{
//      assert(msh.tet2tag(0,iele1) != msh.tag[0]);
//      msh.tet2tag(0,iele1) = msh.tag[0];
//      int iele2 = vect[iele1];
//
//      msh.tet2ref[iele1] = msh.tet2ref[iele2];
//      for(int inode = 0; inode < nnode; inode++) 
//        msh.tet2poi(iele1,inode) = msh.tet2poi(iele2,inode);
//
//      iele1 = iele2;
//      if(nchain > msh.nelem) METRIS_THROW_MSG(AlgoExcept(), "Infinite chain (Hilbert elt sort)")
//    }while(iele1 != iele0);
//
//    msh.tet2ref[ielep] = iref0;
//    for(int inode = 0; inode < nnode; inode++) msh.tet2poi(ielep,inode) = lnode[inode];
//  }



  if(ialloc) delete[] lorder;
  if(iallo2) delete[] invord;


  StopParallel(LPlibIdx);
  return;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void reoderHilbert< n > (MeshBase &msh);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

#endif