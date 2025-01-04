//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_LOW_SUBDIVIDE__
#define __METRIS_LOW_SUBDIVIDE__


#if 0

#include "common_includes.hxx"

namespace Metris{

#if 0
template <int ideg>
int subdivtet(int nelem, intAr2 &tet2poi, 
              int npoin, dblAr2 &coord, 
              int ielem, int iprt);
#endif

constexpr std::array<int,METRIS_MAX_DEG_JACOBIAN+1> nsubfac{[]() constexpr{
  std::array<int,METRIS_MAX_DEG_JACOBIAN+1> ret{};
  for(int ideg=0;ideg<METRIS_MAX_DEG_JACOBIAN+1;ideg++){
    ret[ideg] = 1;
    for(int k = 1; k < ideg; k++){
      ret[ideg] += 2*k+1;
    }
  }
  return ret;
}()};

// See subref
constexpr std::array<int,METRIS_MAX_DEG_JACOBIAN+1> msubtet{[]() constexpr{
  std::array<int,METRIS_MAX_DEG_JACOBIAN+1> ret{};
  for(int ideg=0;ideg<METRIS_MAX_DEG_JACOBIAN+1;ideg++){

    ret[0] = 1;

    for(int ideg = 1; ideg <= METRIS_MAX_DEG_JACOBIAN; ideg++){
      const int nneigh = 81;
      ret[ideg] = 0;
      for(int i = 0; i <= ideg; i++){
        for(int j = 0; j <= ideg - i ; j++){
          for(int k = 0; k <= ideg - i - j; k++){
            int l = ideg - i - j - k;
            int irnk = mul2nod(i,j,k,l);
            int irns[nneigh] = {0};
            for(int ii = 0; ii < nneigh; ii++) irns[ii] = -1;

            int inei = 0;
            int idx[4] = {0}; 
            for(int sg1 = -1; sg1 <= 1; sg1++){
              for(int sg2 = -1; sg2 <= 1; sg2++){
                for(int sg3 = -1; sg3 <= 1; sg3++){
                  int sg4 = - (sg1 + sg2 + sg3); 
                  ////if(iprt > 0) printf("Step %d %d %d %d \n",sg1,sg2,sg3,sg4);
                  if((sg4 == -1 || sg4 == 1 || sg4 == 0)
                    && (sg1 != 0 || sg2 != 0 || sg3 != 0 || sg4 != 0)){
                    idx[0] = i + sg1;
                    idx[1] = j + sg2;
                    idx[2] = k + sg3;
                    idx[3] = l + sg4;
                    ////if(iprt > 0) printf("Admissible idx = %d %d %d %d \n",idx[0],idx[1],idx[2],idx[3]);
                    if(idx[0] >= 0 && idx[0] <= ideg 
                    && idx[1] >= 0 && idx[1] <= ideg
                    && idx[2] >= 0 && idx[2] <= ideg
                    && idx[3] >= 0 && idx[3] <= ideg){
                      irns[inei] = mul2nod(idx[0],idx[1],idx[2],idx[3]);
                      inei++; 
                    }
                  }
                }
              }
            }

            for(int ii = 0; ii < inei; ii++){
              for(int jj = ii+1; jj < inei; jj++){
                for(int kk = jj+1; kk < inei; kk++){
                  ret[ideg]++;
                }
              }
            }


          }
        }
      }
    }
    
  }
  return ret;
}()};



/*
  - tdim: topological dimension
  - ref2poi: size neref x (edg|fac|tet)npps[ideg]
*/
struct subref_constructor{
//  constexpr subref_constructor() : fac(), tet(){ 
   constexpr subref_constructor() : fac(), tet(), ntet(){
    // --------------------------------------------------------------------------------------------------
    // TRIANGLES
    // --------------------------------------------------------------------------------------------------
    fac[0][0][0] = 0;
    fac[0][0][1] = 1;
    fac[0][0][2] = 2;

    
    for(int ideg = 1; ideg <= METRIS_MAX_DEG_JACOBIAN; ideg++){
      
      const int nneigh = 6;

      int ntri = 0;
      // Loop over possible multi-indices. 
      for(int i = 0; i <= ideg; i++){
        for(int j = 0; j <= ideg - i ; j++){
          int k = ideg - i - j;
          // Rank will be used to verify we do not duplicate triangles. 
          int irnk = mul2nod(i,j,k);

          ////if(iprt > 0) printf("ideg %d i j k %d %d %d \n",ideg,i,j,k);

          // The neighbours of the node are: 
          //  e^1 - e^3, e^2 - e^3, e^2 - e^1, e^3 - e^1, e^3 - e^2, e^1 - e^2
          int irns[6] = {-1, -1, -1, -1, -1, -1};

          // We can loop over these systematically by doing alternating incrementing first index, then second
          // e.g. (1,3) -> (2,3) -> (2,1) -> (3,1) -> (3,2) -> (1,2)

            //// We can loop over these systematically by picking 2 indices to increment / decrement 
            //int inei = 0;
            //for(int ii = 0; ii < 3; ii++){
            //  for(int jj = ii + 1; jj < 4; jj++){
            //    // Increment ii, decrement jj
            //    int idx[4] = {i,j,k,l};
            //    idx[ii] += 1;
            //    idx[jj] -= 1; 
            //    if(idx[ii] >= 0 && idx[ii] <= ideg 
            //    && idx[jj] >= 0 && idx[jj] <= ideg){
            //      irns[inei] = mul2nod(idx[0],idx[1],idx[2],idx[3]);
            //      inei++; 
            //    }
            //    idx[ii] -= 2;
            //    idx[jj] += 2;
            //    if(idx[ii] >= 0 && idx[ii] <= ideg 
            //    && idx[jj] >= 0 && idx[jj] <= ideg){
            //      irns[inei] = mul2nod(idx[0],idx[1],idx[2],idx[3]);
            //      inei++; 
            //    }
            //  }
            //}

          int ip = 0, jp = 0, kp = 0;

          // + e^1 - e^3
          ip = i + 1;
          jp = j;
          kp = k - 1;
          if(ip >= 0 && ip <= ideg && jp >= 0 && jp <= ideg && kp >= 0 && kp <= ideg){
            irns[0] = mul2nod(ip,jp,kp);
          }
          // + e^2 - e^3
          ip = i;
          jp = j + 1;
          kp = k - 1;
          if(ip >= 0 && ip <= ideg && jp >= 0 && jp <= ideg && kp >= 0 && kp <= ideg){
            irns[1] = mul2nod(ip,jp,kp);
          }
          // + e^2 - e^1
          ip = i - 1;
          jp = j + 1;
          kp = k;
          if(ip >= 0 && ip <= ideg && jp >= 0 && jp <= ideg && kp >= 0 && kp <= ideg){
            irns[2] = mul2nod(ip,jp,kp);
          }
          // + e^3 - e^1
          ip = i - 1;
          jp = j;
          kp = k + 1;
          if(ip >= 0 && ip <= ideg && jp >= 0 && jp <= ideg && kp >= 0 && kp <= ideg){
            irns[3] = mul2nod(ip,jp,kp);
          }
          // + e^3 - e^2
          ip = i;
          jp = j - 1;
          kp = k + 1;
          if(ip >= 0 && ip <= ideg && jp >= 0 && jp <= ideg && kp >= 0 && kp <= ideg){
            irns[4] = mul2nod(ip,jp,kp);
          }
          // + e^1 - e^2
          ip = i + 1;
          jp = j - 1;
          kp = k;
          if(ip >= 0 && ip <= ideg && jp >= 0 && jp <= ideg && kp >= 0 && kp <= ideg){
            irns[5] = mul2nod(ip,jp,kp);
          }
          
          ////if(iprt > 0) printf("neighbours %d %d %d %d %d %d \n",irns[0],irns[1],irns[2],irns[3],irns[4],irns[5]);

          // Go over all two consecutive ranks and form triangles with those yet to be seen. 
          for(int i = 0; i < 6; i++){
            int irn1 = irns[i];
            int irn2 = irns[(i+1)%nneigh];
            if(irn1 > irnk 
            && irn2 > irnk){
              // All three neighbours are defined, and their pair to current will only have been seen once. 
              fac[ideg][ntri][0] = irnk; 
              fac[ideg][ntri][1] = irn1; 
              fac[ideg][ntri][2] = irn2; 
              ntri++;
            }
          }

          //int idx0[3] = {i,j,k};
          //// We now create triangles with i+1, j-1, k; i-1, j+1, k; i+1, j, k-1 and so on. 
          //// excluding those with an index inferior to ours. Those have already been seen. 
          //for(int iidx1 = 0; iidx1 < 2; iidx1++){
          //  for(int jjdx1 = iidx1 + 1; jjdx1 < 3; jjdx1++){
          //    int idx1[3] = {i,j,k}; 
          //    // Loop over +- 1
          //    for(int sg = -1; sg <= 1; sg +=2){
          //      idx1[iidx1] = idx0[iidx1] + sg;
          //      idx1[jjdx1] = idx0[jjdx1] - sg;
          //      // Ensure multi-index valid 
          //      if(idx1[iidx1] < 0 || idx1[iidx1] > ideg) continue;
          //      if(idx1[jjdx1] < 0 || idx1[jjdx1] > ideg) continue;
          //      // Get neighbouring node rank to test for repeats
          //      int irn1 = mul2nod(idx1[0],idx1[1],idx1[2]);
          //      if(irn1 < irnk) continue;
          //      simp[ideg][ntri][0] = irnk;
          //      ntri++;

          //    }
          //  }
          //}

        }
      }
    }


    // --------------------------------------------------------------------------------------------------
    // TETRAHEDRA
    // --------------------------------------------------------------------------------------------------
    tet[0][0][0] = 0;
    tet[0][0][1] = 1;
    tet[0][0][2] = 2;
    tet[0][0][3] = 3;
    ntet[0] = 1;

    
    constexpr int mmhead = MAX(1,msubtet[METRIS_MAX_DEG_JACOBIAN]/3);
    constexpr int mmlist = 4*msubtet[METRIS_MAX_DEG_JACOBIAN];
    int head[mmhead] = {0};
    int list[mmlist][4] = {0};

    for(int ideg = 1; ideg <= METRIS_MAX_DEG_JACOBIAN; ideg++){
      //int iprt = 0;
      //if(ideg <= 2) iprt = 1;

      // A gnerous overestimation
      constexpr int nneigh = 81;

      ntet[ideg] = 0;

      int mhead = MAX(1,msubtet[ideg]/3), mlist = 4*msubtet[ideg];
      //printf("\n\n start ideg = %d mlist = %d mhead %d expect %d \n\n\n",ideg,mlist,mhead,msubtet[ideg]);
      int nlist = 0;
      // Store ORIENTED face idx then next in line

      for(int ii = 0; ii < mhead; ii++) head[ii] = -1;
      for(int ii = 0; ii < mlist; ii++) list[ii][3] = -1;

      // Loop over possible multi-indices. 
      for(int irnk = 0; irnk < tetnpps[ideg]; irnk++){

            //if(iprt > 0) printf("Start loop ideg %d irnk %d \n",ideg,irnk);
            int i = ordtet.s[ideg][irnk][0];
            int j = ordtet.s[ideg][irnk][1];
            int k = ordtet.s[ideg][irnk][2];
            int l = ordtet.s[ideg][irnk][3];

            // The neighbours of the node are: 
            //  e^1 - e^3, e^2 - e^3, e^2 - e^1, e^3 - e^1, e^3 - e^2, e^1 - e^2
            int irns[nneigh] = {-1};
            bool ttag[nneigh] = {false};
            //for(int ii = 0; ii < nneigh; ii++) irns[ii] = -1;
            //for(int ii = 0; ii < nneigh; ii++) ttag[ii] = false;

            // We actually accept e.g. e^1 - e^2 - e^3 + e^4
            // So any move by steps of 1 all summing to 0

            int inei = 0;
            int idx[4] = {0}; 
            for(int sg1 = -1; sg1 <= 1; sg1++){
              for(int sg2 = -1; sg2 <= 1; sg2++){
                for(int sg3 = -1; sg3 <= 1; sg3++){
                  int sg4 = - (sg1 + sg2 + sg3); 
                  ////if(iprt > 0) printf("Step %d %d %d %d \n",sg1,sg2,sg3,sg4);
                  if((sg4 == -1 || sg4 == 1 || sg4 == 0)
                    && (sg1 != 0 || sg2 != 0 || sg3 != 0 || sg4 != 0)){
                    idx[0] = i + sg1;
                    idx[1] = j + sg2;
                    idx[2] = k + sg3;
                    idx[3] = l + sg4;
                    ////if(iprt > 0) printf("Admissible idx = %d %d %d %d \n",idx[0],idx[1],idx[2],idx[3]);
                    if(idx[0] >= 0 && idx[0] <= ideg 
                    && idx[1] >= 0 && idx[1] <= ideg
                    && idx[2] >= 0 && idx[2] <= ideg
                    && idx[3] >= 0 && idx[3] <= ideg){
                      irns[inei] = mul2nod(idx[0],idx[1],idx[2],idx[3]);
                      inei++; 
                    }
                  }
                }
              }
            }

            //if(iprt > 0){
              //printf("Rank %d (F%d) nneigh %d \n",irnk,irnk+1,inei);
              //for(int i = 0; i < inei; i++){
              //  printf(" %d ",irns[i]);
              //}
              //printf("\n");
            //}

            // Pick any three neighbours and form tetrahedra if not done already. 
            for(int ii = 0; ii < inei; ii++){
              for(int jj = ii+1; jj < inei; jj++){
                for(int kk = jj+1; kk < inei; kk++){
                  bool maketet = true;

                  int irn1 = irns[ii];
                  int irn2 = irns[jj];
                  int irn3 = irns[kk]; 

                  //if(iprt > 0) printf("Pick neighbours %d %d %d \n",irn1,irn2,irn3);

                  // Check they are allowed to form a tetrahedron
                  for(int i = 0; i < 4; i++){
                    int a1 = ordtet.s[ideg][irn1][i] - ordtet.s[ideg][irn2][i];
                    if(a1 < 0) a1 *= -1;

                    int a2 = ordtet.s[ideg][irn1][i] - ordtet.s[ideg][irn3][i];
                    if(a2 < 0) a2 *= -1;

                    int a3 = ordtet.s[ideg][irn2][i] - ordtet.s[ideg][irn3][i];
                    if(a3 < 0) a3 *= -1;

                    if(a1 > 1){
                      maketet = false;
                      break; 
                    }else if(a2 > 1){
                      maketet = false;
                      break; 
                    }else if(a3 > 1){
                      maketet = false;
                      break; 
                    }
                  }

                  if(!maketet) continue;


                  // Orient the face for positive volume. 
                  int v1[3] = {ordtet.s[ideg][irn1][0] - ordtet.s[ideg][irnk][0],
                               ordtet.s[ideg][irn1][1] - ordtet.s[ideg][irnk][1],
                               ordtet.s[ideg][irn1][2] - ordtet.s[ideg][irnk][2]};
                  int v2[3] = {ordtet.s[ideg][irn2][0] - ordtet.s[ideg][irnk][0],
                               ordtet.s[ideg][irn2][1] - ordtet.s[ideg][irnk][1],
                               ordtet.s[ideg][irn2][2] - ordtet.s[ideg][irnk][2]};

                  int vec[3] = {v1[1]*v2[2] - v1[2]*v2[1],
                                v1[2]*v2[0] - v1[0]*v2[2],
                                v1[0]*v2[1] - v1[1]*v2[0]};

                  int v3[3] = {ordtet.s[ideg][irn3][0] - ordtet.s[ideg][irnk][0],
                               ordtet.s[ideg][irn3][1] - ordtet.s[ideg][irnk][1],
                               ordtet.s[ideg][irn3][2] - ordtet.s[ideg][irnk][2]};

                  int dot = vec[0]*v3[0] + vec[1]*v3[1] + vec[2]*v3[2];


                  if(dot == 0) continue;


                  //if(iprt > 0) printf("Volume sign %d \n",dot);


                  if(dot < 0){
                    int tmp = irn1;
                    irn1 = irn2;
                    irn2 = tmp;
                  }


                  // Check if any of the faces already exists in the hash table. 
                  // If so, abort. 
                  int irnks[4] = {irn1, irn2, irn3, irnk};
                  for(int ifa = 0; ifa < 4; ifa++){
                    int irl1 = irnks[(ifa+1)%4];
                    int irl2 = 0, irl3 = 0;

                    // Say tetrahedron 1 2 3 4 is positive. 
                    // Then            2 1 4 3 is also positive. 

                    // In K1, the first face is 2 3 4 
                    //    K2      2nd           4 3 2
                    // The faces run opposite one to the other. 
                    // This alternates 1st, 2nd, 3rd, 4th face. 
                    // This is because a 4-cycle is a -1 signature permutation. 
                    // Therefore the i-th face is made the (i+1)-th face of a positive tet
                    // by a 4-cycle (which would preserve the face order) FOLLOWED by a transposition
                    // which inverts the face order. 
                    if(ifa % 2 == 0){
                      irl2 = irnks[(ifa+2)%4];
                      irl3 = irnks[(ifa+3)%4];
                    }else{
                      irl3 = irnks[(ifa+2)%4];
                      irl2 = irnks[(ifa+3)%4];
                    }
                    int ikey = (irl1 + irl2 + irl3) % mhead;
                    int inext = head[ikey];
                    //if(iprt > 0) printf("  Pick new tet face %d %d %d ikey = %d -> %d \n",irl1,irl2,irl3,ikey,inext);
                    while(inext != -1){
                      int jrn1 = list[inext][0];
                      int jrn2 = list[inext][1];
                      int jrn3 = list[inext][2];

                      //if(iprt > 0) printf("    In list check face %d %d %d inext %d \n",jrn1,jrn2,jrn3,list[inext][3]);

                      // Compare up to cycle (preserves sign)
                      if(jrn1 == irl1 && jrn2 == irl2 && jrn3 && irl3 ||
                         jrn2 == irl1 && jrn3 == irl2 && jrn1 && irl3 ||
                         jrn3 == irl1 && jrn1 == irl2 && jrn2 && irl3 ){
                        //if(iprt > 0) printf("      -> Face found in hashtable!\n");
                        maketet = false;
                        break;
                      }
                      inext = list[inext][3];
                    }
                    if(!maketet) break; 
                  }

                  if(!maketet) continue;

                  // None of the faces have been found. 
                  //if(iprt > 0) printf("Create tetra %d %d %d %d N = %d \n",irnk,irn1,irn2,irn3,ntet[ideg]+1);
                  tet[ideg][ntet[ideg]][0] = irnk;
                  tet[ideg][ntet[ideg]][1] = irn1;
                  tet[ideg][ntet[ideg]][2] = irn2;
                  tet[ideg][ntet[ideg]][3] = irn3;
                  ntet[ideg]++; 

                  for(int ifa = 0; ifa < 4; ifa++){
                    int irl1 = irnks[(ifa+1)%4];
                    int irl2 = 0, irl3 = 0;
                    // See explanation above
                    if(ifa % 2 == 0){
                      irl2 = irnks[(ifa+2)%4];
                      irl3 = irnks[(ifa+3)%4];
                    }else{
                      irl3 = irnks[(ifa+2)%4];
                      irl2 = irnks[(ifa+3)%4];
                    }
                    int ikey = (irl1 + irl2 + irl3) % mhead;
                    if(head[ikey] == -1){
                      head[ikey] = nlist;
                      //if(iprt > 0) printf("  Add %d %d %d to head %d \n ",irl1,irl2,irl3,ikey);
                    }else{
                      int inext = head[ikey];
                      int iprev = inext; 
                      while(inext != -1){
                        iprev = inext;
                        inext = list[inext][3];
                      }
                      list[iprev][3] = nlist;
                      //if(iprt > 0) printf("  Add %d %d %d after list %d \n",irl1,irl2,irl3,iprev);
                    }

                    list[nlist][0] = irl1; 
                    list[nlist][1] = irl2; 
                    list[nlist][2] = irl3; 
                    list[nlist][3] = -1;
                    //if(ifa < 3) list[nlist][3] = nlist+1;
                    nlist++;

                  }
                }
              }
            }



            //// This time, we need to pick faces, i.e. triplets of neighbours. 
            //for(int ii = 0; ii < inei; ii++){
            //  for(int jj = ii+1; jj < inei; jj++){
            //    for(int kk = jj+1; kk < inei; kk++){
            //      int irn1 = irns[ii];
            //      int irn2 = irns[jj];
            //      int irn3 = irns[kk]; 
            //      bool dotet1 = irn1 > irnk && irn2 > irnk && irn3 > irnk;
            //      int igp1 = 0, igp2 = 0, igp3 = 0;
            //      for(int i = 0; i < 4; i++){
            //        igp1 += abs(ordtet.s[ideg][irn1][i] - ordtet.s[ideg][irn2][i]);
            //        igp2 += abs(ordtet.s[ideg][irn1][i] - ordtet.s[ideg][irn3][i]);
            //        igp3 += abs(ordtet.s[ideg][irn2][i] - ordtet.s[ideg][irn3][i]);
            //      }
            //      // If the vertices are not all neighbours, they cannot have created a tet. 
            //      bool dotet2 = igp1 > 2 || igp2 > 2 || igp3 > 2;
            //      dotet2 = false;
            //      bool dotet = (dotet1 || dotet2) && (ttag[ii]+ttag[jj]+ttag[kk] < 2);
            //      if(dotet){
            //        tet[ideg][ntet][0] = irnk;
            //        tet[ideg][ntet][1] = irn1;
            //        tet[ideg][ntet][2] = irn2;
            //        tet[ideg][ntet][3] = irn3;
            //        ttag[ii] = true;
            //        ttag[jj] = true;
            //        ttag[kk] = true;
            //        ntet++;
            //        if(iprt == 1){
            //          printf("ntet = F%d (C%d) done with irnk = %d (%d%d%d%d)\n",ntet,ntet-1,irnk,i,j,k,l);
            //          printf("ii, jj, kk  %d %d %d irn1 irn2 irn3 %d %d %d\n",ii,jj,kk,irn1,irn2,irn3);
            //          printf("%d%d%d%d, %d%d%d%d, %d%d%d%d\n",
            //             ordtet.s[ideg][irn1][0],ordtet.s[ideg][irn1][1]
            //            ,ordtet.s[ideg][irn1][2],ordtet.s[ideg][irn1][3]
            //            ,ordtet.s[ideg][irn2][0],ordtet.s[ideg][irn2][1]
            //            ,ordtet.s[ideg][irn2][2],ordtet.s[ideg][irn2][3]
            //            ,ordtet.s[ideg][irn3][0],ordtet.s[ideg][irn3][1]
            //            ,ordtet.s[ideg][irn3][2],ordtet.s[ideg][irn3][3]
            //            );
            //        }
            //      }
            //    }
            //  }
            //}


      } 

    } // End for ideg 
  } // End const
  int fac[1+METRIS_MAX_DEG_JACOBIAN][nsubfac[METRIS_MAX_DEG_JACOBIAN]][3];
  int tet[1+METRIS_MAX_DEG_JACOBIAN][msubtet[METRIS_MAX_DEG_JACOBIAN]][4];
  int ntet[1+METRIS_MAX_DEG_JACOBIAN];
};

constexpr subref_constructor subref;





} // End namespace

#endif // if 0

#endif