//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "aux_utils.hxx"

#include "ho_constants.hxx"
#include "aux_exceptions.hxx"

#include "../SANS/Surreal/SurrealS.h"
#include "../SANS/tools/minmax.h"

#include <boost/preprocessor/iteration/local.hpp>

#include <cmath>
#include <vector>
#include <numeric>


namespace Metris{

//template <typename ftype>
//void fini_diff(int ndim, std::function< ftype(double*) >& fun,
//               double *x0){
//
//}

//     iinter = 1 : linear progression between min and max value
//     iinter = 2 ; geometric
//     nlist: input size
//     rlist: values
//     llist: indices ; for instance, if rlist pertains to edges, perhaps
//                      min and max should print llist(imin) and llist(imax)
//                      instead
//     if llist = NULL, simply index 0 ... nlist - 1
// T is any hash table type with an iterator that has a ->first and ->second double field
void generic_hist(double *rlist,int nlist, int *llist,double vlo,double vhi,
                  const char *vname,const char *title,int iinter){
  int nslot;
  const int mslot = 20;
  int cslot[mslot];
  double vslot[mslot];
  double vmin,vmax,vavg,dslot,v0;
  int imin,imax,nval,nlo,nhi;

  if(nlist <= 0) METRIS_THROW_MSG(WArgExcept(),"Empty list in histo");

  vmin = rlist[0];
  vmax = rlist[0];
  imin = 0;
  imax = 0;
  vavg = 0.0;
  nval = 0;
  nlo = 0;
  nhi = 0;
  for(int ilist = 1; ilist < nlist; ilist++){
    double rlcur = rlist[ilist];
    if(rlcur < vlo) nlo++;
    else if(rlcur > vhi) nhi++;
    else nval++;
    if(rlcur < vmin){
      imin = ilist;
      vmin = rlcur;
    }
    if(rlcur > vmax){
      imax = ilist;
      vmax = rlcur;
    }
    vavg += rlcur;
  }
  vavg /= nlist;

  nslot = nval < mslot ? nval : mslot;
  if(iinter == 1){
    dslot = ( (vmax < vhi ? vmax : vhi) - (vmin > vlo ? vmin : vlo) ) / nslot;
  }else{
    dslot = log( (vmax < vhi ? vmax : vhi) - (vmin > vlo ? vmin : vlo) ) / nslot;
  }

  for(int i = 0; i < nslot; i++){
    if(iinter == 1){
      vslot[i] = (vmin > vlo ? vmin : vlo) + i*dslot;
    }else{
      vslot[i] = exp(log((vmin > vlo ? vmin : vlo))+i*dslot);
    }
    cslot[i] = 0;
  }
  if(nhi > 0) vslot[nslot-1] = vhi;
  for(int ilist = 0; ilist < nlist; ilist++){
    double rlcur = rlist[ilist];
    if(rlcur < vlo || rlcur > vhi) continue;
    int iskp = 0;
    for(int j=0;j<nslot;j++){
      if(rlcur > vslot[j]) continue;
      cslot[j] ++;
      iskp = 1;
      break;
    }
    if(iskp == 0){
      cslot[nslot-1] ++; 
    }
  }
  printf("-- %s histogram\n",title);
  if(llist != NULL){
    imin = llist[imin];
    imax = llist[imax];
  }
  printf("  - min = %9.2e at %d\n",vmin,imin);
  printf("  - avg = %9.2e \n",vavg);
  printf("  - max = %9.2e at %d\n",vmax,imax);
  //if(nval < nlist) printf(" - outside bounds n = %d\n",nlist-nval);

    if(iinter == 1){
      if(nlo > 0){
        printf("                  %s < %9.3f %15d       %6.2f %%\n",vname,vlo,nlo,(100.0*nlo)/nlist);
      }
      if(cslot[0] > 0){
        v0 = vmin;
        if(nlo > 0) v0 = vlo;
        printf("      %9.3f < %s < %9.3f %15d       %6.2f %%\n",v0,vname,vslot[0],cslot[0],(100.0*cslot[0])/nlist);
      }
      for(int i = 0; i < nslot-2; i++){
        if(cslot[i+1]>0)
          printf("      %9.3f < %s < %9.3f %15d       %6.2f %%\n",vslot[i],vname,vslot[i+1],cslot[i+1],(100.0*cslot[i+1])/nlist);
      }
      if(nhi > 0){
        printf("      %9.3f < %s            %d       %6.2f %%\n",vhi,vname,nhi,(100.0*nhi)/nlist);
      }
    }else{
      if(nlo > 0){
        printf("                  %s < %9.2e %15d       %6.2f %%\n",vname,vlo,nlo,(100.0*nlo)/nlist);
      }
      if(cslot[0] > 0){
        v0 = vmin;
        if(nlo > 0) v0 = vlo;
        printf("      %9.2e < %s < %9.2e %15d       %6.2f %%\n",v0,vname,vslot[0],cslot[0],(100.0*cslot[0])/nlist);
      }
      for(int i = 0; i < nslot-2; i++){
        if(cslot[i+1]>0)
          printf("      %9.2e < %s < %9.2e %15d       %6.2f %%\n",vslot[i],vname,vslot[i+1],cslot[i+1],(100.0*cslot[i+1])/nlist);
      }
      if(nhi > 0){
        printf("      %9.2e < %s            %d       %6.2f %%\n",vhi,vname,nhi,(100.0*nhi)/nlist);
      }
    }
}




// Specialization for termination
//template<>  double idpow<0>(double x){return 1.0;}
//template<>  constexpr int ifact<0>(){return 1;}

/**
 * C++ version 0.4 char* style "itoa":
 * Written by Lukás Chmela
 * Released under GPLv3.

 */
char* itoa(int value, char* result, int base) {
  // check that the base if valid
  if (base < 2 || base > 36) { *result = '\0'; return result; }

  char* ptr = result, *ptr1 = result, tmp_char;
  int tmp_value;

  do {
    tmp_value = value;
    value /= base;
    *ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
  } while ( value );

  // Apply negative sign
  if (tmp_value < 0) *ptr++ = '-';
  *ptr-- = '\0';
  while(ptr1 < ptr) {
    tmp_char = *ptr;
    *ptr--= *ptr1;
    *ptr1++ = tmp_char;
  }
  return result;
}

void wait(){
  fflush(stdout);
  char tmp;
  printf("Press key to continue ");
  scanf("%c",&tmp);
}

// Gets 2,3 ints, returns sorted tuple (for hash keys)
std::tuple<int,int,int> stup3(int i1,int i2,int i3){
  int key_[3] = {i1,i2,i3};
  sortupto8_dec(key_,3);
  return std::tuple<int,int,int>({key_[0],key_[1],key_[2]});
}
std::tuple<int,int>     stup2(int i1,int i2){
  int key_[2] = {i1,i2};
  sortupto8_dec(key_,2);
  return std::tuple<int,int>({key_[0],key_[1]});
}
//template<int n>
//typename std::conditional<n == 2, std::tuple<int,int>, std::tuple<int,int,int>>::type 
//stupn(const int *ii){
//  static_assert(n == 2 || n == 3);
//  if constexpr(n == 2){
//    int key_[2] = {ii[0],ii[1]};
//    sortupto8_dec(key_,n);
//    return std::tuple<int,int>({key_[0],key_[1]});
//  }else{
//    int key_[3] = {ii[0],ii[1],ii[2]};
//    sortupto8_dec(key_,n);
//    return std::tuple<int,int,int>({key_[0],key_[1],key_[2]});
//  }
//}
//template std::tuple<int,int> stupn<2>(const int *ii);
//template std::tuple<int,int,int> stupn<3>(const int *ii);


double linearRegression(int n, double *x, double *y) {
  std::vector<double> v_x(n), v_y(n);
  for(int i = 0; i < n; i++){
    v_x[i] = x[i];
    v_y[i] = y[i];
  }
  const auto s_x  = std::accumulate(v_x.begin(), v_x.end(), 0.0);
  const auto s_y  = std::accumulate(v_y.begin(), v_y.end(), 0.0);
  const auto s_xx = std::inner_product(v_x.begin(), v_x.end(), v_x.begin(), 0.0);
  const auto s_xy = std::inner_product(v_x.begin(), v_x.end(), v_y.begin(), 0.0);
  const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  return a;
}



template <int ideg, int idim>
void gen_ordering_Vizir(int *ord){
  int n = 0;
//    std::cout<<"Enter vizir deg, dim "<<deg<<" "<<dim<<std::endl;
  if constexpr(idim==1){
    if(ideg == 0){
      ord[0] = 0;
      ord[1] = 0;
      return;
    }
    ord[0] = ideg;
    ord[1] = 0;
    ord[2] = 0;
    ord[3] = ideg;
    n = 2;
    for(int i1 = 1; i1 < ideg; i1++){
      ord[2*n+0]=ideg-i1;
      ord[2*n+1]=i1;
    }
  }else if(idim==2){
    if(ideg == 0){
      ord[0] = 0;
      ord[1] = 0;
      ord[2] = 0;
      return;
    }
    ord[0] = ideg;
    ord[1] = 0;
    ord[2] = 0;
    ord[3] = 0;
    ord[4] = ideg;
    ord[5] = 0;
    ord[6] = 0;
    ord[7] = 0;
    ord[8] = ideg;
    n = 3;
    for(int i1 = 1;i1 < ideg;i1++){
      ord[3*n+0] = ideg-i1;
      ord[3*n+1] =     i1;
      ord[3*n+2] = 0;
      n++;
    }
    for(int i1 = 1;i1 < ideg;i1++){
      ord[3*n+0] = 0;
      ord[3*n+1] = ideg-i1;
      ord[3*n+2] =     i1;
      n++;
    }
    for(int i1 = 1;i1 < ideg;i1++){
      ord[3*n+0] =     i1;
      ord[3*n+1] = 0;
      ord[3*n+2] = ideg-i1;
      n++;
    }
    if constexpr(ideg >= 3){
      constexpr int nppSub = facnpps[ideg - 3];
      int subOrd[nppSub];
      gen_ordering_Vizir<ideg-3,idim>(subOrd);
//            std::cout<<"first entry in subOrd"<<subOrd[0]<<" "<<subOrd[1]<<" "<<subOrd[2]<<" "<<std::endl;
//            std::cout<<"subOrd size "<<nppSub<<std::endl;
      for(int i = 0;i < nppSub;i++){
        ord[3*n+0] =1+subOrd[3*i+0];
        ord[3*n+1] =1+subOrd[3*i+1];
        ord[3*n+2] =1+subOrd[3*i+2];
        n++;
      }
    }
  }else if(idim==3){
    if(ideg == 0){
      ord[0] = 0;
      ord[1] = 0;
      ord[2] = 0;
      ord[3] = 0;
      return;
    }
    ord[0 ] = ideg;
    ord[1 ] = 0;
    ord[2 ] = 0;
    ord[3 ] = 0;

    ord[4 ] = 0;
    ord[5 ] = ideg;
    ord[6 ] = 0;
    ord[7 ] = 0;

    ord[8 ] = 0;
    ord[9 ] = 0;
    ord[10] = ideg;
    ord[11] = 0;

    ord[12] = 0;
    ord[13] = 0;
    ord[14] = 0;
    ord[15] = ideg;

    n = 4;

    for(int i = 1; i < ideg; i++){
      ord[4*n+0]=ideg-i ;
      ord[4*n+1]=i     ; 
      ord[4*n+2]=0     ;
      ord[4*n+3]=0     ;
      n++ ;
    }
    
    for(int i = 1; i < ideg; i++){
      ord[4*n+0]=0     ;
      ord[4*n+1]=ideg-i   ; 
      ord[4*n+2]=i     ;
      ord[4*n+3]=0     ;
      n++ ;
    }
    
    for(int i = 1; i < ideg; i++){
      ord[4*n+0]=i     ;
      ord[4*n+1]=0     ; 
      ord[4*n+2]=ideg-i  ;
      ord[4*n+3]=0     ;
      n++ ;
    }
    
    for(int i = 1; i < ideg; i++){
      ord[4*n+0]=ideg-i ;
      ord[4*n+1]=0     ; 
      ord[4*n+2]=0     ;
      ord[4*n+3]=i     ;
      n++ ;
    }
    
    for(int i = 1; i < ideg; i++){
      ord[4*n+0]=0     ;
      ord[4*n+1]=ideg-i   ; 
      ord[4*n+2]=0     ;
      ord[4*n+3]=i     ;
      n++ ;
    }
    
    for(int i = 1; i < ideg; i++){
      ord[4*n+0]=0     ;
      ord[4*n+1]=0     ; 
      ord[4*n+2]=ideg-i   ;
      ord[4*n+3]=i     ;
      n++ ;
    }

    if constexpr(ideg >= 3){

      constexpr int nppSub = facnpps[ideg - 3];
      int subOrd[nppSub];
      gen_ordering_Vizir<ideg-3,2>(subOrd);

      for (int i = 0 ; i < nppSub ; i++){
        ord[4*n+0] = 1+subOrd[3*i+0];
        ord[4*n+1] = 1+subOrd[3*i+1]; 
        ord[4*n+2] = 1+subOrd[3*i+2]; 
        ord[4*n+3] = 0;
        n++;
      }
      
      for (int i = 0 ; i < nppSub ; i++){
        ord[4*n+0] = 1+subOrd[3*i+0];
        ord[4*n+1] = 1+subOrd[3*i+1]; 
        ord[4*n+2] = 0              ; 
        ord[4*n+3] = 1+subOrd[3*i+2];
        n++;
      }
      
      for (int i = 0 ; i < nppSub ; i++){
        ord[4*n+0] = 0              ;
        ord[4*n+1] = 1+subOrd[3*i+0]; 
        ord[4*n+2] = 1+subOrd[3*i+1]; 
        ord[4*n+3] = 1+subOrd[3*i+2];
        n++;
      }
      
      for (int i = 0 ; i < nppSub ; i++){
        ord[4*n+0] = 1+subOrd[3*i+1];
        ord[4*n+1] = 0              ; 
        ord[4*n+2] = 1+subOrd[3*i+0]; 
        ord[4*n+3] = 1+subOrd[3*i+2];
        n++;
      }
    }
    if constexpr(ideg >= 4){
      constexpr int nppSub = tetnpps[ideg - 4];
      int subOrd[nppSub];
      gen_ordering_Vizir<ideg-4,idim>(subOrd);
      for (int i = 0 ; i < nppSub ; i++){
        ord[4*n+0]=1+subOrd[4*i+0];
        ord[4*n+1]=1+subOrd[4*i+1];
        ord[4*n+2]=1+subOrd[4*i+2];
        ord[4*n+3]=1+subOrd[4*i+3];
        n++;
      }
    }
  }
}

#define BOOST_PP_LOCAL_MACRO(n)\
template void gen_ordering_Vizir< n , 1 >(int *ord);\
template void gen_ordering_Vizir< n , 2 >(int *ord);\
template void gen_ordering_Vizir< n , 3 >(int *ord);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

///*
//This is a modified version of Vizir's recursive ordering based on P.-L.'s books
//This one is coherent with neighbour ordering. 
//Indeed, the i-th neighbour is opposite the i-th vertex. 
//Naturally, the i-th face of a tet, or the i-th edge of a triangle, should not contain the i-th vertex...
//*/
//void gen_ordering_Vizir_coherentHO(int deg,int dim, int *ord){
//    int n = 0;
////    std::cout<<"Enter vizir deg, dim "<<deg<<" "<<dim<<std::endl;
//    if(dim==1){
//        if(deg == 0){
//            ord[0] = 0;
//            ord[1] = 0;
//            return;
//        }
//        ord[0] = deg;
//        ord[1] = 0;
//        ord[2] = 0;
//        ord[3] = deg;
//        n = 2;
//        for(int i1 = 1; i1 < deg; i1++){
//            ord[2*n+0]=deg-i1;
//            ord[2*n+1]=i1;
//        }
//    }else if(dim==2){
//        if(deg == 0){
//            ord[0] = 0;
//            ord[1] = 0;
//            ord[2] = 0;
//            return;
//        }
//        ord[0] = deg;
//        ord[1] = 0;
//        ord[2] = 0;
//        ord[3] = 0;
//        ord[4] = deg;
//        ord[5] = 0;
//        ord[6] = 0;
//        ord[7] = 0;
//        ord[8] = deg;
//        n = 3;
//        for(int i1 = 1;i1 < deg;i1++){
//            ord[3*n+0] = 0;
//            ord[3*n+1] = deg-i1;
//            ord[3*n+2] =     i1;
//            n++;
//        }
//        for(int i1 = 1;i1 < deg;i1++){
//            ord[3*n+0] =     i1;
//            ord[3*n+1] = 0;
//            ord[3*n+2] = deg-i1;
//            n++;
//        }
//        for(int i1 = 1;i1 < deg;i1++){
//            ord[3*n+0] = deg-i1;
//            ord[3*n+1] =     i1;
//            ord[3*n+2] = 0;
//            n++;
//        }
//        if(deg >= 3){
//            int nppSub = facnpps[deg - 3];
//            int subOrd[nppSub];
//            gen_ordering_Vizir_coherentHO(deg-3,dim,subOrd);
////            std::cout<<"first entry in subOrd"<<subOrd[0]<<" "<<subOrd[1]<<" "<<subOrd[2]<<" "<<std::endl;
////            std::cout<<"subOrd size "<<nppSub<<std::endl;
//            for(int i = 0;i < nppSub;i++){
//                ord[3*n+0] =1+subOrd[3*i+0];
//                ord[3*n+1] =1+subOrd[3*i+1];
//                ord[3*n+2] =1+subOrd[3*i+2];
//                n++;
//            }
//        }
//    }else if(dim==3){
//        if(deg == 0){
//            ord[0] = 0;
//            ord[1] = 0;
//            ord[2] = 0;
//            ord[3] = 0;
//            return;
//        }
//        ord[0 ] = deg;
//        ord[1 ] = 0;
//        ord[2 ] = 0;
//        ord[3 ] = 0;
//
//        ord[4 ] = 0;
//        ord[5 ] = deg;
//        ord[6 ] = 0;
//        ord[7 ] = 0;
//
//        ord[8 ] = 0;
//        ord[9 ] = 0;
//        ord[10] = deg;
//        ord[11] = 0;
//
//        ord[12] = 0;
//        ord[13] = 0;
//        ord[14] = 0;
//        ord[15] = deg;
//
//        n = 4;
//
//        for(int i = 1; i < deg; i++){
//            ord[4*n+0]=deg-i ;
//            ord[4*n+1]=i     ; 
//            ord[4*n+2]=0     ;
//            ord[4*n+3]=0     ;
//            n++ ;
//        }
//  
//        for(int i = 1; i < deg; i++){
//            ord[4*n+0]=0     ;
//            ord[4*n+1]=deg-i   ; 
//            ord[4*n+2]=i     ;
//            ord[4*n+3]=0     ;
//            n++ ;
//        }
//  
//        for(int i = 1; i < deg; i++){
//            ord[4*n+0]=i     ;
//            ord[4*n+1]=0     ; 
//            ord[4*n+2]=deg-i  ;
//            ord[4*n+3]=0     ;
//            n++ ;
//        }
//  
//        for(int i = 1; i < deg; i++){
//            ord[4*n+0]=deg-i ;
//            ord[4*n+1]=0     ; 
//            ord[4*n+2]=0     ;
//            ord[4*n+3]=i     ;
//            n++ ;
//        }
//  
//        for(int i = 1; i < deg; i++){
//            ord[4*n+0]=0     ;
//            ord[4*n+1]=deg-i   ; 
//            ord[4*n+2]=0     ;
//            ord[4*n+3]=i     ;
//            n++ ;
//        }
//  
//        for(int i = 1; i < deg; i++){
//            ord[4*n+0]=0     ;
//            ord[4*n+1]=0     ; 
//            ord[4*n+2]=deg-i   ;
//            ord[4*n+3]=i     ;
//            n++ ;
//        }
//
//        if(deg >= 3){
//
//            int nppSub = facnpps[deg - 3];
//            int subOrd[nppSub];
//            gen_ordering_Vizir_coherentHO(deg-3,2,subOrd);
//  
//            for (int i = 0 ; i < nppSub ; i++){
//                ord[4*n+0] = 0              ;
//                ord[4*n+1] = 1+subOrd[3*i+0]; 
//                ord[4*n+2] = 1+subOrd[3*i+1]; 
//                ord[4*n+3] = 1+subOrd[3*i+2];
//                n++;
//            }
//
//            for (int i = 0 ; i < nppSub ; i++){
//                ord[4*n+0] = 1+subOrd[3*i+1];
//                ord[4*n+1] = 0              ; 
//                ord[4*n+2] = 1+subOrd[3*i+0]; 
//                ord[4*n+3] = 1+subOrd[3*i+2];
//                n++;
//            }
//
//            for (int i = 0 ; i < nppSub ; i++){
//                ord[4*n+0] = 1+subOrd[3*i+0];
//                ord[4*n+1] = 1+subOrd[3*i+1]; 
//                ord[4*n+2] = 0              ; 
//                ord[4*n+3] = 1+subOrd[3*i+2];
//                n++;
//            }
//            
//            for (int i = 0 ; i < nppSub ; i++){
//                ord[4*n+0] = 1+subOrd[3*i+0];
//                ord[4*n+1] = 1+subOrd[3*i+1]; 
//                ord[4*n+2] = 1+subOrd[3*i+2]; 
//                ord[4*n+3] = 0;
//                n++;
//            }
//  
//  
//        }
//        if(deg >= 4){
//            int nppSub = tetnpps[deg - 4];
//            int subOrd[nppSub];
//            gen_ordering_Vizir_coherentHO(deg-4,dim,subOrd);
//            for (int i = 0 ; i < nppSub ; i++){
//                ord[4*n+0]=1+subOrd[4*i+0];
//                ord[4*n+1]=1+subOrd[4*i+1];
//                ord[4*n+2]=1+subOrd[4*i+2];
//                ord[4*n+3]=1+subOrd[4*i+3];
//                n++;
//            }
//        }
//    }
//}
//

// Optimal sorting networks for up to 8 elements
// Sort as decreasing
template<typename T>
void sortupto8_dec(T *tab,int n){

  if(n == 8){
    if(tab[1-1]<tab[2-1])  swi(tab[1-1],tab[2-1]);
    if(tab[3-1]<tab[4-1])  swi(tab[3-1],tab[4-1]);
    if(tab[5-1]<tab[6-1])  swi(tab[5-1],tab[6-1]);
    if(tab[7-1]<tab[8-1])  swi(tab[7-1],tab[8-1]);
    
    if(tab[1-1]<tab[3-1])  swi(tab[1-1],tab[3-1]);
    if(tab[2-1]<tab[4-1])  swi(tab[2-1],tab[4-1]) ;
    if(tab[5-1]<tab[7-1])  swi(tab[5-1],tab[7-1]) ;
    if(tab[6-1]<tab[8-1])  swi(tab[6-1],tab[8-1]) ;
    
    if(tab[2-1]<tab[3-1])  swi(tab[2-1],tab[3-1]);
    if(tab[6-1]<tab[7-1])  swi(tab[6-1],tab[7-1]);
    if(tab[1-1]<tab[5-1])  swi(tab[1-1],tab[5-1]);
    if(tab[4-1]<tab[8-1])  swi(tab[4-1],tab[8-1]);
    
    if(tab[2-1]<tab[6-1])  swi(tab[2-1],tab[6-1]);
    if(tab[3-1]<tab[7-1])  swi(tab[3-1],tab[7-1]);
    
    if(tab[2-1]<tab[5-1])  swi(tab[2-1],tab[5-1]);
    if(tab[4-1]<tab[7-1])  swi(tab[4-1],tab[7-1]);
    
    if(tab[3-1]<tab[5-1])  swi(tab[3-1],tab[5-1]);
    if(tab[4-1]<tab[6-1])  swi(tab[4-1],tab[6-1]);
    
    if(tab[4-1]<tab[5-1])  swi(tab[4-1],tab[5-1]);

  }else if(n==7){

    if(tab[2-1]<tab[3-1] )  swi(tab[2-1],tab[3-1]);
    if(tab[4-1]<tab[5-1] )  swi(tab[4-1],tab[5-1]);
    if(tab[6-1]<tab[7-1] )  swi(tab[6-1],tab[7-1]);
    
    if(tab[1-1]<tab[3-1] )  swi(tab[1-1],tab[3-1]);
    if(tab[4-1]<tab[6-1] )  swi(tab[4-1],tab[6-1]);
    if(tab[5-1]<tab[7-1] )  swi(tab[5-1],tab[7-1]);
    
    if(tab[1-1]<tab[2-1] )  swi(tab[1-1],tab[2-1]);
    if(tab[5-1]<tab[6-1] )  swi(tab[5-1],tab[6-1]);
    if(tab[3-1]<tab[7-1] )  swi(tab[3-1],tab[7-1]);
    
    if(tab[1-1]<tab[5-1] )  swi(tab[1-1],tab[5-1]);
    if(tab[2-1]<tab[6-1] )  swi(tab[2-1],tab[6-1]);
    
    if(tab[1-1]<tab[4-1] )  swi(tab[1-1],tab[4-1]);
    if(tab[3-1]<tab[6-1] )  swi(tab[3-1],tab[6-1]);
    
    if(tab[2-1]<tab[4-1] )  swi(tab[2-1],tab[4-1]);
    if(tab[3-1]<tab[5-1] )  swi(tab[3-1],tab[5-1]);
    
    if(tab[3-1]<tab[4-1] )  swi(tab[3-1],tab[4-1]);
    
    
  }else if(n==6){
    
    
    if(tab[2-1]<tab[3-1] )  swi(tab[2-1],tab[3-1]);
    if(tab[5-1]<tab[6-1] )  swi(tab[5-1],tab[6-1]);
    
    if(tab[1-1]<tab[3-1] )  swi(tab[1-1],tab[3-1]);
    if(tab[4-1]<tab[6-1] )  swi(tab[4-1],tab[6-1]);
    
    if(tab[1-1]<tab[2-1] )  swi(tab[1-1],tab[2-1]);
    if(tab[4-1]<tab[5-1] )  swi(tab[4-1],tab[5-1]);
    if(tab[3-1]<tab[6-1] )  swi(tab[3-1],tab[6-1]);
    
    if(tab[1-1]<tab[4-1] )  swi(tab[1-1],tab[4-1]);
    if(tab[2-1]<tab[5-1] )  swi(tab[2-1],tab[5-1]);
    
    if(tab[3-1]<tab[5-1] )  swi(tab[3-1],tab[5-1]);
    if(tab[2-1]<tab[4-1] )  swi(tab[2-1],tab[4-1]);
    
    if(tab[3-1]<tab[4-1] )  swi(tab[3-1],tab[4-1]);
    
    
  }else if(n==5){
    
   
    if(tab[1-1]<tab[2-1] )  swi(tab[1-1],tab[2-1]);
    if(tab[4-1]<tab[5-1] )  swi(tab[4-1],tab[5-1]);
    
    if(tab[3-1]<tab[5-1] )  swi(tab[3-1],tab[5-1]);
    
    if(tab[3-1]<tab[4-1] )  swi(tab[3-1],tab[4-1]);
    if(tab[2-1]<tab[5-1] )  swi(tab[2-1],tab[5-1]);
    
    if(tab[1-1]<tab[4-1] )  swi(tab[1-1],tab[4-1]);
    
    if(tab[1-1]<tab[3-1] )  swi(tab[1-1],tab[3-1]);
    if(tab[2-1]<tab[4-1] )  swi(tab[2-1],tab[4-1]);
    
    if(tab[2-1]<tab[3-1] )  swi(tab[2-1],tab[3-1]);
    
    
  }else if(n==4){
    
   
    if(tab[1-1]<tab[2-1] )  swi(tab[1-1],tab[2-1]);
    if(tab[3-1]<tab[4-1] )  swi(tab[3-1],tab[4-1]);
    
    if(tab[1-1]<tab[3-1] )  swi(tab[1-1],tab[3-1]);
    if(tab[2-1]<tab[4-1] )  swi(tab[2-1],tab[4-1]);
    
    if(tab[2-1]<tab[3-1] )  swi(tab[2-1],tab[3-1]);
    
    
  }else if(n==3){
    
   
    if(tab[2-1]<tab[3-1] )  swi(tab[2-1],tab[3-1]);
    
    if(tab[1-1]<tab[3-1] )  swi(tab[1-1],tab[3-1]);
    
    if(tab[1-1]<tab[2-1] )  swi(tab[1-1],tab[2-1]);
    

  }else if(n==2){

    if(tab[1-1]<tab[2-1] )  swi(tab[1-1],tab[2-1]);

  }else if(n<1){
    METRIS_THROW(WArgExcept());
  }
}

// Optimal sorting networks for up to 8 elements
// Sort as decreasing
template<typename T>
void sortupto8_dec(T *tab, int *idx, int n){

  if(n == 8){
    if(tab[idx[1-1]]<tab[idx[2-1]])  swi(idx[1-1],idx[2-1]);
    if(tab[idx[3-1]]<tab[idx[4-1]])  swi(idx[3-1],idx[4-1]);
    if(tab[idx[5-1]]<tab[idx[6-1]])  swi(idx[5-1],idx[6-1]);
    if(tab[idx[7-1]]<tab[idx[8-1]])  swi(idx[7-1],idx[8-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]])  swi(idx[1-1],idx[3-1]);
    if(tab[idx[2-1]]<tab[idx[4-1]])  swi(idx[2-1],idx[4-1]) ;
    if(tab[idx[5-1]]<tab[idx[7-1]])  swi(idx[5-1],idx[7-1]) ;
    if(tab[idx[6-1]]<tab[idx[8-1]])  swi(idx[6-1],idx[8-1]) ;
    
    if(tab[idx[2-1]]<tab[idx[3-1]])  swi(idx[2-1],idx[3-1]);
    if(tab[idx[6-1]]<tab[idx[7-1]])  swi(idx[6-1],idx[7-1]);
    if(tab[idx[1-1]]<tab[idx[5-1]])  swi(idx[1-1],idx[5-1]);
    if(tab[idx[4-1]]<tab[idx[8-1]])  swi(idx[4-1],idx[8-1]);
    
    if(tab[idx[2-1]]<tab[idx[6-1]])  swi(idx[2-1],idx[6-1]);
    if(tab[idx[3-1]]<tab[idx[7-1]])  swi(idx[3-1],idx[7-1]);
    
    if(tab[idx[2-1]]<tab[idx[5-1]])  swi(idx[2-1],idx[5-1]);
    if(tab[idx[4-1]]<tab[idx[7-1]])  swi(idx[4-1],idx[7-1]);
    
    if(tab[idx[3-1]]<tab[idx[5-1]])  swi(idx[3-1],idx[5-1]);
    if(tab[idx[4-1]]<tab[idx[6-1]])  swi(idx[4-1],idx[6-1]);
    
    if(tab[idx[4-1]]<tab[idx[5-1]])  swi(idx[4-1],idx[5-1]);

  }else if(n==7){

    if(tab[idx[2-1]]<tab[idx[3-1]] )  swi(idx[2-1],idx[3-1]);
    if(tab[idx[4-1]]<tab[idx[5-1]] )  swi(idx[4-1],idx[5-1]);
    if(tab[idx[6-1]]<tab[idx[7-1]] )  swi(idx[6-1],idx[7-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]] )  swi(idx[1-1],idx[3-1]);
    if(tab[idx[4-1]]<tab[idx[6-1]] )  swi(idx[4-1],idx[6-1]);
    if(tab[idx[5-1]]<tab[idx[7-1]] )  swi(idx[5-1],idx[7-1]);
    
    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);
    if(tab[idx[5-1]]<tab[idx[6-1]] )  swi(idx[5-1],idx[6-1]);
    if(tab[idx[3-1]]<tab[idx[7-1]] )  swi(idx[3-1],idx[7-1]);
    
    if(tab[idx[1-1]]<tab[idx[5-1]] )  swi(idx[1-1],idx[5-1]);
    if(tab[idx[2-1]]<tab[idx[6-1]] )  swi(idx[2-1],idx[6-1]);
    
    if(tab[idx[1-1]]<tab[idx[4-1]] )  swi(idx[1-1],idx[4-1]);
    if(tab[idx[3-1]]<tab[idx[6-1]] )  swi(idx[3-1],idx[6-1]);
    
    if(tab[idx[2-1]]<tab[idx[4-1]] )  swi(idx[2-1],idx[4-1]);
    if(tab[idx[3-1]]<tab[idx[5-1]] )  swi(idx[3-1],idx[5-1]);
    
    if(tab[idx[3-1]]<tab[idx[4-1]] )  swi(idx[3-1],idx[4-1]);
    
    
  }else if(n==6){
    
    
    if(tab[idx[2-1]]<tab[idx[3-1]] )  swi(idx[2-1],idx[3-1]);
    if(tab[idx[5-1]]<tab[idx[6-1]] )  swi(idx[5-1],idx[6-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]] )  swi(idx[1-1],idx[3-1]);
    if(tab[idx[4-1]]<tab[idx[6-1]] )  swi(idx[4-1],idx[6-1]);
    
    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);
    if(tab[idx[4-1]]<tab[idx[5-1]] )  swi(idx[4-1],idx[5-1]);
    if(tab[idx[3-1]]<tab[idx[6-1]] )  swi(idx[3-1],idx[6-1]);
    
    if(tab[idx[1-1]]<tab[idx[4-1]] )  swi(idx[1-1],idx[4-1]);
    if(tab[idx[2-1]]<tab[idx[5-1]] )  swi(idx[2-1],idx[5-1]);
    
    if(tab[idx[3-1]]<tab[idx[5-1]] )  swi(idx[3-1],idx[5-1]);
    if(tab[idx[2-1]]<tab[idx[4-1]] )  swi(idx[2-1],idx[4-1]);
    
    if(tab[idx[3-1]]<tab[idx[4-1]] )  swi(idx[3-1],idx[4-1]);
    
    
  }else if(n==5){
    
   
    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);
    if(tab[idx[4-1]]<tab[idx[5-1]] )  swi(idx[4-1],idx[5-1]);
    
    if(tab[idx[3-1]]<tab[idx[5-1]] )  swi(idx[3-1],idx[5-1]);
    
    if(tab[idx[3-1]]<tab[idx[4-1]] )  swi(idx[3-1],idx[4-1]);
    if(tab[idx[2-1]]<tab[idx[5-1]] )  swi(idx[2-1],idx[5-1]);
    
    if(tab[idx[1-1]]<tab[idx[4-1]] )  swi(idx[1-1],idx[4-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]] )  swi(idx[1-1],idx[3-1]);
    if(tab[idx[2-1]]<tab[idx[4-1]] )  swi(idx[2-1],idx[4-1]);
    
    if(tab[idx[2-1]]<tab[idx[3-1]] )  swi(idx[2-1],idx[3-1]);
    
    
  }else if(n==4){
    
   
    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);
    if(tab[idx[3-1]]<tab[idx[4-1]] )  swi(idx[3-1],idx[4-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]] )  swi(idx[1-1],idx[3-1]);
    if(tab[idx[2-1]]<tab[idx[4-1]] )  swi(idx[2-1],idx[4-1]);
    
    if(tab[idx[2-1]]<tab[idx[3-1]] )  swi(idx[2-1],idx[3-1]);
    
    
  }else if(n==3){
    
   
    if(tab[idx[2-1]]<tab[idx[3-1]] )  swi(idx[2-1],idx[3-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]] )  swi(idx[1-1],idx[3-1]);
    
    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);
    

  }else if(n==2){

    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);

  }else if(n<1){
    METRIS_THROW(WArgExcept());
  }
}
#if 0
// Optimal sorting networks for up to 8 elements
// Sort as decreasing
template<typename T, int n>
void sortupto8_dec(T tab[n], int idx[n]){

  if constexpr (n == 8){
    if(tab[idx[1-1]]<tab[idx[2-1]])  swi(idx[1-1],idx[2-1]);
    if(tab[idx[3-1]]<tab[idx[4-1]])  swi(idx[3-1],idx[4-1]);
    if(tab[idx[5-1]]<tab[idx[6-1]])  swi(idx[5-1],idx[6-1]);
    if(tab[idx[7-1]]<tab[idx[8-1]])  swi(idx[7-1],idx[8-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]])  swi(idx[1-1],idx[3-1]);
    if(tab[idx[2-1]]<tab[idx[4-1]])  swi(idx[2-1],idx[4-1]) ;
    if(tab[idx[5-1]]<tab[idx[7-1]])  swi(idx[5-1],idx[7-1]) ;
    if(tab[idx[6-1]]<tab[idx[8-1]])  swi(idx[6-1],idx[8-1]) ;
    
    if(tab[idx[2-1]]<tab[idx[3-1]])  swi(idx[2-1],idx[3-1]);
    if(tab[idx[6-1]]<tab[idx[7-1]])  swi(idx[6-1],idx[7-1]);
    if(tab[idx[1-1]]<tab[idx[5-1]])  swi(idx[1-1],idx[5-1]);
    if(tab[idx[4-1]]<tab[idx[8-1]])  swi(idx[4-1],idx[8-1]);
    
    if(tab[idx[2-1]]<tab[idx[6-1]])  swi(idx[2-1],idx[6-1]);
    if(tab[idx[3-1]]<tab[idx[7-1]])  swi(idx[3-1],idx[7-1]);
    
    if(tab[idx[2-1]]<tab[idx[5-1]])  swi(idx[2-1],idx[5-1]);
    if(tab[idx[4-1]]<tab[idx[7-1]])  swi(idx[4-1],idx[7-1]);
    
    if(tab[idx[3-1]]<tab[idx[5-1]])  swi(idx[3-1],idx[5-1]);
    if(tab[idx[4-1]]<tab[idx[6-1]])  swi(idx[4-1],idx[6-1]);
    
    if(tab[idx[4-1]]<tab[idx[5-1]])  swi(idx[4-1],idx[5-1]);

  }else if(n==7){

    if(tab[idx[2-1]]<tab[idx[3-1]] )  swi(idx[2-1],idx[3-1]);
    if(tab[idx[4-1]]<tab[idx[5-1]] )  swi(idx[4-1],idx[5-1]);
    if(tab[idx[6-1]]<tab[idx[7-1]] )  swi(idx[6-1],idx[7-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]] )  swi(idx[1-1],idx[3-1]);
    if(tab[idx[4-1]]<tab[idx[6-1]] )  swi(idx[4-1],idx[6-1]);
    if(tab[idx[5-1]]<tab[idx[7-1]] )  swi(idx[5-1],idx[7-1]);
    
    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);
    if(tab[idx[5-1]]<tab[idx[6-1]] )  swi(idx[5-1],idx[6-1]);
    if(tab[idx[3-1]]<tab[idx[7-1]] )  swi(idx[3-1],idx[7-1]);
    
    if(tab[idx[1-1]]<tab[idx[5-1]] )  swi(idx[1-1],idx[5-1]);
    if(tab[idx[2-1]]<tab[idx[6-1]] )  swi(idx[2-1],idx[6-1]);
    
    if(tab[idx[1-1]]<tab[idx[4-1]] )  swi(idx[1-1],idx[4-1]);
    if(tab[idx[3-1]]<tab[idx[6-1]] )  swi(idx[3-1],idx[6-1]);
    
    if(tab[idx[2-1]]<tab[idx[4-1]] )  swi(idx[2-1],idx[4-1]);
    if(tab[idx[3-1]]<tab[idx[5-1]] )  swi(idx[3-1],idx[5-1]);
    
    if(tab[idx[3-1]]<tab[idx[4-1]] )  swi(idx[3-1],idx[4-1]);
    
    
  }else if(n==6){
    
    
    if(tab[idx[2-1]]<tab[idx[3-1]] )  swi(idx[2-1],idx[3-1]);
    if(tab[idx[5-1]]<tab[idx[6-1]] )  swi(idx[5-1],idx[6-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]] )  swi(idx[1-1],idx[3-1]);
    if(tab[idx[4-1]]<tab[idx[6-1]] )  swi(idx[4-1],idx[6-1]);
    
    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);
    if(tab[idx[4-1]]<tab[idx[5-1]] )  swi(idx[4-1],idx[5-1]);
    if(tab[idx[3-1]]<tab[idx[6-1]] )  swi(idx[3-1],idx[6-1]);
    
    if(tab[idx[1-1]]<tab[idx[4-1]] )  swi(idx[1-1],idx[4-1]);
    if(tab[idx[2-1]]<tab[idx[5-1]] )  swi(idx[2-1],idx[5-1]);
    
    if(tab[idx[3-1]]<tab[idx[5-1]] )  swi(idx[3-1],idx[5-1]);
    if(tab[idx[2-1]]<tab[idx[4-1]] )  swi(idx[2-1],idx[4-1]);
    
    if(tab[idx[3-1]]<tab[idx[4-1]] )  swi(idx[3-1],idx[4-1]);
    
    
  }else if(n==5){
    
   
    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);
    if(tab[idx[4-1]]<tab[idx[5-1]] )  swi(idx[4-1],idx[5-1]);
    
    if(tab[idx[3-1]]<tab[idx[5-1]] )  swi(idx[3-1],idx[5-1]);
    
    if(tab[idx[3-1]]<tab[idx[4-1]] )  swi(idx[3-1],idx[4-1]);
    if(tab[idx[2-1]]<tab[idx[5-1]] )  swi(idx[2-1],idx[5-1]);
    
    if(tab[idx[1-1]]<tab[idx[4-1]] )  swi(idx[1-1],idx[4-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]] )  swi(idx[1-1],idx[3-1]);
    if(tab[idx[2-1]]<tab[idx[4-1]] )  swi(idx[2-1],idx[4-1]);
    
    if(tab[idx[2-1]]<tab[idx[3-1]] )  swi(idx[2-1],idx[3-1]);
    
    
  }else if(n==4){
    
   
    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);
    if(tab[idx[3-1]]<tab[idx[4-1]] )  swi(idx[3-1],idx[4-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]] )  swi(idx[1-1],idx[3-1]);
    if(tab[idx[2-1]]<tab[idx[4-1]] )  swi(idx[2-1],idx[4-1]);
    
    if(tab[idx[2-1]]<tab[idx[3-1]] )  swi(idx[2-1],idx[3-1]);
    
    
  }else if(n==3){
    
   
    if(tab[idx[2-1]]<tab[idx[3-1]] )  swi(idx[2-1],idx[3-1]);
    
    if(tab[idx[1-1]]<tab[idx[3-1]] )  swi(idx[1-1],idx[3-1]);
    
    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);
    

  }else if(n==2){

    if(tab[idx[1-1]]<tab[idx[2-1]] )  swi(idx[1-1],idx[2-1]);

  }else if(n<1){
    METRIS_THROW(WArgExcept());
  }
}
#endif

template<typename T>
void sortupto8_inc(T *tab,int n){

  if(n == 8){
    if(tab[1-1]>tab[2-1])  swi(tab[1-1],tab[2-1]);
    if(tab[3-1]>tab[4-1])  swi(tab[3-1],tab[4-1]);
    if(tab[5-1]>tab[6-1])  swi(tab[5-1],tab[6-1]);
    if(tab[7-1]>tab[8-1])  swi(tab[7-1],tab[8-1]);
    
    if(tab[1-1]>tab[3-1])  swi(tab[1-1],tab[3-1]);
    if(tab[2-1]>tab[4-1])  swi(tab[2-1],tab[4-1]) ;
    if(tab[5-1]>tab[7-1])  swi(tab[5-1],tab[7-1]) ;
    if(tab[6-1]>tab[8-1])  swi(tab[6-1],tab[8-1]) ;
    
    if(tab[2-1]>tab[3-1])  swi(tab[2-1],tab[3-1]);
    if(tab[6-1]>tab[7-1])  swi(tab[6-1],tab[7-1]);
    if(tab[1-1]>tab[5-1])  swi(tab[1-1],tab[5-1]);
    if(tab[4-1]>tab[8-1])  swi(tab[4-1],tab[8-1]);
    
    if(tab[2-1]>tab[6-1])  swi(tab[2-1],tab[6-1]);
    if(tab[3-1]>tab[7-1])  swi(tab[3-1],tab[7-1]);
    
    if(tab[2-1]>tab[5-1])  swi(tab[2-1],tab[5-1]);
    if(tab[4-1]>tab[7-1])  swi(tab[4-1],tab[7-1]);
    
    if(tab[3-1]>tab[5-1])  swi(tab[3-1],tab[5-1]);
    if(tab[4-1]>tab[6-1])  swi(tab[4-1],tab[6-1]);
    
    if(tab[4-1]>tab[5-1])  swi(tab[4-1],tab[5-1]);

  }else if(n==7){

    if(tab[2-1]>tab[3-1] )  swi(tab[2-1],tab[3-1]);
    if(tab[4-1]>tab[5-1] )  swi(tab[4-1],tab[5-1]);
    if(tab[6-1]>tab[7-1] )  swi(tab[6-1],tab[7-1]);
    
    if(tab[1-1]>tab[3-1] )  swi(tab[1-1],tab[3-1]);
    if(tab[4-1]>tab[6-1] )  swi(tab[4-1],tab[6-1]);
    if(tab[5-1]>tab[7-1] )  swi(tab[5-1],tab[7-1]);
    
    if(tab[1-1]>tab[2-1] )  swi(tab[1-1],tab[2-1]);
    if(tab[5-1]>tab[6-1] )  swi(tab[5-1],tab[6-1]);
    if(tab[3-1]>tab[7-1] )  swi(tab[3-1],tab[7-1]);
    
    if(tab[1-1]>tab[5-1] )  swi(tab[1-1],tab[5-1]);
    if(tab[2-1]>tab[6-1] )  swi(tab[2-1],tab[6-1]);
    
    if(tab[1-1]>tab[4-1] )  swi(tab[1-1],tab[4-1]);
    if(tab[3-1]>tab[6-1] )  swi(tab[3-1],tab[6-1]);
    
    if(tab[2-1]>tab[4-1] )  swi(tab[2-1],tab[4-1]);
    if(tab[3-1]>tab[5-1] )  swi(tab[3-1],tab[5-1]);
    
    if(tab[3-1]>tab[4-1] )  swi(tab[3-1],tab[4-1]);
    
    
  }else if(n==6){
    
    
    if(tab[2-1]>tab[3-1] )  swi(tab[2-1],tab[3-1]);
    if(tab[5-1]>tab[6-1] )  swi(tab[5-1],tab[6-1]);
    
    if(tab[1-1]>tab[3-1] )  swi(tab[1-1],tab[3-1]);
    if(tab[4-1]>tab[6-1] )  swi(tab[4-1],tab[6-1]);
    
    if(tab[1-1]>tab[2-1] )  swi(tab[1-1],tab[2-1]);
    if(tab[4-1]>tab[5-1] )  swi(tab[4-1],tab[5-1]);
    if(tab[3-1]>tab[6-1] )  swi(tab[3-1],tab[6-1]);
    
    if(tab[1-1]>tab[4-1] )  swi(tab[1-1],tab[4-1]);
    if(tab[2-1]>tab[5-1] )  swi(tab[2-1],tab[5-1]);
    
    if(tab[3-1]>tab[5-1] )  swi(tab[3-1],tab[5-1]);
    if(tab[2-1]>tab[4-1] )  swi(tab[2-1],tab[4-1]);
    
    if(tab[3-1]>tab[4-1] )  swi(tab[3-1],tab[4-1]);
    
    
  }else if(n==5){
    
   
    if(tab[1-1]>tab[2-1] )  swi(tab[1-1],tab[2-1]);
    if(tab[4-1]>tab[5-1] )  swi(tab[4-1],tab[5-1]);
    
    if(tab[3-1]>tab[5-1] )  swi(tab[3-1],tab[5-1]);
    
    if(tab[3-1]>tab[4-1] )  swi(tab[3-1],tab[4-1]);
    if(tab[2-1]>tab[5-1] )  swi(tab[2-1],tab[5-1]);
    
    if(tab[1-1]>tab[4-1] )  swi(tab[1-1],tab[4-1]);
    
    if(tab[1-1]>tab[3-1] )  swi(tab[1-1],tab[3-1]);
    if(tab[2-1]>tab[4-1] )  swi(tab[2-1],tab[4-1]);
    
    if(tab[2-1]>tab[3-1] )  swi(tab[2-1],tab[3-1]);
    
    
  }else if(n==4){
    
   
    if(tab[1-1]>tab[2-1] )  swi(tab[1-1],tab[2-1]);
    if(tab[3-1]>tab[4-1] )  swi(tab[3-1],tab[4-1]);
    
    if(tab[1-1]>tab[3-1] )  swi(tab[1-1],tab[3-1]);
    if(tab[2-1]>tab[4-1] )  swi(tab[2-1],tab[4-1]);
    
    if(tab[2-1]>tab[3-1] )  swi(tab[2-1],tab[3-1]);
    
    
  }else if(n==3){
    
   
    if(tab[2-1]>tab[3-1] )  swi(tab[2-1],tab[3-1]);
    
    if(tab[1-1]>tab[3-1] )  swi(tab[1-1],tab[3-1]);
    
    if(tab[1-1]>tab[2-1] )  swi(tab[1-1],tab[2-1]);
    

  }else if(n==2){

    if(tab[1-1]>tab[2-1] )  swi(tab[1-1],tab[2-1]);

  }else if(n<1){
    METRIS_THROW(WArgExcept());
  }
}

template void sortupto8_dec<int>(int *tab,int n);
template void sortupto8_inc<int>(int *tab,int n);

template void sortupto8_dec<double>(double *tab,int n);
template void sortupto8_dec<double>(double *tab,int * idx, int n);
//template void sortupto8_dec<double,2>(double tab[2],int idx[2]);
//template void sortupto8_dec<double,3>(double tab[3],int idx[3]);
//template void sortupto8_dec<double,4>(double tab[4],int idx[4]);
template void sortupto8_inc<double>(double *tab,int n);

template void sortupto8_dec<SANS::SurrealS<2,double>>(SANS::SurrealS<2,double> *tab,int n);
template void sortupto8_inc<SANS::SurrealS<2,double>>(SANS::SurrealS<2,double> *tab,int n);

template void sortupto8_dec<SANS::SurrealS<3,double>>(SANS::SurrealS<3,double> *tab,int n);
template void sortupto8_inc<SANS::SurrealS<3,double>>(SANS::SurrealS<3,double> *tab,int n);



} // End namespace

