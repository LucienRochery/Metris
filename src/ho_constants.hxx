//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __HO_CONSTANTS__
#define __HO_CONSTANTS__


#include "metris_constants.hxx"

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/iteration/local.hpp>

#include <array>

namespace Metris{

/*
These arrays started their lives as dynamically allocated intAr2/3 etc type objects, constructed using
a move constructor like so:
const foo(std::move(init_foo()))
which, unlike the copy constructor, cedes control of the pointer to the receptor. 

This was not optimal. These arrays are basically compile time constants and I only included them
here algorithmically rather than as a bunch of hard coded numbers for two reasons:
- compactness; higher degree tets would have been painful. For instance, P10 tetrahedra will have 286 
nodes, times 4 multi-indices (we could compact by one), times how many considered degrees (we could 
also compress here). This is possibly dozens of thousands of entries just for one array. 
- flexibility: this will more easily allow us to switch ordering schemes if we want to.
- on a related note: debuggability.  

The current solution, apart from the simplest arrays, is to (ab)use constexpr statements in order
to have the compiler pre-compute everything. 
This seems to be very powerful as we have the advantages of "dynamically" constructed arrays as well as 
those of statically sized arrays. 
This + a couple of insignificant optimizations made cycle count drop from 160k to 134k for main, 
putting the C++ implementation of degree elevation slightly ahead in speed of the Fortran one. 

TODO: Can we further propagate "constness" into other areas of the code? Some loops?
The functions to generate new HO nodes from a given degree element? 
Would this be advantageous?
*/

/*
These are basic mesh constants
*/

constexpr void ini_ordedg();
constexpr void ini_ordfac();
constexpr void ini_ordtet();
constexpr void ini_invordedg();
constexpr void ini_invordfac();
constexpr void ini_invordtet();
constexpr void ini_cbzedg();
constexpr void ini_cbzfac();
constexpr void ini_cbztet();


// Degree of the maximum pnorm (<= 2) interpolation error integrand
#ifndef METRIS_MAX_DEG_INTPINTGD
#define METRIS_MAX_DEG_INTPINTGD METRIS_MAX_DEG_JACOBIAN + 2*(METRIS_MAX_DEG+1)*METRIS_MAX_DEG
#endif

#ifndef METRIS_MAX_DEG_ORDERING
#define METRIS_MAX_DEG_ORDERING METRIS_MAX_DEG_INTPINTGD
#endif 

// Maximum we need is to support pointwise error of degree (MAX_DEG + 1)*MAX_DEG
// Squaring (or any power) is handled using Bernstein multiplication, no Lag2Bez
// This is for debugging purposes (exact error). Hopefully we can cut it closer 
// to the bone at MAX_DEG + 1.
#ifndef METRIS_MAX_DEG_LAG2BEZ
#define METRIS_MAX_DEG_LAG2BEZ (METRIS_MAX_DEG+1)*METRIS_MAX_DEG
#endif

#ifndef METRIS_MAX_DEG_EVAL
#define METRIS_MAX_DEG_EVAL MAX(METRIS_MAX_DEG_LAG2BEZ, METRIS_MAX_DEG_JACOBIAN)
#endif
//constexpr int METRIS_MAX_DEG = 4;

#define SMOO_DEGJ(x) x - 1

#define ORDELT(n) []() -> auto {\
         if constexpr(n == 1){return ordedg.s;}\
    else if constexpr(n == 2){return ordfac.s;}\
    else                     {return ordtet.s;}\
  }();



constexpr int getnnod1(int i){
  return i+1;
}

constexpr int getnnod2(int i){
  return ((i+2)*(i+1))/2;
}

constexpr int getnnod3(int i){
  return ((i+3)*(i+2)*(i+1))/6;
}

constexpr int getnnode(int idim, int i){
  switch(idim){
    case 1:
      return getnnod1(i);
    case 2:
      return getnnod2(i);
    case 3:
      return getnnod3(i);
  }
  return -1;
}



// ord(edg|fac|tet).s[ideg][inode][0:tdim+1]

// replace ord(edg|fac|tet) by ord(edg|fac|tet).s
// array is not constexpr compatible in C++17
struct ordedg_constructor{
  constexpr ordedg_constructor() : s() {

    s[0][0][0] = 0;
    s[0][0][1] = 0;

    for(int ideg = 1; ideg <= METRIS_MAX_DEG_ORDERING; ideg++){
      s[ideg][0][0] = ideg;
      s[ideg][0][1] = 0;

      s[ideg][1][0] = 0;
      s[ideg][1][1] = ideg;

      if(ideg == 1) continue;

      int n = 2;
      for(int i = 1; i <= ideg-1; i++){
        n = n + 1;
        s[ideg][n-1][0] = ideg - i;
        s[ideg][n-1][1] =        i;

      }
    }  

  }
  int s[1+METRIS_MAX_DEG_ORDERING][getnnod1(METRIS_MAX_DEG_ORDERING)][2];

  // These guys do not work.
  template<int ideg>
  constexpr std::array<std::array<int,2>,getnnod1(ideg)> get() const{
    std::array<std::array<int,2>,getnnod1(ideg)> arr = 0;

    arr[0][0] = ideg;
    arr[0][1] = 0;

    arr[1][0] = 0;
    arr[1][1] = ideg;

    if(ideg == 1) return arr;

    int n = 2;
    for(int i = 1; i <= ideg-1; i++){
      n = n + 1;
      arr[n-1][0] = ideg - i;
      arr[n-1][1] =        i;
    }

    return arr;
  }
};

constexpr const ordedg_constructor ordedg;

struct ordfac_constructor{
  constexpr ordfac_constructor() : s(){

    s[0][0][0] = 0;
    s[0][0][1] = 0;
    s[0][0][2] = 0;

    for(int ideg = 1; ideg <= METRIS_MAX_DEG_ORDERING; ideg++){
      s[ideg][0][0] = ideg;
      s[ideg][0][1] = 0;
      s[ideg][0][2] = 0;

      s[ideg][1][0] = 0;
      s[ideg][1][1] = ideg;
      s[ideg][1][2] = 0;

      s[ideg][2][0] = 0;
      s[ideg][2][1] = 0;
      s[ideg][2][2] = ideg;
      
      if(ideg == 1) continue;

      int n = 3;
      for(int i = 1; i <= ideg-1; i++){
        n = n + 1;
        s[ideg][n-1][0] = 0;
        s[ideg][n-1][1] = ideg - i;
        s[ideg][n-1][2] =        i;
      }
      for(int i = 1; i <= ideg-1; i++){
        n = n + 1;
        s[ideg][n-1][0] =        i;
        s[ideg][n-1][1] = 0;
        s[ideg][n-1][2] = ideg - i;
      } 
      for(int i = 1; i <= ideg-1; i++){
        n = n + 1;
        s[ideg][n-1][0] = ideg - i;
        s[ideg][n-1][1] =        i;
        s[ideg][n-1][2] = 0;
      }

      if(ideg == 2) continue;

      for(int i=0;i<getnnod2(ideg-3);i++){ 
        n = n + 1;
        s[ideg][n-1][0] = 1 + s[ideg-3][i][0];
        s[ideg][n-1][1] = 1 + s[ideg-3][i][1];
        s[ideg][n-1][2] = 1 + s[ideg-3][i][2];
      }
    }  
  }
  int s[1+METRIS_MAX_DEG_ORDERING][getnnod2(METRIS_MAX_DEG_ORDERING)][3];


  template<int ideg>
  constexpr std::array<std::array<int,3>,getnnod2(ideg)> get() const{
    std::array<std::array<int,3>,getnnod2(ideg)> arr;

    arr[0][0] = ideg;
    arr[0][1] = 0;
    arr[0][2] = 0;

    arr[1][0] = 0;
    arr[1][1] = ideg;
    arr[1][2] = 0;

    arr[2][0] = 0;
    arr[2][1] = 0;
    arr[2][2] = ideg;
    
    if(ideg == 1) return arr;

    int n = 3;
    for(int i = 1; i <= ideg-1; i++){
      n = n + 1;
      arr[n-1][0] = 0;
      arr[n-1][1] = ideg - i;
      arr[n-1][2] =        i;
    }
    for(int i = 1; i <= ideg-1; i++){
      n = n + 1;
      arr[n-1][0] =        i;
      arr[n-1][1] = 0;
      arr[n-1][2] = ideg - i;
    } 
    for(int i = 1; i <= ideg-1; i++){
      n = n + 1;
      arr[n-1][0] = ideg - i;
      arr[n-1][1] =        i;
      arr[n-1][2] = 0;
    }

    if(ideg == 2) return arr;

    for(int i=0;i<getnnod2(ideg-3);i++){ 
      n = n + 1;
      arr[n-1][0] = 1 + s[ideg-3][i][0];
      arr[n-1][1] = 1 + s[ideg-3][i][1];
      arr[n-1][2] = 1 + s[ideg-3][i][2];
    }
    return arr;
  }
};

constexpr const ordfac_constructor ordfac;

struct ordtet_constructor{
  constexpr ordtet_constructor() : s(), mask(){

    // Degree 0: 1 node
    s[0][0][0] = 0;
    s[0][0][1] = 0;
    s[0][0][2] = 0;
    s[0][0][3] = 0;

    for(int ideg = 1; ideg <= METRIS_MAX_DEG_ORDERING; ideg++){
      // Start with vertices: always the first 4 
      s[ideg][0][0] = ideg;
      s[ideg][0][1] = 0;
      s[ideg][0][2] = 0;
      s[ideg][0][3] = 0;

      s[ideg][1][0] = 0;
      s[ideg][1][1] = ideg;
      s[ideg][1][2] = 0;
      s[ideg][1][3] = 0;

      s[ideg][2][0] = 0;
      s[ideg][2][1] = 0;
      s[ideg][2][2] = ideg;
      s[ideg][2][3] = 0;

      s[ideg][3][0] = 0;
      s[ideg][3][1] = 0;
      s[ideg][3][2] = 0;
      s[ideg][3][3] = ideg;

      if(ideg == 1) continue;

      // Edges 
      // Current writable index
      int n = 4;
      for(int i = 1; i <= ideg-1; i++){
        n = n + 1;
        s[ideg][n-1][0] = ideg - i;
        s[ideg][n-1][1] =        i;
        s[ideg][n-1][2] = 0;
        s[ideg][n-1][3] = 0;
      }

      for(int i = 1; i <= ideg-1; i++){
        n = n + 1;
        s[ideg][n-1][0] = 0;
        s[ideg][n-1][1] = ideg - i;
        s[ideg][n-1][2] =        i;
        s[ideg][n-1][3] = 0;
      }

      for(int i = 1; i <= ideg-1; i++){
        n = n + 1;
        s[ideg][n-1][0] =        i;
        s[ideg][n-1][1] = 0;
        s[ideg][n-1][2] = ideg - i;
        s[ideg][n-1][3] = 0;
      } 

      for(int i = 1; i <= ideg-1; i++){
        n = n + 1;
        s[ideg][n-1][0] = ideg - i;
        s[ideg][n-1][1] = 0;
        s[ideg][n-1][2] = 0;
        s[ideg][n-1][3] = i;
      } 

      for(int i = 1; i <= ideg-1; i++){
        n = n + 1;
        s[ideg][n-1][0] = 0;
        s[ideg][n-1][1] = ideg - i;
        s[ideg][n-1][2] = 0;
        s[ideg][n-1][3] = i;
      } 

      for(int i = 1; i <= ideg-1; i++){
        n = n + 1;
        s[ideg][n-1][0] = 0;
        s[ideg][n-1][1] = 0;
        s[ideg][n-1][2] = ideg - i;
        s[ideg][n-1][3] = i;
      } 


      if(ideg <= 2) continue;

      // Faces: inherited from deg - 3 faces. 
      for(int i=0;i<getnnod2(ideg-3);i++){
        n = n + 1;
        s[ideg][n-1][0] = 0              ;
        s[ideg][n-1][1] = 1+ordfac.s[ideg-3][i][0] ;
        s[ideg][n-1][2] = 1+ordfac.s[ideg-3][i][1] ;
        s[ideg][n-1][3] = 1+ordfac.s[ideg-3][i][2];
      }

      for(int i=0;i<getnnod2(ideg-3);i++){
        n = n + 1;
        s[ideg][n-1][0] = 1+ordfac.s[ideg-3][i][1];
        s[ideg][n-1][1] = 0      ;
        s[ideg][n-1][2] = 1+ordfac.s[ideg-3][i][0] ;
        s[ideg][n-1][3] = 1+ordfac.s[ideg-3][i][2];
      }

      for(int i=0;i<getnnod2(ideg-3);i++){
        n = n + 1;
        s[ideg][n-1][0] = 1+ordfac.s[ideg-3][i][0];
        s[ideg][n-1][1] = 1+ordfac.s[ideg-3][i][1] ;
        s[ideg][n-1][2] = 0      ;
        s[ideg][n-1][3] = 1+ordfac.s[ideg-3][i][2];
      }

      for(int i=0;i<getnnod2(ideg-3);i++){
        n = n + 1;
        s[ideg][n-1][0] = 1+ordfac.s[ideg-3][i][0];
        s[ideg][n-1][1] = 1+ordfac.s[ideg-3][i][1] ;
        s[ideg][n-1][2] = 1+ordfac.s[ideg-3][i][2] ;
        s[ideg][n-1][3] = 0;
      }

      if(ideg <= 3) continue; 

      for(int i=0;i<getnnod3(ideg-4);i++){
        n = n + 1;
        s[ideg][n-1][0] = 1+s[ideg-4][i][0];
        s[ideg][n-1][1] = 1+s[ideg-4][i][1];
        s[ideg][n-1][2] = 1+s[ideg-4][i][2];
        s[ideg][n-1][3] = 1+s[ideg-4][i][3];
      }
    }

    mask[0][0] = 0b1111;
    for(int ideg = 1; ideg <= METRIS_MAX_DEG_ORDERING; ideg++){
      for(int i = 0; i < getnnod3(ideg); i++){
        mask[ideg][i] = (1 << 3)*(s[ideg][i][0] > 0)
                      + (1 << 2)*(s[ideg][i][1] > 0)
                      + (1 << 1)*(s[ideg][i][2] > 0)
                      + (1 << 0)*(s[ideg][i][3] > 0);
      }
    }

  }

  int s[1+METRIS_MAX_DEG_ORDERING][getnnod3(METRIS_MAX_DEG_ORDERING)][4];
  unsigned char mask[1+METRIS_MAX_DEG_ORDERING][getnnod3(METRIS_MAX_DEG_ORDERING)];

  template<int ideg>
  constexpr std::array<std::array<int,4>,getnnod3(ideg)> get() const{
    std::array<std::array<int,4>,getnnod3(ideg)> arr;

    // Start with vertices: always the first 4 
    arr[0][0] = ideg;
    arr[0][1] = 0;
    arr[0][2] = 0;
    arr[0][3] = 0;

    arr[1][0] = 0;
    arr[1][1] = ideg;
    arr[1][2] = 0;
    arr[1][3] = 0;

    arr[2][0] = 0;
    arr[2][1] = 0;
    arr[2][2] = ideg;
    arr[2][3] = 0;

    arr[3][0] = 0;
    arr[3][1] = 0;
    arr[3][2] = 0;
    arr[3][3] = ideg;

    if(ideg == 1) return arr;

    // Edges 
    // Current writable index
    int n = 4;
    for(int i = 1; i <= ideg-1; i++){
      n = n + 1;
      arr[n-1][0] = ideg - i;
      arr[n-1][1] =        i;
      arr[n-1][2] = 0;
      arr[n-1][3] = 0;
    }

    for(int i = 1; i <= ideg-1; i++){
      n = n + 1;
      arr[n-1][0] = 0;
      arr[n-1][1] = ideg - i;
      arr[n-1][2] =        i;
      arr[n-1][3] = 0;
    }

    for(int i = 1; i <= ideg-1; i++){
      n = n + 1;
      arr[n-1][0] =        i;
      arr[n-1][1] = 0;
      arr[n-1][2] = ideg - i;
      arr[n-1][3] = 0;
    } 

    for(int i = 1; i <= ideg-1; i++){
      n = n + 1;
      arr[n-1][0] = ideg - i;
      arr[n-1][1] = 0;
      arr[n-1][2] = 0;
      arr[n-1][3] = i;
    } 

    for(int i = 1; i <= ideg-1; i++){
      n = n + 1;
      arr[n-1][0] = 0;
      arr[n-1][1] = ideg - i;
      arr[n-1][2] = 0;
      arr[n-1][3] = i;
    } 

    for(int i = 1; i <= ideg-1; i++){
      n = n + 1;
      arr[n-1][0] = 0;
      arr[n-1][1] = 0;
      arr[n-1][2] = ideg - i;
      arr[n-1][3] = i;
    } 


    if(ideg <= 2) return arr;

    // Faces: inherited from deg - 3 faces. 
    for(int i=0;i<getnnod2(ideg-3);i++){
      n = n + 1;
      arr[n-1][0] = 0              ;
      arr[n-1][1] = 1+ordfac.s[ideg-3][i][0] ;
      arr[n-1][2] = 1+ordfac.s[ideg-3][i][1] ;
      arr[n-1][3] = 1+ordfac.s[ideg-3][i][2];
    }

    for(int i=0;i<getnnod2(ideg-3);i++){
      n = n + 1;
      arr[n-1][0] = 1+ordfac.s[ideg-3][i][1];
      arr[n-1][1] = 0      ;
      arr[n-1][2] = 1+ordfac.s[ideg-3][i][0] ;
      arr[n-1][3] = 1+ordfac.s[ideg-3][i][2];
    }

    for(int i=0;i<getnnod2(ideg-3);i++){
      n = n + 1;
      arr[n-1][0] = 1+ordfac.s[ideg-3][i][0];
      arr[n-1][1] = 1+ordfac.s[ideg-3][i][1] ;
      arr[n-1][2] = 0      ;
      arr[n-1][3] = 1+ordfac.s[ideg-3][i][2];
    }

    for(int i=0;i<getnnod2(ideg-3);i++){
      n = n + 1;
      arr[n-1][0] = 1+ordfac.s[ideg-3][i][0];
      arr[n-1][1] = 1+ordfac.s[ideg-3][i][1] ;
      arr[n-1][2] = 1+ordfac.s[ideg-3][i][2] ;
      arr[n-1][3] = 0;
    }

    if(ideg <= 3) return arr; 

    for(int i=0;i<getnnod3(ideg-4);i++){
      n = n + 1;
      arr[n-1][0] = 1+s[ideg-4][i][0];
      arr[n-1][1] = 1+s[ideg-4][i][1];
      arr[n-1][2] = 1+s[ideg-4][i][2];
      arr[n-1][3] = 1+s[ideg-4][i][3];
    }
    return arr;
  }
};


constexpr const ordtet_constructor ordtet;


//constexpr int idx_invord_tab_fac(int d,int i,int j,int k);
//constexpr int idx_invord_tab_tet(int d,int i,int j,int k,int l);

constexpr  int idx_invord_tab_fac(int i,int j,int k){
  if(i==(i+j+k)){
    return 0;
  }
  else{
    return j + getnnod2(j+k-1);
  }
}
constexpr  int idx_invord_tab_tet(int i,int j,int k,int l){
  if(i == (i+j+k+l)){
    return 0;
  }
  else{
//        printf("Debug d-i={} j = {} k = {} l = {} fac = {} \n",d-i,j,k,l,idx_invord_tab_fac(d-i,j,k,l));
    return idx_invord_tab_fac(j,k,l) + getnnod3(j+k+l-1);
  }

}

struct invordedg_constructor{
  constexpr invordedg_constructor() : s(){
    s[0][0] = 0;
    for (int ideg = 1;ideg <= METRIS_MAX_DEG_ORDERING; ideg++){
      int npp = getnnod1(ideg);
      for(int ip = 0; ip <npp; ip++){
        int i = ordedg.s[ideg][ip][0];
        //int j = ordedg.s[ideg][ip][1];
        s[ideg][i] = ip;
      }
    }
  }
  int s[1+METRIS_MAX_DEG_ORDERING][getnnod1(METRIS_MAX_DEG_ORDERING)];
};

constexpr const invordedg_constructor invordedg;
constexpr int idx_invord_tab_edg(int d,int i,int j);

struct invordfac_constructor{
  constexpr invordfac_constructor() : s(){
    s[0][0] = 0;
    for (int ideg = 1;ideg <= METRIS_MAX_DEG_ORDERING; ideg++){
      int npp = getnnod2(ideg);
      for(int ip = 0; ip <npp; ip++){
        int i = ordfac.s[ideg][ip][0];
        int j = ordfac.s[ideg][ip][1];
        int k = ordfac.s[ideg][ip][2];
        s[ideg][idx_invord_tab_fac(i,j,k)] = ip;
      }
    }
  }
  int s[1+METRIS_MAX_DEG_ORDERING][getnnod2(METRIS_MAX_DEG_ORDERING)] = {};
};

constexpr const invordfac_constructor invordfac;

struct invordtet_constructor{
  constexpr invordtet_constructor() : s(){
    s[0][0] = 0;
    for (int ideg = 1;ideg <= METRIS_MAX_DEG_ORDERING; ideg++){
      int npp = getnnod3(ideg);
      for(int ip = 0; ip < npp; ip++){
        int i = ordtet.s[ideg][ip][0];
        int j = ordtet.s[ideg][ip][1];
        int k = ordtet.s[ideg][ip][2];
        int l = ordtet.s[ideg][ip][3];
        ;
        s[ideg][idx_invord_tab_tet(i,j,k,l)] = ip;
      }
    }
  }
  int s[1+METRIS_MAX_DEG_ORDERING][getnnod3(METRIS_MAX_DEG_ORDERING)];
};
constexpr const invordtet_constructor invordtet;

//constexpr const int ordedg[1+METRIS_MAX_DEG_JACOBIAN][getnnod1(METRIS_MAX_DEG_JACOBIAN)][2];
//constexpr const int ordfac[1+METRIS_MAX_DEG_JACOBIAN][getnnod2(METRIS_MAX_DEG_JACOBIAN)][3] = {};
//constexpr const int ordtet[1+METRIS_MAX_DEG_JACOBIAN][getnnod3(METRIS_MAX_DEG_JACOBIAN)][4] = {};

// These are not accessed directly. Call instead mul2nod, mul2nod and mul2nod.
//const int invordedg[METRIS_MAX_DEG_JACOBIAN][getnnod1(METRIS_MAX_DEG_JACOBIAN)] = {};
//const int invordfac[METRIS_MAX_DEG_JACOBIAN][getnnod2(METRIS_MAX_DEG_JACOBIAN)] = {};
//const int invordtet[METRIS_MAX_DEG_JACOBIAN][getnnod3(METRIS_MAX_DEG_JACOBIAN)] = {};

// Bézier coefficients to evaluate isolated Bernstein polynomials.
// These are mathematically speaking integers, but they are destined to be multiplied to other doubles.
// Might as well save a type conversion or slower int x dbl instructions. 
// Computations are carried out as integers and conversion is carried out in the end. 


constexpr int mul2nod(int i,int j){
  return invordedg.s[i+j][i];
}
constexpr int mul2nod(int i,int j,int k){
  return invordfac.s[i+j+k][idx_invord_tab_fac(i,j,k)];
}
constexpr int mul2nod(int i,int j,int k,int l){
  return invordtet.s[i+j+k+l][idx_invord_tab_tet(i,j,k,l)];
}
inline int mul2nod(int n, const int* i){
  return n == 1 ? mul2nod(i[0],i[1]) :
         n == 2 ? mul2nod(i[0],i[1],i[2]) :
                  mul2nod(i[0],i[1],i[2],i[3]) ;
}

struct cbzedg_constructor{
  constexpr cbzedg_constructor() : s() {
    s[0][0] = 1;
    for(int ideg = 1; ideg <= METRIS_MAX_DEG_EVAL; ideg++){
      int npp = getnnod1(ideg-1);
      // Loop over previous degree multiply each by B_{e_1}+ .. and accumulate
//      printf("\nDebug ideg = {} \n",ideg);
      for(int i=0;i<npp;i++){
        int i1 = ordedg.s[ideg-1][i][0];
        int i2 = ordedg.s[ideg-1][i][1];
//        printf("Debug i = {} i1, i2 = {} {}, ",i,i1,i2);
        int irnk = mul2nod(1+i1,i2);
//        printf("irnk1 = {} ",irnk);
        s[ideg][irnk] += s[ideg-1][i];

        irnk = mul2nod(i1,1+i2);
//        printf("irnk2 = {} \n",irnk);
        s[ideg][irnk] += s[ideg-1][i];
      }
    }
  }
  double s[1+METRIS_MAX_DEG_EVAL][getnnod1(METRIS_MAX_DEG_EVAL)];
};
struct cbzfac_constructor{
  constexpr cbzfac_constructor() : s() {
    s[0][0] = 1;
    for(int ideg = 1; ideg <= METRIS_MAX_DEG_EVAL; ideg++){
      int npp = getnnod2(ideg-1);
      // Loop over previous degree multiply each by B_{e_1}+ .. and accumulate
      for(int i=0;i<npp;i++){
        int i1 = ordfac.s[ideg-1][i][0];
        int i2 = ordfac.s[ideg-1][i][1];
        int i3 = ordfac.s[ideg-1][i][2];

        int irnk = mul2nod(1+i1,i2,i3);
        s[ideg][irnk] += s[ideg-1][i];

        irnk = mul2nod(i1,1+i2,i3);
        s[ideg][irnk] += s[ideg-1][i];

        irnk = mul2nod(i1,i2,1+i3);
        s[ideg][irnk] += s[ideg-1][i];
      }
    }
  }
  double s[1+METRIS_MAX_DEG_EVAL][getnnod2(METRIS_MAX_DEG_EVAL)];
};
struct cbztet_constructor{
  constexpr cbztet_constructor() : s() {

    s[0][0] = 1;
    for(int ideg = 1; ideg <= METRIS_MAX_DEG_EVAL; ideg++){
      int npp = getnnod3(ideg-1);
      // Loop over previous degree multiply each by B_{e_1}+ .. and accumulate
      for(int i=0;i<npp;i++){
        int i1 = ordtet.s[ideg-1][i][0];
        int i2 = ordtet.s[ideg-1][i][1];
        int i3 = ordtet.s[ideg-1][i][2];
        int i4 = ordtet.s[ideg-1][i][3];

        int irnk = mul2nod(1+i1,i2,i3,i4);
        s[ideg][irnk] += s[ideg-1][i];

        irnk = mul2nod(i1,1+i2,i3,i4);
        s[ideg][irnk] += s[ideg-1][i];

        irnk = mul2nod(i1,i2,1+i3,i4);
        s[ideg][irnk] += s[ideg-1][i];

        irnk = mul2nod(i1,i2,i3,1+i4);
        s[ideg][irnk] += s[ideg-1][i];
      }
    }

  }
  double s[1+METRIS_MAX_DEG_EVAL][getnnod3(METRIS_MAX_DEG_EVAL)];
};


constexpr const cbzedg_constructor cbzedg;//[1+METRIS_MAX_DEG_JACOBIAN][getnnod1(METRIS_MAX_DEG_JACOBIAN)] = {};
constexpr const cbzfac_constructor cbzfac;//[1+METRIS_MAX_DEG_JACOBIAN][getnnod2(METRIS_MAX_DEG_JACOBIAN)] = {};
constexpr const cbztet_constructor cbztet;//[1+METRIS_MAX_DEG_JACOBIAN][getnnod3(METRIS_MAX_DEG_JACOBIAN)] = {};

//const double cbzedg[1+METRIS_MAX_DEG_JACOBIAN][getnnod1(METRIS_MAX_DEG_JACOBIAN)] = {};
//const double cbzfac[1+METRIS_MAX_DEG_JACOBIAN][getnnod2(METRIS_MAX_DEG_JACOBIAN)] = {};
//const double cbztet[1+METRIS_MAX_DEG_JACOBIAN][getnnod3(METRIS_MAX_DEG_JACOBIAN)] = {};

// Give it a multi-index, get rank in connectivity arrays. 
// This is necessary because inverse ordering goes through composition by a "trivial inverse order"
constexpr int mul2nod(int ideg,int i,int j);
constexpr int mul2nod(int ideg,int i,int j,int k);
constexpr int mul2nod(int ideg,int i,int j,int k,int l);

} // End namespace

#endif