//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_AUX_UTILS__
#define __METRIS_AUX_UTILS__


#include "inc_hana.hxx"

#include "aux_exceptions.hxx"


#include <tuple>
#include <cstdlib>
#include <utility>
#include <string>
#include <sstream>
#include <unistd.h>
#include <fcntl.h>



#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif


namespace Metris{


inline unsigned long long int getSysMem(){
  return sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGE_SIZE);
}

char* itoa(int value, char* result, int base);

void wait();

void print_ordering_tet(int ideg);
void print_ordering_fac(int ideg);

void print_EGADS_error(std::string fname, int ierro);



// sortupto8_dec implemented using optimal sorting networks
// to minimize comparisons. Swaps could be reduced by not calling the helpers. 

template<typename T>
inline void swi(T &x,T &y){
  T tmp = x;
  x = y;
  y = tmp;
}

template<typename T>
void sortupto8_dec(T *tab,int n);
template<typename T>
void sortupto8_dec(T *tab,int * idx,int n);
//template<typename T, int n>
//void sortupto8_dec(T tab[n],int idx[n]);

// Optimal sorting networks for up to 8 elements
// Sort as decreasing
template<typename T, int n>
inline void sortupto8_dec(T tab[n], int idx[n]){

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


template<typename T>
void sortupto8_inc(T *tab,int n);

// Gets 2,3 ints, returns sorted tuple (for hash keys)
std::tuple<int,int,int> stup3(int i1,int i2,int i3);
std::tuple<int,int>     stup2(int i1,int i2);
//template<int n>
//typename std::conditional<n == 2, std::tuple<int,int>, std::tuple<int,int,int>>::type
//stupn(const int *ii);

template<int n>
typename std::conditional<n == 2, std::tuple<int,int>, std::tuple<int,int,int>>::type 
stupn(const int *ii){
  static_assert(n == 2 || n == 3);
  if constexpr(n == 2){
    int key_[2] = {ii[0],ii[1]};
    sortupto8_dec(key_,n);
    return std::tuple<int,int>({key_[0],key_[1]});
  }else{
    int key_[3] = {ii[0],ii[1],ii[2]};
    sortupto8_dec(key_,n);
    return std::tuple<int,int,int>({key_[0],key_[1],key_[2]});
  }
}
//template std::tuple<int,int> stupn<2>(const int *ii);
//template std::tuple<int,int,int> stupn<3>(const int *ii);


//void generic_hist(int nlist,double *rlist,int *llist,double vlo,double vhi,
//                  char *vname,char *title,int iinter);
/* 
Experiments on compiler hand-holding
*/
template<int n>
inline double idpow(double x){
  if constexpr(n%2 == 0){
    return idpow<n/2>(x) * idpow<n/2>(x);
  }else{
    return idpow<n/2>(x) * idpow<n/2>(x) * x;
  }
}
template<> inline double idpow<0>(double){return 1.0;}

template<int n, class ftype>
inline ftype irpow(ftype x){
  if constexpr(n%2 == 0){
    return irpow<n/2,ftype>(x) * irpow<n/2,ftype>(x);
  }else{
    return irpow<n/2,ftype>(x) * irpow<n/2,ftype>(x) * x;
  }
}
template<> inline double irpow<0,double>(double){return 1.0;}


#ifdef USE_MULTIPRECISION
  template<> inline float4 irpow<0,float4>(float4){return 1.0;}
  template<> inline float8 irpow<0,float8>(float8){return 1.0;}

  template<int n>
  inline float4 id4pow(float4 x){
    if(n%2 == 0){
      return id4pow<n/2>(x) * id4pow<n/2>(x);
    }else{
      return id4pow<n/2>(x) * id4pow<n/2>(x) * x;
    }
  }
  template<> inline float4 id4pow<0>(float4){return 1.0;}


  template<int n>
  inline float8 id8pow(float8 x){
    if(n%2 == 0){
      return id8pow<n/2>(x) * id8pow<n/2>(x);
    }else{
      return id8pow<n/2>(x) * id8pow<n/2>(x) * x;
    }
  }
  template<> inline float8 id8pow<0>(float8){return 1.0;}
#endif


template<int n>
inline int iipow(int x){
  if(n%2 == 0){
    return iipow<n/2>(x) * iipow<n/2>(x);
  }else{
    return iipow<n/2>(x) * iipow<n/2>(x) * x;
  }
}
template<> inline int iipow<0>(int ){return 1;}


template<int n> 
inline constexpr int ifact(){
	return n*ifact<n-1>();
}
template<> inline constexpr int ifact<0>(){return 1;}

template<int ideg, int idim>
void gen_ordering_Vizir(int *ord);


inline void gen_argv(int *argc, char **argv, std::string cmd){
  *argc = 1;
  // cxxopts ignores first argument.
  argv[0] = (char *) malloc(4*sizeof(char));
  strncpy(argv[0],std::string("dum").c_str(),4);
  std::string str;
  std::stringstream cmd_(cmd);
  while(std::getline(cmd_, str, ' ')){
    int n = str.length();
    argv[*argc] = (char *) malloc((n+1)*sizeof(char));
    strncpy(argv[*argc],str.c_str(),n+1);
    (*argc)++;
  }
}


double linearRegression(int n, double *x, double *y);




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
                  const char *vname,const char *title,int iinter);




} // End namespace

#endif