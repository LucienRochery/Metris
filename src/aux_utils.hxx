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

inline int stdoutSilence(){
  fflush(stdout);
  int stdout_fd = dup(STDOUT_FILENO);
  int redir_fd = open("/dev/null", O_WRONLY);
  dup2(redir_fd, STDOUT_FILENO);
  close(redir_fd);
  return stdout_fd;
}

inline void stdoutRestore(int stdout_fd){
  fflush(stdout);
  dup2(stdout_fd, STDOUT_FILENO);
  close(stdout_fd);
}

#if 0
// https://stackoverflow.com/questions/21917529/is-it-possible-to-initialize-stdvector-over-already-allocated-memory
// To pass pre-allocated memory to std::vector. 
template <typename T>
class PreAllocator
{
    private:
        T* memory_ptr;
        std::size_t memory_size;

    public:
        typedef std::size_t     size_type;
        typedef T*              pointer;
        typedef T               value_type;

        PreAllocator(T* memory_ptr, std::size_t memory_size) : memory_ptr(memory_ptr), memory_size(memory_size) {}

        PreAllocator(const PreAllocator& other) throw() : memory_ptr(other.memory_ptr), memory_size(other.memory_size) {}

        template<typename U>
        PreAllocator(const PreAllocator<U>& other) throw() : memory_ptr(other.memory_ptr), memory_size(other.memory_size) {}

        template<typename U>
        PreAllocator& operator = (const PreAllocator<U>& other) { return *this; }
        PreAllocator<T>& operator = (const PreAllocator& other) { return *this; }
        ~PreAllocator() {}


        pointer allocate(size_type , const void* hint = 0) {return memory_ptr;}
        void deallocate(T* , size_type ) {}

        size_type max_size() const {return memory_size;}
};
#endif

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

//template<int n, typename T>
//inline int ipow(T x){
//  return ipow<n-1,T>(x)*x;
//}
//template<> inline int ipow<1,double>(double x){return x;}
//template<> inline int ipow<1,int>(int x){return x;}
//template<> inline int ipow<1,SurrealS<3,double>>(SurrealS<3,double> x){
//  return x;
//}
//

template<int n> 
inline constexpr int ifact(){
	return n*ifact<n-1>();
}
template<> inline constexpr int ifact<0>(){return 1;}

template<int ideg, int idim>
void gen_ordering_Vizir(int *ord);


//void gen_argv(int argc, char **argv, std::string s...){
//  va_list args;
//  va_start(args, s);
//  for(int i=0; i < argc; i++){
//    std::string str = va_arg(args, std::string);
//    argv[i] = (char *) malloc((str.length()+1)*sizeof(char));
//    strncpy(argv[i],str.c_str(),str.length()+1);
//  }
//  va_end(args);
//}


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




/*
// This version is intended to be used with hash tables
//     iinter = 1 : linear progression between min and max value
//     iinter = 2 ; geometric
//     nlist: input size
//     rlist: values
//     llist: indices ; for instance, if rlist pertains to edges, perhaps
//                      min and max should print llist(imin) and llist(imax)
//                      instead
//     if llist = NULL, simply index 0 ... nlist - 1
// T is any hash table type with an iterator that has a ->first and ->second double field
template<typename T>
void generic_hist(T &rlist,int *llist,double vlo,double vhi,
                  const char *vname,const char *title,int iinter){
  int nslot,r;
  const int mslot = 20;
  int cslot[mslot];
  double vslot[mslot];
  double vmin,vmax,vavg,dslot,v0;
  int imin,imax,nval,nlo,nhi;

  auto rlcur = rlist.begin();
  if(rlcur == rlist.end()){
    printf("## generic_hist: EMPTY HASH TABLE\n");
    return;
  }

  vmin = rlcur->second;
  vmax = rlcur->second;
  imin = 0;
  imax = 0;
  vavg = 0.0;
  nval = 0;
  nlo = 0;
  nhi = 0;
  int nlist = 0;
  int i=0;
  while(rlcur != rlist.end()){
    if(rlcur->second < vlo) nlo++;
    else if(rlcur->second > vhi) nhi++;
    else nval++;
    if(rlcur->second < vmin){
      imin = i;
      vmin = rlcur->second;
    }
    if(rlcur->second > vmax){
      imax = i;
      vmax = rlcur->second;
    }
    vavg += rlcur->second;
    nlist++;
  	rlcur++;
  	i++;
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
  for(rlcur = rlist.begin(); rlcur != rlist.end(); rlcur++){
    if(rlcur->second < vlo || rlcur->second > vhi) continue;
    int iskp = 0;
    for(int j=0;j<nslot;j++){
      if(rlcur->second > vslot[j]) continue;
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
  if(nval < nlist) printf(" - outside bounds n = %d\n",nlist-nval);

//  if( abs(vhi) > 1.0e6 || abs(vlo) > 1.0e6){
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
*/


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

/*
Credit for replace_at_helper and replate_at_c goes to 
https://stackoverflow.com/questions/61225367/boosthana-tuple-best-way-to-modify-a-value
*/
template <typename Xs, typename X, std::size_t ...before, std::size_t ...after>
constexpr auto replace_at_helper(Xs&& xs, X&&x, std::index_sequence<before...>,
                                      std::index_sequence<after...>) {
  return hana::make_tuple(
      hana::at_c<before>(std::forward<Xs>(xs))...,
      std::forward<X>(x),
      hana::at_c<after + sizeof...(before) + 1>(std::forward<Xs>(xs))...);
}
template <int  n>
constexpr auto replace_at_c = [](auto&& xs, auto&& x) {
    constexpr auto len = decltype(hana::length(xs))::value;
    return replace_at_helper(static_cast<decltype(xs)>(xs),
                             static_cast<decltype(x)>(x),
                             std::make_index_sequence<n>{},
                             std::make_index_sequence<len - n - 1>{});
};

// Convert hana::tuple<stuff> into an std::tuple<stuff>
// by using C++ template parameter deduction to get "stuff"
template<typename... Ts>
auto to_std_tuple(hana::tuple<Ts...> in){
  return std::tuple<Ts...>{};
}

template <typename T>
struct tuple_wrapper{
  tuple_wrapper(T t):tup(t){};

  // Integral constant, or how to pass a vlue by type
  // so we can bend template parameter deduction to our will
  // The lengths we go to to write rfld[i(_c)] !
  template<typename S>
  auto& operator[](S v){
    return std::get<(int)v>(tup);
  }
  template<typename S>
  const auto& operator[](S v) const{
    return std::get<(int)v>(tup);
  }
  T tup;
};


} // End namespace

#endif