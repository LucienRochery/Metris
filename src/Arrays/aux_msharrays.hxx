//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __SRC_BASIC_TYPES__
#define __SRC_BASIC_TYPES__


#ifndef ALWAYS_INLINE

// ALWAYS_INLINE is a macro to further encourage the compiler to inline a function

#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#else
#warning Not forcing inline with this compiler... (Please add this compiler to tools/always_inline.h)
#define ALWAYS_INLINE inline
#endif

#endif

#include "../aux_exceptions.hxx"

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <type_traits>

#include <climits>
#include <limits>

/*
  Rewriting array classes. See docs/Arrays for preliminary benchmarks and explanations. 
  Classes are templated to allow either an underlying type T*, or shared_ptr<T[]>
*/


namespace Metris{


template <typename T, typename INT1 = int>
class MeshArray1D{
public:

    MeshArray1D();
    MeshArray1D(INT1 m);
    //MeshArray1D(const MeshArray1D &x);
    MeshArray1D(INT1 n, T *a);
    MeshArray1D(INT1 n, const T *a);
    MeshArray1D(const std::initializer_list<T> & list);
    MeshArray1D(MeshArray1D &&cpy);
    void set_buffer(INT1 n, T *a);

    // Set n > 0 to copy first n elements and set n1 = n.
    bool allocate(INT1 m);
    void set_n(INT1 n){
      METRIS_ASSERT(n >= 0 && n <= m1); 
      n1 = n;
    }
    ALWAYS_INLINE INT1 get_n() const {
      return n1;
    }

    void free();

    void fill(INT1 m, T x);
    void fill(T x);

    void copyTo(MeshArray1D<T,INT1> &out, INT1 ncopy = -1) const;

    ~MeshArray1D();


    MeshArray1D<T,INT1>& operator=(const MeshArray1D &cpy); 
    MeshArray1D<T,INT1>& operator=(MeshArray1D &&mve); 
    

    void print(INT1 n) const;
    void print() const;

    operator const T*() const{return array;}
    operator  T*() {return array;}

    inline INT1 size() const {return m1;}
    inline INT1 size1() const {return m1;}

    ALWAYS_INLINE T &operator[](const INT1 &i){
        METRIS_ASSERT_MSG(i >= 0 && i < n1, "Array1D out of bounds (1) i = "
                   <<i<<" >= N = "<<n1<<" addr "<<array<<" m1 = "<<m1);
        return array[i];
    }
    ALWAYS_INLINE const T &operator[](const INT1 &i) const {
        METRIS_ASSERT_MSG(i >= 0 && i < n1, "Array1D out of bounds (2) i = "
                                          <<i<<" >= N = "<<m1<<" addr "<<array);
        return array[i];
    }

    void stack(T val);
    T pop();


    template<class AA> 
    class iterator{
    public:
      iterator(AA &arr_,INT1 idx_) : arr(arr_), idx(idx_){}
      bool operator != (iterator &other) const {
        return idx != other.idx;
      }

      typename std::conditional<std::is_const<AA>::value, const T&, T&>::type
      operator* () { 
        return arr[idx];
      } 
      const T& operator* () const { 
        return arr[idx];
      } 

      const iterator& operator++(){ 
        ++idx; 
        return *this; 
      } 

    private:
      AA &arr;
      INT1 idx; 
    };

    using iterator_mutbl = iterator<MeshArray1D<T,INT1>>;
    using iterator_const = iterator<MeshArray1D<T,INT1> const>;

    iterator_mutbl begin(){
      return iterator_mutbl(*this,0);
    }
    iterator_mutbl end(){
      return iterator_mutbl(*this,n1);
    }
    iterator_const begin() const{
      return iterator_const(*this,0);
    }
    iterator_const end() const{
      return iterator_const(*this,n1);
    }

public:
    const INT1 &n1_ = n1;

protected:
    INT1 m1, n1;
    // Owning or not 
    int ialloc; 
    T   *__restrict__ array;

};



template <typename T, typename INT1 = int, typename INT2 = int>
class MeshArray2D{
public:

    using INTL = typename std::conditional< (std::numeric_limits<INT1>::max()
                                            >std::numeric_limits<INT2>::max()),
                                             INT1, INT2>::type ;

    MeshArray2D();
    MeshArray2D(INT1 m, INT2 s); 
    MeshArray2D(INT2 s, const std::initializer_list<T> & list); 
    MeshArray2D(INT1 n, INT2 s, const T* ar); 
    MeshArray2D(INT1 n, INT2 s, T* ar);
    //MeshArray2D(const MeshArray2D &cpy);
    MeshArray2D(MeshArray2D &&cpy);

    // The additional argument is to ensure reallocation is deliberate
    bool allocate(INT1 m, INT2 s);
    void free();


    void set_n(INT1 n){
      METRIS_ASSERT(n >= 0 && n <= m1); 
      n1 = n;
    }
    ALWAYS_INLINE int get_n() const {
      return n1;
    }
    void inc_n();



    void fill(INT1 n, INT2 s, T x);
    void fill(T x);


    void copyTo(MeshArray2D<T,INT1,INT2> &out, INT1 ncopy = -1) const;

    void print(INT1 n) const;
    void print() const;



    ~MeshArray2D();

    MeshArray2D<T,INT1,INT2>& operator=(const std::initializer_list<T> & list);
    MeshArray2D<T,INT1,INT2>& operator=(const MeshArray2D &cpy); 
    MeshArray2D<T,INT1,INT2>& operator=(MeshArray2D &&cpy); 

    ALWAYS_INLINE MeshArray1D<T,INT1> get_arr1D(INT1 i){
      METRIS_ASSERT_MSG(i >= 0 && i < n1, "i out of bounds i = "<<i<<" n = "<<n1
        <<" m = "<<m1);
      return MeshArray1D<T,INT1>(stride,&array[i*stride]);
    }
    ALWAYS_INLINE const MeshArray1D<T,INT1> get_arr1D(INT1 i) const{
      METRIS_ASSERT_MSG(i >= 0 && i < n1, "i out of bounds i = "<<i<<" n = "<<n1);
      return MeshArray1D<T,INT1>(stride,&array[i*stride]);
    }

    #ifndef ARRAY2D_EXPENSIVE
      ALWAYS_INLINE T* operator[](INT1 i){
        METRIS_ASSERT_MSG(i >= 0 && i < n1, "i out of bounds i = "<<i<<" n = "<<n1
          <<" m = "<<m1);
        return &array[i*stride];
      }
      ALWAYS_INLINE const T* operator[](INT1 i) const{
        METRIS_ASSERT_MSG(i >= 0 && i < n1, "i out of bounds i = "<<i<<" n = "<<n1);
        return &array[i*stride];
      }
    #else
      ALWAYS_INLINE MeshArray1D<T> operator[](INT1 i){
        return get_arr1D(i);
      }
      ALWAYS_INLINE const MeshArray1D<T> operator[](INT1 i) const{
        return get_arr1D(i);
      }
    #endif


    ALWAYS_INLINE T& operator()(INT1 i, INT2 j){
      METRIS_ASSERT_MSG(i >= 0 && i < n1, "i = "<<i<<" out of bounds n = "<<n1);
      METRIS_ASSERT_MSG(j >= 0 && j < stride, "j = "<<j<<" out of bounds n = "<<stride);
      return array[i*stride + j];
    }
    ALWAYS_INLINE const T& operator()(INT1 i, INT2 j) const{
      METRIS_ASSERT_MSG(i >= 0 && i < n1, "i = "<<i<<" out of bounds n = "<<n1);
      METRIS_ASSERT_MSG(j >= 0 && j < stride, "j = "<<j<<" out of bounds n = "<<stride);
      return array[i*stride + j];
    }

    operator const T*() const{return array;}
    operator  T*() {return array;}

    inline INT2 get_stride() const {return stride;}

    void set_stride(INT2 s){
      stride = s;
      // Force reallocation if necessary
      set_n(n1);
    }
    inline INTL size() const {return nmemalc;}
    inline INT1 size1() const {return m1;}
    inline INT2 size2() const {return stride;}

protected:
    INT2 stride; 
    //int64_t nmemalc;
    INT1 m1, n1; 
    // Largest of both types 
    INTL nmemalc;   
    //int64_t nmemalc;
    int ialloc;
    T *__restrict__ array;
};



/* We need 3D arrays too.*/
template <typename T>
class MeshArray3D{
public:

    MeshArray3D(){
        nmemalc = 0;
        s1      = 0;
        s2      = 0;
        array   = NULL;
        ialloc  = 0;
        dbgid   = 0; 
    }
    MeshArray3D(int n,int s1_,int s2_){
        dbgid   = 0; 
        if(n < 0 || s1_ <= 0 || s2_ < 0) METRIS_THROW(WArgExcept());
        s1  = s1_;
        s2  = s2_;
        nmemalc = n*s1*s2;
        array   = new T[nmemalc];
        //printf("Called new 3[]\n");
        if(array == NULL) METRIS_THROW(DMemExcept());
        ialloc  = 1;
    }

    MeshArray3D(int n,int s1_,int s2_, T* buff){
        dbgid   = 0; 
        if(n < 0 || s1_ <= 0 || s2_ < 0) METRIS_THROW(WArgExcept());
        s1  = s1_;
        s2  = s2_;
        nmemalc = n*s1*s2;
        array   = buff;
        ialloc = 0;
        if(array == NULL) METRIS_THROW(DMemExcept());
    }

    MeshArray3D(int s1_, int s2_,const std::initializer_list<T> & list){
        dbgid   = 0; 
        if(s1_ <= 0 || s2_ <= 0) METRIS_THROW(WArgExcept());

        if(list.size()%(s2_*s1_) != 0) METRIS_THROW(WArgExcept());

        s1  = s1_;
        s2  = s2_;
        nmemalc = list.size();
        array   = new T[nmemalc];
        //printf("Called new 3[]\n");
        if(array == NULL) METRIS_THROW(DMemExcept());

        std::copy_n(list.begin(),list.size(),array);
        ialloc  = 1;
    }

    MeshArray3D(const MeshArray3D &cpy){
        dbgid   = 0; 
        ialloc = 0;
        array  = cpy.array;
        nmemalc= cpy.nmemalc;
        s1  = cpy.s1;
        s2  = cpy.s2;
    }

//  The only difference is we take control of the array.
    MeshArray3D(MeshArray3D &&cpy){
        dbgid   = 0; 
        ialloc = 1;
        cpy.ialloc = 0;
        array  = cpy.array;
        nmemalc= cpy.nmemalc;
        s1  = cpy.s1;
        s2  = cpy.s2;
    }
//    MeshArray3D(int m, int s, const T* a){
//        stride  = s;
//        nmemalc = m*stride;
//        array   = new T[nmemalc];
//        for(int i =0; i< m;i++){
//            for(int j =0; j<s; j++){
//                array[i*stride+j] = a[i*stride+j];
//            }
//        }
//    }
//

    //void set_dbgid(int dbgid){
    //  this->dbgid = dbgid;
    //}

    void allocate(int m,int s1_, int s2_){
        if(nmemalc > 0) METRIS_THROW(WArgExcept());
        if(m < 0 || s1_ <= 0 || s2_ <= 0) METRIS_THROW(WArgExcept());

        s1  = s1_;
        s2  = s2_;
        nmemalc = m*s1*s2;
        ialloc = 1;
        array = new T[nmemalc];
        //printf("Called new 3[]\n");
        if(array == NULL) METRIS_THROW(DMemExcept());
    }

    void fill(int m, T x){
        for(int i = 0;i < m; i++) array[i] = x;
    }

    void fill(T x){
        for(int i = 0;i < nmemalc; i++) array[i] = x;
    }
    

    void print(int n) const {
        int m = n < (nmemalc/s1/s2) ? n : (nmemalc/s1/s2);
        for(int i = 0; i < m; i++){
            std::printf("%4d: ",i);
            for(int j = 0; j < s1; j++){
                std::printf("%4d:",j);
                for(int k = 0; k < s2; k++){
                    std::printf(" %d",array[i*s1*s2+j*s2+k]);
                }
                std::printf("\n      ");
            }
            std::printf("\n");
        }
    }
    void print() const{
        this->print(nmemalc);
    }




    ~MeshArray3D(){
        nmemalc = 0;
        if(array != NULL && ialloc > 0) delete[] array;
        array = NULL;
    }

//    inline       T &operator[](int i)       {return array[i*stride];}
//    inline const T &operator[](int i) const {return array[i*stride];}
    
    MeshArray3D<T>& operator=(const std::initializer_list<T> & list){
        if(list.size()%(s1*s2) != 0) METRIS_THROW(WArgExcept());
        if(list.size() > nmemalc)    METRIS_THROW(WArgExcept());
        std::copy_n(list.begin(),list.size(),array);
        nmemalc = list.size();
        return *this;
    }

//    operator T*() const{return array;}
    inline  T &operator()(int i, int j, int k) {
        #ifdef ARRAY_CHECKS
            if(k >= s2) METRIS_THROW(WArgExcept());
            if(j >= s1) METRIS_THROW(WArgExcept());
            if(i*s1*s2 + j*s2 + k >= nmemalc) METRIS_THROW(WArgExcept());
        #endif
        return array[i*s1*s2 + j*s2 + k]; 
    }
    inline const T &operator()(int i, int j, int k) const {
        #ifdef ARRAY_CHECKS
            if(k >= s2) METRIS_THROW(WArgExcept());
            if(j >= s1) METRIS_THROW(WArgExcept());
            if(i*s1*s2 + j*s2 + k >= nmemalc) METRIS_THROW(WArgExcept());
        #endif
        return array[i*s1*s2 + j*s2 + k]; 
    }
//    inline const T &operator()(int i,int j) {
//        #ifdef OURDEBUG
//            if(j >= stride) throw std::invalid_argument("INDEX LARGER THAN STRIDE !");
//            if(i*stride + j >= nmemalc) throw std::invalid_argument("INDEX LARGER THAN SIZE !");
//        #endif
//        return array[i*stride + j]; 
//    }

//    inline int get_s1() const {return s1;}

//    inline T  *cptr(){return array;}//&array[(ielem-1)*stride + i - 1];}

//    inline int fitsFor(int i) const {return (stride*i < nmemalc);}

    void set_stride(int s1_, int s2_){
        if(s1_ <= 0 || s2_ <= 0) METRIS_THROW(WArgExcept());
        if(nmemalc%(s1_*s2_) != 0) METRIS_THROW(WArgExcept());
        s1 = s1_;
        s2 = s2_;
    }
    inline int size() const {return nmemalc;}

protected:
    int s1,s2;
    int64_t nmemalc;
    int ialloc;
    T *__restrict__ array;
    int dbgid;
};


} // End namespace

#endif