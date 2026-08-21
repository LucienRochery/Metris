// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_DLA_MATRIX_UTIL_H
#define METRIS_DLA_MATRIX_UTIL_H

#if defined(METRIS_DLA_BLAS_ATLAS) || \
    defined(METRIS_DLA_BLAS_GOTO) || \
    defined(METRIS_DLA_BLAS_MKL) || \
    defined(METRIS_DLA_BLAS_ACCELERATE) || \
    defined(METRIS_DLA_BLAS_GENERIC)

#define METRIS_DLA_BLAS

#endif

#include "../../../tools/SANSnumerics.h" // Metris::Real

#include <cmath>   // abs
#include <utility> // std::move

//Matrix utility functions that take advantage of BLAS routines if they are available

namespace Metris
{
namespace DLA
{

template<class T, class Ta>
struct MatrixUtil_Native_Scal
{
  //Scale an array
  inline static void scal(T* __restrict x, const Ta& a, const int size)
  {
    //This must be multiplied from the left in case T and Ta are matrices
    for (int i = 0; i < size; ++i)
      x[i] = std::move(T(a*x[i]));
  }
};

template<class T>
struct MatrixUtil_Native_Scal<T, Metris::Real>
{
  //Scale an array
  inline static void scal(T* __restrict x, const Metris::Real& a, const int size)
  {
    for (int i = 0; i < size; ++i)
      x[i] *= a;
  }
};

template<class T, class Ta>
struct MatrixUtil_Native : public MatrixUtil_Native_Scal<T, Ta>
{
  //Swap two arrays
  inline static void swap(T* __restrict x, T* __restrict y, const int size)
  {
    for (int j = 0; j < size; ++j)
    {
      T tmp = x[j];
      x[j]  = y[j];
      y[j]  = std::move(tmp);
    }
  }

  // y = a*x + y
  inline static void axpy(const T* __restrict x, T* __restrict y, const Ta& a, const int size)
  {
    for (int i = 0; i < size; ++i)
      y[i] += a*x[i];
  }

  //Find the maximum value in a column
  inline static int max_row_in_col(const T* __restrict x, const int size, const int stride, const int start = 0)
  {
    //vera++ does not allow for "using namespace" in header files. However, it is ok here as it's in a local scope.
    //vera++ has been suppressed with scripts/vera/exclusions.txt
    using namespace std;
    T tmp = fabs(x[start*stride]);
    T abs_x;
    int iloc = start;
    for (int i = start+1; i < size; ++i)
    {
      abs_x = fabs(x[i*stride]);
      if ( tmp < abs_x )
      {
        tmp = abs_x;
        iloc = i;
      }
    }
    return iloc;
  }
};


//=============================================================================
// BLAS implementation
#ifdef METRIS_DLA_BLAS

template<class T>
struct MatrixUtil_BLAS;

template<>
struct MatrixUtil_BLAS<float>
{
  static void swap(float* __restrict x, float* __restrict y, const int size);

  static void scal(float* __restrict x, const float& a, const int size);

  static void axpy(const float* __restrict x, float* __restrict y, const float& a, const int size);

  static int max_row_in_col(const float* __restrict x, const int size, const int stride, const int start = 0);
};

template<>
struct MatrixUtil_BLAS<double>
{
  static void swap(double* __restrict x, double* __restrict y, const int size);

  static void scal(double* __restrict x, const double& a, const int size);

  static void axpy(const double* __restrict x, double* __restrict y, const double& a, const int size);

  static int max_row_in_col(const double* __restrict x, const int size, const int stride, const int start = 0);
};

#endif

template<class T, class Ta>
struct MatrixUtil : public MatrixUtil_Native<T,Ta> {};

#ifdef METRIS_DLA_BLAS

template<>
struct MatrixUtil<float, float> : public MatrixUtil_BLAS<float> {};
template<>
struct MatrixUtil<double, double> : public MatrixUtil_BLAS<double> {};

#endif

} //namespace DLA
} //namespace Metris

#endif // METRIS_DLA_MATRIX_UTIL_H
