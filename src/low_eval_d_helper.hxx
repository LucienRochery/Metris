//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __LOW_EVAL_D_HELPER__
#define __LOW_EVAL_D_HELPER__

#include "types.hxx"
#include "ho_constants.hxx"
#include "utils/aux_misc.hxx"
#include "low_eval_d_SurrealS.hxx"
#include "low_eval_d_bezier.hxx"

#include "../libs/SANS/Surreal/SurrealS.h"
#include "../libs/SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
#include <boost/hana.hpp>


namespace Metris{


// Convert hana::tuple<stuff> into an std::tuple<stuff>
// by using C++ template parameter deduction to get "stuff"
template<typename... Ts>
auto to_std_tuple(hana::tuple<Ts...>){
  return std::tuple<Ts...>{};
}
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

template <typename T>
struct tuple_wrapper{
  tuple_wrapper(T t):tup(t){};

  // Integral constant, or how to pass a value by type
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




}//namespace
#endif