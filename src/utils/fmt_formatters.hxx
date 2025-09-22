//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_FMT_FORMATTERS__
#define __METRIS_FMT_FORMATTERS__

#include "fmt/format.h"

#include <egads.h>
#include "../types_arrays.hxx"
#include "../MetrisRunner/MeshStat.hxx"
#include "Surreal/SurrealS.h"

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& _os, const std::pair<T1, T2>& _p) {
  _os << '(' << _p.first << ',' << _p.second << ')';
  return _os;
}

//std::ostream& operator<<(std::ostream& _os, const ego& _p) {
//  _os << "egObject(mtype = " << _p->mtype << ", oclass = " << _p->oclass << ")";
//  return _os;
//}


#include "fmt/ostream.h"

//template<>
//struct fmt::formatter<ego> {
//
//  // Default, we disregard the format. See https://hackingcpp.com/cpp/libs/fmt.html#panel-fold64a section "Fully Custom Formatter"
//  constexpr auto parse(fmt::format_parse_context& ctx){
//    return ctx.begin();
//  }
//
//  template <typename FmtContext>
//  auto format(const ego &obj, FmtContext& ctx){
//    std::ostringstream oss;
//    oss << obj; 
//    return fmt::format_to(ctx.out(), "{}", oss.str());
//  }
//};

template<typename T1, typename T2>
struct fmt::formatter<std::pair<T1, T2>> : fmt::ostream_formatter {};

template<int N, typename T>
struct fmt::formatter<SANS::SurrealS<N, T>> : fmt::ostream_formatter {};

template<typename T, typename INT1>
struct fmt::formatter<Metris::MeshArray1D<T, INT1>> : fmt::ostream_formatter {};

template<typename T, typename INT1, typename INT2>
struct fmt::formatter<Metris::MeshArray2D<T, INT1, INT2>> : fmt::ostream_formatter {};

template<>
struct fmt::formatter<Metris::MinMaxAvg> : fmt::ostream_formatter {};

template<>
struct fmt::formatter<Metris::MeshStat> : fmt::ostream_formatter {};

#endif