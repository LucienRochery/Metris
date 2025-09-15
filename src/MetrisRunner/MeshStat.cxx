//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "MeshStat.hxx"

#include "nlohmann/json.hpp"

#include "fmt/ostream.h"

#include <cmath>
namespace Metris{

MeshStat::MeshStat(const MeshStat& rhs){
  *this = rhs;
}

void MeshStat::reset(){
  pctunit.set_n(0);
  len.set_n(0);
  quaP1.set_n(0);
  quaPk.set_n(0);

  tdim_max = 0;
}

void MeshStat::setLength(int tdim, const dblAr1& rlened){
  tdim_max = MAX(tdim_max,tdim);
  len.set_n(tdim_max);
  pctunit.set_n(tdim_max);

  int nunit = 0;
  for(double lened : rlened){
    if(lened >= 1.0/sqrt(2.0) && lened <= sqrt(2.0)) nunit++;
    len[tdim-1] += lened;
  }
  pctunit[tdim-1] = 100*((double)nunit) / rlened.get_n();
}

void MeshStat::setQuality(int tdim, const dblAr1& rqualel, AsDeg asdeg){
  tdim_max = MAX(tdim_max,tdim);
  MeshArray1D<MinMaxAvg>& qua = asdeg == AsDeg::P1 ? quaP1 : quaPk;
  qua.allocate(tdim_max);
  qua.set_n(tdim_max);
  for(double qualel : rqualel) qua[tdim-1] += qualel;
}


MeshStat& MeshStat::operator=(const MeshStat& rhs){

  rhs.pctunit.copyTo(this->pctunit);
  rhs.len.copyTo(this->len);
  rhs.quaP1.copyTo(this->quaP1);
  rhs.quaPk.copyTo(this->quaPk);
  this->tdim_max = rhs.tdim_max;

  return *this;
}

bool MeshStat::operator==(const MeshStat& rhs){
  for(int tdim = 1; tdim <= pctunit.get_n(); tdim++){
    if(std::abs(this->pctunit[tdim-1] - rhs.pctunit[tdim-1]) >= 1.0e-3) return false;
  }
  for(int tdim = 1; tdim <= len.get_n(); tdim++){
    if( std::abs(this->len[tdim-1].min() - rhs.len[tdim-1].min()) >= 1.0e-1
     || std::abs(this->len[tdim-1].avg() - rhs.len[tdim-1].avg()) >= 1.0e-2
     || std::abs(this->len[tdim-1].max() - rhs.len[tdim-1].max()) >= 1.0e-1 ) return false;
  }
  for(int tdim = 1; tdim <= quaP1.get_n(); tdim++){
    if( std::abs(this->quaP1[tdim-1].min() - rhs.quaP1[tdim-1].min()) >= 1.0e-3
     || std::abs(this->quaP1[tdim-1].avg() - rhs.quaP1[tdim-1].avg()) >= 1.0e-3
     || std::abs(this->quaP1[tdim-1].max() - rhs.quaP1[tdim-1].max()) >= 1.0e-3 ) return false;
  }
  for(int tdim = 1; tdim <= quaPk.get_n(); tdim++){
    if( std::abs(this->quaPk[tdim-1].min() - rhs.quaPk[tdim-1].min()) >= 1.0e-3
     || std::abs(this->quaPk[tdim-1].avg() - rhs.quaPk[tdim-1].avg()) >= 1.0e-3
     || std::abs(this->quaPk[tdim-1].max() - rhs.quaPk[tdim-1].max()) >= 1.0e-3 ) return false;
  }

  return true;
}

//void MeshStat::print(std::string name, FILE* logfile) const {
//  fmt::print(logfile, "-- Mesh stat summary {}:\n",name.c_str());
//  for(int tdim = 1; tdim <= len.get_n(); tdim++)
//    fmt::print(logfile, " - Length ({}D) : {:.2f}% unit w/ {:.2f} < l ~= {:.2f} < {:.2f} \n",tdim,pctunit[tdim-1],len[tdim-1].min(),len[tdim-1].avg(),len[tdim-1].max());
//  for(int tdim = 1; tdim <= quaP1.get_n(); tdim++)
//    fmt::print(logfile, " - Quality P1 ({}D) : {:.2e} < q ~= {:.2e} < {:.2e} \n",tdim,quaP1[tdim-1].min(),quaP1[tdim-1].avg(),quaP1[tdim-1].max());
//  for(int tdim = 1; tdim <= quaPk.get_n(); tdim++)
//    fmt::print(logfile, " - Quality Pk ({}D) : {:.2e} < q ~= {:.2e} < {:.2e} \n",tdim,quaPk[tdim-1].min(),quaPk[tdim-1].avg(),quaPk[tdim-1].max());
//}
template<class STREAM_T>
void MeshStat::print(std::string name, STREAM_T logfile) const {
  fmt::print(logfile, "-- Mesh stat summary {}:\n",name.c_str());
  for(int tdim = 1; tdim <= len.get_n(); tdim++)
    fmt::print(logfile, " - Length ({}D) : {:.2f}% unit w/ {:.2f} < l ~= {:.2f} < {:.2f} \n",tdim,pctunit[tdim-1],len[tdim-1].min(),len[tdim-1].avg(),len[tdim-1].max());
  for(int tdim = 1; tdim <= quaP1.get_n(); tdim++)
    fmt::print(logfile, " - Quality P1 ({}D) : {:.2e} < q ~= {:.2e} < {:.2e} \n",tdim,quaP1[tdim-1].min(),quaP1[tdim-1].avg(),quaP1[tdim-1].max());
  for(int tdim = 1; tdim <= quaPk.get_n(); tdim++)
    fmt::print(logfile, " - Quality Pk ({}D) : {:.2e} < q ~= {:.2e} < {:.2e} \n",tdim,quaPk[tdim-1].min(),quaPk[tdim-1].avg(),quaPk[tdim-1].max());
}

template void MeshStat::print<std::FILE*>(std::string name, std::FILE* logfile) const;
template void MeshStat::print<std::ostream&>(std::string name, std::ostream& logfile) const;
//std::ostream& operator<<(std::ostream& _os, const MeshStat& stat){
//  stat.print<std::ostream&>("", _os);
//  return _os;   
//}


void to_json(nlohmann::json& jj, const MeshStat& stat) {
  jj = nlohmann::json::object();
  for(int tdim = 1; tdim <= stat.pctunit.get_n(); tdim++){
    jj["pctunit" + std::to_string(tdim)] = stat.pctunit[tdim-1];
  }
  for(int tdim = 1; tdim <= stat.len.get_n(); tdim++){
    jj["minlen" + std::to_string(tdim)] = stat.len[tdim-1].min();
    jj["avglen" + std::to_string(tdim)] = stat.len[tdim-1].avg();
    jj["maxlen" + std::to_string(tdim)] = stat.len[tdim-1].max();
  }
  for(int tdim = 1; tdim <= stat.quaP1.get_n(); tdim++){
    jj["minquaP1_" + std::to_string(tdim)] = stat.quaP1[tdim-1].min();
    jj["avgquaP1_" + std::to_string(tdim)] = stat.quaP1[tdim-1].avg();
    jj["maxquaP1_" + std::to_string(tdim)] = stat.quaP1[tdim-1].max();
  }
  for(int tdim = 1; tdim <= stat.quaPk.get_n(); tdim++){
    jj["minquaPk_" + std::to_string(tdim)] = stat.quaPk[tdim-1].min();
    jj["avgquaPk_" + std::to_string(tdim)] = stat.quaPk[tdim-1].avg();
    jj["maxquaPk_" + std::to_string(tdim)] = stat.quaPk[tdim-1].max();
  }
}

void from_json(const nlohmann::json& jj, MeshStat& stat){
  stat.reset();
  for(int tdim = 1; tdim <= 3; tdim++){
    bool itdim = jj.contains("pctunit" + std::to_string(tdim))
               || jj.contains("minlen" + std::to_string(tdim))
               || jj.contains("avglen" + std::to_string(tdim))
               || jj.contains("maxlen" + std::to_string(tdim))
               || jj.contains("minquaP1_" + std::to_string(tdim))
               || jj.contains("avgquaP1_" + std::to_string(tdim))
               || jj.contains("maxquaP1_" + std::to_string(tdim))
               || jj.contains("minquaPk_" + std::to_string(tdim))
               || jj.contains("avgquaPk_" + std::to_string(tdim))
               || jj.contains("maxquaPk_" + std::to_string(tdim));
    if(!itdim) continue;

    bool iHO = jj.contains("minquaPk_" + std::to_string(tdim))
            || jj.contains("avgquaPk_" + std::to_string(tdim))
            || jj.contains("maxquaPk_" + std::to_string(tdim));
    // If we have one tdim entry, we should have them all
    // We're ready for json to throw
    stat.pctunit.set_n(tdim);
    stat.len.set_n(tdim);
    stat.quaP1.set_n(tdim);
    if(iHO) stat.quaPk.set_n(tdim);

    stat.pctunit[tdim-1] = jj.at("pctunit" + std::to_string(tdim)).get<double>();
    stat.len.allocate(tdim);
    stat.len[tdim-1].force_min(jj.at("minlen" + std::to_string(tdim)).get<double>());
    stat.len[tdim-1].force_avg(jj.at("avglen" + std::to_string(tdim)).get<double>());
    stat.len[tdim-1].force_max(jj.at("maxlen" + std::to_string(tdim)).get<double>());

    stat.quaP1[tdim-1].force_min(jj.at("minquaP1_" + std::to_string(tdim)).get<double>());
    stat.quaP1[tdim-1].force_avg(jj.at("avgquaP1_" + std::to_string(tdim)).get<double>());
    stat.quaP1[tdim-1].force_max(jj.at("maxquaP1_" + std::to_string(tdim)).get<double>());

    if(!iHO) continue;
    stat.quaPk[tdim-1].force_min(jj.at("minquaPk_" + std::to_string(tdim)).get<double>());
    stat.quaPk[tdim-1].force_avg(jj.at("avgquaPk_" + std::to_string(tdim)).get<double>());
    stat.quaPk[tdim-1].force_max(jj.at("maxquaPk_" + std::to_string(tdim)).get<double>());
  }
}


} // namespace Metris