//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "MetrisRunner.hxx"

#include "nlohmann/json.hpp"

namespace Metris{

bool MeshStat::operator==(const MeshStat& rhs){
  return abs(this->pctunit - rhs.pctunit) < 1.0e-3
      && abs(this->minlen  - rhs.minlen ) < 1.0e-1
      && abs(this->maxlen  - rhs.maxlen ) < 1.0e-1
      && abs(this->avglen  - rhs.avglen ) < 1.0e-2

      && abs(this->pctunit_bdry - rhs.pctunit_bdry) < 1.0e-3
      && abs(this->minlen_bdry  - rhs.minlen_bdry ) < 1.0e-1
      && abs(this->maxlen_bdry  - rhs.maxlen_bdry ) < 1.0e-1
      && abs(this->avglen_bdry  - rhs.avglen_bdry ) < 1.0e-2


      && abs(this->minqua - rhs.minqua ) < 1.0e-3
      && abs(this->maxqua - rhs.maxqua ) < 1.0e-3
      && abs(this->avgqua - rhs.avgqua ) < 1.0e-3 

      && abs(this->minqua_bdry - rhs.minqua_bdry ) < 1.0e-3
      && abs(this->maxqua_bdry - rhs.maxqua_bdry ) < 1.0e-3
      && abs(this->avgqua_bdry - rhs.avgqua_bdry ) < 1.0e-3 ;
}
void MeshStat::print(std::string name, FILE* logfile){
  fmt::print(logfile, "-- Mesh stat summary {}:\n",name.c_str());
  fmt::print(logfile, " - Length       : {:.2f}% unit w/ {:.2f} < l ~= {:.2f} < {:.2f} \n",pctunit,minlen,avglen,maxlen);
  fmt::print(logfile, " - Length (bdry): {:.2f}% unit w/ {:.2f} < l ~= {:.2f} < {:.2f} \n",pctunit_bdry,minlen_bdry,avglen_bdry,maxlen_bdry);
  fmt::print(logfile, " - Conf. err.       : {:.2e} < q ~= {:.2e} < {:.2e} \n",minqua,avgqua,maxqua);
  fmt::print(logfile, " - Conf. err. (bdry): {:.2e} < q ~= {:.2e} < {:.2e} \n",minqua_bdry,avgqua_bdry,maxqua_bdry);
}

void to_json(nlohmann::json& jj, const MeshStat& stat) {
  jj = nlohmann::json{{"pctunit", stat.pctunit}
                      ,{"minlen", stat.minlen}
                      ,{"maxlen", stat.maxlen}
                      ,{"avglen", stat.avglen}
                      ,{"pctunit_bdry", stat.pctunit_bdry}
                      ,{"minlen_bdry", stat.minlen_bdry}
                      ,{"maxlen_bdry", stat.maxlen_bdry}
                      ,{"avglen_bdry", stat.avglen_bdry}
                      ,{"minqua", stat.minqua}
                      ,{"avgqua", stat.avgqua}
                      ,{"maxqua", stat.maxqua}
                      ,{"minqua_bdry", stat.minqua_bdry}
                      ,{"avgqua_bdry", stat.avgqua_bdry}
                      ,{"maxqua_bdry", stat.maxqua_bdry}
                     };
}

void from_json(const nlohmann::json& jj, MeshStat& stat) {
  jj.at("pctunit").get_to(stat.pctunit);
  jj.at("minlen").get_to(stat.minlen);
  jj.at("maxlen").get_to(stat.maxlen);
  jj.at("avglen").get_to(stat.avglen);
  jj.at("pctunit_bdry").get_to(stat.pctunit_bdry);
  jj.at("minlen_bdry").get_to(stat.minlen_bdry);
  jj.at("maxlen_bdry").get_to(stat.maxlen_bdry);
  jj.at("avglen_bdry").get_to(stat.avglen_bdry);
  jj.at("minqua").get_to(stat.minqua);
  jj.at("avgqua").get_to(stat.avgqua);
  jj.at("maxqua").get_to(stat.maxqua);
  jj.at("minqua_bdry").get_to(stat.minqua_bdry);
  jj.at("avgqua_bdry").get_to(stat.avgqua_bdry);
  jj.at("maxqua_bdry").get_to(stat.maxqua_bdry);
}

} // namespace Metris