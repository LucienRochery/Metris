//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "MetrisRunner.hxx"

#include "nlohmann/json.hpp"

namespace Metris{

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