//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __GEN_CCOEF__
#define __GEN_CCOEF__

#include <string>
#include <map>

#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"

void simpfrac(int x, int y, int *xs, int *ys);
void insert_map2(std::map<std::pair<int, int>, std::vector<SANS::DLA::VectorS<4, int>>> &d_ccoef_map,
            std::pair<int,int> &key, int up, int lo, int irnk1, int irnk2);

void gen_ccoef();
//void eval_lag_func(const int* idx,
//                   StrExpr<'*'> &eval, StrExpr<'*'> dlag[]);
void gen_lageval();
void gen_lageval_alldim();

void gen_ccoef3();
void gen_ccoef2();

void gen_ccoef2_d_coord();
void gen_ccoef2_d_coord0(int ideg);

void gen_ccoef3_d_coord();
void gen_ccoef3_d_coord0(int ideg);

void gen_ccoef2_d_pt();
void gen_ccoef2_d_pt0(int ideg);



#endif