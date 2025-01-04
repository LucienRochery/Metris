//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __GEN_CCOEF__
#define __GEN_CCOEF__

#include <string>

void simpfrac(int x, int y, int *xs, int *ys);

void gen_ccoef();
//void eval_lag_func(const int* idx,
//                   StrExpr<'*'> &eval, StrExpr<'*'> dlag[]);
void gen_lageval();
void gen_lageval_alldim();

void gen_ccoeff3();
void gen_ccoeff2();

void gen_ccoeff2_d();
void getccoef2_map_coord(int ideg);


void gen_ccoef2_d_pt();
void get_point_derivatives(int ideg);

void gen_ccoef3_d();
void get_ccoeff3d(int ideg);


#endif