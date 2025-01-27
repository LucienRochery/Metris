//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "mprintf.hxx"

namespace Metris{
  
int DepthCounter::depth = 0;

#define SPAC4 ' ', ' ', ' ', ' '
#define SPAC8 SPAC4, SPAC4 
#define SPAC16 SPAC8, SPAC8
#define SPAC32 SPAC16, SPAC16
#define SPAC64 SPAC32, SPAC32
char  DepthCounter::spaces[65] = {SPAC64, '\0'};

}