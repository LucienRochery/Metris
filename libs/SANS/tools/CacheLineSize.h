// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_CACHE_LINE_SIZE_H
#define METRIS_CACHE_LINE_SIZE_H

#ifndef METRIS_CACHE_LINE_SIZE
#define METRIS_CACHE_LINE_SIZE 64 //This should be fetched from the system with getconf LEVEL1_DCACHE_LINESIZE
#endif

#endif // METRIS_CACHE_LINE_SIZE_H
