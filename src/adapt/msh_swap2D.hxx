//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_SWAP2D__
#define __METRIS_MSH_SWAP2D__

#include "../Mesh/MeshFwd.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"


namespace Metris{


struct swapOptions{

  // By default, accept any swap that decreases maximum quality 
  swapOptions() = delete;

  swapOptions(int max_iter, int swap_norm, double swap_thres){
    this->max_iter   = max_iter;
    this->swap_norm  = swap_norm;
    this->swap_thres = swap_thres;  
  }

  swapOptions(const MetrisParameters &param){
    max_iter   = param.opt_swap_niter;
    swap_norm  = param.opt_swap_pnorm;
    swap_thres = param.opt_swap_thres;
  }

  // Maximum global swap iterations. Stops before if nothing to do. 
  int max_iter; 

  // lp norm to improve: governed by option -qswap-norm in module optim
  // - 0 is infinity (max). 
  // - 1+ is Lp 
  // - -1 is length-based: not accessible by option (adapt exclusive)
  int swap_norm;

  // Minimum increase in quality function to go through with a swap. Governed by
  // option -qswap-thres in optimization module (not adaptation). 
  double swap_thres;
};


namespace Defaults{
  //const swapOptions swapOptAdapt(100, 0, 0.0);
  const swapOptions swapOptAdapt(100, -1, 0.0);
}


template<class MetricFieldType, int gdim, int ideg>
double swap2D(Mesh<MetricFieldType> &msh, swapOptions opt, int *nswap, 
              int ithrd1 = 0, int ithrd2 = 1);


} // end namespace

#endif
