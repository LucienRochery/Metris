
//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __FACT_MULIDX__
#define __FACT_MULIDX__


#include "ho_constants.hxx"
#include "utils/fact_pow.hxx"

#include <cstdint>
#include <array>

namespace Metris{



// Fastest types supporting at least required integer limit. Alternative is smallest: uint_leastX_t
typedef std::conditional<(METRIS_MAX_DEG_ORDERING<13),uint_fast32_t,uint_fast64_t>::type INT_FACT;

class fact_mulidx{

private:
  // Empty initializer list construction is value initialized (to 0). 
  // THIS IS NOT TRUE OF THE DEFAULT CONSTRUCTOR!
  fact_mulidx(){
    for(int idim = 1; idim <= 3; idim++){
      for(int ideg = 0; ideg <= METRIS_MAX_DEG_ORDERING; ideg++){
        for(int inode = 0; inode < getnnode(idim,ideg); inode++){
          cache[idim-1][ideg][inode] = 0;
        }
      }
    }
  }

  // Delete copy constructor and assignment operator to prevent copying
  fact_mulidx(const fact_mulidx&) = delete;
  fact_mulidx& operator=(const fact_mulidx&) = delete;

  std::array<std::array<std::array<INT_FACT,getnnod3(METRIS_MAX_DEG_ORDERING)>,
                                            METRIS_MAX_DEG_ORDERING+1>,
                                            3> cache;

public:
  //// Static method to get the singleton instance
  //static fact_mulidx& getInstance() {
  //  static fact_mulidx instance; // thread-safe in C++17
  //  return instance;
  //}

  static INT_FACT get(int idim, int ideg, int inode){
    METRIS_ASSERT(idim >= 1 && idim <= 3);
    METRIS_ASSERT(ideg >= 0 && ideg <= METRIS_MAX_DEG_ORDERING);
    METRIS_ASSERT_MSG(inode >= 0 && inode < getnnode(idim,ideg),
      "inode = "<<inode<<" >= 0? <? "<<getnnode(idim,ideg));

    static fact_mulidx instance; // thread-safe in C++17


    if(instance.cache[idim-1][ideg][inode] == 0){
      if(idim == 1){
        instance.cache[idim-1][ideg][inode] = ifact(ordedg.s[ideg][inode][0])
                                            * ifact(ordedg.s[ideg][inode][1]);
      }else if(idim == 2){
        instance.cache[idim-1][ideg][inode] = ifact(ordfac.s[ideg][inode][0])
                                            * ifact(ordfac.s[ideg][inode][1])
                                            * ifact(ordfac.s[ideg][inode][2]);
      }else{
        instance.cache[idim-1][ideg][inode] = ifact(ordtet.s[ideg][inode][0])
                                            * ifact(ordtet.s[ideg][inode][1])
                                            * ifact(ordtet.s[ideg][inode][2])
                                            * ifact(ordtet.s[ideg][inode][3]);
      }
      //printf("## fact mulidxcache miss {} {} {} value %llu \n",idim-1,ideg,inode,
      //  instance.cache[idim-1][ideg][inode]);
    }

    return instance.cache[idim-1][ideg][inode];
  }


};



} // namespace
#endif