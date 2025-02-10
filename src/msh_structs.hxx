//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __MSH_STRUCTS__
#define __MSH_STRUCTS__

#include "types.hxx"
#include "metris_constants.hxx"

#include <boost/pool/poolfwd.hpp>
#ifdef USE_ABSL
#include <absl/container/flat_hash_map.h>
#endif
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>

#include <egads.h>

#include <tuple>
#include <vector>

//#define __MAX_ALLOC_REGIONS__ 64


namespace Metris{

/*
This class wraps around a boost::fast_pool_allocator to easily manage 
allocations and frees. The idea is that we have a boost::fast_pool_allocator
for the duration of the program, which is then passed to routines to quickly
"allocate" work arrays. These routines should then call free() on the memory pool
on exit. The destructor of RoutineWorkMemory handles these calls to free by 
keeping track of allocations. 
64 allocations per routine were put as an upper bound, this should be largely sufficient !
*/
// TODO: not threadsafe. We don't know if someone else has been requesting the 
// from the same chunk of memory, unless fast_pool_allocators are thread safe. 

// This can be a ressource hog for strange reasons. 
// In one case, the 2D swap was taking 98% of runtime because of a single allocate(cav.ncfac (=2)) call.. 
// On the other hand, size 100 allocates were perfectly fine. 
// Replacing that with a max(100,ncfac) made the run go from 2.7s to about 0.280s...
template<typename T>
class RoutineWorkMemory{
public:
	RoutineWorkMemory() = delete; 
	RoutineWorkMemory(boost::fast_pool_allocator<T> &pool_) : pool(pool_) {}
	// Free chunks from memory pool: this does not entail actual deallocation.
	~RoutineWorkMemory(){
    for(std::size_t i = 0; i < ptralc.size(); i++)
			pool.deallocate(ptralc[i],sizalc[i]);
	}

	T* allocate(int n){
		//if(nalloc >= __MAX_ALLOC_REGIONS__){
		//	printf("## INCREASE __MAX_ALLOC_REGIONS__ CONSTANT\n");
		//	exit(1);
		//}
    if(n == 0) return NULL;

    METRIS_ASSERT(n>=0);

		T* ptr = pool.allocate(n);
    ptralc.push_back(ptr);
    sizalc.push_back(n);
		//ptralc[nalloc] = ptr;
		//sizalc[nalloc] = n;
		//nalloc ++;
		return ptr;
	}
private:
	boost::fast_pool_allocator<T> &pool; 
	//int nalloc; // How many locations have been allocated
  std::vector<T*>  ptralc;
  std::vector<int> sizalc;

	//T*  ptralc[__MAX_ALLOC_REGIONS__]; // The pointers to the allocated regions
	//int sizalc[__MAX_ALLOC_REGIONS__]; // The sizes of the alloc'd regions
};



} // End namespace

#endif