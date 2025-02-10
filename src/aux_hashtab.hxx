//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __SRC_AUX_HASHTAB__
#define __SRC_AUX_HASHTAB__
#include <cstdlib>

#ifndef USE_ABSL
  #include <boost/functional/hash.hpp>
#endif
/*
Neighbour initialization takes 46s with -O3 using stl std::unordered_map on a 19.3M elt mesh
Let's see if we can speed it up by ditching pointers and tailoring to our case
Conclusion: yes, 45->29s. 
*/

namespace Metris{
	
class FaceHashTable{
public:
	FaceHashTable(int mhead_,int mlist_){
		mhead = mhead_ > 0 ? mhead_ : 1;
		mlist = mlist_;
		head = (int*) calloc (mlist,sizeof(int)); // This guy just points to first element in chained list
		wdth = 5; // this could be templated if indeed faster
		list = new int[wdth*mlist]; // i1 i2 i3 ielem inext
		nlist = 0;
	}
	~FaceHashTable(){
		free(head);
		delete[] list;
	}
	// User assumed to have called find already to check whether this is in the table or not
	int insert(int i1, int i2, int i3, int ielem, int iprev);

//   Return 0 if first with this hash
//        < 0 if not found, -index of where to put it
//        > 0 if     fount @index
	int find(int i1, int i2, int i3, int *ilist);


//   Return 0 if first with this hash
//        < 0 if not found, -index of where to put it
//        > 0 if     fount @index
//   If found, delete the element.
	int find_and_rem(int i1, int i2, int i3, int *ilist);


	int hashkey(int i1,int i2,int i3);


private:
	int *head, *list;
	int mhead,mlist,nlist;
	int wdth;
};

#ifndef USE_ABSL
namespace tup2_hash{
  struct hash
  {
    std::size_t operator()(const std::tuple<int,int>& key) const {
      std::size_t seed = 0;
      boost::hash_combine(seed, std::get<0>(key));
      boost::hash_combine(seed, std::get<1>(key));
      return seed;
    }
  };
}
namespace tup3_hash{
  struct hash
  {
    std::size_t operator()(const std::tuple<int,int,int>& key) const {
      std::size_t seed = 0;
      boost::hash_combine(seed, std::get<0>(key));
      boost::hash_combine(seed, std::get<1>(key));
      boost::hash_combine(seed, std::get<2>(key));
      return seed;
    }
  };
}
#endif


} // End namespace



#endif