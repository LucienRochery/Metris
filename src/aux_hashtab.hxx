//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __SRC_AUX_HASHTAB__
#define __SRC_AUX_HASHTAB__
#include <cstdlib>

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


//namespace std {
//
//	template <>
//	struct  hash<std::array<int, 2>>
//	{
//	  std::size_t operator()(const std::array<int, 2>& k) const
//	  {
//	    using std::size_t;
//     	using std::hash;
//	
//	    return ((k[0] ^ (k[1] << 1)) >> 1);
//	  }
//	};
//	
//	template <>
//	struct hash<std::array<int, 3>>
//	{
//	  std::size_t operator()(const std::array<int, 3>& k) const
//	  {
//	    using std::size_t;
//	    using std::hash;
//	
//	    return ((k[0] ^ (k[1] << 1)) >> 1) ^ (k[2] << 1);
//	  }
//	};
//
//
//	template <>
//	struct  hash<std::tuple<int, int>>
//	{
//	  std::size_t operator()(const std::tuple<int, int>& k) const
//	  {
//	    using std::size_t;
//     	using std::hash;
//	
//	    return ((std::get<0>(k) ^ (std::get<1>(k) << 1)) >> 1);
//	  }
//	};
//	
//	template <>
//	struct hash< std::tuple<int, int, int>>
//	{
//	  std::size_t operator()(const std::tuple<int, int, int>& k) const
//	  {
//	    using std::size_t;
//	    using std::hash;
//	
//	    return ((std::get<0>(k) ^ (std::get<1>(k) << 1)) >> 1) ^ (std::get<2>(k) << 1);
//	  }
//	};
//
//
//}
} // End namespace



#endif