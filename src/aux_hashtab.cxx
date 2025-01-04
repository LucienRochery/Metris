//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "aux_hashtab.hxx"

#include <tuple>
#include <absl/hash/hash.h>

namespace Metris{

	// User assumed to have called find already to check whether this is in the table or not
int FaceHashTable::insert(int i1, int i2, int i3, int ielem, int iprev){
	if(nlist >= mlist) return 1;

	if(iprev > 0){
		list[wdth*(iprev-1)+4] = nlist + 1;
	} else{
		int key = hashkey(i1,i2,i3);
		int inext = head[key]; 
		if(!inext){
			head[key] = nlist + 1;
		}else{
			int iprev = inext;
			while(inext){
				iprev = inext;
				inext = list[wdth*(inext-1)+4];
			}
	
			list[wdth*(iprev-1)+4] = nlist + 1;
		}
	}

	list[wdth*nlist + 0] = i1;
	list[wdth*nlist + 1] = i2;
	list[wdth*nlist + 2] = i3;

	list[wdth*nlist + 3] = ielem;

	list[wdth*nlist + 4] = 0;

	nlist++;

	return 0;
}

//   Return 0 if first with this hash
//        < 0 if not found, -index of where to put it
//        > 0 if     fount @index
int FaceHashTable::find(int i1, int i2, int i3, int *ilist){
	int key = hashkey(i1,i2,i3);
	int inext = head[key]; 
	*ilist = 0;
	if(!inext) return 0;

	int ifnd = 0;
	while(inext){
		int j1 = list[wdth*(inext-1)+0];
		int j2 = list[wdth*(inext-1)+1];
		int j3 = list[wdth*(inext-1)+2];
		if(i1 == j1 && i2 == j2 && i3 == j3){
			ifnd = inext;
			break;
		}
		ifnd  = -inext;
		inext = list[wdth*(inext-1)+4];
	}
	*ilist = ifnd;
	if(ifnd <= 0) return 0;
	return list[wdth*(ifnd-1)+3];
}


//   Return 0 if first with this hash
//        < 0 if not found, -index of where to put it
//        > 0 if     fount @index
//   If found, delete the element.
int FaceHashTable::find_and_rem(int i1, int i2, int i3, int *ilist){
	int key = hashkey(i1,i2,i3);
	int inext = head[key]; 
	*ilist = 0;
	if(!inext) return 0;

	int ifnd  = 0;
	int iprev = 0; 
	while(inext){
		int j1 = list[wdth*(inext-1)+0];
		int j2 = list[wdth*(inext-1)+1];
		int j3 = list[wdth*(inext-1)+2];
		if(i1 == j1 && i2 == j2 && i3 == j3){
			ifnd = inext;
			break;
		}
		iprev = inext; 
		ifnd  = -inext;
		inext = list[wdth*(inext-1)+4];
	}
	*ilist = ifnd;
	if(ifnd <= 0) return 0;
	inext = list[wdth*(ifnd-1)+4]; 
	if(iprev == 0){
		head[key] = inext;
	}else{
		list[wdth*(iprev-1)+4] = inext; 
	}
	return list[wdth*(ifnd-1)+3];
}


int FaceHashTable::hashkey(int i1,int i2,int i3){
//	int key = ((i1 ^ (i2 << 1)) >> 1) ^ (i3 << 1);
	int key = absl::Hash<std::tuple<int,int,int>>{} ({i1,i2,i3});
	return key%mhead;
}


} // End namespace
