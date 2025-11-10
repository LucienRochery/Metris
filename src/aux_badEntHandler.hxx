//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AUX_BADENTHANDLER
#define AUX_BADENTHANDLER

#include <vector>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <unordered_map>
#include <functional>

#include "common_setup.hxx"

/*

Handler for keeping track of the top X% worst entities in mesh, call it K, and the remainder of entities, call it R.

--- K ---

  1- We need quick access to the worst entity in K, to fetch it and process it.
  2- We need quick access to the quality of the best entity in K, to decide if a new/updated entity belongs in K.

  The adopted solution is to use a std::set, with a sort rule from best to worst.
  This keeps fast access to both ends, and avoids O(n) shifts that would happen in a contiguous sorted array like std::vector.


  In particular, access/deleting either end entries is O(1).

  Additionally, whenever we erase or insert something,
  the sorting rule makes sure we end up with with worst and best at the ends.

  When inserting, the cost is O(log |K|).


---------

--- R ---

  1- We need quick access to the worst entity in R, to potentially promote it to K.
  2- We don't care about the rest

  The adopted solution is to use a max-heap by quality, i.e. worst at front, with access O(1).

  A max-heap is a container (e.g. std::vector) together with a heap rule, which sort it in a desired way.
  The heap rule for max-heap is arr[ii] >= arr[2*ii+1] and arr[ii] >= arr[2*ii+2].

  Example: the following array satifies max-heap rule:
           [10, 7, 9, 5, 6]

  This is weaker than a fully sorted array, and has always the largest entry at front.

  We use it with std::vector. And the following functions are used to enforce the heap-rule.

    - std::make_heap: Called once after creation of R, to sort it according to heap rule -- O(|R|)

    - std::pop_heap: This is used when we want to fetch the worst entity; the only one we care about.
                     What it does is to move the worst entity the end of the array, and apply heap rule to the rest. -- O(log |R|)
                     Then, we do R.pop_back() to actually fetch it.

    - std::push_heap: This enforces the heap rule -- O(log |R|)
                      It is used after we add an entity to R by calling R.push_back().
*/

namespace Metris{

struct BadEntHandler{

  BadEntHandler(double topXpctg) : topX(topXpctg), sizeK(0), nentt(0), currentGen(1) {}

  struct EntQual{
    int ient;            // entity id
    double qent;         // quality
    std::uint64_t genStamp; // generation when recorded
  };

  // multiset comparator for K
  struct CmpAscendingQua{
    bool operator()(const EntQual& entA, const EntQual& entB) const {
      if (entA.qent != entB.qent) return entA.qent < entB.qent;
      return entA.ient < entB.ient;
    }
  };

  // max-heap comparator for R (largest quality at front)
  struct CmpMaxQua { bool operator()(const EntQual& entA, const EntQual& entB) const { return entA.qent < entB.qent; } };

  // multiset for K
  std::multiset<EntQual, CmpAscendingQua> K; // top X% worst entities

  // max-heap by q for R
  std::vector<EntQual> R; // remainder

  // state info
  double topX;              // percentage of entities we want in K
  int sizeK;                // current size of K (adapted as mesh grows)
  int nentt;                // current number of entities
  std::uint64_t currentGen; // current generation (incremented after successful operation)

  // latest known generation for entity id
  std::unordered_map<int, std::uint64_t> latestGenStamp;

  // placeholder for getting quality of ent k
  std::function<double(int)> qua_at; // syntax std::function<RType(AType1,AType2,...,ATypeN)
                                     // means a callable returning type RType and taking N arguments of type AType1,... etc

  // placeholder for knowing if ent is dead
  std::function<bool(int)> is_dead;

  // to set the placeholders from the exterior of the struct
  void setCallbacks(std::function<double(int)> qua_atExt, std::function<bool(int)> is_deadExt){
    qua_at = std::move(qua_atExt);
    is_dead = std::move(is_deadExt);
  }

  // to update K after a successful operation
  void updateK(const int nenttNew){
    nentt = nenttNew;
    sizeK = std::max(1, (int)std::round(nentt * topX/100.));
    rebalance();
  }

  // construct K and R from sorted index list (worst quality first)
  void seedFromSortedIDs(const std::vector<int>& sortedIDs){
    METRIS_ASSERT(qua_at && is_dead); // check these have been set

    nentt = sortedIDs.size();
    sizeK = std::max(1, (int)std::round(nentt * topX/100.));

    K.clear(); R.clear(); latestGenStamp.clear();

    int placedInK = 0;
    for (int ii = 0; ii < nentt; ii++){

      int ient = sortedIDs[ii];
      if (is_dead(ient)) continue;

      double q = qua_at(ient);
      EntQual ent{ient, q, currentGen};
      latestGenStamp[ient] = currentGen;

      if (placedInK < sizeK){
        K.insert(ent);
        placedInK++;
      }
      else{
        R.push_back(ent);
      }
    }
    std::make_heap(R.begin(), R.end(), CmpMaxQua{});

    // if dead filtering left K smaller than sizeK, fill from R
    rebalance();
  }

  // fetch worst entity in K -- O(1)
  EntQual worstInK() const {
    auto it = std::prev(K.end());
    return *it;
  }

  // fetch best quality in K, which represent the cutoff -- O(1)
  double cutoff() const {
    return K.begin()->qent;
  }

  // to purge front of R
  void purgeFrontR(){
    while (!R.empty()){

      // check state of worst ent in R (its front)
      const EntQual& ent = R.front();
      const auto it = latestGenStamp.find(ent.ient);
      bool stale = it->second != ent.genStamp;
      bool dead = is_dead(ent.ient);

      if (!stale && !dead) return; // front is valid, done

      // else

      // move invalid worst ent to last slot and re-heapifies the rest -- O(log n)
      std::pop_heap(R.begin(), R.end(), CmpMaxQua{});
      // then remove invalid worst ent
      R.pop_back();
    }
  }

  // to fetch the front of R and re-heapify it
  EntQual getFrontR(){
    purgeFrontR();

    std::pop_heap(R.begin(), R.end(), CmpMaxQua{});
    EntQual ent = R.back();
    R.pop_back();

    return ent;
  }

  // to rebalance K and R, i.e. keep the top X% worst in K as the mesh grows, and the rest in R
  void rebalance(){

    // if K not big enough, fetch worst from R until desired size
    while ((int)K.size() < sizeK && !R.empty()) K.insert(getFrontR());

    // if K too big, send best in K to R
    while ((int)K.size() > sizeK){

      // remove best in K
      auto itBest = K.begin();
      EntQual bestInK = *itBest;
      K.erase(itBest);

      // put it in R
      R.push_back(bestInK);
      // reheapify
      std::push_heap(R.begin(), R.end(), CmpMaxQua{});
    }
  }

  // to fetch the front of K
  EntQual getFrontK(){
    if (K.empty()){
      rebalance(); // try to fill from R once
      if (K.empty()) return EntQual{-1, 0., 0};
    }

    auto itWorst = std::prev(K.end());
    EntQual ent = *itWorst;
    K.erase(itWorst);

    rebalance();
    return ent;
  }

  // to update quality in entity or insert new entity
  void insertEnt(int ient, double qent){

    // this works the same if ient exists already or not
    latestGenStamp[ient] = currentGen;

    // get the quality gate to decide if entity goes in K or R
    const double quaGate = cutoff();
    if (qent > quaGate){

      K.insert(EntQual{ient, qent, currentGen});
    }
    else{
      R.push_back(EntQual{ient, qent, currentGen});
      std::push_heap(R.begin(), R.end(), CmpMaxQua{});
    }

    rebalance();
  }
};

} // namespace Metris

#endif // AUX_BADENTHANDLER
