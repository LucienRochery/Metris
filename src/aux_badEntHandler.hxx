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
#include <set>
#include <cmath>
#include <utility>


/*

Handler for keeping track of the top X% worst entities in mesh, call it K, and the remainder of entities, call it R.

--- K ---

  1- We need quick access to the worst entity in K, to fetch it and process it.
  2- We need quick access to the quality of the best entity in K, to decide if a new/updated entity belongs in K.
  3- Must have only fresh and alive entities, because its size represent the true X% worst elements.
     So, we need to keep track of where the entities are in K, for deletion of them when they die during an operation.

  The adopted solution is to use a std::set, with a sort rule from worst to best.
  And, we use an std::unordered_map<int,set::iterator> to map entity to position in K.

  This keeps fast access to both ends, and relatively fast insertion/deletions, and avoids O(N) shifts that would happen in a contiguous sorted array like std::vector.

  Access/deleting either end entries is O(1).
  Eerase/insertion has cost O(log|K|).

---------

--- R ---

  1- We need quick access to the worst entity in R, to potentially promote it to K.
  2- We don't care about the rest. So, don't do deletions in R and instead let the elements there get outdated.
     When inserting an entity we mark it with a time stamp and set its latest time stamp in a different array.
     Then, we can check later if an entity in R is stale by comparing is time stamp with the latest stamp for that entity.
     We call this lazy insertion.

  The adopted solution is to use a max-heap by quality, i.e. worst at front, with access O(1).

  A max-heap is a container (e.g. std::vector) together with a heap rule, which sort it in a desired way.
  The heap rule for max-heap is arr[ii] >= arr[2*ii+1] and arr[ii] >= arr[2*ii+2].

  Example: the following array satisfies max-heap rule:
           [10, 7, 9, 5, 6]

  This is weaker than a fully sorted array, and has always the largest entry at front.

  We use it with std::vector. And the following functions are used to enforce the heap-rule.

    - std::make_heap: Called once after creation of R, to sort it according to heap rule -- O(|R|)

    - std::pop_heap: This is used when we want to fetch the worst entity; the only one we care about.
                     What it does is to move the worst entity to the end of the array, and apply heap rule to the rest. -- O(log |R|)
                     Then, we do R.pop_back() to actually fetch it.

    - std::push_heap: This enforces the heap rule -- O(log |R|)
                      It is used after we add an entity to R by calling R.push_back().

  If we want to promote the worst element in R to K, we purge the front of R until a fresh (time stamp == latest stamp) entity appears.

*/

namespace Metris{

class BadEntHandler{

public:

  BadEntHandler(double topXpctg = 50., double alphaImpr = 0.5)
  : topX(topXpctg), alpha(alphaImpr), sizeK(0), nentt(0), currentGen(1) {

      affectedEnttsAlive.reserve(100);
      deadEntts.reserve(100);
  }

  struct EntQual{
    int ientt;            // entity id
    double qentt;         // quality
    std::uint64_t genStamp; // generation when recorded
  };

private:

  // set comparator for K (set)
  struct CmpWorstFirst{
    bool operator()(const EntQual& entA, const EntQual& entB) const {
      if (entA.qentt != entB.qentt) return entA.qentt > entB.qentt;
      return entA.ientt < entB.ientt; // we need a tie breaker, otherwise two entities with same quality would be considered the same
    }                               // and in a set we have uniqueness
  };

  // max-heap comparator for R (largest quality at front; confusing because the comparator works in the opposite way here)
  struct CmpMaxQua { bool operator()(const EntQual& entA, const EntQual& entB) const { return entA.qentt < entB.qentt; } };

  // latest known generation for entity id
  std::unordered_map<int, std::uint64_t> latestGenStamp;

  // placeholder for getting quality of entt k
  std::function<double(int)> qua_at; // syntax std::function<RType(AType1,AType2,...,ATypeN)
                                     // means a callable returning type RType and taking N arguments of type AType1,... etc

  // placeholder for knowing if entt is dead
  std::function<bool(int)> is_dead;

public:

  // to set the placeholders from the exterior of the struct
  void setCallbacks(std::function<double(int)> qua_atExt, std::function<bool(int)> is_deadExt){
    qua_at = std::move(qua_atExt);
    is_dead = std::move(is_deadExt);
  }

  // construct K and R from sorted index list (worst quality first)
  void seedFromSortedIDs(const std::vector<int>& sortedIDs){
    METRIS_ASSERT(qua_at && is_dead); // check these have been set

    nentt = sortedIDs.size();
    sizeK = std::max(1, (int)std::round(nentt * topX/100.));

    K.clear(); R.clear(); latestGenStamp.clear(); inK.clear();
    R.reserve(nentt); latestGenStamp.reserve(nentt); inK.reserve(nentt);

    int placedInK = 0;
    for (int ii = 0; ii < nentt; ii++){

      int ientt = sortedIDs[ii];
      if (is_dead(ientt)) continue;

      double q = qua_at(ientt);
      EntQual entt{ientt, q, currentGen};
      latestGenStamp[ientt] = currentGen;

      if (placedInK < sizeK){

        // goes to K

        // need to save iterator to later enforce uniqueness
        auto [it, inserted] = K.insert(entt);
        inK[ientt] = it;
        placedInK++;
      }
      else{

        // goes to R

        R.push_back(entt);
      }
    }
    // create max-heap in R
    std::make_heap(R.begin(), R.end(), CmpMaxQua{});

  }

  // to update K after a successful operation
  void updateK(const int nenttNew){

    // std::cout << "Base number of alive entities = " << nentt << std::endl;
    currentGen++;

    if (!deadEntts.empty()) nentt += affectedEnttsAlive.size() - deadEntts.size();
    sizeK = std::max(1, (int)std::round(nentt * topX/100.));

    // std::cout << "Modyfing/Inserting " << affectedEnttsAlive.size() << " entities" << std::endl;
    for (const auto& entt : affectedEnttsAlive) insertEntt(entt.first, entt.second);
    affectedEnttsAlive.clear();

    // std::cout << "Killing " << deadEntts.size() << " entities" << std::endl;
    for (const auto entt : deadEntts) killEnttK(entt);
    deadEntts.clear();

    rebalance();
  }

private:

  // to insert modified/created entity (general)
  void insertEntt(int ientt, double qentt){

    // first thing is to remove the entity if in K
    killEnttK(ientt);

    // then decide if fresh entity goes in K or R
    if (qentt > qualCutoff())  insertInK(ientt,qentt);
    else                       insertInR(ientt,qentt);
  }

  // fetch best quality in K, which represent the cutoff
  double qualCutoff() const {
    return K.empty() ? -std::numeric_limits<double>::infinity() :  std::prev(K.end())->qentt;
  }

  void insertInK(int ientt, double qentt){

    latestGenStamp[ientt] = currentGen;

    EntQual entt{ientt, qentt, currentGen};

    // insert fresh entity
    auto [itNew, inserted] = K.insert(entt);
    inK[ientt] = itNew;
  }

  void insertInR(int ientt, double qentt){

    // lazy insertion in R, we don't care if it is in R already
    // if so, the old version would just become stale (via genStamp)
    // the same logic does not apply in K because the size of K
    // must represent the top X% worst elements, all of them fresh and alive

    latestGenStamp[ientt] = currentGen;

    EntQual entt{ientt, qentt, currentGen};

    R.push_back(entt);
    std::push_heap(R.begin(), R.end(), CmpMaxQua{});
  }

  void killEnttK(int ientt){

    if (auto it = inK.find(ientt); it != inK.end()){

      K.erase(it->second);
      inK.erase(it);
    }
  }

  // to rebalance K and R, i.e. keep the top X% (fresh) worst entities in K as the mesh grows, and the rest in R
  void rebalance(){

    // if K not big enough, fetch worst from R until desired size
    while ((int)K.size() < sizeK && !R.empty()){

      purgeFrontR();
      EntQual frontR = getFrontR();

      // I think we are 100% sure the entity is not in K already
      // but might want to call killEnttK here to be safer

      insertInK(frontR.ientt, frontR.qentt);
    }

    // if K too big, send best in K to R
    while ((int)K.size() > sizeK){

      // remove best in K
      auto itBest = std::prev(K.end());
      EntQual bestInK = *itBest;
      killEnttK(bestInK.ientt);

      // put it in R
      insertInR(bestInK.ientt, bestInK.qentt);
    }

    // std::cout << "New number total alive nentt = " << nentt << std::endl;
    // std::cout << "Size of K = " << K.size() << ". (topX = " << topX  << ")" << std::endl;
    // std::cout << "Size of R = " << R.size() << std::endl;
  }

  // to fetch (and pop) the front of R and re-heapify it
  EntQual getFrontR(){

    if (!R.empty()){

      std::pop_heap(R.begin(), R.end(), CmpMaxQua{});
      EntQual entt = R.back();
      R.pop_back();

      return entt;
    }

    return EntQual{-1, 0., 0};
  }

  // to make the front of R fresh
  void purgeFrontR(){
    while (!R.empty()){

      // check state of worst entt in R (its front)
      const EntQual& entt = R.front();
      const auto it = latestGenStamp.find(entt.ientt);
      bool stale = it->second != entt.genStamp;
      bool dead = is_dead(entt.ientt);

      if (!stale && !dead) return; // front is valid, done

      // else

      // move invalid worst entt to last slot and re-heapifies the rest -- O(log |R|)
      std::pop_heap(R.begin(), R.end(), CmpMaxQua{});
      // then remove invalid worst entt
      R.pop_back();
    }
  }

public:

  bool checkSuccess(const double quaNew, const double quaOld) const {

    std::cout << "quaNew = " << quaNew << std::endl;
    std::cout << "quaOld = " << quaOld << std::endl;
    return quaNew <= ((1. - alpha/100.) * quaOld );
  }

  // useful if we want to iterate over several alpha values,
  // to not need a new handler from scratch
  void setAlpha(double alphaImprv) { alpha = alphaImprv; }

  // Members
  //-----------------------------------------------------------//

public:

  // set for K
  std::set<EntQual, CmpWorstFirst> K;                       // top X% worst entities, all fresh and alive (no lazy insertions)
  using SetIt = std::set<EntQual, CmpWorstFirst>::iterator; // to keep track of position of entities in K
  std::unordered_map<int, SetIt> inK;                       // map entt ID to iterator in K to keep track of "where" the entities are in K

private:

  // max-heap by quality for R
  std::vector<EntQual> R;                                   // remainder of elements, max-heap with lazy insertions

  // state info
  const double topX;                                        // percentage of entities we want in K
  int sizeK;                                                // current size of K (adapted as mesh grows)
  int nentt;                                                // current number of entities
  std::uint64_t currentGen;                                 // current generation (incremented after successful operation)

  double alpha;                                             // percentage of improvement to consider success

public:

  std::unordered_map<int, double> affectedEnttsAlive;       // keep track of entities modified/created during successful operation
                                                            // feed externally during the operation itself

  std::vector<int> deadEntts;                               // keep track of entities killed during successful operation
                                                            // feed externally during the operation itself

};

} // namespace Metris

#endif // AUX_BADENTHANDLER
