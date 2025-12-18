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

  BadEntHandler(const int tdim_ = 2, double topXpctg = 50., double alphaImpr = 0.5)
  : tdim(tdim_), topX(topXpctg), alpha(alphaImpr), sizeK(0), nentt(0) {

      affectedEnttsAlive.reserve(100);
      neighbToAffected.reserve(100);
      seenEntts.reserve(100);
      deadEntts.reserve(100);
  }

  struct EntQual{
    int ientt;            // entity id
    double qentt;         // quality
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

    K.clear(); inK.clear(); inK.reserve(nentt);
    R.clear(); inR.clear(); inK.reserve(nentt);

    int placedInK = 0;
    for (int ii = 0; ii < nentt; ii++){

      int ientt = sortedIDs[ii];
      if (is_dead(ientt)) continue;

      double q = qua_at(ientt);
      EntQual entt{ientt, q};

      if (placedInK < sizeK){

        // goes to K

        // need to save iterator to later enforce uniqueness
        auto [it, inserted] = K.insert(entt);
        inK[ientt] = it;
        placedInK++;
      }
      else{

        // goes to R

        auto [it, inserted] = R.insert(entt);
        inR[ientt] = it;
        // R.push_back(entt);
      }
    }
    // create max-heap in R
    // std::make_heap(R.begin(), R.end(), CmpMaxQua{});

  }

  // to update K after a successful operation
  void updateK(const int icurrentEntt, const intAr2& ent2ent, intAr2& ent2tag, const int tag, const int ithrd){

    // std::cout << "Base number of alive entities = " << nentt << std::endl;

    if (!deadEntts.empty()) nentt += affectedEnttsAlive.size() - deadEntts.size();
    sizeK = std::max(1, (int)std::round(nentt * topX/100.));

    // std::cout << "Modyfing/Inserting " << affectedEnttsAlive.size() << " entities" << std::endl;

    purgeFrontKUpTo(icurrentEntt);

    for (const auto& entt: affectedEnttsAlive) ent2tag(ithrd,entt.first) = tag;

    for (const auto& entt : affectedEnttsAlive){

      const int ientt = entt.first;
      const double quaentt = entt.second;

      for (int ii = 0; ii < tdim+1; ii++){

        const int ienei = ent2ent(ientt,ii);

        if (ienei < 0 || is_dead(ienei) || ent2tag(ithrd,ienei) == tag) continue;

        if (auto it = seenEntts.find(ienei); it != seenEntts.end()){
          neighbToAffected[ienei] = seenEntts[ienei];
          seenEntts.erase(ienei);
          deadEntts.push_back(ienei); // being extra safe, I need to go back through all the logic again because I changed a few things
        }
      }
      insertEntt(ientt,quaentt);
    }
    affectedEnttsAlive.clear();

    // std::cout << "Killing " << deadEntts.size() << " entities" << std::endl;
    for (const auto entt : deadEntts) killEntt(entt);
    deadEntts.clear();

    for (const auto& entt : neighbToAffected) insertEntt(entt.first,entt.second);
    neighbToAffected.clear();

    rebalance();
  }

private:

  void purgeFrontKUpTo(const int ientt){

    auto findInK = inK.find(ientt);
    METRIS_ASSERT(findInK != inK.end());

    auto itKientt = findInK->second;

    std::vector<int> ids;
    ids.reserve(std::distance(K.begin(),itKientt));
    for (auto it = K.begin(); it != itKientt; it++){

      ids.push_back(it->ientt);
      seenEntts[it->ientt] = it->qentt;
    }

    for (const int id : ids) inK.erase(id);

    K.erase(K.begin(),itKientt);
  }

  // to insert modified/created entity (general)
  void insertEntt(int ientt, double qentt, const int insertionFlag = 0){

    METRIS_ASSERT(!is_dead(ientt));

    EntQual entt{ientt, qentt};

    if (insertionFlag > 0){
      auto [itNew, inserted] = K.insert(entt);
      inK[ientt] = itNew;
      return;
    }

    if (insertionFlag < 0){
      auto [itNew, inserted] = R.insert(entt);
      inR[ientt] = itNew;
      return;
    }

    // first thing is to remove the entity if already in our arrays
    killEntt(ientt);

    // insert fresh entity
    if (qentt > qualCutoff()){
      auto [itNew, inserted] = K.insert(entt);
      inK[ientt] = itNew;
    }
    else{
      auto [itNew, inserted] = R.insert(entt);
      inR[ientt] = itNew;
    }
  }

  // fetch best quality in K, which represent the cutoff
  double qualCutoff() const {
    return K.empty() ? -std::numeric_limits<double>::infinity() :  std::prev(K.end())->qentt;
  }

  void killEntt(int ientt){

    if (auto it = inK.find(ientt); it != inK.end()){

      K.erase(it->second);
      inK.erase(it);
    }

    if (auto it = inR.find(ientt); it != inR.end()){

      R.erase(it->second);
      inR.erase(it);
    }
  }

  // to rebalance K and R, i.e. keep the top X% (fresh) worst entities in K as the mesh grows, and the rest in R
  void rebalance(){

    // if K not big enough, fetch worst from R until desired size
    while ((int)K.size() < sizeK && !R.empty()){

      purgeFrontR();
      EntQual frontR = getFrontR();

      // I think we are 100% sure the entity is not in K already
      // but might want to call killEntt here to be safer

      insertEntt(frontR.ientt, frontR.qentt, 1);
    }

    // if K too big, send best in K to R
    while ((int)K.size() > sizeK){

      // remove best in K
      auto itBest = std::prev(K.end());
      EntQual bestInK = *itBest;
      killEntt(bestInK.ientt);

      // put it in R
      insertEntt(bestInK.ientt, bestInK.qentt,-1);
    }

    // std::cout << "New number total alive nentt = " << nentt << std::endl;
    // std::cout << "Size of K = " << K.size() << ". (topX = " << topX  << ")" << std::endl;
    // std::cout << "Size of R = " << R.size() << std::endl;
  }

  // to fetch (and pop) the front of R and re-heapify it
  EntQual getFrontR(){

    if (!R.empty()){

      auto it = R.begin();
      EntQual entt = *it;

      R.erase(it);
      inR.erase(entt.ientt);

      return entt;
    }

    return EntQual{-1, 0.};
  }

  // to make the front of R fresh
  void purgeFrontR(){
    while (!R.empty()){

      // check state of worst entt in R (its front)
      auto it = R.begin();
      const EntQual& entt = *it;
      bool dead = is_dead(entt.ientt);

      if (!dead) return; // front is valid, done

      R.erase(it);
      inR.erase(entt.ientt);
    }
  }

public:

  bool checkSuccess(const double quaNew, const double quaOld) const {

    // std::cout << "quaNew = " << quaNew << std::endl;
    // std::cout << "quaOld = " << quaOld << std::endl;
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
  // std::vector<EntQual> R;                                   // remainder of elements, max-heap with lazy insertions

  std::set<EntQual, CmpWorstFirst> R;
  std::unordered_map<int, SetIt> inR;

  std::unordered_map<int, double> seenEntts;

  // state info
  const double topX;                                        // percentage of entities we want in K
  int sizeK;                                                // current size of K (adapted as mesh grows)
  int nentt;

private:

  double alpha;                                             // percentage of improvement to consider success
  const int tdim;

public:

  std::unordered_map<int, double> affectedEnttsAlive;       // keep track of entities modified/created during successful operation
                                                            // feed externally during the operation itself
  std::unordered_map<int,double> neighbToAffected;

  std::vector<int> deadEntts;                               // keep track of entities killed during successful operation
                                                            // feed externally during the operation itself

};

} // namespace Metris

#endif // AUX_BADENTHANDLER
