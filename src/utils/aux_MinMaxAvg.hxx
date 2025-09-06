//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_UTILS_MINMAXAVG__
#define __METRIS_UTILS_MINMAXAVG__

#include <iostream>

namespace Metris{

class MinMaxAvg {
public:
  MinMaxAvg(){
    reset();
  }

  void reset(){
    min_ = 1.0e30;
    avg_ = 0.0;
    max_ = -1.0e30;
    navg = 0.0;
  }

  MinMaxAvg& operator+=(double val){
    if (val < min_) min_ = val;
    if (val > max_) max_ = val;
    avg_ += val;
    navg++;
    return *this;
  }

  double min() const {return min_;}
  double avg() const {return navg > 0 ? avg_/navg : 0;}
  double max() const {return max_;}

  bool operator<=(double value) const {
    return max_ <= value;
  }

  friend bool operator<=(double lhs, const MinMaxAvg& rhs) {
    return lhs <= rhs.max_;
  }

  bool operator<(double value) const {
    return max_ < value;
  }

  friend bool operator<(double lhs, const MinMaxAvg& rhs) {
    return lhs < rhs.max_;
  }


  bool operator>=(double value) const {
    return max_ >= value;
  }

  friend bool operator>=(double lhs, const MinMaxAvg& rhs) {
    return lhs >= rhs.max_;
  }

  bool operator>(double value) const {
    return max_ > value;
  }

  friend bool operator>(double lhs, const MinMaxAvg& rhs) {
    return lhs > rhs.max_;
  }

  friend std::ostream& operator<<(std::ostream& os, const MinMaxAvg& mma) {
    os << "min: " << mma.min_ << ", "
      << "avg: " << mma.avg() << ", "
      << "max: " << mma.max_ << ", "
      << "count: " << mma.navg;
    return os;
  }
  

  // Shouldn't be used, they're a hack.
  void force_min(double v) { min_ = v; }
  void force_avg(double v) { navg = 1; avg_ = v;}
  void force_max(double v) { max_ = v; }

private:
  double min_, avg_, max_;
  unsigned long long int navg;
};


} // namespace Metris


#endif