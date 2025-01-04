// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef DENSEMATRIX_IDENTITY_H
#define DENSEMATRIX_IDENTITY_H

namespace SANS
{
namespace DLA
{

//This a datatype used to indicate that a matrix should be set to Identity
class Identity
{
public:
  Identity() : val_(1) {}
  operator int() const { return val_; }

  Identity operator-() { return Identity(-val_); }
  Identity operator+() { return *this; }

protected:
  explicit Identity(int val) : val_(val) {}

  int val_;
};

} //namespace DLA
} //namespace SANS


#endif //DENSEMATRIX_IDENTITY_H
