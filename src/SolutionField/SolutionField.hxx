//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_SOLUTION_FIELD__
#define __METRIS_SOLUTION_FIELD__

#include "anasol.hxx"

#include "../metris_constants.hxx"
#include "../aux_exceptions.hxx"
#include "../linalg/explogmet.hxx"
#include "../low_eval.hxx"
#include "../msh_anamet.hxx"

#include <boost/pool/pool_alloc.hpp>


// Scalar solution fields: discrete (FE) and analytical

namespace Metris{


class MeshBase;
class SolutionFieldAnalytical;

enum class SolutionClass{None, SolutionFieldFE, SolutionFieldAnalytical};


class SolutionFieldBase{
public:
  virtual SolutionClass solutionClass() const { return SolutionClass::None; }

  SolutionFieldBase() = delete;
  SolutionFieldBase(const MeshBase &msh_) : msh(&msh_) {}

  // virtual incurs costs
  inline const double& operator[]([[maybe_unused]] int i) const {
    METRIS_THROW_MSG(WArgExcept(),"SolutionFieldBase::operator[] should not be called");
  }
  inline double& operator[]([[maybe_unused]] int i){
    METRIS_THROW_MSG(WArgExcept(),"SolutionFieldBase::operator[] should not be called");
  }

  SolutionFieldBase &operator=(const SolutionFieldBase& inp){
    msh = inp.msh;
    return *this;
  };

  // Derivatives of order 1, then 2, then ...
  // Maximum order 3 for now. Use sym2idx (i,j) and sym3idx (i,j,k) for ordering.
  double getSolBary([[maybe_unused]] int tdim, [[maybe_unused]] int ielem, 
                    [[maybe_unused]] const double* __restrict__ bary, 
                    [[maybe_unused]] std::initializer_list<double*> dfun = {}) const {
    METRIS_THROW_MSG(WArgExcept(),"SolutionFieldBase::getSolBary() should not be called");
  }

  const MeshBase *msh;
};



class SolutionFieldAnalytical : public SolutionFieldBase {

public:
  SolutionClass solutionClass() const override { return SolutionClass::SolutionFieldAnalytical; }

  SolutionFieldAnalytical() = delete;
  SolutionFieldAnalytical(const MeshBase &msh_) : 
    SolutionFieldBase::SolutionFieldBase(msh_), ianasol(-1), anasol(NULL){}

  void setAnalyticalSolution(int ianamet_);
  void setAnalyticalSolution(anasol_proto fptr);

  SolutionFieldAnalytical &operator=(const SolutionFieldAnalytical& inp);

  // AsDeg just to ensure compatible prototypes, unused option
  double getSolBary(int tdim, int ielem, 
                    const double* __restrict__ bary, 
                    std::initializer_list<double*> dfun = {}) const;

protected:
  int ianasol;
public:
  anasol_proto anasol;
  //double(*anasol)(void*,const double*__restrict__,std::initializer_list<double*>);
};






#if 0
class SolutionFieldFE : public SolutionFieldBase{
public:
  virtual SolutionClass SolutionClass() const { return SolutionClass::SolutionFieldFE; }

  SolutionFieldFE() = delete;
  SolutionFieldFE(MeshBase &msh_);

  void allocate();

  inline const double& operator[](int i) const {return rfld[i];}
  inline double& operator[](int i){return rfld[i];}

  SolutionFieldFE &operator=(const SolutionFieldFE& inp);

  FEBasis getBasis()const{return ibasis;}

  void setBasis(FEBasis ibasn){
    if(ibasn == getBasis()) return;
    METRIS_ASSERT(ibasn == FEBasis::Lagrange || ibasn == FEBasis::Bezier);

    if(ibasn == FEBasis::Lagrange) setLagrange();
    else                           setBezier();
  }

  // To be used in initialization
  void forceBasisFlag(FEBasis ibasn){ibasis = ibasn;}

  void readSolutionFile(std::string outname);
  void writeSolutionFile(std::string outname, bool iprefix = true);

  // Derivatives of order 1, then 2, then ...
  // Maximum order 3 for now. Use sym2idx (i,j) and sym3idx (i,j,k) for ordering.
  double getSolBary(AsDeg asdmet, DifVar idiff,
                    const int*__restrict__ ent2pol, 
                    int tdim,  const double* __restrict__ bary, 
                    double*__restrict__ dfun...);

public: 
  dblAr1 rfld;

protected:
  FEBasis ibasis;
  
  const MeshBase &msh;

  void setLagrange();
  void setBezier();

private:
  template<int ideg>
  double getSolBary0(DifVar idiff, MetSpace tarspac, const int*__restrict__ ent2pol, 
                     int tdim, const double*__restrict__ bary, 
                     double*__restrict__ dfun...);
};
#endif














} // End namespace





































#endif