//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __COMMON_SETUP_MSH__
#define __COMMON_SETUP_MSH__

#include <string>
#include "src/types.hxx"
#include "src/io_libmeshb.hxx"
#include "src/ho_constants.hxx"
#include "src/aux_topo.hxx"
#include "src/utils/aux_misc.hxx"
#include "src/low_topo.hxx"

#include "src/low_geo.hxx"
#include "src/linalg/matprods.hxx"
#include "src/msh_degelev.hxx"
#include "src/aux_exceptions.hxx"
#include "src/MetrisRunner/MetrisRunner.hxx"

#include "src/main_adap.hxx"
#include <fcntl.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp> 
#include <boost/timer/progress_display.hpp>

#include <sstream>

namespace Metris{

class scriptArrayString {
public:
  scriptArrayString() : ifirst(true), iinit(false) {}

  scriptArrayString(const std::string &name__) : ifirst(true){
    iinit = true;
    name_ = name__;
    setName(name__);
  }

  void setName(const std::string &name__){
    iinit = true;
    name_ = name__;
    oss << name__ << " = [";
  }

  void finish(){
    oss << "]";
  }

  std::string str() const {
    return oss.str();
  }

  std::string name() const {
    return name_;
  }

  template<typename T>
  scriptArrayString& operator+=(T val){
    METRIS_ASSERT(iinit);
    if(!ifirst){
      oss << ",";
    }else{
      ifirst = false;
    }
    oss << val;
    return *this;
  }
private:
  std::ostringstream oss;
  bool ifirst, iinit;
  std::string name_;
};

class MinMaxAvg {
public:
  MinMaxAvg() : min_(1.0e30), avg_(0), max_(-1.0e30), navg(0){}

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

private:
  double min_, avg_, max_;
  unsigned long long int navg;
  bool iinit;
};


template<typename T = double>
struct LinReg{
public:
  LinReg(int nn, T *x, T *y){
    T s_x = 0, s_y = 0, s_xx = 0, s_xy = 0;
    for(int ii = 0; ii < nn; ii++){
      s_x += x[ii];
      s_y += y[ii];
      s_xx += x[ii]*x[ii];
      s_xy += x[ii]*y[ii];
    }
    T det = s_xx*nn - s_x*s_x;
    if(abs((double) det) < 1.0e-30){
      printf("## linearRegression nan coeff a using n = %d sx = %e sy = %e\n",
        nn,(double)s_x,(double)s_y);
      printf(" s_xx %e s_xy %e\n",(double)s_xx,(double)s_xy);
      printf("x : ");
      MeshArray1D<T>(nn, x).print();
      printf("y : ");
      MeshArray1D<T>(nn, y).print();
      METRIS_THROW(GeomExcept());
    }

    slope = (nn*s_xy - s_x*s_y)/det;
    origin = (-s_x*s_xy + s_xx*s_y)/det;
  }
  LinReg(MeshArray1D<T> &x,MeshArray1D<T> &y) : LinReg(x.get_n(), &x[0], &y[0]) {}

  T slope;
  T origin;
};

double linearRegression(int nn, double *x, double *y){
  LinReg linreg(nn, x, y);
  return linreg.slope;
}

#if 0
template<class MetricFieldType>
struct 
MeshTestSetup{
	MeshTestSetup() = delete;
//	MeshTestSetup(std::string meshFileName, int tardeg = 1){
//
//  	int usrMaxDeg = METRIS_MAX_DEG, usrMinDeg = tardeg; // Empty values
//
//  // usrMaxDeg is the very maximum the user is allowing for storage. It is hard bounded by the constant METRIS_MAX_DEG
//  // usrMinDeg is the minimum degree the user wants. 
//  if(usrMaxDeg > METRIS_MAX_DEG){
//    std::cout<<"!! MAX DEG "<<usrMaxDeg<<" SET ABOVE HARD LIMIT OF "<<METRIS_MAX_DEG<<std::endl;
//    usrMaxDeg = METRIS_MAX_DEG;
//  }
//  if(usrMinDeg > METRIS_MAX_DEG) std::cout<<"!! TARGET MIN DEG "<<usrMinDeg<<" SET ABOVE HARD LIMIT OF "<<METRIS_MAX_DEG<<std::endl;
//  if(usrMinDeg > usrMaxDeg){
//    std::cout<<"!! TARGET MIN DEG "<<usrMinDeg<<" ABOVE MAX ADMISSIBLE "<<usrMaxDeg<<std::endl;
//    usrMinDeg = usrMaxDeg;
//  }
//
//  nt npoin  = 0, nedge = 0, nface = 0, nelem = 0;
//  nt maxDeg = 0; // This will be determined by the reader 
//
//  nt ilagMsh;
//  td::string CADName = "";
//  if( iniMesh(meshFileName , CADName, usrMinDeg , usrMaxDeg , msh) ) {
//    std::cout<<"## FAILED TO READ MESH"<<std::endl;
//    exit(1);
//  }
//
//
//
//  int ierro;
//     // This is a compile-time for loop
//  hana::while_(hana::less_equal.than(hana::int_c<METRIS_MAX_DEG>), 1_c, [&](auto ideg_c){
//  constexpr int ideg = ideg_c;
//  hana::while_(hana::less_equal.than(hana::int_c<METRIS_MAX_DEG>), ideg_c+1_c, [&](auto tdeg_c){
//    constexpr int tdeg = tdeg_c;
//    if(ideg == msh.curdeg && tdeg == usrMinDeg){
//      if(msh.ilag == 1){
//        deg_elevate<ideg,tdeg,1>(msh,  &ierro);
//      }else{
//        deg_elevate<ideg,tdeg,0>(msh,  &ierro);
//      }
//    }
//    return tdeg_c+1_c;});
//  return ideg_c+1_c;});
//
//
//  }

  // This one allows us to control the initializer more finely with the full range of options
  MeshTestSetup(int argc, char** argv) : run(argc,argv){
    iniFromArgs(argc,argv);
  }

  //MeshTestSetup(std::string meshFileName, int tardeg = 1){
  //  char **argv = (char**) malloc(256*sizeof(char*));
  //  int argc;
  //  
  //  gen_argv(&argc,argv,"-i "+meshFileName+" -tardeg "+std::to_string(tardeg)+" -nosort");

  //  iniFromArgs(argc,argv);

  //  for(int i = 0; i < argc; i++) free(argv[i]);
  //  free(argv);
  //}


  void iniFromArgs(int argc, char** argv){
    // In Release mode, get rid of all initialization prints.
    run.degElevate();

    //METRIS_ENFORCE( (run.runnerMetricFE         && std::is_same<MetricFieldType, MetricFieldFE        >::value
    //            || run.runnerMetricAnalytical && std::is_same<MetricFieldType, MetricFieldAnalytical>::value));



    msh = (Mesh<MetricFieldType>*) run.msh_g;


    //try{
    //  int ret = main_metris(argc,argv,msh);
    //}catch(const MetrisExcept &e){
    //  printf("## EXCEPTION CAUGHT BY iniFromArgs\n");

    //  std::cout<<"## Type: "<<e.what()<<std::endl;
  
    //  if(std::string const * ms=boost::get_error_info<excMessage>(e) )
    //    std::cout<<"## Message: "<<*ms; 
    //  if(boost::stacktrace::stacktrace const * tr=boost::get_error_info<excStackTrace>(e) )
    //    std::cerr << "## Call stack: \n" << *tr;
    //}
    //#ifdef NDEBUG
    //  fflush(stdout);
    //  dup2(stdout_fd, STDOUT_FILENO);
    //  close(stdout_fd);
    //#endif
  }

	~MeshTestSetup(){} // Mesh and MshConComps already handle the deallocation

	Mesh<MetricFieldType> *msh;
  MetrisRunner run;
};
#endif



} // end namespace

#endif