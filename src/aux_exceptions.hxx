//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_EXCEPTIONS__
#define __METRIS_EXCEPTIONS__


#include <iostream>
#include <string>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>

#include <exception>

/*
Exceptions powered by Boost. Use with METRIS_THROW for backtrace printing. 
Additional data can be passed by streaming a boost::error_info type. 
METRIS_THROW(something() << excStackTrace(boost::stacktrace::stacktrace()))
This is redundant with the macro, but other types can be used. 
See below for other types, or expand typedefs. 
Using the BOOST_THROW_EXCEPTION macro already guarantees what() is printed 
as well as the throw location (but not call stack). 

~~ Throw data:
	- excCallerName. Name of function throwing exception. 
~~ Exception types:
  - AlgoExcept: algorithm failed by no fault of the mesh, e.g. optimization. 
  Intended to be recoverable. 
  - WArgExcept: wrong arguments passed, e.g. negative value where > 0 expected. 
	- TopoExcept: Failure due to wrong integer or hash tables.
	- RealExcept: Real op failed to converge, be distinct from zero, etc. 
	- GeomExcept: mesh invalid in some sense. Negative Jacobian, etc. 
	- SMemExcept: some temp array is too small. Could be fixed by increasing a constant.
	- DMemExcept: dynamic allocation failed.


~~ Boost docs:
 - boost::exception https://www.boost.org/doc/libs/1_75_0/libs/exception/doc/boost-exception.html
for boost::exception docs. 
 - boost::throw_exception https://www.boost.org/doc/libs/1_83_0/libs/throw_exception/doc/html/throw_exception.html 
 - boost::stacktrace https://www.boost.org/doc/libs/1_65_0/doc/html/stacktrace/getting_started.html#stacktrace.getting_started.how_to_print_current_call_stack
 Print stacktrace using std::cout<<boost::stacktrace::stacktrace().
*/

//#define METRIS_STACKTRACE (excStackTrace(boost::stacktrace::stacktrace()))
// Disabled stacktrace useage
#ifndef NO_BOOST_EXCEPT
#define NO_BOOST_EXCEPT
#endif


#ifndef NO_BOOST_EXCEPT
  #include <boost/stacktrace.hpp>
  #define METRIS_THROW(x) \
   (boost::throw_exception(x << excStackTrace(boost::stacktrace::stacktrace()),\
    BOOST_CURRENT_LOCATION));
  #define METRIS_THROW_MSG(x,m) \
   (boost::throw_exception(x << excStackTrace(boost::stacktrace::stacktrace())\
                             << excMessage((std::stringstream()<<m<<"\n").str()),\
                           BOOST_CURRENT_LOCATION));
  #define METRIS_TRY0(x) try{x;}catch(const MetrisExcept &__e__){}
  #define METRIS_TRYFULL(x) try{x}catch(const Metris::MetrisExcept &e){\
        printf("\n\n## MAIN_METRIS THROWS EXCEPTION:\n");\
        std::cout<<"## Type: "<<e.what()<<std::endl;\
        if(std::string const * ms=boost::get_error_info<Metris::excMessage>(e) )\
          std::cout<<"## Message: "<<*ms; \
        if(boost::stacktrace::stacktrace const * tr=boost::get_error_info<Metris::excStackTrace>(e) )\
          std::cerr << "## Call stack: \n" << *tr;\
      }
#else
  #define METRIS_THROW(x) throw((x));
  #define METRIS_THROW_MSG(x,m) {std::cout<<m<<"\n";throw((x));}
  #define METRIS_TRY0(x) try{x;}catch(const MetrisExcept &__e__){}
  #define METRIS_TRYFULL(x) try{x}catch(const Metris::MetrisExcept &e){\
        printf("\n\n## MAIN_METRIS THROWS EXCEPTION:\n");\
        std::cout<<"## Type: "<<e.what()<<std::endl;\
      }
#endif


#define METRIS_ENFORCE(x) if(!(x)) {METRIS_THROW_MSG(MetrisEnforce(), #x)};
#define METRIS_ENFORCE_MSG(x,m) if(!(x)) {METRIS_THROW_MSG(MetrisEnforce(),"assert failed "<< #x<< "\n message: "<< m<<"\n")};

#ifndef NDEBUG
 //#define METRIS_ASSERT(x) if(!(x)) {METRIS_THROW(MetrisAssert())};
 //#define METRIS_ASSERT_MSG(x,m) if(!(x)) {METRIS_THROW_MSG(MetrisAssert(),m)};
  #define METRIS_ASSERT(x) METRIS_ENFORCE(x) ;
  #define METRIS_ASSERT_MSG(x,m) METRIS_ENFORCE_MSG(x,m) ;
#else
 #define METRIS_ASSERT(x) ;
 #define METRIS_ASSERT_MSG(x,m) ;
#endif



namespace Metris{

#ifndef NO_BOOST_EXCEPT
  typedef boost::error_info<struct tag_stackTrace,boost::stacktrace::stacktrace> 
  excStackTrace;
  typedef boost::error_info<struct tag_errorMessage,std::string> 
  excMessage;
#endif
//typedef boost::error_info<struct tag_int, int> excDbgInt;

//struct my_error: virtual boost::exception, virtual std::exception { }; //(2)

#ifndef NO_BOOST_EXCEPT
  struct MetrisExcept: virtual boost::exception, virtual std::exception{};
#else
  struct MetrisExcept: virtual std::exception{};
#endif

// For use with assert
struct MetrisAssert: public MetrisExcept
{
  virtual const char* what() const throw(){
    return "MetrisAssert: METRIS_ASSERT called.";
  }
};
// For use with assert
struct MetrisEnforce: public MetrisExcept
{
  virtual const char* what() const throw(){
    return "MetrisEnforce: METRIS_ASSERT called.";
  }
};

struct AlgoExcept: public MetrisExcept
{
  virtual const char* what() const throw(){
    return "AlgoExcept: Algorithm failure despite valid data.";
  }
};

struct WArgExcept: public MetrisExcept
{
  virtual const char* what() const throw(){
    return "WArgExcept: Incompatible arguments passed or invalid values.";
  }
};

struct TopoExcept: public MetrisExcept
{
  virtual const char* what() const throw(){
    return "TopoExcept: Corrupted topology or topologically invalid input mesh.";
  }
};

struct RealExcept: public MetrisExcept
{
  virtual const char* what() const throw(){
    return "RealExcept: Real computation failed (division by zero, failure to converge etc).";
  }
};

struct GeomExcept: public MetrisExcept
{
  virtual const char* what() const throw(){
    return "GeomExcept: Mesh geometrically invalid (negative Jacobian etc).";
  }
};

// Lacking stack memory. Typically mshell, etc, with constant sizes. 
// Change constant size. 
struct SMemExcept: public MetrisExcept
{
  virtual const char* what() const throw(){
    return "SMemExcept: A stack array is too small.";
  }
};

// Malloc fail
struct DMemExcept: public MetrisExcept
{
  virtual const char* what() const throw(){
    return "DMemExcept: Dynamic allocation failed.";
  }
};

struct TODOExcept: public MetrisExcept
{
  virtual const char* what() const throw(){
    return "TODOExcept: TODO: implement this case.";
  }
};

} // End namespace

#endif