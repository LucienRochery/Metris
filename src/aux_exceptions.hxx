//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_EXCEPTIONS__
#define __METRIS_EXCEPTIONS__


#include <iostream>
#include <string>

#include <exception>

#include "fmt/format.h"

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

#if 0
#define METRIS_THROW() throw(MetrisExcept());
#define METRIS_THROW_MSG(FMT__,...) {throw(MetrisExcept(fmt::format(FMT__, ##__VA_ARGS__)));}
#define METRIS_TRY0(x) try{x;}catch(const MetrisExcept &__e__){}

#define METRIS_ENFORCE(x) if(!(x)) {throw(MetrisAssert(#x,""));};
#define METRIS_ENFORCE_MSG(x,FMT__,...) if(!(x)) {throw(MetrisAssert(#x, fmt::format(FMT__, ##__VA_ARGS__)));};

#ifndef NDEBUG
  #define METRIS_ASSERT(x) METRIS_ENFORCE(x) ;
  #define METRIS_ASSERT_MSG(x,FMT__,...) METRIS_ENFORCE_MSG(x,FMT__, ##__VA_ARGS__) ;
#else
 #define METRIS_ASSERT(x) ;
 #define METRIS_ASSERT_MSG(x,FMT__,...) ;
#endif

#endif


#ifndef NO_EXCEPT_MESSAGES
  #define METRIS_THROW_MSG(...) throw Metris::MetrisExcept(__VA_ARGS__);
#else
  #define METRIS_THROW_MSG(...) throw Metris::MetrisExcept();
#endif
#define METRIS_THROW() throw Metris::MetrisExcept();
#define METRIS_ENFORCE(x) \
    if(!(x)) { throw Metris::MetrisAssert(#x, ""); }

// Using fmt::format is deadly when compiling with clang -fsanitize=address
#ifndef NO_EXCEPT_MESSAGES
  #define METRIS_ENFORCE_MSG(x, ...) \
      if(!(x)) { throw Metris::MetrisAssert(#x, fmt::format(__VA_ARGS__)); }
#else
  #define METRIS_ENFORCE_MSG(x, ...) METRIS_ENFORCE(x)
#endif
#ifndef NDEBUG
    #define METRIS_ASSERT(x) METRIS_ENFORCE(x)
    #define METRIS_ASSERT_MSG(x, ...) METRIS_ENFORCE_MSG(x, __VA_ARGS__)
#else
    #define METRIS_ASSERT(x) {};
    #define METRIS_ASSERT_MSG(x, ...) {};
#endif
#define METRIS_TRY0(x) try{x;}catch(const Metris::MetrisExcept &){}




namespace Metris{


//struct my_error: virtual boost::exception, virtual std::exception { }; //(2)

struct MetrisExcept: virtual std::exception{

  // Using mutable allows us to modify this within the const function what()
  mutable std::string formatted_message; // Cache survives function call
  std::string message;

  // Compile-time string litterals:
  template<size_t N>
  MetrisExcept(const char (&literal)[N]) : message(literal) {}

  // For runtime strings, avoid conversions between char* and std::string
  // using explicit. 
  explicit MetrisExcept(const char* runtime_msg) : message(runtime_msg == NULL ? "" : runtime_msg) {}
  explicit MetrisExcept(const std::string& runtime_msg) : message(runtime_msg) {}

  #ifndef NO_EXCEPT_MESSAGES
    // For formatted strings with compile-time validation
    template<typename... Args>
    MetrisExcept(const char* fmt, const Args&... args) 
        : message(fmt::format(fmt, args...)) {}
  #endif

  MetrisExcept() : MetrisExcept("") {}

  virtual const char* what() const noexcept override {
    // We can't build a string here as .c_str() will point to 
    // dead stack memory. Instead, we need to use our member variable.
    if(formatted_message.empty()){
      formatted_message = "\n## Metris exception.\n## Message: " + message;
    }
    return formatted_message.c_str(); 
  }
  
  
};

// For use with METRIS_ASSERT
struct MetrisAssert: public MetrisExcept
{
  MetrisAssert(std::string LoC, std::string message){
    this->message = "\n## Assertion failed: " + LoC
                  + "\n## Message: " + message;
  }
  MetrisAssert() = delete;

  virtual const char* what() const throw() {
    return message.c_str();
  }
};

} // End namespace

#endif