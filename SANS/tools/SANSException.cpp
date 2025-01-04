// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "SANSException.h"
#include <iostream>
#include <cstdarg>
#include <cstdio>

//=============================================================================
SANSException::SANSException() : errString(""), barrier("\n#=================================================#\n")
{
}

//=============================================================================
SANSException::SANSException( const SANSException& e ) :
  boost::exception(e), std::exception(e),
  errString(e.errString), barrier(e.barrier), outString(e.outString)
{
}

//=============================================================================
char const* SANSException::what() const noexcept
{
  outString  = barrier;
  outString += "#                   SANS ERROR                    #";
  outString += barrier;
  outString += "\n";

  outString += errString;

  outString += "\n";
  outString += barrier;
  outString += "\n";

  //boost test has a limit on the string size due to the use of vsprintf
  //Backtrace information can be much longer than the buffer, so just dump
  //it to the screen
  std::cerr << outString << std::flush;

  return "";
}

//=============================================================================
SANSException::~SANSException() noexcept {}

//=============================================================================
AssertionException::AssertionException(const std::string& assertion )
  : AssertionException(assertion.c_str())
{
}

AssertionException::AssertionException(  const char *assertion )
{
  errString += "Developer Error!!!\n";
  errString += "Assertion \"" + std::string(assertion) + "\" failed ";
}

//=============================================================================
AssertionException::AssertionException( const std::string& assertion, const std::string& message )
{
  errString += "Developer Error!!!\n";
  errString += "Assertion \"" + assertion + "\" failed ";
  errString += message;
}

//=============================================================================
AssertionException::AssertionException(const std::string& assertion, const char *fmt, ...)
 : AssertionException(assertion.c_str(), fmt)
{
}

//=============================================================================
AssertionException::AssertionException( const AssertionException& e )
 : BackTraceException(e)
{
}

//=============================================================================
AssertionException::AssertionException(const char *assertion, const char *fmt, ...)
{
  va_list args, args_copy ;

  try
  {
    // Dynamically allocate the size of the buffer
    va_start( args, fmt ) ;
    va_copy( args_copy, args ) ;

    const int sz = vsnprintf( nullptr, 0, fmt, args ) + 1;

    std::string buffer( sz, ' ' ) ;
    vsnprintf( &buffer.front(), sz, fmt, args_copy ) ;

    va_end(args_copy) ;
    va_end(args) ;

    errString += "Developer Error!!!\n";
      errString += "Assertion \"" + std::string(assertion) + "\" failed.\n";
    errString += buffer;
  }
  catch ( const std::bad_alloc& )
  {
    va_end(args_copy) ;
    va_end(args) ;

    errString += "Developer Error!!!\n";
    errString += "Failed to allocate memory for the error message";
  }
}

//=============================================================================
AssertionException::~AssertionException() noexcept {}


//=============================================================================
DeveloperException::DeveloperException( const std::string& message )
{
  errString += "Developer Error!!!\n";
  errString += message;
}

//=============================================================================
DeveloperException::DeveloperException(const char *fmt, ...)
{
  /*
  // Static sized buffer version found in most examples online
  char buffer[4096];
  va_list argp;

  va_start(argp, fmt);
  vsnprintf(buffer, 4096, fmt, argp);
  va_end(argp);
  */
  va_list args, args_copy ;

  try
  {
    // Dynamically allocate the size of the buffer
    va_start( args, fmt ) ;
    va_copy( args_copy, args ) ;

    const int sz = vsnprintf( nullptr, 0, fmt, args ) + 1;

    std::string buffer( sz, ' ' ) ;
    vsnprintf( &buffer.front(), sz, fmt, args_copy ) ;

    va_end(args_copy) ;
    va_end(args) ;

    errString += "Developer Error!!!\n";
    errString += buffer;
  }
  catch ( const std::bad_alloc& )
  {
    va_end(args_copy) ;
    va_end(args) ;

    errString += "Developer Error!!!\n";
    errString += "Failed to allocate memory for the error message";
  }
}

//=============================================================================
DeveloperException::~DeveloperException() noexcept {}

//=============================================================================
RuntimeException::RuntimeException( const std::string& message )
{
  errString += message;
}

//=============================================================================
RuntimeException::RuntimeException(const char *fmt, ...)
{
  /*
  // Static sized buffer version found in most examples online
  char buffer[4096];
  va_list argp;

  va_start(argp, fmt);
  vsnprintf(buffer, 4096, fmt, argp);
  va_end(argp);
  */
  va_list args, args_copy ;

  try
  {
    // Dynamically allocate the size of the buffer
    va_start( args, fmt ) ;
    va_copy( args_copy, args ) ;

    const int sz = vsnprintf( nullptr, 0, fmt, args ) + 1;

    std::string buffer( sz, ' ' ) ;
    vsnprintf( &buffer.front(), sz, fmt, args_copy ) ;

    va_end(args_copy) ;
    va_end(args) ;

    errString += buffer;
  }
  catch ( const std::bad_alloc& )
  {
    va_end(args_copy) ;
    va_end(args) ;

    errString += "Runtime Error!!!\n";
    errString += "Failed to allocate memory for the error message";
  }
}

//=============================================================================
RuntimeException::~RuntimeException() noexcept {}

//=============================================================================
SANS_BOOST_EXCEPTION_INSTANTIATE(AssertionException)
SANS_BOOST_EXCEPTION_INSTANTIATE(DeveloperException)
SANS_BOOST_EXCEPTION_INSTANTIATE(RuntimeException)
