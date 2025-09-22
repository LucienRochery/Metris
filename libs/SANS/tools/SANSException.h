// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANSEXCEPTION_H
#define SANSEXCEPTION_H

#include <exception>
#include <string>

#include <boost/exception/exception.hpp>
#include <boost/throw_exception.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/version.hpp>

#ifndef unlikely
#ifdef __GNUC__
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define unlikely(x) x
#endif
#endif

//=============================================================================
class SANSException : public boost::exception, public std::exception
{
protected:
  //Derived classes fill the errString with the error message
  std::string errString;
  std::string barrier;

private:
  mutable std::string outString;

public:
  virtual char const* what() const noexcept;

  SANSException();
  SANSException( const SANSException& e );
  SANSException& operator=( const SANSException& ) = delete;

  virtual ~SANSException() noexcept;
};

//=============================================================================
struct BackTraceException : public SANSException
{
  BackTraceException();
  BackTraceException( const BackTraceException& e) : SANSException(e) {}

  virtual ~BackTraceException() noexcept;
};

//=============================================================================
//Exception for assertions
struct AssertionException : public BackTraceException
{
  explicit AssertionException( const char *assertion );
  // cppcheck-suppress noExplicitConstructor
  AssertionException( const std::string& assertion );
  AssertionException( const std::string& assertion, const std::string& message );
  AssertionException( const char *assertion, const char *fmt, ...);
  AssertionException( const std::string& assertion, const char *fmt, ...);
  AssertionException( const AssertionException& e );

  virtual ~AssertionException() noexcept;
};

#define SANS_ASSERT( assertion ) \
  if ( unlikely(!(assertion)) ) \
     BOOST_THROW_EXCEPTION( AssertionException( BOOST_PP_STRINGIZE( assertion ) ) )

#define SANS_ASSERT_MSG( assertion, fmt... ) \
  if ( unlikely(!(assertion)) ) \
    BOOST_THROW_EXCEPTION( AssertionException( BOOST_PP_STRINGIZE( assertion ), fmt ) )

//=============================================================================
//Exception for development errors
struct DeveloperException : public BackTraceException
{
  // cppcheck-suppress noExplicitConstructor
  DeveloperException( const std::string& message );
  DeveloperException( const char *fmt, ...);

  virtual ~DeveloperException() noexcept;
};

#define SANS_DEVELOPER_EXCEPTION( msg... ) \
  BOOST_THROW_EXCEPTION( DeveloperException( msg ) )


//=============================================================================
//Exception for generic runtime errors that contain a simple message
struct RuntimeException : public SANSException
{
  // cppcheck-suppress noExplicitConstructor
  RuntimeException( const std::string& message );
  RuntimeException( const char *fmt, ...);

  virtual ~RuntimeException() noexcept;
};

#define SANS_RUNTIME_EXCEPTION( msg... ) \
  BOOST_THROW_EXCEPTION( RuntimeException( msg ) )

#if BOOST_VERSION >= 107300
#define SANS_BOOST_EXCEPTION_EXTERN( E ) \
  extern template void boost::throw_exception<E>(E const&, boost::source_location const &);

#define SANS_BOOST_EXCEPTION_INSTANTIATE( E ) \
  template void boost::throw_exception<E>(E const&, boost::source_location const &);

#else

#define SANS_BOOST_EXCEPTION_EXTERN( E ) \
extern template class boost::exception_detail::clone_impl<E>; \
extern template void boost::throw_exception<E>(E const&); \
extern template void boost::exception_detail::throw_exception_<E>(E const&, char const*, char const*, int);

#define SANS_BOOST_EXCEPTION_INSTANTIATE( E ) \
template class boost::exception_detail::clone_impl<E>; \
template void boost::throw_exception<E>(E const&); \
template void boost::exception_detail::throw_exception_<E>(E const&, char const*, char const*, int);
#endif

// Reduce compile time by explicitly instantiating these
SANS_BOOST_EXCEPTION_EXTERN(AssertionException)
SANS_BOOST_EXCEPTION_EXTERN(DeveloperException)
SANS_BOOST_EXCEPTION_EXTERN(RuntimeException)

#endif
