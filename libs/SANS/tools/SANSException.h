// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_EXCEPTION_H
#define METRIS_EXCEPTION_H

#include <exception>
#include <string>

#include <boost/exception/exception.hpp>
#include <boost/throw_exception.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/version.hpp>

#ifndef METRIS_SUPPORT_UNLIKELY
#ifdef __GNUC__
#define METRIS_SUPPORT_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define METRIS_SUPPORT_UNLIKELY(x) x
#endif
#endif

namespace Metris
{

//=============================================================================
class MetrisException : public boost::exception, public std::exception
{
protected:
  //Derived classes fill the errString with the error message
  std::string errString;
  std::string barrier;

private:
  mutable std::string outString;

public:
  virtual char const* what() const noexcept;

  MetrisException();
  MetrisException( const MetrisException& e );
  MetrisException& operator=( const MetrisException& ) = delete;

  virtual ~MetrisException() noexcept;
};

//=============================================================================
struct BackTraceException : public MetrisException
{
  BackTraceException();
  BackTraceException( const BackTraceException& e) : MetrisException(e) {}

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

#define METRIS_SUPPORT_ASSERT( assertion ) \
  if ( METRIS_SUPPORT_UNLIKELY(!(assertion)) ) \
     BOOST_THROW_EXCEPTION( ::Metris::AssertionException( BOOST_PP_STRINGIZE( assertion ) ) )

#define METRIS_SUPPORT_ASSERT_MSG( assertion, fmt... ) \
  if ( METRIS_SUPPORT_UNLIKELY(!(assertion)) ) \
    BOOST_THROW_EXCEPTION( ::Metris::AssertionException( BOOST_PP_STRINGIZE( assertion ), fmt ) )

//=============================================================================
//Exception for development errors
struct DeveloperException : public BackTraceException
{
  // cppcheck-suppress noExplicitConstructor
  DeveloperException( const std::string& message );
  DeveloperException( const char *fmt, ...);

  virtual ~DeveloperException() noexcept;
};

#define METRIS_SUPPORT_DEVELOPER_EXCEPTION( msg... ) \
  BOOST_THROW_EXCEPTION( ::Metris::DeveloperException( msg ) )


//=============================================================================
//Exception for generic runtime errors that contain a simple message
struct RuntimeException : public MetrisException
{
  // cppcheck-suppress noExplicitConstructor
  RuntimeException( const std::string& message );
  RuntimeException( const char *fmt, ...);

  virtual ~RuntimeException() noexcept;
};

} // namespace Metris

#define METRIS_SUPPORT_RUNTIME_EXCEPTION( msg... ) \
  BOOST_THROW_EXCEPTION( ::Metris::RuntimeException( msg ) )

#if BOOST_VERSION >= 107300
#define METRIS_SUPPORT_BOOST_EXCEPTION_EXTERN( E ) \
  extern template void boost::throw_exception<E>(E const&, boost::source_location const &);

#define METRIS_SUPPORT_BOOST_EXCEPTION_INSTANTIATE( E ) \
  template void boost::throw_exception<E>(E const&, boost::source_location const &);

#else

#define METRIS_SUPPORT_BOOST_EXCEPTION_EXTERN( E ) \
extern template class boost::exception_detail::clone_impl<E>; \
extern template void boost::throw_exception<E>(E const&); \
extern template void boost::exception_detail::throw_exception_<E>(E const&, char const*, char const*, int);

#define METRIS_SUPPORT_BOOST_EXCEPTION_INSTANTIATE( E ) \
template class boost::exception_detail::clone_impl<E>; \
template void boost::throw_exception<E>(E const&); \
template void boost::exception_detail::throw_exception_<E>(E const&, char const*, char const*, int);
#endif

// Reduce compile time by explicitly instantiating these
METRIS_SUPPORT_BOOST_EXCEPTION_EXTERN(Metris::AssertionException)
METRIS_SUPPORT_BOOST_EXCEPTION_EXTERN(Metris::DeveloperException)
METRIS_SUPPORT_BOOST_EXCEPTION_EXTERN(Metris::RuntimeException)

#endif
