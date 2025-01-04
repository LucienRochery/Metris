// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "SANSException.h"

#include <iostream>
#include <sstream>
#include <algorithm> // std::all_of
#include <cstdlib>
#include <cstdarg>
#include <cctype>  // isspace

#include "demangle.h"

//#include <boost/stacktrace.hpp>

//#include <backtrace-supported.h> libbacktrace leaks memory... can't use it until it doesn't leak

#if BACKTRACE_SUPPORTED && defined(DEBUG_BACKTRACE) && !defined(__INTEL_COMPILER) && !defined(__clang__)
#include <backtrace.h>
#define PRINT_BACKTRACE


/* Passed to callbacks.  */

struct print_data
{
  typedef std::list< std::string > string_list;
  string_list& btsymbols;

  print_data(string_list& btsymbols) : btsymbols(btsymbols) {};
};


/* Print one level of a backtrace.  */

static int
print_callback (void *data, uintptr_t pc, const char *filename, int lineno,
                const char *function)
{
  struct print_data *pdata = (struct print_data *) data;
  char str[512] = {'\0'};

  static size_t funcnamesize = 512;
  char* funcname = (char*)malloc(funcnamesize);
  memset(funcname, '\0', funcnamesize);
  int status;
  char* ret = abi::__cxa_demangle(function, funcname, &funcnamesize, &status);

  if (status != 0)
  {
    free( funcname );
    funcname = const_cast<char*>(function);
  }
  else
  {
    //free( funcname );
    funcname = ret;
  }

  if ( funcname )
  {
    if ( "SANSException" != std::string(funcname) &&
         "DeveloperException" != std::string(funcname) &&
         "AssertionException" != std::string(funcname) )
    {
    //snprintf (str, 511, "0x%lx %s\n    %s:%d",
      snprintf (str, 511, "%s\n    %s:%d",
  //           (unsigned long) pc,
               funcname == NULL ? "???" : funcname,
               filename == NULL ? "???" : filename,
               lineno);
      pdata->btsymbols.push_back(str);
    }
  }

  if (status == 0)
    free( funcname );


  return 0;
}

/* Prmint errors to stderr.  */

static void
error_callback (void *data, const char *msg, int errnum)
{
  //struct print_data *pdata = (struct print_data *) data;

  //fprintf (stderr, "libbacktrace: %s", msg);
  //if (errnum > 0)
  //  fprintf (stderr, ": %s", strerror (errnum));
  //fputc ('\n', stderr);
}


#elif !defined(__CYGWIN__) && !defined(__MINGW32__) && !defined(_MSC_VER)
#define PRINT_BACKTRACE
#include <execinfo.h>
#endif

BackTraceException::BackTraceException()
{
  // if valgrind complains about this function the funcnamesize is likely too small

#if 0
  errString += "Backtrace Information:\n";

  boost::stacktrace::stacktrace st;

  errString += boost::stacktrace::to_string(st);

  errString += "\n";
  errString += barrier;
  errString += "\n";
#endif

#if BACKTRACE_SUPPORTED && defined(DEBUG_BACKTRACE) && !defined(__INTEL_COMPILER) && !defined(__clang__)
  errString += "Backtrace Information:\n";

  struct print_data data(btsymbols);
  int i = 0;

  struct backtrace_state *state = backtrace_create_state (NULL, BACKTRACE_SUPPORTS_THREADS,
                                        error_callback, NULL);

  i = backtrace_full (state, 0, print_callback, error_callback, &data);

  if (i != 0)
  {
    errString += "Failed to retreive backtrace information.";
    return;
  }

  errString += "\n";
  errString += barrier;
  errString += "\n";

#elif !defined(__CYGWIN__) && !defined(__MINGW32__) && !defined(_MSC_VER) && !defined(__APPLE__)

  // stacktrace.h (c) 2008, Timo Bingmann from http://idlebox.net/
  // published under the WTFPL v2.0
#define max_frames 30
  // storage array for stack trace address data
  void* addrlist[max_frames] = {NULL};

  // retrieve current stack addresses
  int addrlen = backtrace(addrlist, max_frames);

  if (addrlen == 0)
    return;

  errString += "Backtrace Information:\n";

  // resolve addresses into strings containing "filename(function+address)",
  // this array must be free()-ed
  char** symbollist = backtrace_symbols(addrlist, addrlen);

  // iterate over the returned symbol lines. skip the first, it is the
  // address of this function.
  int line = 0;
  for (int i = 1; i < addrlen; i++)
  {
    char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

    // find parentheses and +address offset surrounding the mangled name:
    // ./module(function+0x15c) [0x8048a6d]
    for (char *p = symbollist[i]; *p; ++p)
    {
      if (*p == '(')
        begin_name = p;
      else if (*p == '+')
        begin_offset = p;
      else if (*p == ')' && begin_offset)
      {
        end_offset = p;
        break;
      }
    }

    if (begin_name && begin_offset && end_offset
        && begin_name < begin_offset)
    {
      *begin_name++ = '\0';
      *begin_offset++ = '\0';
      *end_offset = '\0';

      // mangled name is now in [begin_name, begin_offset) and caller
      // offset in [begin_offset, end_offset). now apply
      // __cxa_demangle():
      std::string name = SANS::demangle(begin_name);

      bool whiteSpacesOnly = std::all_of(name.begin(),name.end(),isspace);
      if (whiteSpacesOnly) continue;

      if ( name.find("SANSException") == std::string::npos &&
           name.find("BackTraceException") == std::string::npos &&
           name.find("DeveloperException") == std::string::npos &&
           name.find("AssertionException") == std::string::npos )
      {
        errString += std::to_string(line++) + ") ";
        errString += " " + name + "\n";
      }
    }
    else
    {
      // couldn't parse the line? print the whole line.
      //fprintf(out, "  %s\n", symbollist[i]);
      std::string symbol( symbollist[i] );
      // Check if s consists only of whitespaces
      bool whiteSpacesOnly = std::all_of(symbol.begin(),symbol.end(),isspace);
      if (!whiteSpacesOnly)
      {
        errString += std::to_string(line++) + ") ";
        errString += " " + symbol + "\n";
      }
    }
  }

  free(symbollist);

  errString += "\n";
  errString += barrier;
  errString += "\n";
#endif
}

BackTraceException::~BackTraceException() noexcept {}
