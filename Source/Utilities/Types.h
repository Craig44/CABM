/**
 * @file Types.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 22/10/2012
 * @section LICENSE
 *
 * Copyright NIWA Science ©2012 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This file contains conditional types that we'll use in CASAL2 instead
 * of the default types to give us more flexibility.
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef TYPES_H_
#define TYPES_H_

#include <cstdlib>
#include <memory>
#include <cxxabi.h>
#include <string>
#include <vector>


// Namespaces
namespace niwa {
namespace utilities {

/**
 * double conditional depending on if we're using auto differentiation or not
 */

#ifndef USE_AUTODIFF

#define AS_DOUBLE(x) x
typedef float Double;
#endif

typedef std::vector<std::vector<std::vector<std::vector<Double>>>> Vector4;

/**
 * This code is used to demangle the typeid(x).name information
 */
inline std::string demangle(const char* name) {
  int status = -1; // some arbitrary value to eliminate the compiler warning

  char   *realname;
  std::string val = "";
  realname = abi::__cxa_demangle(name, 0, 0, &status);
  val = realname;
  free(realname);

  val = val == "std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >" ? "string" : val;
  val = val.length() > 38 && val.substr(0, 38) == "std::vector<std::__cxx11::basic_string" ? "vector<string>" : val;
  val = val == "std::vector<double, std::allocator<double> >" ? "vector<double>" : val;
  val = val == "std::vector<unsigned int, std::allocator<unsigned int> >" ? "vector<unsigned int>" : val;
  return (status==0) ? val : name ;
}


} /* namespace utilities */
} /* namespace niwa */



#endif /* TYPES_H_ */
