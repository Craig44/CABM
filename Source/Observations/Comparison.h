/**
 * @file Comparison.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 4/02/2014
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * << Add Description >>
 */
#ifndef OBSERVATIONS_COMPARISON_H_
#define OBSERVATIONS_COMPARISON_H_

// headers
#include <string>

#include "Utilities/Types.h"

// namespaces
namespace niwa {
namespace observations {

using std::string;
using utilities::Double;

struct Comparison {
  unsigned  age_ = 0;
  float    length_ = 0;
  float    expected_ = 0;
  float    simulated_ = 0;
  float    error_value_ = 0;
};

} /* namespace observations */
} /* namespace niwa */
#endif /* OBSERVATIONS_COMPARISON_H_ */
