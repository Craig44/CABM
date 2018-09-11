/**
 * @file CallbackBaranov.h
 * @author  C.Marsh
 * @version 1.0
 * @date 10/9/2018
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This is a callback function that i have tailed to solve the baranov catch equation for the process EfforBasedBaranov. Basically the call () will evaluate the catch equation.
 */
#ifndef MINIMISERS_GAMMADIFF_CALLBACK_BARANOV_H_
#define MINIMISERS_GAMMADIFF_CALLBACK_BARANOV_H_

// Headers
#include <vector>

#include "Processes/Children/Mortality.h"

// namespaces
namespace niwa {

namespace minimisers {
namespace gammadiff {

using std::vector;
using processes::Mortality;

/**
 * Class definition
 */
class CallbackBaranov {
public:
  CallbackBaranov(Mortality* mortality);
  virtual                     ~CallbackBaranov() = default;
  double                      operator()(const vector<double>& Parameters);

private:
  Mortality*                    mortality_;
};

} /* namespace gammadiff */
} /* namespace minimiser */
} /* namespace niwa */

#endif /* MINIMISERS_GAMMADIFF_CALLBACK_BARANOV_H_ */
