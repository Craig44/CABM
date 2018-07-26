/*
 * Constant.cpp
 *
 *  Created on: 9/01/2013
 *      Author: Admin
 */

#include "Constant.h"

#include <boost/math/distributions/lognormal.hpp>
#include <cmath>

#include "Model/Model.h"

namespace niwa {
namespace selectivities {

/**
 * Default Constructor
 */
Constant::Constant(Model* model)
: Selectivity(model) {

  parameters_.Bind<float>(PARAM_C, &c_, "C", "");

 // RegisterAsAddressable(PARAM_C, &c_);
}

/**
 * Return the constant result regardless of the
 * age or length specified
 *
 * @param age_or_length unsused in this selectivity
 * @return the constant value
 */
float Constant::GetResult(unsigned age_or_length) {
  return c_;
}


} /* namespace selectivities */
} /* namespace niwa */
