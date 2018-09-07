/**
 * @file DoubleNormal.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 26/7/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */
// Local Headers
#include "DoubleNormal.h"

#include "Utilities/DoubleCompare.h"

namespace niwa {
namespace preference_functions {

DoubleNormal::DoubleNormal(Model* model) : PreferenceFunction(model) {
  parameters_.Bind<float>(PARAM_MU, &mu_, "", "");
  parameters_.Bind<float>(PARAM_SIGMA_L, &sigma_l_, "", "");
  parameters_.Bind<float>(PARAM_SIGMA_R, &sigma_r_, "", "");

}

/*
 * Validate
*/
void DoubleNormal::DoValidate() {

}


/*
 * Build
*/
void DoubleNormal::DoBuild() {


}


/*
 * Get result
*/
float DoubleNormal::get_result(float value) {
  if (value < mu_) {
    return_value_ = pow(2.0,-((value - mu_) / sigma_l_ * (value - mu_) / sigma_l_ ));
  } else {
    return_value_ = pow(2.0,-((value - mu_) / sigma_r_ * (value - mu_) / sigma_r_ ));
  }
  LOG_FINEST() << "return value = " << return_value_ << " value = " << value << " sig l " << sigma_l_ << " sigma R = " << sigma_r_;
  return utilities::doublecompare::ZeroFun(pow(return_value_,alpha_), ZERO);
}

} /* namespace niwa */
} /* namespace preference_functions */

