/**
 * @file Normal.cpp
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
#include "Normal.h"

#include "Utilities/DoubleCompare.h"
#include "Layers/Children/NumericLayer.h"


namespace niwa {
namespace preference_functions {

Normal::Normal(Model* model) : PreferenceFunction(model) {
  parameters_.Bind<float>(PARAM_MU, &mu_, "", "");
  parameters_.Bind<float>(PARAM_SIGMA, &sigma_, "", "");
}

/*
 * Validate
*/
void Normal::DoValidate() {

}


/*
 * Build
*/
void Normal::DoBuild() {


}


/*
 * Get result
*/
float Normal::get_result(float value) {
  return_value_ = pow(2.0,-((value - mu_) / sigma_ * (value - mu_)/ sigma_));
  return utilities::doublecompare::ZeroFun(pow(return_value_,alpha_), ZERO);
}
} /* namespace niwa */
} /* namespace preference_functions */
