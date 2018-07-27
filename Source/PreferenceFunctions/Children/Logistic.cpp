/**
 * @file Logistic.cpp
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
#include "Logistic.h"

#include "Utilities/DoubleCompare.h"

namespace niwa {
namespace preference_functions {

Logistic::Logistic(Model* model) : PreferenceFunction(model) {
  parameters_.Bind<float>(PARAM_A50, &a50_, "the a50 parameter of the logistic function", "");
  parameters_.Bind<float>(PARAM_ATO95, &ato95_, "the ato95 parameter of the logistic function", "");
}

/*
 * Validate
*/
void Logistic::DoValidate() {

}


/*
 * Build
*/
void Logistic::DoBuild() {


}


/*
 * Get result
*/
float Logistic::get_result(float value) {
  float temp = (a50_ - value) / ato95_;
  return_value_  = 1.0 / (1.0 + pow(19.0, temp));
  return utilities::doublecompare::ZeroFun(pow(return_value_,alpha_), ZERO);
}

} /* namespace niwa */
} /* namespace preference_functions */

