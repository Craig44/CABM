/**
 * @file DoubleNormal.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 26/7/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */
#ifndef DOUBLE_NORMAL_PREFERENCE_FUNCTION_H_
#define DOUBLE_NORMAL_PREFERENCE_FUNCTION_H_

// Local Headers
#include "PreferenceFunctions/PreferenceFunction.h"

namespace niwa {
namespace preference_functions {


class DoubleNormal : public PreferenceFunction {
public:
  // Functions
  DoubleNormal(Model* model);
  virtual                    ~DoubleNormal() = default;
  void                       DoValidate() override final;
  void                       DoBuild() override final;

  // accessor
  float                     get_result(float value) override final;

protected:

  // Variables
  float                     mu_;
  float                     sigma_l_;
  float                     sigma_r_;

};

} /* namespace niwa */
} /* namespace preference_functions */


#endif /*DOUBLE_NORMAL_PREFERENCE_FUNCTION_H_*/
