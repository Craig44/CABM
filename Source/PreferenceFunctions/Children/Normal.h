/**
 * @file Normal.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 26/7/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */
#ifndef NORMAL_PREFERENCE_FUNCTION_H_
#define NORMAL_PREFERENCE_FUNCTION_H_

// Local Headers
#include "PreferenceFunctions/PreferenceFunction.h"


namespace niwa {
namespace preference_functions {

class Normal : public PreferenceFunction {
public:
  // Functions
  Normal(Model* model);
  virtual                    ~Normal() = default;
  void                       DoValidate() override final;
  void                       DoBuild() override final;

  // accessor
  float                     get_result(float value) override final;

protected:

  // Variables
  float                     mu_;
  float                     sigma_;


};

} /* namespace niwa */
} /* namespace preference_functions */

#endif /*NORMAL_PREFERENCE_FUNCTION_H_*/
