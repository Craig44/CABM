/**
 * @file Logistic.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 26/7/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */
#ifndef LOGISTIC_PREFERENCE_FUNCTION_H_
#define LOGISTIC_PREFERENCE_FUNCTION_H_

// Local Headers
#include "PreferenceFunctions/PreferenceFunction.h"

namespace niwa {
namespace preference_functions {

//**********************************************************************
//
//
//**********************************************************************
class Logistic : public PreferenceFunction {
public:
  // Functions
  Logistic(Model* model);
  virtual                    ~Logistic() = default;
  void                       DoValidate() override final;
  void                       DoBuild() override final;

  // accessor
  float                     get_result(float value) override final;

protected:

  // Variables
  float                     ato95_;
  float                     a50_;

};

} /* namespace niwa */
} /* namespace preference_functions */


#endif /*LOGISTIC_PREFERENCE_FUNCTION_H_*/
