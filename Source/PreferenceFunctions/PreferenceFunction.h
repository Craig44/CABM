/*
 * PreferenceFunction.h
 *
 *  Created on: 21/12/2012
 *      Author: Admin
 */

#ifndef PREFERENCE_FUNCTION_H_
#define PREFERENCE_FUNCTION_H_

// Headers

#include "BaseClasses/Object.h"
#include "Model/Model.h"
#include "Utilities/Types.h"

// Namespaces
namespace niwa {

// Using
using std::map;

/**
 * Class Definition
 */
class PreferenceFunction : public niwa::base::Object {
public:
  // Methods
  PreferenceFunction() = delete;
  explicit PreferenceFunction(Model* model);
  virtual                     ~PreferenceFunction() = default;
  void                        Validate();
  virtual void                Build();
  void                        Reset();

protected:
  // pure methods
  virtual void                DoValidate() = 0;
  virtual void                DoBuild() = 0;
  // Members
  Model*                      model_ = nullptr;
};
} /* namespace niwa */
#endif /* PREFERENCE_FUNCTION_H_ */
