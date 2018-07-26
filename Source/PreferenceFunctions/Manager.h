/**
 * @file Manager.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 21/12/2012
 * @section LICENSE
 *
 * Copyright NIWA Science ©2012 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * The time class represents a moment of time.
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef PREFERENCE_FUNCTION_MANAGER_H_
#define PREFERENCE_FUNCTION_MANAGER_H_

// Headers
#include "BaseClasses/Manager.h"
#include "PreferenceFunctions/PreferenceFunction.h"

// Namespaces
namespace niwa {
namespace preference_functions {

/**
 * Class defintiion
 */
class Manager : public niwa::base::Manager<niwa::preference_functions::Manager, niwa::PreferenceFunction> {
  friend class niwa::base::Manager<niwa::preference_functions::Manager, niwa::PreferenceFunction>;
  friend class niwa::Managers;
public:
  // methods
  virtual                     ~Manager() noexcept(true) {};
  PreferenceFunction*                GetPreferenceFunction(const string& label);

protected:
  // methods
  Manager();
};

} /* namespace selectivities */
} /* namespace niwa */
#endif /* PREFERENCE_FUNCTION_MANAGER_H_ */
