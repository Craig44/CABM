/**
 * @file Manager.cpp
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 3/02/2015
 * @section LICENSE
 *
 * Copyright NIWA Science �2015 - www.niwa.co.nz
 *
 */

// headers
#include "Manager.h"

// namespaces
namespace niwa {
namespace preference_functions {

/**
 * Default constructor
 */
Manager::Manager() {

}

/**
 *
 */
PreferenceFunction* Manager::GetPreferenceFunction(const string& label) {
  for(auto preference_function : objects_) {
    if (preference_function->label() == label) {
      return preference_function;
    }
  }

  return nullptr;
}

} /* namespace selectivities */
} /* namespace niwa */
