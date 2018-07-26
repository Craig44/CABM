/**
 * @file Factory.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 23/01/2013
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * The time class represents a moment of time.
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers
#include "Factory.h"

#include "Model/Model.h"
#include "Model/Managers.h"
#include "PreferenceFunctions/Manager.h"


// Namespaces
namespace niwa {
namespace preference_functions {

/**
 * Create the instance of our object as defined by the two parameters
 * object_type and sub_type.
 *
 * @param object_type The type of object to create (e.g age_size, process)
 * @param sub_type The child type of the object to create (e.g ageing, schnute)
 * @return shared_ptr to the object we've created
 */
PreferenceFunction* Factory::Create(Model* model, const string& object_type, const string& sub_type) {
  PreferenceFunction* result = nullptr;

  if (object_type == PARAM_PREFERENCE_FUNCTION || object_type == PARAM_PREFERENCE_FUNCTIONS) {

    if (result)
      model->managers().preference_function()->AddObject(result);
  }

  return result;
}

} /* namespace preference_functions */
} /* namespace niwa */
