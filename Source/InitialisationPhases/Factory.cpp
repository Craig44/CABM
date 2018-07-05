/**
 * @file Factory.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 13/12/2012
 * @section LICENSE
 *
 * Copyright NIWA Science ©2017 - www.niwa.co.nz
 *
 */

// Headers
#include "Factory.h"

#include "Model/Model.h"
#include "Model/Managers.h"
#include "InitialisationPhases/Manager.h"

// Namespaces
namespace niwa {
namespace initialisationphases {

/**
 * Create the instance of our object as defined by the two parameters
 * object_type and sub_type.
 *
 * @param object_type The type of object to create (e.g age_size, process)
 * @param sub_type The child type of the object to create (e.g ageing, schnute)
 * @return shared_ptr to the object we've created
 */
InitialisationPhase* Factory::Create(Model* model, const string& object_type, const string& sub_type) {
  InitialisationPhase* result = nullptr;



  if (result)
    model->managers().initialisation_phase()->AddObject(result);

  return result;
}

} /* namespace processes */
} /* namespace niwa */
