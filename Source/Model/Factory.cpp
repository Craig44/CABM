/**
 * @file Factory.cpp
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 8/10/2015
 * @section LICENSE
 *
 * Copyright NIWA Science ©2015 - www.niwa.co.nz
 *
 */

// headers
#include "Factory.h"

// Common Factories
#include "BaseClasses/Object.h"
#include "AgeingErrors/Factory.h"
#include "Asserts/Factory.h"
#include "Layers/Factory.h"
#include "DerivedQuantities/Factory.h"
#include "InitialisationPhases/Factory.h"
#include "Likelihoods/Factory.h"
#include "Model/Model.h"
#include "Observations/Factory.h"
#include "Processes/Factory.h"
#include "Reports/Factory.h"
#include "Selectivities/Factory.h"
#include "TimeSteps/Factory.h"
#include "TimeVarying/Factory.h"
#include "Utilities/To.h"


// Length Factories
// namespaces
namespace niwa {

/**
 * Default constructor
 */
Factory::Factory(Model* model) : model_(model) { }

/**
 * Create an ObjectPtr for a specific class type in our system. This method
 * will check the object_type and find the appropriate child factory to call
 * so the object can be created properly.
 *
 * This design was picked to simplify how the factories work because while having
 * a single super-factory would be the simplest solution we need objects to be created
 * in their child form (e.g AgeSizePtr) to be used properly in the code without
 * having to do dynamic casting (slow).
 *
 * @param object_type The type of object to create (e.g process, selectivity)
 * @param sub_type The specialisation/sub_type of the object to create
 * @param partition_type The specifically defined partition type for this object
 * @return A shared_ptr to the object we've created
 */
base::Object* Factory::CreateObject(const string& object_type, const string& sub_type) {
  string lwr_object_type    = utilities::ToLowercase(object_type);
  string lwr_sub_type       = utilities::ToLowercase(sub_type);

  if (lwr_object_type == PARAM_ASSERT)
    return asserts::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_LAYER)
    return layers::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_AGEING_ERROR || lwr_object_type == PARAM_AGEING_ERRORS)
    return ageingerrors::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_DERIVED_QUANTITY || lwr_object_type == PARAM_DERIVED_QUANTITIES)
    return derivedquantities::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_INITIALISATION_PHASE || lwr_object_type == PARAM_INITIALISATION_PHASES)
    return initialisationphases::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_LIKELIHOOD)
    return likelihoods::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_MODEL)
    return model_;
  else if (lwr_object_type == PARAM_OBSERVATION)
    return observations::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_PROCESS || lwr_object_type == PARAM_PROCESSES)
    return processes::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_STATE || lwr_object_type == PARAM_TAG || lwr_object_type == PARAM_TRANSITION) // @process specialisation
    return processes::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_REPORT)
    return reports::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_SELECTIVITY || lwr_object_type == PARAM_SELECTIVITIES)
    return selectivities::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_TIME_STEP || lwr_object_type == PARAM_TIME_STEPS)
    return timesteps::Factory::Create(model_, lwr_object_type, lwr_sub_type);
  else if (lwr_object_type == PARAM_TIME_VARYING)
    return timevarying::Factory::Create(model_, lwr_object_type, lwr_sub_type);

  return nullptr;
}

} /* namespace niwa */
