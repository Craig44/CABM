/**
 * @file Factory.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 13/12/2012
 * @section LICENSE
 *
 * Copyright NIWA Science ©2012 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers

#include "Factory.h"

#include "Model/Model.h"
#include "Model/Managers.h"
#include "Reports/Manager.h"
#include "Reports/Children/SummariseAgents.h"
#include "Reports/Children/DerivedQuantity.h"
#include "Reports/Children/AgeFrequencyByCell.h"
#include "Reports/Children/Process.h"
#include "Reports/Children/InitialisationPartition.h"
#include "Reports/Children/Observation.h"
#include "Reports/Children/Selectivity.h"
#include "Reports/Children/TimeVarying.h"
#include "Reports/Children/NumericLayer.h"
#include "Reports/Children/WorldAgeFrequency.h"
#include "Reports/Children/ModelAttributes.h"

// Namespaces
namespace niwa {
namespace reports {

/**
 * Create the instance of our object as defined by the two parameters
 * object_type and sub_type.
 *
 * @param object_type The type of object to create (e.g age_size, process)
 * @param sub_type The child type of the object to create (e.g ageing, schnute)
 * @return shared_ptr to the object we've created
 */
Report* Factory::Create(Model* model, const string& object_type, const string& sub_type) {
  Report* result = nullptr;

  if (object_type == PARAM_REPORT) {
    if (sub_type == PARAM_PROCESS)
      result = new Process(model);
    else if (sub_type == PARAM_AGE_FREQUENCY_BY_CELL)
      result = new AgeFrequencyByCell(model);
    else if (sub_type == PARAM_INITIALISATION_PARTITION)
      result = new InitialisationPartition(model);
    else if (sub_type == PARAM_SUMMARISE_AGENTS)
      result = new SummariseAgents(model);
    else if (sub_type == PARAM_WORLD_AGE_FREQUENCY)
      result = new WorldAgeFrequency(model);
    else if (sub_type == PARAM_DERIVED_QUANTITY)
      result = new DerivedQuantity(model);
    else if (sub_type == PARAM_OBSERVATION)
      result = new Observation(model);
    else if (sub_type == PARAM_MODEL_ATTRIBUTES)
      result = new ModelAttributes(model);
    else if (sub_type == PARAM_NUMERIC_LAYER)
      result = new NumericLayer(model);
    else if (sub_type == PARAM_SELECTIVITY)
      result = new Selectivity(model);
    else if (sub_type == PARAM_TIME_VARYING)
      result = new TimeVarying(model);
    if (result) {
      LOG_FINE() << "Creating report " << sub_type;
      model->managers().report()->AddObject(result);
    }
  }

  return result;
}

} /* namespace reports */
} /* namespace niwa */
