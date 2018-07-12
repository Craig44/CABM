/**
 * @file Agent.cpp
 * @author  C.Marsh
 * @version 1.0
 * @date 11/06/2018
 * @section LICENSE
 *
 * Copyright NIWA Science ©2012 - www.niwa.co.nz
 *
 */

// Headers
#include "Agent.h"

#include "Model/Model.h"
#include "Utilities/String.h"
#include "Utilities/To.h"

// Namespaces
namespace niwa {

/**
 * Default constructor
 *
 * Bind any parameters that are allowed to be loaded from the configuration files.
 * Set bounds on registered parameters
 * Register any parameters that can be an estimated or utilised in other run modes (e.g profiling, yields, projections etc)
 * Set some initial values
 *
 * Note: The constructor is parsed to generate Latex for the documentation.
 */
Agent::Agent(Model* model) : model_(model) {
	// Perhaps think if Agent attributes would ever differ from Model attributes.
}

/**
 * Populate any parameters,
 * Validate values are within expected ranges when we cannot use bind<>() overloads
 *
 * Note: all parameters are populated from configuration files
 */
void Agent::Validate() {
  // Check that we actually had a categories block
  parameters_.Populate(model_);

}

/**
 * Build any objects that will need to be utilised by this object.
 * Obtain smart_pointers to any objects that will be used by this object.
 */

void Agent::Build() {

}

} /* namespace niwa */
