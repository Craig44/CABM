/**
 * @file Partition.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 7/11/2012
 * @section LICENSE
 *
 * Copyright NIWA Science �2012 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// defines
#define _USE_MATH_DEFINES

// Headers
#include "Partition.h"

#include <cmath>

#include "Partition/Agent.h"
#include "Model/Model.h"
#include "Logging/Logging.h"

// Namespaces
namespace niwa {

/**
 *
 */
void Partition::Validate() {
}

/**
 * Build our partition structure now. This involves getting
 *
 * We're not interested in the range of years that each
 * category has because this will be addressed with the
 * accessor objects.
 */
void Partition::Build() {
  // Get a pointer to the environment from the model and give each Agent a pointer to that environment.

  //
  for(unsigned i = 0; i < model_->number_of_agents_to_seed(); ++i) {
    LOG_FINEST() << "Adding Agent to the partition";
    agents::Agent* new_agent = new agents::Agent(model_);
    // Initialise Partition here.

    partition_.push_back(*new_agent);
  }
}

/**
 * Reset our partition so all data values are 0.0
 */
void Partition::Reset() {

}

} /* namespace niwa */
