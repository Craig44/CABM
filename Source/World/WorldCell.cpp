/**
 * @file WorldCell.cpp
 * @author  C.Marsh
 * @version 1.0
 * @date 12/7/2018
 * @section LICENSE
 *
 *
 */

// Headers
#include "WorldCell.h"

#include "Agents/Agent.h"

// Namespaces
namespace niwa {

/**
 *
 */
void WorldCell::Validate() {
  LOG_TRACE();

}

/**
 * Build our partition structure now. This involves getting
 *
 * We're not interested in the range of years that each
 * category has because this will be addressed with the
 * accessor objects.
 */
void WorldCell::Build() {
  LOG_TRACE();

  // Get a pointer to the environment from the model and give each Agent a pointer to that environment.

}

/**
 * Reset our partition so all data values are 0.0
 */
void WorldCell::Reset() {
  LOG_TRACE();
}

/**
 * This method is called in Initialisation where we seed agents over the spatial domain before we start iterating
 * each agent that is created will be call seed() this will seed an agent with an equilibrium age structure
 */
void WorldCell::seed_agents(unsigned number_agents_to_seed) {
  for (unsigned agent = 0; agent < number_agents_to_seed; ++agent) {
    Agent new_agent; // seed it with lat long, K, L_inf
    new_agent.seed();
    agents_.push_back(new_agent);
  }
}

} /* namespace niwa */
