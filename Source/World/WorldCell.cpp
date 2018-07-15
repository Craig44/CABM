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
#include "Utilities/RandomNumberGenerator.h"

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
void WorldCell::Build(unsigned row, unsigned col, double lat, double lon, unsigned min_age, unsigned max_age) {
  LOG_TRACE();
  row_ = row;
  col_ = col;
  lat_ = lat;
  lon_ = lon;
  min_age_ = min_age;
  max_age_ = max_age;
  age_spread_ = max_age - min_age + 1;
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
void WorldCell::seed_agents(unsigned number_agents_to_seed, const vector<double>&  mort_par, const vector<vector<double>>&  growth_par, const double& seed_z) {
  LOG_TRACE();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  if (number_agents_to_seed != mort_par.size()) {
    LOG_CODE_ERROR() << "number_agents_to_seed != mort_par.size(), this must be a code error as these should always be true";
  }
  if (number_agents_to_seed != growth_par.size()) {
    LOG_CODE_ERROR() << "number_agents_to_seed != growth_par.size(), this must be a code error as these should always be true";
  }


  unsigned age;
  for (unsigned agent = 0; agent < number_agents_to_seed; ++agent) {
    age = std::max(std::min((unsigned)rng.exponential(seed_z), max_age_),min_age_); // truncate age to between min_age and max_age
    //LOG_FINEST() << age;

    Agent new_agent(growth_par[agent][0], growth_par[agent][1], mort_par[agent], age); // seed it with lat long, K, L_inf
    agents_.push_back(new_agent);
  }
}

/*
 * Returns the age frequency of agents in this cell
*/
void  WorldCell::get_age_frequency(vector<unsigned>& age_freq) {
  age_freq.clear();
  age_freq.resize(age_spread_,0);
  for (auto iter = agents_.begin(); iter != agents_.end(); ++iter) {
    age_freq[(*iter).age() - min_age_]++;
  }
}

} /* namespace niwa */
