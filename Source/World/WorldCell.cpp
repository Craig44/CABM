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
#include "Processes/Manager.h"

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
void WorldCell::Build(unsigned row, unsigned col, float lat, float lon, Model* model) {
  LOG_TRACE();
  row_ = row;
  col_ = col;
  lat_ = lat;
  lon_ = lon;
  model_ = model;
  // Get a pointer to the environment from the model and give each Agent a pointer to that environment.

  // Build Growth and mortality functons
  growth_ = model_->managers().process()->GetGrowthProcess(model_->get_growth_process());
  if (!growth_) {
    LOG_CODE_ERROR() << "!growth_ this should have been checked in the model class";
  }

  mortality_ = model_->managers().process()->GetMortalityProcess(model_->get_mortality_process());
  if (!mortality_) {
    LOG_CODE_ERROR() << "!mortality this should have been checked in the model class";
  }


}

/**
 * Reset our partition so all data values are 0.0
 */
void WorldCell::Reset() {
  LOG_TRACE();
  agents_.clear();
}

/**
 * This method is called in Initialisation where we seed agents over the spatial domain before we start iterating
 * each agent that is created will be call seed() this will seed an agent with an equilibrium age structure
 */
void WorldCell::seed_agents(unsigned number_agents_to_seed, const float& seed_z) {
  LOG_TRACE();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  vector<float> mort_par;
  vector<vector<float>> growth_pars;
  mortality_->draw_rate_param(row_, col_, number_agents_to_seed, mort_par);
  growth_->draw_growth_param(row_, col_, number_agents_to_seed, growth_pars);

  if (number_agents_to_seed != mort_par.size()) {
    LOG_CODE_ERROR() << "number_agents_to_seed != mort_par.size(), this must be a code error as these should always be true";
  }
  if (number_agents_to_seed != growth_pars.size()) {
    LOG_CODE_ERROR() << "number_agents_to_seed != growth_pars.size(), this must be a code error as these should always be true";
  }

  unsigned age;
  for (unsigned agent = 0; agent < number_agents_to_seed; ++agent) {
    age = std::max(std::min((unsigned)rng.exponential(seed_z), model_->max_age()),model_->min_age()); // truncate age to between min_age and max_age
    //LOG_FINEST() << age;
    Agent new_agent(lat_, lon_, growth_pars[agent][0], growth_pars[agent][1], mort_par[agent], (model_->current_year() - age), growth_pars[agent][2], growth_pars[agent][3], model_); // seed it with lat long, L_inf, K
    agents_.push_back(new_agent);
  }
}

/**
 * This method is called in Recruitment processes where we create new agents.
 */
void WorldCell::birth_agents(unsigned birth_agents) {
  LOG_TRACE();
  vector<float> mort_par;
  vector<vector<float>> growth_pars;
  mortality_->draw_rate_param(row_, col_, birth_agents, mort_par);
  growth_->draw_growth_param(row_, col_, birth_agents, growth_pars);

  for (unsigned agent = 0; agent < birth_agents; ++agent) {
    Agent new_agent(lat_, lon_, growth_pars[agent][0], growth_pars[agent][1], mort_par[agent], model_->current_year(), growth_pars[agent][2], growth_pars[agent][3], model_); // seed it with lat long, L_inf, K
    agents_.push_back(new_agent);
  }
}

/*
 * Returns the scaled up abundance for this cell
*/
float  WorldCell::get_abundance() {
  float abundance = 0.0;
  for (auto& agent : agents_)
    abundance += agent.get_scalar();
  return abundance;
}

/*
 * Returns the scaled up biomass for this cell
*/
float  WorldCell::get_biomass() {
  float biomass = 0.0;
  for (auto& agent : agents_)
    biomass += agent.get_weight() * agent.get_scalar();
  return biomass;
}

/*
 * Returns the age frequency of agents in this cell
*/
void  WorldCell::get_age_frequency(vector<unsigned>& age_freq) {
  age_freq.clear();
  age_freq.resize(model_->age_spread(),0);
  for (auto iter = agents_.begin(); iter != agents_.end(); ++iter) {
    age_freq[(*iter).age() - model_->min_age()]++;
  }
}

} /* namespace niwa */
