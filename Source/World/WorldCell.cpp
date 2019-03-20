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
#include "Selectivities/Manager.h"
#include <sstream>

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

  vector<string> mature_labels =  model_->get_maturity_ogive();
  for (auto label : mature_labels) {
    Selectivity* temp_selectivity = nullptr;
    temp_selectivity = model_->managers().selectivity()->GetSelectivity(label);
    if (!temp_selectivity)
      LOG_CODE_ERROR()<< "this should have been checked on the ModelDoBuild please check out, issue with " << label << " selectivity";
    selectivity_.push_back(temp_selectivity);
  }

  // Give each cell a label for reporting
  std::stringstream label;
  label << row_ << "-" << col_;
  cell_label_ = label.str();
  LOG_FINEST() << "cell_label = " << cell_label_;
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
  LOG_FINE();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  vector<float> mort_par;
  mortality_->draw_rate_param(row_, col_, number_agents_to_seed, mort_par);

  if (number_agents_to_seed != mort_par.size()) {
    LOG_CODE_ERROR() << "number_agents_to_seed != mort_par.size(), this must be a code error as these should always be true";
  }


  unsigned age, sex;
  bool mature = false, sexed;
  float male_prop = 1.0, probability_mature_at_age;
  sexed = model_->get_sexed();

  vector<vector<float>> growth_pars;
  vector<vector<float>> female_growth_pars;
  growth_->draw_growth_param(row_, col_, number_agents_to_seed, growth_pars, 0);
  if (sexed) {
    // This creates growth vectors that are larger than we need but it shouldn't be to expensive
    male_prop = model_->get_male_proportions();
    growth_->draw_growth_param(row_, col_, number_agents_to_seed, female_growth_pars, 1);
  }
  if (number_agents_to_seed != growth_pars.size()) {
    LOG_CODE_ERROR() << "number_agents_to_seed != growth_pars.size(), this must be a code error as these should always be true";
  }
  LOG_FINE() << "start seeding";
  for (unsigned agent = 0; agent < number_agents_to_seed; ++agent) {
    sex = 0;
    mature = false;
    age = std::max(std::min((unsigned)rng.exponential(seed_z), model_->max_age()), model_->min_age()); // truncate age to between min_age and max_age
    // Need to add Maturity and sex into this
    if (sexed) {
      if (rng.chance() >= male_prop)
        sex = 1;
    }
    probability_mature_at_age = selectivity_[sex]->GetResult(age);
    if (rng.chance() <= probability_mature_at_age)
      mature = true;

    if (sex == 0) {
      Agent new_agent(lat_, lon_, growth_pars[agent][0], growth_pars[agent][1], growth_pars[agent][2], mort_par[agent], (model_->current_year() - age),
          growth_pars[agent][3], growth_pars[agent][4], model_, mature, sex, 1.0, row_, col_, 0);
      agents_.push_back(new_agent);
    } else {
      Agent new_agent(lat_, lon_, female_growth_pars[agent][0], female_growth_pars[agent][1], female_growth_pars[agent][2], mort_par[agent], (model_->current_year() - age),
          female_growth_pars[agent][3], female_growth_pars[agent][4], model_, mature, sex, 1.0, row_, col_, 0);
      agents_.push_back(new_agent);
    }

/*
   if (agent == 0) {
      LOG_MEDIUM() << "number of bytes of an agent class " << sizeof(new_agent);
    }*/
  }
}

/**
 * This method is called in Recruitment processes where we create new agents.
 */
void WorldCell::birth_agents(unsigned birth_agents,float scalar) {
  LOG_TRACE();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  vector<float> mort_par;
  mortality_->draw_rate_param(row_, col_, birth_agents, mort_par);
  bool sexed = model_->get_sexed();
  float male_prop = model_->get_male_proportions();

  vector<vector<float>> growth_pars;
  vector<vector<float>> female_growth_pars;
  growth_->draw_growth_param(row_, col_, birth_agents, growth_pars, 0);
  if (sexed) {
    // This creates growth vectors that are larger than we need but it shouldn't be to expensive
    growth_->draw_growth_param(row_, col_, birth_agents, female_growth_pars, 1);
  }

  unsigned sex;
  unsigned agent_bin = 0;
  for (unsigned agent = 0; agent < birth_agents; ++agent) {
    sex = 0;
    if (sexed) {
      if (rng.chance() >= male_prop)
        sex = 1;
    }
    if (sex == 0) {
      Agent new_agent(lat_, lon_, growth_pars[agent][0], growth_pars[agent][1],growth_pars[agent][2], mort_par[agent], model_->current_year() - model_->min_age(),
          growth_pars[agent][3], growth_pars[agent][4], model_, false, sex, scalar, row_, col_, 0);
      // Check to see if
      while (agent_bin < agents_.size()) {
        if (not agents_[agent_bin].is_alive()) {
          agents_[agent_bin] = new_agent;
          break;
        } else {
          agent_bin++;
        }
      }
      if (agent_bin == agents_.size()) {
        agents_.push_back(new_agent);
      }
    } else {
      Agent new_agent(lat_, lon_, female_growth_pars[agent][0], female_growth_pars[agent][1],female_growth_pars[agent][2], mort_par[agent], model_->current_year() - model_->min_age(),
          female_growth_pars[agent][3], female_growth_pars[agent][4], model_, false, sex, scalar, row_, col_, 0);
      // Check to see if
      while (agent_bin < agents_.size()) {
        if (not agents_[agent_bin].is_alive()) {
          agents_[agent_bin] = new_agent;
          break;
        } else {
          agent_bin++;
        }
      }
      if (agent_bin == agents_.size()) {
        agents_.push_back(new_agent);
      }
    }
  }
}

/*
 * apply_mortality_time_varying applied when we update M parameters using @time_varying block
 */
void  WorldCell::apply_mortality_time_varying() {
  LOG_FINE() << " ";
  unsigned counter = 0;
  vector<float> mort_par;
  mortality_->draw_rate_param(row_, col_, agents_.size(), mort_par);
  LOG_FINE() << " ";
  for (auto iter = agents_.begin(); iter != agents_.end(); ++iter, ++counter) {
    if ( (*iter).is_alive())
      (*iter).set_m(mort_par[counter]);
  }

}

/*
 * apply_growth_time_varying applied when we update M parameters using @time_varying block
 * These update functions just randomly allocate a new agent a new growth parameter based on a distribution.
 */
void  WorldCell::apply_growth_time_varying() {
  LOG_FINE() << " ";
  unsigned counter = 0;
  vector<vector<float>> growth_pars;
  vector<vector<float>> female_growth_pars;
  growth_->draw_growth_param(row_, col_, agents_.size(), growth_pars, 0);
  if (model_->get_sexed()) {
    growth_->draw_growth_param(row_, col_, agents_.size(), female_growth_pars, 1);
  }

  for (auto iter = agents_.begin(); iter != agents_.end(); ++iter, ++counter) {
    if ( (*iter).is_alive()) {
      if ((*iter).get_sex() == 0) {
        (*iter).set_first_age_length_par(growth_pars[counter][0]);
        (*iter).set_second_age_length_par(growth_pars[counter][1]);
        (*iter).set_third_age_length_par(growth_pars[counter][2]);
        (*iter).set_first_length_weight_par(growth_pars[counter][3]);
        (*iter).set_second_length_weight_par(growth_pars[counter][4]);
      } else {
        (*iter).set_first_age_length_par(female_growth_pars[counter][0]);
        (*iter).set_second_age_length_par(female_growth_pars[counter][1]);
        (*iter).set_third_age_length_par(female_growth_pars[counter][2]);
        (*iter).set_first_length_weight_par(female_growth_pars[counter][3]);
        (*iter).set_second_length_weight_par(female_growth_pars[counter][4]);
      }
    }
  }
}
/*
 * Called in WorldView usually after a movement process, where agents get assigned new parameters if we have spatial Mortality, or growth
 * This is for spatial updating
*/
void  WorldCell::update_agent_parameters() {
  LOG_TRACE();
  unsigned counter = 0;
  if (growth_->update_growth() && mortality_->update_mortality()) {
    vector<float> mort_par;
    mortality_->draw_rate_param(row_, col_, agents_.size(), mort_par);
    vector<vector<float>> growth_pars;
    vector<vector<float>> female_growth_pars;
    growth_->draw_growth_param(row_, col_, agents_.size(), growth_pars, 0);
    if (model_->get_sexed()) {
      growth_->draw_growth_param(row_, col_, agents_.size(), female_growth_pars, 1);
    }

    for (auto iter = agents_.begin(); iter != agents_.end(); ++iter, ++counter) {
      if ( (*iter).is_alive()) {
        if ((*iter).get_sex() == 0) {
          (*iter).set_first_age_length_par(growth_pars[counter][0]);
          (*iter).set_second_age_length_par(growth_pars[counter][1]);
          (*iter).set_third_age_length_par(growth_pars[counter][2]);
          (*iter).set_first_length_weight_par(growth_pars[counter][3]);
          (*iter).set_second_length_weight_par(growth_pars[counter][4]);
        } else {
          (*iter).set_first_age_length_par(female_growth_pars[counter][0]);
          (*iter).set_second_age_length_par(female_growth_pars[counter][1]);
          (*iter).set_third_age_length_par(female_growth_pars[counter][2]);
          (*iter).set_first_length_weight_par(female_growth_pars[counter][3]);
          (*iter).set_second_length_weight_par(female_growth_pars[counter][4]);
        }
        (*iter).set_m(mort_par[counter]);
      }
    }
  } else if (!growth_->update_growth() && mortality_->update_mortality()) {
    vector<float> mort_par;
    mortality_->draw_rate_param(row_, col_, agents_.size(), mort_par);
    for (auto iter = agents_.begin(); iter != agents_.end(); ++iter, ++counter) {
      if ( (*iter).is_alive()) {
        (*iter).set_m(mort_par[counter]);
      }
    }
  } else if (growth_->update_growth() && !mortality_->update_mortality()) {
    vector<vector<float>> growth_pars;
    vector<vector<float>> female_growth_pars;
    growth_->draw_growth_param(row_, col_, agents_.size(), growth_pars, 0);
    if (model_->get_sexed()) {
      growth_->draw_growth_param(row_, col_, agents_.size(), female_growth_pars, 1);
    }
    for (auto iter = agents_.begin(); iter != agents_.end(); ++iter, ++counter) {
      if ( (*iter).is_alive()) {
        if ((*iter).get_sex() == 0) {
          (*iter).set_first_age_length_par(growth_pars[counter][0]);
          (*iter).set_second_age_length_par(growth_pars[counter][1]);
          (*iter).set_third_age_length_par(growth_pars[counter][2]);
          (*iter).set_first_length_weight_par(growth_pars[counter][3]);
          (*iter).set_second_length_weight_par(growth_pars[counter][4]);
        } else {
          (*iter).set_first_age_length_par(female_growth_pars[counter][0]);
          (*iter).set_second_age_length_par(female_growth_pars[counter][1]);
          (*iter).set_third_age_length_par(female_growth_pars[counter][2]);
          (*iter).set_first_length_weight_par(female_growth_pars[counter][3]);
          (*iter).set_second_length_weight_par(female_growth_pars[counter][4]);
        }
      }
    }
  }
}

/*
 * Called in WorldView if agents need to be updated with time varying parameters
*/
void  WorldCell::update_mortality_params() {
  LOG_FINE() << "this method assumes the mean values have been changed in the corresponding class, for example the M value has changed";
  vector<float> mort_par;
  unsigned counter = 0;
  mortality_->draw_rate_param(row_, col_, agents_.size(), mort_par);
  for (auto iter = agents_.begin(); iter != agents_.end(); ++iter, ++counter) {
    if ( (*iter).is_alive()) {
      (*iter).set_m(mort_par[counter]);
    }
  }
}

/*
 * Called in WorldView if agents need to be updated with time varying parameters
*/
void  WorldCell::update_growth_params() {
  LOG_FINE() << "Updating growth params based on time varying parameters";
  vector<vector<float>> growth_pars;
  vector<vector<float>> female_growth_pars;
  growth_->draw_growth_param(row_, col_, agents_.size(), growth_pars, 0);
  if (model_->get_sexed()) {
    growth_->draw_growth_param(row_, col_, agents_.size(), female_growth_pars, 1);
  }
  unsigned counter = 0;
  for (auto iter = agents_.begin(); iter != agents_.end(); ++iter, ++counter) {
    if ( (*iter).is_alive()) {
      if ((*iter).get_sex() == 0) {
        (*iter).set_first_age_length_par(growth_pars[counter][0]);
        (*iter).set_second_age_length_par(growth_pars[counter][1]);
        (*iter).set_third_age_length_par(growth_pars[counter][2]);
        (*iter).set_first_length_weight_par(growth_pars[counter][3]);
        (*iter).set_second_length_weight_par(growth_pars[counter][4]);
      } else {
        (*iter).set_first_age_length_par(female_growth_pars[counter][0]);
        (*iter).set_second_age_length_par(female_growth_pars[counter][1]);
        (*iter).set_third_age_length_par(female_growth_pars[counter][2]);
        (*iter).set_first_length_weight_par(female_growth_pars[counter][3]);
        (*iter).set_second_length_weight_par(female_growth_pars[counter][4]);
      }
    }
  }
}

/*
 * Returns the scaled up abundance for this cell
*/
float  WorldCell::get_abundance() {
  float abundance = 0.0;
  for (auto& agent : agents_) {
    if (agent.is_alive())
      abundance += agent.get_scalar();
  }
  return abundance;
}

/*
 * Returns the scaled up biomass for this cell
*/
float  WorldCell::get_biomass() {
  float biomass = 0.0;
  for (auto& agent : agents_) {
    if (agent.is_alive())
      biomass += agent.get_weight() * agent.get_scalar();
  }
  return biomass;
}

/*
 * Returns the scaled up mature biomass for this cell
*/
float  WorldCell::get_mature_biomass() {
  float biomass = 0.0;
  for (auto& agent : agents_) {
    if (agent.is_alive() & agent.get_maturity()) {
      //LOG_FINEST() << "weight = " << agent.get_weight() << " scalar = " << agent.get_scalar();
      biomass += agent.get_weight() * agent.get_scalar();
    }
  }
  return biomass;
}


/*
 * Returns the age frequency of agents in this cell
*/
void  WorldCell::get_age_frequency(vector<unsigned>& age_freq) {
  age_freq.clear();
  age_freq.resize(model_->age_spread(),0);
  for (auto iter = agents_.begin(); iter != agents_.end(); ++iter) {
    if ((*iter).is_alive()) {
      age_freq[(*iter).get_age_index()]++;
    }
  }
}



/*
 * Return male age frequency
*/
void  WorldCell::get_male_frequency(vector<unsigned>& age_freq) {
  age_freq.clear();
  age_freq.resize(model_->age_spread(),0);
  for (auto iter = agents_.begin(); iter != agents_.end(); ++iter) {
    if ((*iter).is_alive() and (*iter).get_sex() == 0) {
      age_freq[(*iter).get_age_index()]++;
    }
  }
}


/*
 * Returns total female age frequency
*/
void  WorldCell::get_female_frequency(vector<unsigned>& age_freq) {
  age_freq.clear();
  age_freq.resize(model_->age_spread(),0);
  for (auto iter = agents_.begin(); iter != agents_.end(); ++iter) {
    if ((*iter).is_alive() and (*iter).get_sex() == 1) {
      age_freq[(*iter).get_age_index()]++;
    }
  }
}
} /* namespace niwa */
