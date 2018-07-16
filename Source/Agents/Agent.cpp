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
#include "Utilities/RandomNumberGenerator.h"

// Namespaces
namespace niwa {


Agent::Agent(float lat, float lon, float first_age_length_par, float second_age_length_par, float M, unsigned age, float first_length_weight_par, float second_length_weigth_par) :
    lat_(lat),
    lon_(lon),
    first_age_length_par_(first_age_length_par),
    second_age_length_par_(second_age_length_par),
    survival_(1 - exp(-M)),
    age_(age),
    first_length_weight_par_(first_length_weight_par),
    second_length_weight_par_(second_length_weigth_par)

{
  growth_init(); // if age = 0 will set length_ = 0; otherwise will set to what ever the length at age dictates.
}

/*
 * Survival event
*/
void Agent::survival(float& selectivity) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  if (not (rng.chance() > survival_ * selectivity))
    alive_ = false;
}

/*
 * Survival event
*/
void Agent::maturity(float& selectivity) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  if (not mature_) {
    if (rng.chance() < selectivity)
      mature_ = true;
  }
}


/*
 * An internal function to set initial length at age when initially seeding agents in the world,
 * So that we have the equivalent length and weight frequency. Calculate expected length at age assuming von Bert parameters
 *
*/
void Agent::growth_init() {
  length_ = first_age_length_par_ * (1-std::exp(-second_age_length_par_ * (float)age_));
  weight_ = first_length_weight_par_ * pow(length_, second_length_weight_par_); // Just update weight when ever we update length to save executions
}

} /* namespace niwa */
