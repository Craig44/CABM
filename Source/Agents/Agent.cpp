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

// Namespaces
namespace niwa {


Agent::Agent(float lat, float lon, float first_age_length_par, float second_age_length_par, float M, unsigned birth_year, float first_length_weight_par, float second_length_weigth_par, Model* model, bool mature, unsigned sex) :
    lat_(lat),
    lon_(lon),
    first_age_length_par_(first_age_length_par),
    second_age_length_par_(second_age_length_par),
    M_(M),
    birth_year_(birth_year),
    first_length_weight_par_(first_length_weight_par),
    second_length_weight_par_(second_length_weigth_par),
    model_(model),
    mature_(mature),
    sex_(sex)

{
  growth_init(); // if age = 0 will set length_ = 0; otherwise will set to what ever the length at age dictates.
}

// Return Age this allows for implicit ageing, which is handy as it reduces a process out of the dynamics
unsigned Agent::age() {
    unsigned age =  std::min((model_->current_year() - birth_year_), model_->max_age());
    return age;
}

/*
 * An internal function to set initial length at age when initially seeding agents in the world,
 * So that we have the equivalent length and weight frequency. Calculate expected length at age assuming von Bert parameters
 * TODO figure out how to generalise this
 *
*/
void Agent::growth_init() {
  length_ = first_age_length_par_ * (1-std::exp(-second_age_length_par_ * (float)age()));
  weight_ = first_length_weight_par_ * pow(length_, second_length_weight_par_); // Just update weight when ever we update length to save executions
  //LOG_FINEST() << "initialise agent, age = " << age() << " length = " << length_ << " weight = " << weight_;
}

} /* namespace niwa */
