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


Agent::Agent(double first_growth_par, double second_growth_par, double M, unsigned age) : first_growth_par_(first_growth_par), second_growth_par_(second_growth_par), survival_(1 - exp(-M)), age_(age)
{

}

/*
 * Survival event
*/
void Agent::survival() {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  if (not (rng.chance() > survival_))
    alive_ = false;
}


} /* namespace niwa */
