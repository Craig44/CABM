/**
 * @file RandomWalk.cpp
 * @author Craig Marsh
 * @github https://github.com/Zaita
 * @date 2/02/2016
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 */

// headers
#include <TimeVarying/Children/RandomWalk.h>
#include "Utilities/Map.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Model/Objects.h"


// namespaces
namespace niwa {
namespace timevarying {

/**
 * Default constructor
 */
RandomWalk::RandomWalk(Model* model) : TimeVarying(model) {
  parameters_.Bind<float>(PARAM_MEAN, &mu_, "Mean", "", 0);
  parameters_.Bind<float>(PARAM_SIGMA, &sigma_, "Standard deviation", "", 1);
  parameters_.Bind<float>(PARAM_UPPER_BOUND, &upper_bound_, "Upper bound for the random walk", "", 1);
  parameters_.Bind<float>(PARAM_UPPER_BOUND, &lower_bound_, "Lower bound for the random walk", "", 1);
  parameters_.Bind<float>(PARAM_RHO, &rho_, "Auto Correlation parameter", "", 1);
  parameters_.Bind<string>(PARAM_DISTRIBUTION, &distribution_, "distribution", "", PARAM_NORMAL);

  RegisterAsAddressable(PARAM_MEAN, &mu_);
  RegisterAsAddressable(PARAM_SIGMA, &sigma_);
}

/**
 *
 */
void RandomWalk::DoValidate() {
  if (distribution_ != PARAM_NORMAL)
    LOG_ERROR() << "Random Walk can only can draw from a normal distribution currently";
}

/**
 *
 */
void RandomWalk::DoBuild() {
  if(model_->objects().GetAddressableType(parameter_) != addressable::kSingle)
    LOG_ERROR_P(PARAM_TYPE) << "@time_varying blocks of type " << PARAM_RANDOMWALK << " can only be implemented in parameters that are scalars or single values";
}

/**
 *
 */
void RandomWalk::DoUpdate() {
  LOG_FINEST() << "value = " << *addressable_;
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  float value = *addressable_;
  float deviate = rng.normal((mu_), (sigma_));
  value += value * rho_ + deviate;


  LOG_FINEST() << "value after deviate of " << deviate << " = " << value << " for year " << model_->current_year();

  LOG_FINE() << "Setting Value to: " << value;
  parameter_by_year_[model_->current_year()] = value;
  (this->*update_function_)(value);
}

} /* namespace timevarying */
} /* namespace niwa */