/**
 * @file RandomDraw.cpp
 * @author Craig Marsh
 * @github https://github.com/Zaita
 * @date 2/02/2016
 * @section LICENSE
 *
 * Copyright NIWA Science ©2014 - www.niwa.co.nz
 *
 */

// headers
#include "RandomDraw.h"

#include "Utilities/Map.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Model/Objects.h"


// namespaces
namespace niwa {
namespace timevarying {

/**
 * Default constructor
 */
RandomDraw::RandomDraw(Model* model) : TimeVarying(model) {
  parameters_.Bind<float>(PARAM_MEAN, &mu_, "Mean", "", 0);
  parameters_.Bind<float>(PARAM_SIGMA, &sigma_, "Standard deviation", "", 1);
  parameters_.Bind<string>(PARAM_DISTRIBUTION, &distribution_, "distribution", "", PARAM_NORMAL)->set_allowed_values({PARAM_NORMAL,PARAM_LOGNORMAL});
}

/**
 *
 */
void RandomDraw::DoValidate() {

}

/**
 *
 */
void RandomDraw::DoBuild() {

  if(model_->objects().GetAddressableType(parameter_) != addressable::kSingle)
    LOG_ERROR_P(PARAM_TYPE) << "@time_varying blocks of type " << PARAM_RANDOMWALK << " can only be implemented in parameters that are scalars or single values";
  DoReset();
}

/**
 *
 */
void RandomDraw::DoReset() {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  float new_value = 0.0;
  // Draw from the random distribution
  if (distribution_ == PARAM_NORMAL) {
    for (unsigned year : years_) {
    new_value = rng.normal(mu_, sigma_);
    LOG_FINEST() << "with mean = " << mu_ << " and sigma = " << sigma_ << " new value = " << new_value;
    // Set value
    if (has_at_estimate_) {
      if (new_value < lower_bound_) {
        LOG_FINEST() << "hit @estimate lower bound setting value from " << new_value << " to " << lower_bound_;
        new_value = lower_bound_;
      }
      if (new_value > upper_bound_) {
        LOG_FINEST() << "hit @estimate upper bound setting value from " << new_value << " to " << upper_bound_;
        new_value = upper_bound_;
      }
    }
    if (new_value <= 0.0) {
      LOG_WARNING() << "parameter: " << parameter_ << " random draw of value = " << new_value << " a natural lower bound of 0.0 has been forced so resetting the value = 0.01";
      new_value  = 0.01;
    }
    parameter_by_year_[year] = new_value;
    }
  } else if (distribution_ == PARAM_LOGNORMAL)  {
    for (unsigned year : years_) {
      float cv = sigma_ / mu_;
      new_value = rng.lognormal(mu_, cv);
      LOG_FINEST() << "with mean = " << mu_ << " and sigma = " << sigma_ << " new value = " << new_value;
      // Set value
      if (has_at_estimate_) {
        if (new_value < lower_bound_) {
          LOG_FINEST() << "hit @estimate lower bound setting value from " << new_value << " to " << lower_bound_;
          new_value = lower_bound_;
        }
        if (new_value > upper_bound_) {
          LOG_FINEST() << "hit @estimate upper bound setting value from " << new_value << " to " << upper_bound_;
          new_value = upper_bound_;
        }
      }
      if (new_value <= 0.0) {
        LOG_WARNING() << "parameter: " << parameter_ << " random draw of value = " << new_value << " a natural lower bound of 0.0 has been forced so resetting the value = 0.01";
        new_value  = 0.01;
      }
      parameter_by_year_[year] = new_value;
    }
  }
}

/**
 *
 */
void RandomDraw::DoUpdate() {
  LOG_FINE() << "Setting Value to: " << parameter_by_year_[model_->current_year()];
  (this->*update_function_)(parameter_by_year_[model_->current_year()]);
}

} /* namespace timevarying */
} /* namespace niwa */
