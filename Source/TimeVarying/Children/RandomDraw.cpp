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
#include <TimeVarying/Children/RandomDraw.h>
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
  parameters_.Bind<string>(PARAM_DISTRIBUTION, &distribution_label_, "distribution", "", PARAM_NORMAL)->set_allowed_values({PARAM_NORMAL,PARAM_LOGNORMAL});
  parameters_.Bind<float>(PARAM_LOWER_BOUND, &lower_bound_, "Lower bound", "");
  parameters_.Bind<float>(PARAM_UPPER_BOUND, &upper_bound_, "Upper bound", "");

  RegisterAsAddressable(PARAM_MEAN, &mu_);
  RegisterAsAddressable(PARAM_SIGMA, &sigma_);
}

/**
 *
 */
void RandomDraw::DoValidate() {

  if (distribution_label_ == PARAM_NORMAL) {
    distribution_ = Distribution::kNormal;
  } else {
    distribution_ = Distribution::kLogNormal;
  }
}

/**
 *
 */
void RandomDraw::DoBuild() {
  LOG_FINE() << "Dobuild";
  DoReset();
}

/**
 *
 */
void RandomDraw::DoReset() {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  float new_value = 0.0;
  // Draw from the random distribution
  if (distribution_ == Distribution::kNormal) {
    for (unsigned year : years_) {
      new_value = rng.normal((mu_), (sigma_));
      LOG_FINEST() << "with mean = " << mu_ << " and sigma = " << sigma_ << " new value = " << new_value;
      if (new_value < lower_bound_)
        new_value  = lower_bound_;
      if (new_value > upper_bound_)
        new_value  = upper_bound_;
      parameter_by_year_[year] = new_value;
    }
  } else if (distribution_ == Distribution::kLogNormal)  {
    for (unsigned year : years_) {
      float cv = sqrt(exp(sigma_ * sigma_) - 1);
      new_value = rng.lognormal(mu_, cv);
      LOG_FINEST() << "with mean = " << mu_ << " and sigma = " << cv << " new value = " << new_value;
      if (new_value < lower_bound_)
        new_value  = lower_bound_;
      if (new_value > upper_bound_)
        new_value  = upper_bound_;
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
