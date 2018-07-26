/**
 * @file DoubleNormal.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 14/01/2013
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers
#include "DoubleNormal.h"

#include <boost/math/distributions/lognormal.hpp>
#include <cmath>

#include "Model/Model.h"
#include "TimeSteps/Manager.h"

// Namespaces
namespace niwa {
namespace selectivities {

/**
 * Explicit constructor
 */
DoubleNormal::DoubleNormal(Model* model)
: Selectivity(model) {

  parameters_.Bind<float>(PARAM_MU, &mu_, "Mu", "");
  parameters_.Bind<float>(PARAM_SIGMA_L, &sigma_l_, "Sigma L", "");
  parameters_.Bind<float>(PARAM_SIGMA_R, &sigma_r_, "Sigma R", "");
  parameters_.Bind<float>(PARAM_ALPHA, &alpha_, "Alpha", "", 1.0);

  //RegisterAsAddressable(PARAM_MU, &mu_);
  //RegisterAsAddressable(PARAM_SIGMA_L, &sigma_l_);
  //RegisterAsAddressable(PARAM_SIGMA_R, &sigma_r_);
  //RegisterAsAddressable(PARAM_ALPHA, &alpha_);
}

/**
 * Validate this selectivity. This will load the
 * values that were passed in from the configuration
 * file and assign them to the local variables.
 *
 * We'll then do some basic checks on the local
 * variables to ensure they are within the business
 * rules for the model.
 */
void DoubleNormal::DoValidate() {
  if (alpha_ <= 0.0)
    LOG_ERROR_P(PARAM_MU) << ": alpha cannot be less than or equal to 0.0";
  if (sigma_l_ <= 0.0)
    LOG_ERROR_P(PARAM_SIGMA_L) << ": sigma_l cannot be less than or equal to 0.0";
  if (sigma_r_ <= 0.0)
    LOG_ERROR_P(PARAM_SIGMA_R) << ": sigmal_r cannot be less than or equal to 0.0";
}

/**
 * Reset this selectivity so it's ready for the next execution
 * phase in the model.
 *
 * This method will rebuild the cache of selectivity values
 * for each age in the model.
 */
void DoubleNormal::RebuildCache() {
  if (not length_based_) {
    for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age) {
      float temp = (float)age;
      if (temp < mu_)
        values_[age - min_index_] = pow(2.0, -((temp - mu_) / sigma_l_ * (temp - mu_) / sigma_l_)) * alpha_;
      else
        values_[age - min_index_] = pow(2.0, -((temp - mu_) / sigma_r_ * (temp - mu_) / sigma_r_)) * alpha_;
    }
  } else {
    vector<unsigned> length_bins = model_->length_bins();
    for (unsigned length_bin_index = 0; length_bin_index < length_bins.size(); ++length_bin_index) {
      float temp = (float)length_bins[length_bin_index];
      if (temp < mu_)
        length_values_[length_bin_index] = pow(2.0, -((temp - mu_) / sigma_l_ * (temp - mu_) / sigma_l_)) * alpha_;
      else
        length_values_[length_bin_index] = pow(2.0, -((temp - mu_) / sigma_r_ * (temp - mu_) / sigma_r_)) * alpha_;
    }
  }
}

} /* namespace selectivities */
} /* namespace niwa */
