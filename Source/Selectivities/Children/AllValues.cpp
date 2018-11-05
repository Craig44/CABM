/**
 * @file AllValues.cpp
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
#include "AllValues.h"

#include <boost/math/distributions/lognormal.hpp>
#include <cmath>

#include "Model/Model.h"

// Namespaces
namespace niwa {
namespace selectivities {

/**
 * Explicit Constructor
 */
AllValues::AllValues(Model* model)
: Selectivity(model) {

  parameters_.Bind<float>(PARAM_V, &v_, "V", "");

  //RegisterAsAddressable(PARAM_V, &v_);
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
void AllValues::DoValidate() {
  if (not length_based_) {
    if (v_.size() != model_->age_spread()) {
      LOG_ERROR_P(PARAM_V) << ": Number of 'v' values supplied is not the same as the model age spread.\n"
          << "Expected: " << model_->age_spread() << " but got " << v_.size();
    }
  } else {
    if (v_.size() != model_->length_bin_mid_points().size()) {
      LOG_ERROR_P(PARAM_V) << ": Number of 'v' values supplied is not the same as the model length bin count.\n"
          << "Expected: " << model_->length_bin_mid_points().size() << " but got " << v_.size();
    }
  }
}

/**
 * Reset this selectivity so it's ready for the next execution
 * phase in the model.
 *
 * This method will rebuild the cache of selectivity values
 * for each age in the model.
 */
void AllValues::RebuildCache() {
  if (not length_based_) {
    unsigned min_age = model_->min_age();
    for (unsigned i = 0; i < v_.size(); ++i) {
      if (v_[i] < 0.0)
        LOG_FATAL_P(PARAM_V) << "cannot have value < 0.0 in this class. Found value = " << v_[i] << " for age = " << min_age + i;
      values_[i] = v_[i];
    }
  } else {
    for (unsigned i = 0; i < v_.size(); ++i) {
      length_values_[i] = v_[i];
    }
  }
}


} /* namespace selectivities */
} /* namespace niwa */
