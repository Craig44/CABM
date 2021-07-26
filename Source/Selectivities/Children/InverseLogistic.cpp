/**
 * @file InverseLogistic.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 15/01/2013
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers
#include "InverseLogistic.h"

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
InverseLogistic::InverseLogistic(Model* model)
: Selectivity(model) {

  parameters_.Bind<double>(PARAM_A50, &a50_, "A50", "");
  parameters_.Bind<double>(PARAM_ATO95, &ato95_, "aTo95", "");
  parameters_.Bind<double>(PARAM_ALPHA, &alpha_, "Alpha", "", 1.0);

  RegisterAsAddressable(PARAM_A50, &a50_);
  RegisterAsAddressable(PARAM_ATO95, &ato95_);
  RegisterAsAddressable(PARAM_ALPHA, &alpha_);
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
void InverseLogistic::DoValidate() {
  if (alpha_ <= 0.0)
    LOG_ERROR_P(PARAM_ALPHA) << ": alpha (" << alpha_ << ") cannot be less than or equal to 0.0";
  if (ato95_ <= 0.0)
    LOG_ERROR_P(PARAM_ATO95) << ": ato95 (" << ato95_ << ") cannot be less than or equal to 0.0";
}

/**
 * Reset this selectivity so it's ready for the next execution
 * phase in the model.
 *
 * This method will rebuild the cache of selectivity values
 * for each age in the model.
 */
void InverseLogistic::RebuildCache() {
  if (not length_based_) {
    double threshold = 0.0;

    for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age) {
      double temp = (double)age;
      threshold = (double)(a50_ - temp) / ato95_;

      if (not include_zero_age_values_ & (age == 0)) {
        values_[age - min_index_] = 0;
      } else if (threshold > 5.0)
        values_[age - min_index_] = alpha_;
      else if (threshold < -5.0)
        values_[age - min_index_] = 0.0;
      else
        values_[age - min_index_] = alpha_ - (alpha_ / (1.0 + pow(19.0, threshold)));
    }
  } else {
    double threshold = 0.0;
    vector<double> length_bins = model_->length_bin_mid_points();

    for (unsigned length_bin_index = 0; length_bin_index < length_bins.size(); ++length_bin_index) {
      double temp = (double)length_bins[length_bin_index];
      threshold = (double)(a50_ - temp) / ato95_;
      if (threshold > 5.0)
        length_values_[length_bin_index] = alpha_;
      else if (threshold < -5.0)
        length_values_[length_bin_index] = 0.0;
      else
        length_values_[length_bin_index] = alpha_ - (alpha_ / (1.0 + pow(19.0, threshold)));
    }
  }
}

} /* namespace selectivities */
} /* namespace niwa */
