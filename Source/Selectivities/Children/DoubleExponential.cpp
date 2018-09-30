/**
 * @file floatExponential.cpp
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
#include "DoubleExponential.h"

#include <boost/math/distributions/lognormal.hpp>
#include <cmath>

#include "Model/Model.h"
#include "TimeSteps/Manager.h"

// Namespaces
namespace niwa {
namespace selectivities {

/**
 * Explicit Constructor
 */
DoubleExponential::DoubleExponential(Model* model)
: Selectivity(model) {

  parameters_.Bind<float>(PARAM_X0, &x0_, "X0", "");
  parameters_.Bind<float>(PARAM_X1, &x1_, "X1", "");
  parameters_.Bind<float>(PARAM_X2, &x2_, "X2", "");
  parameters_.Bind<float>(PARAM_Y0, &y0_, "Y0", "");
  parameters_.Bind<float>(PARAM_Y1, &y1_, "Y1", "");
  parameters_.Bind<float>(PARAM_Y2, &y2_, "Y2", "");
  parameters_.Bind<float>(PARAM_ALPHA, &alpha_, "Alpha", "", 1.0);

  //RegisterAsAddressable(PARAM_X0, &x0_);
  //RegisterAsAddressable(PARAM_Y0, &y0_);
  //RegisterAsAddressable(PARAM_Y1, &y1_);
  //RegisterAsAddressable(PARAM_Y2, &y2_);
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
void DoubleExponential::DoValidate() {
  // Param: x0, x1, x2 - Check that x1 is between x0 and x2
  if (x0_ < x1_ || x0_ > x2_)
    LOG_ERROR_P(PARAM_X0) << "x0 ( " << x0_ << ") must be between x1 (" << x1_ << ") and x2 (" << x2_ << ")";

  // Param: y0, y1, y2
  if (y0_ < 0.0)
    LOG_ERROR_P(PARAM_Y0) << ": y0 (" << y0_ << ") is less than 0.0";
  if (y1_ < 0.0)
    LOG_ERROR_P(PARAM_Y1) << ": y1 (" << y1_ << ") is less than 0.0";
  if (y2_ < 0.0)
    LOG_ERROR_P(PARAM_Y2) << ": y2 (" << y2_ << ") is less than 0.0";

  // Param: alpha
  if (alpha_ <= 0.0)
    LOG_ERROR_P(PARAM_ALPHA) << ": alpha (" << alpha_ << ") is less than or equal to 0.0";
}

/**
 * Reset this selectivity so it's ready for the next execution
 * phase in the model.
 *
 * This method will rebuild the cache of selectivity values
 * for each age in the model.
 */
void DoubleExponential::RebuildCache() {
  if (not length_based_) {
    for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age) {
      if (not include_zero_age_values_ & (age == 0)) {
        values_[age - min_index_] = 0;
      } else if ((float)age <= x0_) {
        values_[age - min_index_] = alpha_ * y0_ * pow((y1_ / y0_), ((float)age - x0_)/(x1_ - x0_));
      } else if ((float)age > x0_ && (float)age <= x2_) {
        values_[age - min_index_] = alpha_ * y0_ * pow((y2_ / y0_), ((float)age - x0_)/(x2_ - x0_));
      } else {
        values_[age - min_index_] = y2_;
      }
    }
  } else {
    vector<unsigned> length_bins = model_->length_bins();
    for (unsigned length_bin_index = 0; length_bin_index < length_bins.size(); ++length_bin_index) {
      if ((float)length_bins[length_bin_index] <= x0_) {
        length_values_[length_bin_index] = alpha_ * y0_ * pow((y1_ / y0_), ((float)length_bins[length_bin_index] - x0_)/(x1_ - x0_));
      } else if ((float)length_bins[length_bin_index] > x0_ && (float)length_bins[length_bin_index] <= x2_) {
        length_values_[length_bin_index] = alpha_ * y0_ * pow((y2_ / y0_), ((float)length_bins[length_bin_index] - x0_)/(x2_ - x0_));
      } else {
        length_values_[length_bin_index] = y2_;
      }
    }
  }
}

} /* namespace selectivities */
} /* namespace niwa */
