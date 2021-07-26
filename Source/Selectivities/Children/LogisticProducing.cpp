/**
 * @file LogisticProducing.cpp
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
#include "LogisticProducing.h"

#include <boost/math/distributions/lognormal.hpp>
#include <cmath>

#include "Model/Model.h"

// namespaces
namespace niwa {
namespace selectivities {

/**
 * Explicit constructor
 */
LogisticProducing::LogisticProducing(Model* model)
: Selectivity(model) {

  parameters_.Bind<unsigned>(PARAM_L, &low_, "Low", "");
  parameters_.Bind<unsigned>(PARAM_H, &high_, "High", "");
  parameters_.Bind<double>(PARAM_A50, &a50_, "A50", "");
  parameters_.Bind<double>(PARAM_ATO95, &ato95_, "Ato95", "");
  parameters_.Bind<double>(PARAM_ALPHA, &alpha_, "Alpha", "", 1.0);

  //RegisterAsAddressable(PARAM_A50, &a50_);
  //RegisterAsAddressable(PARAM_ATO95, &ato95_);
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
void LogisticProducing::DoValidate() {
  if (low_ > a50_ || high_ < a50_ || low_ >= high_) {
    LOG_ERROR_P(PARAM_A50) << "low (l) can't be greater than a50, or high (h) cant be less than a50 and low can't be greater than high";
  }

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
void LogisticProducing::RebuildCache() {
  if (not length_based_) {
    for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age) {

      if (not include_zero_age_values_ & (age == 0)) {
        values_[age - min_index_] = 0;
      } else if (age < low_)
        values_[age - min_index_] = 0.0;
      else if (age >= high_)
        values_[age - min_index_] = alpha_;
      else if (age == low_)
        values_[age - min_index_] = 1.0 / (1.0 + pow(19.0, (a50_ - (double)age) / ato95_)) * alpha_;
      else {
        double lambda2 = 1.0 / (1.0 + pow(19.0, (a50_- ((double)age - 1)) / ato95_));
        if (lambda2 > 0.9999) {
          values_[age - min_index_] = alpha_;
        } else {
          double lambda1 = 1.0 / (1.0 + pow(19.0, (a50_ - (double)age) / ato95_));
          values_[age - min_index_] = (lambda1 - lambda2) / (1.0 - lambda2) * alpha_;
          LOG_FINEST() << "age = " << age << " lambda1 = " << lambda1 << " lambda2 = " << lambda2 << " value = " <<  values_[age];
        }
      }
    }
  } else {
    vector<double> length_bins = model_->length_bin_mid_points();

    for (unsigned length_bin_index = 0; length_bin_index < length_bins.size(); ++length_bin_index) {
      double temp = (double)length_bins[length_bin_index];
      if (temp < low_)
        length_values_[length_bin_index] = 0.0;
      else if (temp >= high_)
        length_values_[length_bin_index] = alpha_;
      else if (temp == low_)
        length_values_[length_bin_index] = 1.0 / (1.0 + pow(19.0, (a50_ - (double)temp) / ato95_)) * alpha_;
      else {
        double lambda2 = 1.0 / (1.0 + pow(19.0, (a50_- ((double)temp - 1)) / ato95_));
        if (lambda2 > 0.9999) {
          length_values_[length_bin_index] = alpha_;
        } else {
          double lambda1 = 1.0 / (1.0 + pow(19.0, (a50_ - (double)temp) / ato95_));
          length_values_[length_bin_index] = (lambda1 - lambda2) / (1.0 - lambda2) * alpha_;
          LOG_FINEST() << "length = " << temp << " lambda1 = " << lambda1 << " lambda2 = " << lambda2 << " value = " <<  length_values_[length_bin_index];
        }
      }
    }
  }
}


} /* namespace selectivities */
} /* namespace niwa */
