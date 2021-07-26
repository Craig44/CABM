/**
 * @file KnifeEdge.cpp
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
#include "KnifeEdge.h"

#include <boost/math/distributions/lognormal.hpp>
#include <cmath>

#include "Model/Model.h"
#include "TimeSteps/Manager.h"

// namespaces
namespace niwa {
namespace selectivities {

/**
 * Explicit Constructor
 */
KnifeEdge::KnifeEdge(Model* model)
: Selectivity(model) {

  parameters_.Bind<double>(PARAM_E, &edge_, "Edge", "");
  parameters_.Bind<double>(PARAM_ALPHA, &alpha_, "Alpha", "", 1.0);

  //RegisterAsAddressable(PARAM_ALPHA, &alpha_);
 // RegisterAsAddressable(PARAM_E, &edge_);
}

/**
 * Reset this selectivity so it's ready for the next execution
 * phase in the model.
 *
 * This method will rebuild the cache of selectivity values
 * for each age in the model.
 */
void KnifeEdge::RebuildCache() {
  if (not length_based_) {
    for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age) {
      double temp = age * 1.0;
      if (not include_zero_age_values_ & (age == 0)) {
        values_[age - min_index_] = 0;
      } else if (temp >= edge_)
        values_[age - min_index_] = alpha_;
      else
        values_[age - min_index_] = 0.0;
    }
  } else {
    vector<double> length_bins = model_->length_bin_mid_points();

    for (unsigned length_bin_index = 0; length_bin_index < length_bins.size(); ++length_bin_index) {
      double temp = (double)length_bins[length_bin_index];
      if (temp >= edge_)
        length_values_[length_bin_index] = alpha_;
      else
        length_values_[length_bin_index] = 0.0;
    }
  }
}

} /* namespace selectivities */
} /* namespace niwa */
