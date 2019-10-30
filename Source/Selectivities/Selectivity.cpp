/**
 * @file Selectivity.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 11/01/2013
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers
#include "Selectivity.h"


#include "Model/Model.h"
#include <boost/math/distributions/normal.hpp>

// Namesapces
namespace niwa {

/**
 * Explicit Constructor
 */
Selectivity::Selectivity(Model* model)
: model_(model) {

  parameters_.Bind<string>(PARAM_LABEL, &label_, "The label for this selectivity", "");
  parameters_.Bind<string>(PARAM_TYPE, &type_, "The type of selectivity", "");
  parameters_.Bind<bool>(PARAM_LENGTH_BASED, &length_based_, "Is the selectivity length based", "", false);
  parameters_.Bind<bool>(PARAM_INCLUDE_AGE_ZERO_INDIVIDUALS, &include_zero_age_values_, "Include 0 aged fish in selectivity (more for comparing with population models that start modelling fish at age = 1)", "", false);

}

/**
 *
 */
void Selectivity::Validate() {
  parameters_.Populate(model_);
  min_index_ = model_->min_age();

  DoValidate();

  if (not length_based_) {
    values_.assign(model_->age_spread(), 0.0);
  } else {
    vector<float> lengths = model_->length_bin_mid_points();
    length_values_.assign(lengths.size(), 0.0);
    LOG_FINE() << "number of bins = " << length_values_.size();
  }

}


/**
 *
 */
void Selectivity::Reset() {
  RebuildCache();
}

/**
 * Return the cached value for the specified age or length from
 * our internal map
 *
 * @param age_or_length The age or length to get selectivity value for
 * @return The value stored in the map or 0.0 as default
 */

float Selectivity::GetResult(unsigned age_or_length) {
  if (not length_based_)
    return values_[age_or_length];
  else
    return length_values_[age_or_length];
}



} /* namespace niwa */
