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

  RegisterAsAddressable(PARAM_VALUES, &values_, addressable::kLookup);
  RegisterAsAddressable(PARAM_LENGTH_VALUES, &length_values_, addressable::kLookup);
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
    length_values_.assign(model_->length_bins().size(), 0.0);
  }
}


/**
 *
 */
void Selectivity::Reset() {

}

/**
 * Return the cached value for the specified age or length from
 * our internal map
 *
 * @param age_or_length The age or length to get selectivity value for
 * @return The value stored in the map or 0.0 as default
 */

Double Selectivity::GetResult(unsigned age_or_length) {
    return values_[age_or_length - min_index_];
}



} /* namespace niwa */
