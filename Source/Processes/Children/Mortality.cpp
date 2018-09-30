/**
 * @file Mortality.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */

// headers
#include "Mortality.h"

#include "Layers/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
Mortality::Mortality(Model* model) : Process(model) {
  process_type_ = ProcessType::kMortality;

}

void Mortality::DoReset() {
  LOG_FINE() << "clearing containers";
  removals_by_age_and_area_.clear();
  removals_by_length_and_area_.clear();
  removals_census_.clear();
  removals_tag_recapture_.clear();
}


} /* namespace processes */
} /* namespace niwa */
