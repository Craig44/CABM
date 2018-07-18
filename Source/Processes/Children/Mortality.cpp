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


} /* namespace processes */
} /* namespace niwa */
