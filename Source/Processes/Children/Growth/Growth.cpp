/**
 * @file Growth.cpp
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
#include "Growth.h"

#include "Layers/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
Growth::Growth(Model* model) : Process(model) {
  process_type_ = ProcessType::kGrowth;
  //parameters_.Bind<string>(PARAM_L_INF_LAYER_LABEL, &l_inf_layer_label_, "Label for the numeric layer that describes mean L_inf through space", ""); // Maybe add a cv by space?
  parameters_.Bind<float>(PARAM_TIME_STEP_PROPORTIONS, &time_step_proportions_, "A vector of proportions that describe how much of the annual increment to add to the lenght of an agent in each time step this process is applied", "");
  parameters_.Bind<string>(PARAM_DISTRIBUTION, &distribution_, "the distribution to allocate the parameters to the agents", "");
  parameters_.Bind<float>(PARAM_CV, &cv_, "The cv of the distribution", "");
}


} /* namespace processes */
} /* namespace niwa */
