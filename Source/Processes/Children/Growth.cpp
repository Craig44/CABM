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
#include "TimeSteps/Manager.h"

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
  parameters_.Bind<bool>(PARAM_UPDATE_GROWTH_PARAMETERS, &update_growth_parameters_, "If an agent/individual moves do you want it to take on the growth parameters of the new spatial cell", "");
  parameters_.Bind<float>(PARAM_CV, &cv_, "The cv of the distribution", "");
}


void Growth::DoBuild() {
  LOG_FINE();
  vector<TimeStep*> time_steps = model_->managers().time_step()->ordered_time_steps();
  LOG_FINEST() << "time_steps.size(): " << time_steps.size();
  vector<unsigned> active_time_steps;
  for (unsigned i = 0; i < time_steps.size(); ++i) {
    if (time_steps[i]->HasProcess(label_))
      active_time_steps.push_back(i);
  }

  if (time_step_proportions_.size() == 0) {
    for (unsigned i : active_time_steps)
      time_step_proportions_[i] = 1.0;
  } else {
    if (time_step_proportions_.size() != active_time_steps.size())
      LOG_FATAL_P(PARAM_TIME_STEP_PROPORTIONS) << " length (" << time_step_proportions_.size()
          << ") does not match the number of time steps this process has been assigned to (" << active_time_steps.size() << ")";

    for (float value : time_step_proportions_) {
      if (value < 0.0 || value > 1.0)
        LOG_ERROR_P(PARAM_TIME_STEP_PROPORTIONS) << " value (" << value << ") must be between 0.0 (exclusive) and 1.0 (inclusive)";
    }

    for (unsigned i = 0; i < time_step_proportions_.size(); ++i) {
      time_step_ratios_[active_time_steps[i]] = time_step_proportions_[i];
      LOG_FINE() << "setting growth in time step " << active_time_steps[i] << " = " << time_step_proportions_[i];
    }
  }
}

} /* namespace processes */
} /* namespace niwa */
