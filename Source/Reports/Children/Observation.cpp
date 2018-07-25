/**
 * @file Observation.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 21/01/2014
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 */

// headers
#include "Observation.h"

#include "Observations/Manager.h"
#include "Observations/Observation.h"

// namespaces
namespace niwa {
namespace reports {

namespace obs = niwa::observations;

/**
 *
 */
Observation::Observation(Model* model) : Report(model) {
  LOG_TRACE();
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
  model_state_ = (State::Type)(State::kIterationComplete);
  parameters_.Bind<string>(PARAM_OBSERVATION, &observation_label_, "Observation label", "");
}

/**
 *
 */
void Observation::DoBuild() {
  LOG_TRACE();
  observation_ = model_->managers().observation()->GetObservation(observation_label_);
  if (!observation_) {
    auto observations = model_->managers().observation()->objects();
    for (auto observation : observations)
      cout << observation->label() << endl;
    LOG_ERROR_P(PARAM_OBSERVATION) << " (" << observation_label_ << ") could not be found. Have you defined it?";
  }
}

/**
 *	Execute the report
 */
void Observation::DoExecute() {
  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";
  cache_ << "observation_type: " << observation_->type() << "\n";
  cache_ << "likelihood: " << observation_->likelihood() << "\n";
  cache_ << "Values " <<REPORT_R_DATAFRAME <<"\n";
  map<unsigned, vector<obs::Comparison>>& comparisons = observation_->comparisons();
  // report raw residuals
  cache_ << "year age length expected simulated error_value\n";
  for (auto iter = comparisons.begin(); iter != comparisons.end(); ++iter) {
    for (obs::Comparison comparison : iter->second) {
      cache_ << iter->first << " " << comparison.age_ << " " << comparison.length_ << " " << comparison.expected_ << " " << comparison.simulated_ << " " << comparison.error_value_ << "\n";
    }
  }
  ready_for_writing_ = true;
}

} /* namespace reports */
} /* namespace niwa */
