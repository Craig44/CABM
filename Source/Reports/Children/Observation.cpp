/**
 * @file Observation.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 21/01/2014
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
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
  model_state_ = State::kIterationComplete;
  skip_tags_   = false;
  parameters_.Bind<string>(PARAM_OBSERVATION, &observation_label_, "Observation label", "");
  //parameters_.Bind<bool>(PARAM_SIMULATE_IN_NEW_FILE, &simulate_new_file_, "Produce sets of new files that are for each set of simulations", "", false);
}


/**
 * Validate method
 */
void Observation::DoValidate() {
  if ((model_->run_mode() == RunMode::kMSE) || (model_->run_mode() == RunMode::kBasic) || (model_->run_mode() == RunMode::kSimulation)) {
    run_mode_ = model_->run_mode();
    //model_state_ = State::kIterationComplete;
  }
  //if(model_->run_mode() == RunMode::kMSE)
  //  model_state_ = State::kHCR;
}

/**
 *
 */
void Observation::DoBuild() {
  LOG_FINE();
  observation_ = model_->managers().observation()->GetObservation(observation_label_);
  if (!observation_) {
    auto observations = model_->managers().observation()->objects();
    for (auto observation : observations)
      LOG_FINE() << observation->label() << endl;
    LOG_ERROR_P(PARAM_OBSERVATION) << " (" << observation_label_ << ") could not be found. Have you defined it?";
  }
}

/**
 *	Execute the report
 */
void Observation::DoExecute() {
  LOG_FINE();
  map<unsigned,map<string,vector<obs::Comparison> > >& comparisons = observation_->comparisons();
  LOG_FINE() << "generate normal observation report";
  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";
  cache_ << "observation_type: " << observation_->type() << "\n";
  cache_ << "likelihood: " << observation_->likelihood() << "\n";
  cache_ << "Values " <<REPORT_R_DATAFRAME <<"\n";
  // report raw residuals
  cache_ << "year cell sex age length expected simulated error_value cell_biomass\n";
  for (auto iter = comparisons.begin(); iter != comparisons.end(); ++iter) {
    for (auto second_iter = iter->second.begin(); second_iter != iter->second.end(); ++second_iter) {
      for (obs::Comparison comparison : second_iter->second) {
        cache_ << iter->first << " " << second_iter->first << " " << comparison.sex_ << " "  << comparison.age_ << " " << comparison.length_ << " " << comparison.expected_ << " " << comparison.simulated_ << " " << comparison.error_value_ << " " << comparison.cell_biomass_ << "\n";
      }
    }
  }
  observation_->FillReportCache(cache_);

  ready_for_writing_ = true;
}

} /* namespace reports */
} /* namespace niwa */
