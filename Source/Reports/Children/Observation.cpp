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
  //parameters_.Bind<bool>(PARAM_SIMULATE_IN_NEW_FILE, &simulate_new_file_, "Produce sets of new files that are for each set of simulations", "", false);
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
      cout << observation->label() << endl;
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
  cache_ << "year cell sex age length expected simulated error_value\n";
  for (auto iter = comparisons.begin(); iter != comparisons.end(); ++iter) {
    for (auto second_iter = iter->second.begin(); second_iter != iter->second.end(); ++second_iter) {
      for (obs::Comparison comparison : second_iter->second) {
        cache_ << iter->first << " " << second_iter->first << " " << comparison.sex_ << " "  << comparison.age_ << " " << comparison.length_ << " " << comparison.expected_ << " " << comparison.simulated_ << " " << comparison.error_value_ << "\n";
      }
    }
  }

  if (observation_->type() == PARAM_MORTALITY_SCALES_AGE_FREQUENCY || observation_->type() == PARAM_MORTALITY_EVENT_BIOMASS_SCALED_AGE_FREQUENCY) {
    // Print the age length key for curiosity
    map<unsigned,map<string,vector<float>>>& LF = observation_->get_lf();
    map<unsigned,map<string,vector<vector<float>>>>& ALK = observation_->get_alk();
    cache_ << "length_freq_by_year_stratum "  <<REPORT_R_DATAFRAME <<"\n";
    cache_ << "year_stratum ";
    for (auto len : model_->length_bin_mid_points())
      cache_ << len << " ";
    cache_ << "\n";

    for(auto& len_freq_by_year : LF) {
      for(auto len_freq_by_year_strata : len_freq_by_year.second) {
        cache_ << len_freq_by_year.first << "_" << len_freq_by_year_strata.first << " ";
        for (unsigned i = 0; i < len_freq_by_year_strata.second.size(); ++i) {
          cache_ << len_freq_by_year_strata.second[i] << " ";
        }
        cache_ << "\n";
      }
    }




    for(auto& alk_by_year : ALK) {
      for(auto alk_by_strata : alk_by_year.second) {
        cache_ << "ALK_" << alk_by_year.first << "_" << alk_by_strata.first << " " <<REPORT_R_MATRIX <<"\n";
        for (unsigned i = 0; i < alk_by_strata.second.size(); ++i) {
          for (unsigned j = 0; j < alk_by_strata.second[i].size(); ++j) {
            cache_ << alk_by_strata.second[i][j] << " ";
          }
          cache_ << "\n";
        }
      }
    }


  }
  ready_for_writing_ = true;
}

} /* namespace reports */
} /* namespace niwa */
