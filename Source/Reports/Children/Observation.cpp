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
  parameters_.Bind<bool>(PARAM_SIMULATE_IN_NEW_FILE, &simulate_new_file_, "Produce sets of new files that are for each set of simulations", "", false);
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

  if (simulate_new_file_ & !parameters_.Get(PARAM_FILE_NAME)->has_been_defined())
    LOG_ERROR_P(PARAM_SIMULATE_IN_NEW_FILE) << "You need to supply the subcommand " << PARAM_FILE_NAME << " with this option";
}

/**
 *	Execute the report
 */
void Observation::DoExecute() {
  LOG_FINE();
  map<unsigned,map<string,vector<obs::Comparison> > >& comparisons = observation_->comparisons();

  if (not simulate_new_file_) {
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
    ready_for_writing_ = true;
  } else {
    LOG_FINE() << "generate simulated style report";
    cache_ << CONFIG_SECTION_SYMBOL << PARAM_OBSERVATION << " " << label_ << "\n";
    bool biomass_abundance_obs = false;
    ParameterList& parameter_list = observation_->parameters();
    const map<string, Parameter*>& parameters = parameter_list.parameters();
    for (auto iter = parameters.begin(); iter != parameters.end(); ++iter) {
      if (iter->first == PARAM_SIMULATION_LIKELIHOOD) {
        if (iter->second->values()[0] == PARAM_PSEUDO)
          cache_ << PARAM_LIKELIHOOD << " " << parameter_list.Get(PARAM_OBSERVATION)->values()[0] << "\n";
        else
          cache_ << PARAM_LIKELIHOOD << " " << iter->second->values()[0] << "\n";

        continue;
      }

      if (iter->first == PARAM_OBS || iter->first == PARAM_ERROR_VALUE || iter->first == PARAM_LABEL)
        continue;

      if (iter->second->values().size() > 0) {
        cache_ << iter->first << " ";
        const vector<string>& values = iter->second->values();
        for (string value : values) {
          if (iter->first == PARAM_TYPE)
            LOG_FINE() << "Type of observation simulating = " << value;
          if ((iter->first == PARAM_TYPE && value == PARAM_BIOMASS) || (iter->first == PARAM_TYPE && value == PARAM_ABUNDANCE)) {
            biomass_abundance_obs = true;
          }
          cache_ << value << " ";
        }
        cache_ << "\n";
      }
    }
    if (biomass_abundance_obs) {
      // biomass obs
      cache_ << PARAM_OBS << " ";
      for (auto iter = comparison.begin(); iter != comparison.end(); ++iter) {
        for (obs::Comparison comparison : iter->second)
          cache_ << comparison.simualted_ << " ";
      }
      cache_ << "\n";
    } else {
      // proportion at age obs
      cache_ << PARAM_TABLE << " " << PARAM_OBS << "\n";
      for (auto iter = comparison.begin(); iter != comparison.end(); ++iter) {
        cache_ << iter->first << " ";
        for (obs::Comparison comparison : iter->second) {
          cache_ << comparison.simualted_  << " ";
        }
        cache_ << "\n";
      }
      cache_ << PARAM_END_TABLE << "\n";
    }
    // Print Error values
    if (biomass_abundance_obs) {
      // biomass error
      cache_ << PARAM_ERROR_VALUE << " ";
      for (auto iter = comparison.begin(); iter != comparison.end(); ++iter) {
        for (obs::Comparison comparison : iter->second)
          cache_ << comparison.error_value_ << " ";
      }
      cache_ << "\n";
    } else {
      // proportion at age obs
      cache_ << PARAM_TABLE << " " << PARAM_ERROR_VALUES << "\n";
      for (auto iter = comparison.begin(); iter != comparison.end(); ++iter) {
        cache_ << iter->first << " ";
        for (obs::Comparison comparison : iter->second) {
          cache_ << comparison.error_value_ << " ";
        }
        cache_ << "\n";
      }
      cache_ << PARAM_END_TABLE << "\n";
    }
    cache_ << "\n";
    ready_for_writing_ = true;
  }
}

} /* namespace reports */
} /* namespace niwa */
