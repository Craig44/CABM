/**
 * @file SimulatedObservation.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 8/05/2014
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 */

// headers
#include "SimulatedObservation.h"

#include "Observations/Manager.h"
#include "Observations/Observation.h"

// namespaces
namespace niwa {
namespace reports {

/**
 * default constructor
 */
SimulatedObservation::SimulatedObservation(Model* model) : Report(model) {
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
  model_state_ = State::kIterationComplete;
  skip_tags_   = true;

  parameters_.Bind<string>(PARAM_OBSERVATION, &observation_label_, "Observation label", "");
  parameters_.Bind<string>(PARAM_CELLS, &cells_, "Cells to aggregate the observaton over.", "", true);

}

/**
 * build method
 */
void SimulatedObservation::DoBuild() {
  observation_ = model_->managers().observation()->GetObservation(observation_label_);
  if (!observation_)
    LOG_ERROR_P(PARAM_OBSERVATION) << "(" << observation_label_ << ") could not be found. Have you defined it?";
}

/**
 * execute method
 */
void SimulatedObservation::DoExecute() {
  LOG_FINE() << "generate simulated style report";
  cache_ << CONFIG_SECTION_SYMBOL << PARAM_OBSERVATION << " " << label_ << "\n";
  map<unsigned,map<string,vector<obs::Comparison> > >& comparisons = observation_->comparisons();

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
    for (auto iter = comparisons.begin(); iter != comparisons.end(); ++iter) {
      for (auto second_iter = iter->second.begin(); second_iter != iter->second.end(); ++second_iter) {
        for (obs::Comparison comparison : second_iter->second) {
          cache_ << comparison.simulated_ << " ";
        }
      }
    }
    cache_ << "\n";
  } else {
    // proportion at age obs
    cache_ << PARAM_TABLE << " " << PARAM_OBS << "\n";
    for (auto iter = comparisons.begin(); iter != comparisons.end(); ++iter) {
      cache_ << iter->first << " ";
      for (auto second_iter = iter->second.begin(); second_iter != iter->second.end(); ++second_iter) {
        for (obs::Comparison comparison : second_iter->second) {
          cache_ << comparison.simulated_  << " ";
        }
      }
      cache_ << "\n";
    }
    cache_ << PARAM_END_TABLE << "\n";
  }
  // Print Error values
  if (biomass_abundance_obs) {
    // biomass error
    cache_ << PARAM_ERROR_VALUE << " ";
    for (auto iter = comparisons.begin(); iter != comparisons.end(); ++iter) {
      cache_ << " ";
      for (auto second_iter = iter->second.begin(); second_iter != iter->second.end(); ++second_iter) {
        for (obs::Comparison comparison : second_iter->second) {
          cache_ << comparison.error_value_ << " ";
        }
      }
    }
    cache_ << "\n";
  } else {
    // proportion at age obs
    cache_ << PARAM_TABLE << " " << PARAM_ERROR_VALUES << "\n";
    for (auto iter = comparisons.begin(); iter != comparisons.end(); ++iter) {
      cache_ << iter->first << " ";
      for (auto second_iter = iter->second.begin(); second_iter != iter->second.end(); ++second_iter) {
        for (obs::Comparison comparison : second_iter->second) {
          cache_ << comparison.error_value_ << " ";
        }
      }
      cache_ << "\n";
    }
    cache_ << PARAM_END_TABLE << "\n";
  }
  cache_ << "\n";
  ready_for_writing_ = true;


}

} /* namespace reports */
} /* namespace niwa */