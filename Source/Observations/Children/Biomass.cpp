/**
 * @file Biomass.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 7/01/2014
 * @section LICENSE
 *
 * Copyright NIWA Science @2013 - www.niwa.co.nz
 *
 */

// headers
#include "Biomass.h"

//#include "Catchabilities/Manager.h"

#include "TimeSteps/Manager.h"
#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "Likelihoods/Manager.h"

#include "World/WorldView.h"
#include "World/WorldCell.h"

#include <omp.h>

#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/To.h"
#include "Utilities/Map.h"
#include "Utilities/DoubleCompare.h"

// namespaces
namespace niwa {
namespace observations {

namespace utils = niwa::utilities;

/**
 * Default constructor
 */
Biomass::Biomass(Model* model) : Observation(model) {
  parameters_.Bind<float>(PARAM_CATCHABILITY, &catchability_value_, "The Catchability multiplier", ""); // TODO not sure if neccessary
  parameters_.Bind<string>(PARAM_TIME_STEP, &time_step_label_, "The label of time-step that the observation occurs in", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "The years of the observed values", "");
  parameters_.Bind<float>(PARAM_ERROR_VALUE, &error_values_, "The error values of the observed values (note the units depend on the likelihood)", "");
  parameters_.Bind<string>(PARAM_SELECTIVITIES, &selectivity_labels_, "Labels of the selectivities", "", true);
  parameters_.Bind<float>(PARAM_PROPORTION_TRHOUGH_MORTALITY, &time_step_proportion_, "Proportion through the mortality block of the time step to infer observation with", "", float(0.5))->set_range(0.0, 1.0);
  parameters_.Bind<string>(PARAM_LAYER_OF_CELLS, &layer_label_, "The layer that indicates what area to summarise observations over.", "");
  parameters_.Bind<string>(PARAM_CELLS, &cells_, "The cells we want to generate observations for from the layer of cells supplied", "");
  parameters_.Bind<string>(PARAM_SIMULATION_LIKELIHOOD, &simulation_likelihood_label_, "Simulation likelihood to use", "");
  parameters_.Bind<bool>(PARAM_ABUNDANCE, &obs_abundance_, "Is the index biomass (abundance = false) or abundance (abundance = true)", "", false);
  RegisterAsAddressable(PARAM_CATCHABILITY, &catchability_value_);

  allowed_likelihood_types_.push_back(PARAM_NORMAL);
  allowed_likelihood_types_.push_back(PARAM_LOGNORMAL);
  allowed_likelihood_types_.push_back(PARAM_PSEUDO);
}

/**
 *
 */
void Biomass::DoValidate() {
  LOG_TRACE();
  // Error Value
  if (error_values_.size() == 1 && years_.size() > 1)
    error_values_.assign(years_.size(), error_values_[0]);
  if (error_values_.size() != years_.size())
    LOG_ERROR_P(PARAM_ERROR_VALUE) << ": error_value length (" << error_values_.size()
        << ") must be same length as years (" << years_.size() << ")";

  error_values_by_year_ = utils::Map::create(years_, error_values_);
  for (unsigned i = 0; i < years_.size(); ++i) {
    // Check error values
    if (error_values_[i] <= 0.0)
      LOG_ERROR_P(PARAM_ERROR_VALUE) << "for year '" << years_[i] << "' we found an error term that is less than or equal to 0.0, this is not allowed please check it";
  }
}

/**
 *
 */
void Biomass::DoBuild() {
  LOG_MEDIUM();

  likelihood_ = model_->managers().likelihood()->GetOrCreateLikelihood(model_, label_, simulation_likelihood_label_);
  if (!likelihood_) {
    LOG_FATAL_P(PARAM_SIMULATION_LIKELIHOOD) << "(" << simulation_likelihood_label_ << ") could not be found or constructed.";
    return;
  }
  if (std::find(allowed_likelihood_types_.begin(), allowed_likelihood_types_.end(), likelihood_->type()) == allowed_likelihood_types_.end()) {
    string allowed = boost::algorithm::join(allowed_likelihood_types_, ", ");
    LOG_FATAL_P(PARAM_SIMULATION_LIKELIHOOD) << ": likelihood " << likelihood_->type() << " is not supported by the " << type_ << " observation."
        << " Allowed types are: " << allowed;
  }
  //catchability_ = model_->managers().catchability()->GetCatchability(catchability_label_);
  //if (!catchability_)
  //  LOG_FATAL_P(PARAM_CATCHABILITY) << ": catchability " << catchability_label_ << " could not be found. Have you defined it?";


  // Build and validate layers
  layer_ = model_->managers().layer()->GetCategoricalLayer(layer_label_);
  if (!layer_)
    LOG_FATAL_P(PARAM_LAYER_OF_CELLS) << "could not find layer " << layer_label_ << " does it exist?, if it exists is of type categorical?";


  if (selectivity_labels_.size() > 2)
    LOG_ERROR_P(PARAM_SELECTIVITIES) << "you supplied " << selectivity_labels_.size() << " you cannot supply over 2 selectivities one for each sex";

  // Build Selectivity pointers
  bool first = true;
  for (auto label : selectivity_labels_) {
    Selectivity* temp_selectivity = model_->managers().selectivity()->GetSelectivity(label);
    if (!temp_selectivity)
      LOG_ERROR_P(PARAM_SELECTIVITIES) << ": selectivity " << label << " does not exist. Have you defined it?";

    selectivities_.push_back(temp_selectivity);
    if (first) {
      first = false;
      selectivity_length_based_ = temp_selectivity->is_length_based();
    } else {
      if (selectivity_length_based_ != temp_selectivity->is_length_based()) {
        LOG_ERROR_P(PARAM_SELECTIVITIES) << "The selectivity  " << label << " was not the same type (age or length based) as the previous selectivity label";
      }
    }
  }
  if (selectivities_.size() == 1)
    selectivities_.assign(2, selectivities_[0]);

  LOG_FINEST() << "finished building selectivities";
  world_ = model_->world_view();
  if (!world_)
    LOG_CODE_ERROR() << "!world_ could not create pointer to world view model, something is wrong";

  LOG_FINEST() << "finished building pointer to world";

  // Check all the cells supplied are in the layer
  for (auto cell :  cells_) {
    bool cell_found = false;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        if (layer_->get_value(row,col) == cell)
          cell_found = true;
      }
    }
    if (not cell_found)
      LOG_ERROR_P(PARAM_CELLS) << "could not find the cell '" << cell << "' in the layer " << layer_label_ << " please make sure that you supply cell labels that are consistent with the layer.";
  }

  LOG_FINEST() << "finished checking cells";


  for (auto year : years_) {
    if((year < model_->start_year()) || (year > model_->final_year()))
      LOG_FATAL_P(PARAM_YEARS) << "Years can't be less than start_year (" << model_->start_year() << "), or greater than final_year (" << model_->final_year() << "). Please fix this.";
  }

  // Subscribe this observation to the timestep
  auto time_step = model_->managers().time_step()->GetTimeStep(time_step_label_);
  if (!time_step) {
    LOG_ERROR_P(PARAM_TIME_STEP) << time_step_label_ << " could not be found. Have you defined it?";
  } else {
    for (unsigned year : years_)
      time_step->SubscribeToBlock(this, year);
  }

  LOG_FINEST() << "finished subscribing to time step";

}

void Biomass::DoReset() {
 // clear containers
  for (auto& year : years_) {
    for(auto& cell : cells_) {
      pre_obs_values_by_year_[year][cell] = 0.0;
      obs_values_by_year_[year][cell] = 0.0;
    }
  }

}
/**
 *
 */
void Biomass::PreExecute() {
  LOG_FINE() << label_;
  unsigned year = model_->current_year();
  if (find(years_.begin(), years_.end(), year) != years_.end()) {

    utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
    if (utilities::doublecompare::IsOne(time_step_proportion_)) {
      return void();
    } else {
      if (selectivity_length_based_) {
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            string cell_label = layer_->get_value(row, col);
            if (find(cells_.begin(), cells_.end(), cell_label) == cells_.end())
              continue;
            WorldCell* cell = world_->get_base_square(row, col);
            if (cell->is_enabled()) {
              LOG_FINEST() << "about to convert " << cell->agents_.size() << " through the maturity process";
              //unsigned counter = 1;
              float probability;
              if (obs_abundance_) {
                for (Agent& agent : cell->agents_) {
                  if (agent.is_alive()) {
                    probability = selectivities_[agent.get_sex()]->GetResult(agent.get_length_bin_index());
                    if (rng.chance() <= probability) {
                      pre_obs_values_by_year_[year][cell_label] += agent.get_scalar();
                    }
                  }
                }
              } else {
                for (Agent& agent : cell->agents_) {
                  if (agent.is_alive()) {
                    probability = selectivities_[agent.get_sex()]->GetResult(agent.get_length_bin_index());
                    if (rng.chance() <= probability) {
                      pre_obs_values_by_year_[year][cell_label] += agent.get_weight() * agent.get_scalar();
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            string cell_label = layer_->get_value(row, col);
            if (find(cells_.begin(), cells_.end(), cell_label) == cells_.end())
              continue;
            WorldCell* cell = world_->get_base_square(row, col);
            if (cell->is_enabled()) {
              LOG_FINEST() << "about to convert " << cell->agents_.size() << " through the maturity process";
              //unsigned counter = 1;
              float probability;
              if (obs_abundance_) {
                for (Agent& agent : cell->agents_) {
                  if (agent.is_alive()) {
                    probability = selectivities_[agent.get_sex()]->GetResult(agent.get_age_index());
                    if (rng.chance() <= probability) {
                      pre_obs_values_by_year_[year][cell_label] += agent.get_scalar();
                    }
                  }
                }
              } else {
                for (Agent& agent : cell->agents_) {
                  if (agent.is_alive()) {
                    probability = selectivities_[agent.get_sex()]->GetResult(agent.get_age_index());
                    if (rng.chance() <= probability) {
                      pre_obs_values_by_year_[year][cell_label] += agent.get_weight() * agent.get_scalar();
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

/**
 *
 */
void Biomass::Execute() {
  LOG_FINE() << label_;
  unsigned year = model_->current_year();
  if (find(years_.begin(), years_.end(), year) != years_.end()) {
    utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
    if (selectivity_length_based_) {
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          string cell_label = layer_->get_value(row, col);
          if (find(cells_.begin(), cells_.end(), cell_label) == cells_.end())
            continue;
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            LOG_FINEST() << "about to convert " << cell->agents_.size() << " through the maturity process";
            //unsigned counter = 1;
            float probability;
            if (obs_abundance_) {
              for (Agent& agent : cell->agents_) {
                if (agent.is_alive()) {
                  probability = selectivities_[agent.get_sex()]->GetResult(agent.get_length_bin_index());
                  if (rng.chance() <= probability) {
                    obs_values_by_year_[year][cell_label] += agent.get_scalar();
                  }
                }
              }
            } else {
              for (Agent& agent : cell->agents_) {
                if (agent.is_alive()) {
                  probability = selectivities_[agent.get_sex()]->GetResult(agent.get_length_bin_index());
                  if (rng.chance() <= probability) {
                    obs_values_by_year_[year][cell_label] += agent.get_weight() * agent.get_scalar();
                  }
                }
              }
            }
          }
        }
      }
    } else {
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          string cell_label = layer_->get_value(row, col);
          if (find(cells_.begin(), cells_.end(), cell_label) == cells_.end())
            continue;
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            LOG_FINEST() << "about to convert " << cell->agents_.size() << " through the maturity process";
            //unsigned counter = 1;
            float probability;
            if (obs_abundance_) {
              for (Agent& agent : cell->agents_) {
                if (agent.is_alive()) {
                  probability = selectivities_[agent.get_sex()]->GetResult(agent.get_age_index());
                  if (rng.chance() <= probability) {
                    obs_values_by_year_[year][cell_label] += agent.get_scalar();
                  }
                }
              }
            } else {
              for (Agent& agent : cell->agents_) {
                if (agent.is_alive()) {
                  probability = selectivities_[agent.get_sex()]->GetResult(agent.get_age_index());
                  if (rng.chance() <= probability) {
                    obs_values_by_year_[year][cell_label] += agent.get_weight() * agent.get_scalar();
                  }
                }
              }
            }
          }
        }
      }
    }

    for (auto& second_iter : obs_values_by_year_[year]) {
      if (utilities::doublecompare::IsZero(time_step_proportion_)) {
        LOG_FINE() << label_ << " Biomass obs mort = 0.0";
        SaveComparison(0, 0, second_iter.first, pre_obs_values_by_year_[year][second_iter.first] * catchability_value_, 0.0, error_values_by_year_[year], year);
      } else if (utilities::doublecompare::IsOne(time_step_proportion_)) {
        LOG_FINE() << label_ << " Biomass obs mort = 1.0";
        SaveComparison(0, 0, second_iter.first, obs_values_by_year_[year][second_iter.first] * catchability_value_, 0.0, error_values_by_year_[year], year);
      } else {
        LOG_FINE() << label_ << " Biomass obs mort between 0.0 - 1.0";
        float value = 0.0;
        value = pre_obs_values_by_year_[year][second_iter.first] + ((obs_values_by_year_[year][second_iter.first] - pre_obs_values_by_year_[year][second_iter.first]) * time_step_proportion_);
        SaveComparison(0, 0, second_iter.first, value * catchability_value_, 0.0, error_values_by_year_[year], year);
      }
    }
  }
}

/**
 *
 */
void Biomass::Simulate() {
  /**
   * Simulate or generate results
   * During simulation mode we'll simulate results for this observation
   */
	LOG_FINEST() << "Calculating score for observation = " << label_;
  likelihood_->SimulateObserved(comparisons_);
}

/**
 * Extra Reporting
 */
void Biomass::FillReportCache(ostringstream& cache) {
  cache << "catchability: " << catchability_value_ << "\n";
}

} /* namespace observations */
} /* namespace niwa */

