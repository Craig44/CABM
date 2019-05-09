/**
 * @file ProportionsAtAge.cpp
 * @author  C.Marsh
 * @date 7/01/2014
 * @section LICENSE
 *
 * Copyright NIWA Science @2013 - www.niwa.co.nz
 *
 */

// headers
#include "ProportionsAtAge.h"


#include "TimeSteps/Manager.h"
#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "Layers/Manager.h"
#include "AgeingErrors/Manager.h"
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
ProportionsAtAge::ProportionsAtAge(Model* model) : Observation(model) {
  error_values_table_ = new parameters::Table(PARAM_ERROR_VALUES);
  parameters_.Bind<unsigned>(PARAM_MIN_AGE, &min_age_, "Minimum age", "");
  parameters_.Bind<unsigned>(PARAM_MAX_AGE, &max_age_, "Maximum age", "");
  parameters_.Bind<bool>(PARAM_SEXED, &sexed_, "Observation split by sex", "" , false);
  parameters_.Bind<string>(PARAM_SELECTIVITIES, &selectivity_labels_, "Labels of the selectivities", "");
  parameters_.Bind<float>(PARAM_PROPORTION_TRHOUGH_MORTALITY, &time_step_proportion_, "Proportion through the mortality block of the time step to infer observation with", "", float(0.5))->set_range(0.0, 1.0);
  parameters_.Bind<bool>(PARAM_PLUS_GROUP, &plus_group_, "max age is a plus group", "", true);
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years for which to calculate an observation", "");
  parameters_.Bind<string>(PARAM_AGEING_ERROR, &ageing_error_label_, "Label of ageing error to use", "", "");
  parameters_.BindTable(PARAM_ERROR_VALUES, error_values_table_, "Table of error values of the observed values (note the units depend on the likelihood)", "", false);
  parameters_.Bind<string>(PARAM_LAYER_OF_CELLS, &layer_label_, "The layer that indicates what area to summarise observations over.", "");
  parameters_.Bind<string>(PARAM_CELLS, &cells_, "The cells we want to generate observations for from the layer of cells supplied", "");
  parameters_.Bind<string>(PARAM_TIME_STEP, &time_step_label_, "The label of time-step that the observation occurs in", "");
  parameters_.Bind<string>(PARAM_SIMULATION_LIKELIHOOD, &simulation_likelihood_label_, "Simulation likelihood to use", "");

  allowed_likelihood_types_.push_back(PARAM_LOGNORMAL);
  allowed_likelihood_types_.push_back(PARAM_MULTINOMIAL);
  allowed_likelihood_types_.push_back(PARAM_DIRICHLET);
  allowed_likelihood_types_.push_back(PARAM_LOGISTIC_NORMAL);
  allowed_likelihood_types_.push_back(PARAM_PSEUDO);
}

/**
 * Destructor
 */
ProportionsAtAge::~ProportionsAtAge() {
  delete error_values_table_;
}

/**
 *
 */
void ProportionsAtAge::DoValidate() {
  LOG_MEDIUM();
  age_spread_ = (max_age_ - min_age_) + 1;

  if(sexed_) {
    if (not model_->get_sexed())
      LOG_ERROR_P(PARAM_SEXED) << "you have asked for a sexed observation in an unsexed model please set this to sexed false, are create a sexed model. Chairs";
  }
  /**
   * Do some simple checks
   */
  if (min_age_ < model_->min_age())
    LOG_ERROR_P(PARAM_MIN_AGE) << ": min_age (" << min_age_ << ") is less than the model's min_age (" << model_->min_age() << ")";
  if (max_age_ > model_->max_age())
    LOG_ERROR_P(PARAM_MAX_AGE) << ": max_age (" << max_age_ << ") is greater than the model's max_age (" << model_->max_age() << ")";

  for (auto year : years_) {
    LOG_FINEST() << "year : " << year;
    if((year < model_->start_year()) || (year > model_->final_year()))
      LOG_ERROR_P(PARAM_YEARS) << "Years can't be less than start_year (" << model_->start_year() << "), or greater than final_year (" << model_->final_year() << "). Please fix this.";
  }

  /**
   * Validate the number of obs provided matches age spread * category_labels * years
   * This is because we'll have 1 set of obs per category collection provided.
   * categories male+female male = 2 collections
   */
  unsigned obs_expected = age_spread_ + 1;

  /**
   * Build our error value map
   */
  vector<vector<string>>& error_values_data = error_values_table_->data();
  if (error_values_data.size() != years_.size()) {
    LOG_ERROR_P(PARAM_ERROR_VALUES) << " has " << error_values_data.size() << " rows defined, but we expected " << years_.size()
        << " to match the number of years provided";
  }

  for (vector<string>& error_values_data_line : error_values_data) {
    if (error_values_data_line.size() != 2 && error_values_data_line.size() != obs_expected) {
      LOG_FATAL_P(PARAM_ERROR_VALUES) << " has " << error_values_data_line.size() << " values defined, but we expected " << obs_expected
          << " to match the age speard  + 1 (for year)";
    }

    unsigned year = 0;
    if (!utilities::To<unsigned>(error_values_data_line[0], year))
      LOG_ERROR_P(PARAM_ERROR_VALUES) << " value " << error_values_data_line[0] << " could not be converted in to an unsigned integer. It should be the year for this line";
    if (std::find(years_.begin(), years_.end(), year) == years_.end())
      LOG_ERROR_P(PARAM_ERROR_VALUES) << " value " << year << " is not a valid year for this observation";
    for (unsigned i = 1; i < error_values_data_line.size(); ++i) {
      float value = 0;

      if (!utilities::To<float>(error_values_data_line[i], value))
        LOG_FATAL_P(PARAM_ERROR_VALUES) << " value (" << error_values_data_line[i] << ") could not be converted to a float";
      if (simulation_likelihood_label_ == PARAM_LOGNORMAL && value <= 0.0) {
        LOG_ERROR_P(PARAM_ERROR_VALUES) << ": error_value (" << value << ") cannot be equal to or less than 0.0";
      } else if (simulation_likelihood_label_ == PARAM_MULTINOMIAL && value < 0.0) {
        LOG_ERROR_P(PARAM_ERROR_VALUES) << ": error_value (" << value << ") cannot be less than 0.0";
      }

      error_values_by_year_[year].push_back(value);
    }
    if (error_values_by_year_[year].size() == 1) {
      error_values_by_year_[year].assign(obs_expected - 1, error_values_by_year_[year][0]);
    }
    LOG_FINEST() << "number of error values in year " << year << " = " << error_values_by_year_[year].size();
    if (error_values_by_year_[year].size() != obs_expected - 1)
      LOG_FATAL_P(PARAM_ERROR_VALUES) << "We counted " << error_values_by_year_[year].size() << " error values by year but expected " << obs_expected -1 << " based on the obs table";
  }

}

/**
 *
 */
void ProportionsAtAge::DoBuild() {
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

  // Build and validate layers
  layer_ = model_->managers().layer()->GetCategoricalLayer(layer_label_);
  if (!layer_)
    LOG_FATAL_P(PARAM_LAYER_OF_CELLS) << "could not find layer " << layer_label_ << " does it exist?, if it exists is of type categorical?";

  if( ageing_error_label_ != "") {
    ageing_error_ = model_->managers().ageing_error()->GetAgeingError(ageing_error_label_);
    if (!ageing_error_)
      LOG_ERROR_P(PARAM_AGEING_ERROR) << "(" << ageing_error_label_ << ") could not be found. Have you defined it?";
  }
  if (ageing_error_label_ == "") {
    LOG_WARNING() << "You are suppling a an age based observation with no ageing_misclassification";
  }

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

  LOG_FINE() << "finished building selectivities";
  world_ = model_->world_view();
  if (!world_)
    LOG_CODE_ERROR() << "!world_ could not create pointer to world via model, something is wrong";

  LOG_FINE() << "finished building pointer to world";

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
      LOG_FATAL_P(PARAM_CELLS) << "could not find the cell '" << cell << "' in the layer " << layer_label_ << " please make sure that you supply cell labels that are consistent with the layer.";
    unsigned n_bins = model_->age_spread();
    if (sexed_)
      n_bins *= 2;
    pre_age_freq_[cell].resize(n_bins);
    age_freq_[cell].resize(n_bins);
    final_age_freq_[cell].resize(n_bins);
  }

  LOG_FINE() << "finished checking cells";


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

/**
 *
 */
void ProportionsAtAge::PreExecute() {
  LOG_FINE();
  // allocate memory
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
              fill(pre_age_freq_[cell_label].begin(),pre_age_freq_[cell_label].end(),0);
              LOG_FINEST() << "about to convert " << cell->agents_.size() << " through the maturity process";
              //unsigned counter = 1;
              float probability;
              for (Agent& agent : cell->agents_) {
                if (agent.is_alive()) {
                  probability = selectivities_[agent.get_sex()]->GetResult(agent.get_length_bin_index());
                  if (rng.chance() <= probability) {
                    pre_age_freq_[cell_label][agent.get_age_index()] += agent.get_scalar();
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
              fill(pre_age_freq_[cell_label].begin(),pre_age_freq_[cell_label].end(),0);

              LOG_FINEST() << "about to convert " << cell->agents_.size() << " through the maturity process";
              //unsigned counter = 1;
              float probability;
              for (Agent& agent : cell->agents_) {
                if (agent.is_alive()) {
                  probability = selectivities_[agent.get_sex()]->GetResult(agent.get_age_index());
                  if (rng.chance() <= probability) {
                    pre_age_freq_[cell_label][agent.get_age_index()] += agent.get_scalar();
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
void ProportionsAtAge::Execute() {
  LOG_FINE();
  unsigned year = model_->current_year();
  if (find(years_.begin(), years_.end(), year) != years_.end()) {
    unsigned model_age_spread = model_->age_spread();
    utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
    if (!utilities::doublecompare::IsZero(time_step_proportion_)) {
      LOG_FINE() << "need to execute as need to account for some mortality";
      if (selectivity_length_based_) {
        LOG_FINE() << "Selectivity length based";
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            LOG_FINE() << "row = " << row << " col = " << col;
            string cell_label = layer_->get_value(row, col);
            if (find(cells_.begin(), cells_.end(), cell_label) == cells_.end())
              continue;
            WorldCell* cell = world_->get_base_square(row, col);
            if (cell->is_enabled()) {
              fill(age_freq_[cell_label].begin(),age_freq_[cell_label].end(),0);
              LOG_FINE() << "about to convert " << cell->agents_.size() << " through the prop at age process";
              //unsigned counter = 1;
              float probability;
              if (sexed_) {
                for (Agent& agent : cell->agents_) {
                  if (agent.is_alive()) {
                    probability = selectivities_[agent.get_sex()]->GetResult(agent.get_length_bin_index());
                    if (rng.chance() <= probability) {
                      age_freq_[cell_label][agent.get_age_index() + agent.get_sex() * model_age_spread] += agent.get_scalar();
                    }
                  }
                }
              } else {
                for (Agent& agent : cell->agents_) {
                  if (agent.is_alive()) {
                    probability = selectivities_[agent.get_sex()]->GetResult(agent.get_length_bin_index());
                    if (rng.chance() <= probability) {
                      age_freq_[cell_label][agent.get_age_index()] += agent.get_scalar();
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        LOG_FINE() << "Selectivity age based";
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            string cell_label = layer_->get_value(row, col);
            if (find(cells_.begin(), cells_.end(), cell_label) == cells_.end())
              continue;
            WorldCell* cell = world_->get_base_square(row, col);
            if (cell->is_enabled()) {
              LOG_FINE() << "row = " << row << " col = " << col;

              fill(age_freq_[cell_label].begin(),age_freq_[cell_label].end(),0);
              LOG_FINE() << age_freq_[cell_label].size();
              LOG_FINE() << "agents = " << cell->agents_.size() << " number of selectivities " << selectivities_.size();
              //unsigned counter = 1;
              if (sexed_) {
                LOG_FINE() << "observation is sexed specific";
                for (Agent& agent : cell->agents_) {
                  if (agent.is_alive()) {
                    if (rng.chance() <= selectivities_[agent.get_sex()]->GetResult(agent.get_age_index())) {
                      //LOG_FINEST() << "index = " << agent.get_age_index() << " scalar = " <<  agent.get_scalar() << agent.get_sex() << " sex = " << agent.get_sex() << " age = " << agent.get_age();
                      age_freq_[cell_label][agent.get_age_index() + agent.get_sex() * model_age_spread] += agent.get_scalar();
                    }
                  }
                }
              } else {
                for (Agent& agent : cell->agents_) {
                  if (agent.is_alive()) {
                    if (rng.chance() <= selectivities_[agent.get_sex()]->GetResult(agent.get_age_index())) {
                      //LOG_FINEST() << "index = " << agent.get_age_index() << " scalar = " <<  agent.get_scalar() << agent.get_sex() << " sex = " << agent.get_sex() << " age = " << agent.get_age();
                      age_freq_[cell_label][agent.get_age_index()] += agent.get_scalar();
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    LOG_FINE() << "Need to take into account ";
    vector<float> accumulated_age_frequency;
    unsigned age_offset = min_age_ - model_->min_age();
    for (auto& cell_label : cells_) {
      LOG_FINE() << "for cell " << cell_label << " we need to truncate, apply ageing error and save the value";
      if (utilities::doublecompare::IsZero(time_step_proportion_)) {
        // use pre ve
        final_age_freq_[cell_label] = pre_age_freq_[cell_label];
      } else if (utilities::doublecompare::IsOne(time_step_proportion_)) {
        final_age_freq_[cell_label] = age_freq_[cell_label];
      } else {
        float value;
        for (unsigned age_ndx = 0; age_ndx < age_freq_[cell_label].size(); ++age_ndx) {
          value = pre_age_freq_[cell_label][age_ndx] + ((age_freq_[cell_label][age_ndx] - pre_age_freq_[cell_label][age_ndx]) * time_step_proportion_);
          final_age_freq_[cell_label][age_ndx] = value;
        }
      }
      // Truncate and Apply ageing error if it exists
      if (ageing_error_) {
        if (not sexed_) {
          vector<float> temp(model_age_spread, 0.0);
          vector<vector<float>> &mis_matrix = ageing_error_->mis_matrix();
          for (unsigned i = 0; i < mis_matrix.size(); ++i) {
            for (unsigned j = 0; j < mis_matrix[i].size(); ++j) {
              temp[j] += final_age_freq_[cell_label][i] * mis_matrix[i][j];
            }
          }
          accumulated_age_frequency = temp;
        } else {
          vector<float> temp(model_age_spread * 2, 0.0);
          vector<vector<float>> &mis_matrix = ageing_error_->mis_matrix();
          for (unsigned sex = 0; sex <= 1; ++sex) {
            for (unsigned i = 0; i < mis_matrix.size(); ++i) {
              for (unsigned j = 0; j < mis_matrix[i].size(); ++j) {
                temp[j + sex * model_age_spread] += final_age_freq_[cell_label][i  + sex * model_age_spread] * mis_matrix[i][j];
              }
            }
          }
          accumulated_age_frequency = temp;
        }
      } else {
        accumulated_age_frequency = final_age_freq_[cell_label];
      }
      /*
       *  Now collapse the number_age into the expected_values for the observation
       */
      for (unsigned sex = 0; sex <= 1; ++sex) {
        if (not sexed_ and sex > 0)
          break;
        float plus_group = 0.0;
        for (unsigned k = 0; k < model_->age_spread(); ++k) {
          if (k >= age_offset && (k - age_offset + min_age_) < max_age_) {
            SaveComparison(k + model_->min_age(), 0, cell_label, accumulated_age_frequency[k + sex * model_age_spread], 0.0, error_values_by_year_[year][k - age_offset], year);
            LOG_FINE() << "age = " << k + model_->min_age() << " expected = " << accumulated_age_frequency[k + sex * model_age_spread];
          }
          // Deal with the plus group
          if (((k - age_offset + min_age_) >= max_age_) && plus_group_)
            plus_group += accumulated_age_frequency[k + sex * model_age_spread];
          else if (((k - age_offset + min_age_) == max_age_) && !plus_group_)
            plus_group = accumulated_age_frequency[k + sex * model_age_spread]; // no plus group and we are max age
        }
        SaveComparison(max_age_, sex, 0, cell_label, plus_group, 0.0, error_values_by_year_[year][max_age_ - min_age_], year);
        LOG_FINE() << "age = " << max_age_  << " expected = " << plus_group;
      }
    }
  }
}

/**
 *
 */
void ProportionsAtAge::Simulate() {
  /**
   * Simulate or generate results
   * During simulation mode we'll simulate results for this observation
   */
	LOG_FINE() << "Calculating score for observation = " << label_;
  // Convert to propotions before simulating for each year and cell sum = 1
  for (auto& iter : comparisons_) {  // year
    for (auto& second_iter : iter.second) {  // cell
      float total_expec = 0.0;
      for (auto& comparison : second_iter.second)
        total_expec += comparison.expected_;
      for (auto& comparison : second_iter.second)
        comparison.expected_ /= total_expec;
    }
  }
  likelihood_->SimulateObserved(comparisons_);
  // Simualte numbers at age, but we want proportion
  for (auto& iter : comparisons_) {
    for (auto& second_iter : iter.second) {  // cell
      float total = 0.0;
      for (auto& comparison : second_iter.second)
        total += comparison.simulated_;
      for (auto& comparison : second_iter.second)
        comparison.simulated_ /= total;
    }
  }}


} /* namespace observations */
} /* namespace niwa */

