/**
 * @file ProcessRemovalsByAge.cpp
 * @author  C Marsh
 * @version 1.0
 * @date 25/08/15
 * @section LICENSE
 *
 * Copyright NIWA Science 2016 - www.niwa.co.nz
 */

// Headers
#include "ProcessRemovalsByAge.h"

#include <algorithm>

#include "Model/Model.h"
#include "Processes/Manager.h"
#include "AgeingErrors/Manager.h"
#include "Utilities/Map.h"
#include "Utilities/Math.h"
#include "Utilities/To.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
// Namespaces
namespace niwa {
namespace observations {

/**
 * Default constructor
 */
ProcessRemovalsByAge::ProcessRemovalsByAge(Model* model) : Observation(model) {
  error_values_table_ = new parameters::Table(PARAM_ERROR_VALUES);

  parameters_.Bind<unsigned>(PARAM_MIN_AGE, &min_age_, "Minimum age", "");
  parameters_.Bind<unsigned>(PARAM_MAX_AGE, &max_age_, "Maximum age", "");
  parameters_.Bind<bool>(PARAM_PLUS_GROUP, &plus_group_, "Use age plus group", "", true);
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years for which there are observations", "");
  parameters_.Bind<string>(PARAM_AGEING_ERROR, &ageing_error_label_, "Label of ageing error to use", "", "");
  parameters_.BindTable(PARAM_ERROR_VALUES, error_values_table_, "Table of error values of the observed values (note the units depend on the likelihood)", "", false);
  parameters_.Bind<string>(PARAM_PROCESS_LABEL, &process_label_, "The label of the mortality process for the observation", "");

  allowed_likelihood_types_.push_back(PARAM_LOGNORMAL);
  allowed_likelihood_types_.push_back(PARAM_MULTINOMIAL);
}

/**
 * Destructor
 */
ProcessRemovalsByAge::~ProcessRemovalsByAge() {
  delete error_values_table_;
}

/**
 * Validate configuration file parameters
 */
void ProcessRemovalsByAge::DoValidate() {
  age_spread_ = (max_age_ - min_age_) + 1;
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
 * Build any runtime relationships we may have and ensure
 * the labels for other objects are valid.
 */
void ProcessRemovalsByAge::DoBuild() {
  // Create a pointer to misclassification matrix
    if( ageing_error_label_ != "") {
      ageing_error_ = model_->managers().ageing_error()->GetAgeingError(ageing_error_label_);
      if (!ageing_error_)
        LOG_ERROR_P(PARAM_AGEING_ERROR) << "(" << ageing_error_label_ << ") could not be found. Have you defined it?";
    }
    if (ageing_error_label_ == "") {
      LOG_WARNING() << "You are suppling a an age based observation with no ageing_misclassification";
    }

    mortality_process_ = model_->managers().process()->GetMortalityProcess(process_label_);
    if (!mortality_process_)
      LOG_FATAL_P(PARAM_PROCESS_LABEL) << "could not find the process " << process_label_ << ", please make sure it exists";

}

/**
 * We have to have a pre-execute and execute but all the information is stored on the process so I am just going to do
 * all the calculations in teh Simulate call
 */
void ProcessRemovalsByAge::PreExecute() {

}

/**
 *
 */
void ProcessRemovalsByAge::Execute() {

}

/**
 * This method is called at the end of a model iteration
 * to generate simulated data
 */
void ProcessRemovalsByAge::Simulate() {
  /**
   * Simulate or generate results
   * During simulation mode we'll simulate results for this observation
   */
  LOG_MEDIUM() << "Simulating data for observation = " << label_;
  map<unsigned, vector<unsigned>>& age_frequency = mortality_process_->get_removals_by_age();
  LOG_FINEST() << "number of years for this observation = " << age_frequency.size();
  unsigned age_offset = min_age_ - model_->min_age();
  // Apply ageing error
  float plus_group = 0.0;
  if (ageing_error_) {
    vector<vector<float>> &mis_matrix = ageing_error_->mis_matrix();
    for (auto year_age_freq : age_frequency) {
      vector<float> temp(model_->age_spread(), 0.0);
      if (find(years_.begin(), years_.end(), year_age_freq.first) != years_.end()) {
        for (unsigned i = 0; i < mis_matrix.size(); ++i) {
          for (unsigned j = 0; j < mis_matrix[i].size(); ++j) {
            temp[j] += (float)year_age_freq.second[i] * mis_matrix[i][j];
          }
        }
        /*
         *  Now collapse the number_age into the expected_values for the observation
         */
        for (unsigned k = 0; k < model_->age_spread(); ++k) {
          if (k >= age_offset && (k - age_offset + min_age_) < max_age_)
            SaveComparison(k + model_->min_age(), 0, temp[k], 0.0, error_values_by_year_[year_age_freq.first][k - age_offset], year_age_freq.first);
          // Deal with the plus group
          if (((k - age_offset + min_age_) >= max_age_) && plus_group_)
            plus_group += temp[k];
          else if (((k - age_offset + min_age_) == max_age_) && !plus_group_)
            plus_group = temp[k]; // no plus group and we are max age

        }
        SaveComparison(max_age_, 0, plus_group, 0.0, error_values_by_year_[year_age_freq.first][max_age_ - min_age_], year_age_freq.first);
      }
    }
  } else {
    for (auto year_age_freq : age_frequency) {
      LOG_FINEST() << "year = " << year_age_freq.first << " age offset = " << age_offset;
      if (find(years_.begin(), years_.end(), year_age_freq.first) != years_.end()) {
        for (unsigned k = 0; k < model_->age_spread(); ++k) {
          if (k >= age_offset && (k - age_offset + min_age_) < max_age_)
            SaveComparison(k + model_->min_age(), 0, (float)year_age_freq.second[k], 0.0, error_values_by_year_[year_age_freq.first][k - age_offset], year_age_freq.first);
          if (((k - age_offset + min_age_) >= max_age_) && plus_group_)
            plus_group += (float)year_age_freq.second[k];
          else if (((k - age_offset + min_age_) == max_age_) && !plus_group_)
            plus_group = (float)year_age_freq.second[k]; // no plus group and we are max age
        }
        SaveComparison(max_age_, 0, plus_group, 0.0, error_values_by_year_[year_age_freq.first][max_age_ - min_age_], year_age_freq.first);
      }
    }
  }

  // Convert to propotions before simulating
  for (auto& iter : comparisons_) {
    float total_expec = 0.0;
    for (auto& comparison : iter.second)
      total_expec += comparison.expected_;
    for (auto& comparison : iter.second)
      comparison.expected_ /= total_expec;
  }
  likelihood_->SimulateObserved(comparisons_);
  // Simualte numbers at age, but we want proportion
  for (auto& iter : comparisons_) {
    float total = 0.0;
    for (auto& comparison : iter.second)
      total += comparison.simulated_;
    for (auto& comparison : iter.second)
      comparison.simulated_ /= total;
  }
}

} /* namespace observations */
} /* namespace niwa */
