/**
 * @file MortalityScaledAgeFrequency.cpp
 * @author  C.Marsh github.com/Craig44
 * @date 4/11/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "MortalityScaledAgeFrequency.h"


#include "Processes/Manager.h"
#include "Layers/Manager.h"
#include "AgeingErrors/Manager.h"

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
MortalityScaledAgeFrequency::MortalityScaledAgeFrequency(Model* model) : Observation(model) {
  sample_table_ = new parameters::Table(PARAM_SAMPLES);
  parameters_.Bind<string>(PARAM_TIME_STEP, &time_step_label_, "The label of time-step that the observation occurs in", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "The years of the observed values", "");
  parameters_.Bind<string>(PARAM_AGEING_ERROR, &ageing_error_label_, "Label of ageing error to use", "", PARAM_NONE);
  parameters_.Bind<string>(PARAM_PROCESS_LABEL, &process_label_, "Label of of removal process", "", "");

  parameters_.Bind<string>(PARAM_STRATUM_WEIGHT_METHOD, &stratum_weight_method_, "Method to weight stratum estimates by", "", PARAM_BIOMASS)->set_allowed_values({PARAM_BIOMASS,PARAM_AREA});
  parameters_.Bind<string>(PARAM_LAYER_OF_CELLS, &layer_label_, "The layer that indicates what the stratum boundaries are.", "");
  parameters_.Bind<string>(PARAM_CELLS, &cells_, "The cells which represent individual stratum to be included in the analysis, default is all cells are used from the layer", "", true);

  allowed_likelihood_types_.push_back(PARAM_LOGNORMAL);
  allowed_likelihood_types_.push_back(PARAM_MULTINOMIAL);
  allowed_likelihood_types_.push_back(PARAM_DIRICHLET);
  allowed_likelihood_types_.push_back(PARAM_LOGISTIC_NORMAL);
}
/**
 * Destructor
 */
MortalityScaledAgeFrequency::~MortalityScaledAgeFrequency() {
  delete sample_table_;
}
/**
 *
 */
void MortalityScaledAgeFrequency::DoValidate() {
  LOG_TRACE();
  for (auto year : years_) {
    LOG_FINEST() << "year : " << year;
    if ((year < model_->start_year()) || (year > model_->final_year()))
      LOG_ERROR_P(PARAM_YEARS) << "Years can't be less than start_year (" << model_->start_year() << "), or greater than final_year (" << model_->final_year()
          << "). Please fix this.";
  }

}

/**
 *
 */
void MortalityScaledAgeFrequency::DoBuild() {
  LOG_TRACE();
  // Create a pointer to misclassification matrix
  if (ageing_error_label_ != PARAM_NONE) {
    ageing_error_ = model_->managers().ageing_error()->GetAgeingError(ageing_error_label_);
    if (!ageing_error_)
      LOG_ERROR_P(PARAM_AGEING_ERROR) << "(" << ageing_error_label_ << ") could not be found. Have you defined it?";
  }
  if (ageing_error_label_ == PARAM_NONE) {
    LOG_WARNING() << "You are suppling a an age based observation with no ageing_misclassification";
  }

  mortality_process_ = model_->managers().process()->GetMortalityProcess(process_label_);
  if (!mortality_process_)
    LOG_FATAL_P(PARAM_PROCESS_LABEL)<< "could not find the process " << process_label_ << ", please make sure it exists";

    // Build and validate layers
  layer_ = model_->managers().layer()->GetCategoricalLayer(layer_label_);
  if (!layer_)
    LOG_FATAL_P(PARAM_LAYER_OF_CELLS)<< "could not find layer " << layer_label_ << " does it exist?, if it exists is of type categorical?";

  if (parameters_.Get(PARAM_CELLS)->has_been_defined()) {
    // Check all the cells supplied are in the layer
    for (auto cell : cells_) {
      LOG_FINEST() << "checking cell " << cell << " exists in layer";
      bool cell_found = false;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          LOG_FINEST() << "checking row = " << row << " col = " << col << " value = " << layer_->get_value(row, col);
          if (layer_->get_value(row, col) == cell) {
            cell_found = true;
            stratum_rows_[cell].push_back(row);
            stratum_cols_[cell].push_back(col);
            break;
          }
        }
      }
      if (not cell_found)
        LOG_ERROR_P(PARAM_CELLS) << "could not find the cell '" << cell << "' in the layer " << layer_label_
            << " please make sure that you supply cell labels that are consistent with the layer.";
    }
  } else {
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        string temp_cell = layer_->get_value(row, col);
        stratum_rows_[temp_cell].push_back(row);
        stratum_cols_[temp_cell].push_back(col);
        if (find(cells_.begin(), cells_.end(), temp_cell) == cells_.end())
          cells_.push_back(temp_cell);
      }
    }
  }
  if (not mortality_process_->check_years(years_)) {
    LOG_ERROR_P(PARAM_YEARS)
        << "there was a year that the mortality process doesn't not execute in, can you please check that the years you have supplied for this observation are years that the mortality process occurs in cheers.";
  }

  /*
   * Build sample table
  */
  LOG_FINE() << "Build sample table";
  vector<vector<string>>& sample_data = sample_table_->data();
  if (sample_data.size() != years_.size()) {
    LOG_ERROR_P(PARAM_SAMPLES) << " has " << sample_data.size() << " rows defined, but we expected " << years_.size()
        << " to match the number of years provided";
  }
  unsigned row_counter = 0;
  vector<string>  stratum_order_from_sample_table = sample_data[0];
  for (vector<string>& sample_values_data_line : sample_data) {
    if (sample_values_data_line.size() != cells_.size()) {
      LOG_FATAL_P(PARAM_SAMPLES) << " has " << sample_values_data_line.size() << " values defined, but we expected " << cells_.size()
          << " to match a number for each stratum";
    }
    if (row_counter == 0) {
      //load stratum look ups
      continue;
    }
    unsigned year = 0;
    if (!utilities::To<unsigned>(sample_values_data_line[0], year))
      LOG_ERROR_P(PARAM_SAMPLES) << " value " << sample_values_data_line[0] << " could not be converted in to an unsigned integer. It should be the year for this line";
    if (std::find(years_.begin(), years_.end(), year) == years_.end())
      LOG_ERROR_P(PARAM_SAMPLES) << " value " << year << " is not a valid year for this observation";
    for (unsigned i = 1; i < sample_values_data_line.size(); ++i) {
      unsigned value = 0;

      if (!utilities::To<unsigned>(sample_values_data_line[i], value))
        LOG_FATAL_P(PARAM_SAMPLES) << " value (" << sample_values_data_line[i] << ") could not be converted to a unsigned";
      if (value <= 0) {
        LOG_ERROR_P(PARAM_SAMPLES) << "at row = " << row_counter << " and column " << i << " the vaue given = " << value << " this needs to be a positive integer.";
      }

      samples_by_year_and_stratum_[year][stratum_order_from_sample_table[i - 1]].push_back(value);
    }
    ++row_counter;
  }

}

/**
 *
 */
void MortalityScaledAgeFrequency::PreExecute() {

}

/**
 *
 */
void MortalityScaledAgeFrequency::Execute() {

}

/**
 *
 */
void MortalityScaledAgeFrequency::Simulate() {

  LOG_MEDIUM() << "Simulating data for observation = " << label_;
  /**
   * Simulate or generate results
   * During simulation mode we'll simulate results for this observation
   */

}


} /* namespace observations */
} /* namespace niwa */

