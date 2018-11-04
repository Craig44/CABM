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
  parameters_.BindTable(PARAM_SAMPLES, sample_table_, "Table of sample sizes used to generate age length key for each stratum and each year.", "", false);
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "The years of the observed values", "");
  parameters_.Bind<string>(PARAM_AGEING_ERROR, &ageing_error_label_, "Label of ageing error to use", "", PARAM_NONE);
  parameters_.Bind<string>(PARAM_PROCESS_LABEL, &process_label_, "Label of of removal process", "", "");

  parameters_.Bind<string>(PARAM_STRATUM_WEIGHT_METHOD, &stratum_weight_method_, "Method to weight stratum estimates by", "", PARAM_BIOMASS)->set_allowed_values({PARAM_BIOMASS,PARAM_AREA});
  parameters_.Bind<string>(PARAM_LAYER_OF_STRATUM_DEFINITIONS, &layer_label_, "The layer that indicates what the stratum boundaries are.", "");
  parameters_.Bind<string>(PARAM_STRATUMS_TO_INCLUDE, &cells_, "The cells which represent individual stratum to be included in the analysis, default is all cells are used from the layer", "", true);

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
    LOG_FATAL_P(PARAM_LAYER_OF_STRATUM_DEFINITIONS)<< "could not find layer " << layer_label_ << " does it exist?, if it exists is of type categorical?";

  LOG_FINE() << "Check stratum are consistent";
  if (parameters_.Get(PARAM_STRATUMS_TO_INCLUDE)->has_been_defined()) {
    // Check all the cells supplied are in the layer
    for (auto cell : cells_) {
      LOG_FINE() << "checking cell " << cell << " exists in layer";
      bool cell_found = false;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          LOG_FINE() << "checking row = " << row << " col = " << col;
          LOG_FINE() << "value = " << layer_->get_value(row, col);
          if (layer_->get_value(row, col) == cell) {
            cell_found = true;
            stratum_rows_[cell].push_back(row);
            stratum_cols_[cell].push_back(col);
          }
        }
      }
      if (not cell_found)
        LOG_ERROR_P(PARAM_STRATUMS_TO_INCLUDE) << "could not find the cell '" << cell << "' in the layer " << layer_label_
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
  LOG_FINE() << "Check mortality process is consistent with observation";
  if (not mortality_process_->check_years(years_)) {
    LOG_ERROR_P(PARAM_YEARS)
        << "there was a year that the mortality process doesn't not execute in, can you please check that the years you have supplied for this observation are years that the mortality process occurs in cheers.";
  }

  for (auto& row_map : stratum_rows_) {
    LOG_FINE() << "rows in cell " << row_map.first;
    for (auto val : row_map.second)
      LOG_FINE() << val;
  }
  for (auto& col_map : stratum_cols_) {
    LOG_FINE() << "cols in cell " << col_map.first;
    for (auto val : col_map.second)
      LOG_FINE() << val;
  }
  /*
   * Build sample table
  */
  LOG_FINE() << "Build sample table";
  vector<vector<string>>& sample_data = sample_table_->data();
  if ((sample_data.size() - 1) != years_.size()) {
    LOG_ERROR_P(PARAM_SAMPLES) << " has " << (sample_data.size() - 1) << " rows defined, but we expected " << years_.size()
        << " to match the number of years provided";
  }
  vector<string>  stratum_order_from_sample_table = sample_data[0];

  if (stratum_order_from_sample_table[0] != PARAM_YEAR)
    LOG_FATAL_P(PARAM_SAMPLES) << "Expected the first column to have the following header '" << PARAM_YEAR << "'";

  for (unsigned col_ndx = 1; col_ndx < stratum_order_from_sample_table.size(); ++col_ndx) {
    if (find(cells_.begin(),cells_.end(),stratum_order_from_sample_table[col_ndx]) == cells_.end())
      LOG_FATAL_P(PARAM_SAMPLES) << "Could not find the stratum '" << stratum_order_from_sample_table[col_ndx] << "' (colum header '" << col_ndx + 1 << "') in the parameter " << PARAM_STRATUMS_TO_INCLUDE << " can you please check that the column headers are consistent with this parameter, chairs";

  }



  for (unsigned row_counter = 1; row_counter < sample_data.size();++row_counter) {
    if (sample_data[row_counter].size() != (cells_.size() + 1)) {
      LOG_FATAL_P(PARAM_SAMPLES) << " has " << sample_data[row_counter].size() << " values defined, but we expected " << cells_.size() + 1
          << " to match a number for each stratum";
    }
    unsigned year = 0;
    if (!utilities::To<unsigned>(sample_data[row_counter][0], year))
      LOG_ERROR_P(PARAM_SAMPLES) << " value " << sample_data[row_counter][0] << " could not be converted in to an unsigned integer. It should be the year for this line";
    if (std::find(years_.begin(), years_.end(), year) == years_.end())
      LOG_ERROR_P(PARAM_SAMPLES) << " value " << year << " is not a valid year for this observation";

    LOG_FINE() << "Year = " << year;
    for (unsigned i = 1; i < sample_data[row_counter].size(); ++i) {
      unsigned value = 0;

      if (!utilities::To<unsigned>(sample_data[row_counter][i], value))
        LOG_FATAL_P(PARAM_SAMPLES) << " value (" << sample_data[row_counter][i] << ") could not be converted to a unsigned";
      if (value <= 0) {
        LOG_ERROR_P(PARAM_SAMPLES) << "at row = " << row_counter << " and column " << i << " the vaue given = " << value << " this needs to be a positive integer.";
      }
      LOG_FINE() << "stratum = " << stratum_order_from_sample_table[i] << " samples = " << value;

      samples_by_year_and_stratum_[year][stratum_order_from_sample_table[i]] = value;
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
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  vector<vector<float>> age_length_key;
  age_length_key.resize(model_->age_spread());
  for (unsigned i = 0; i < model_->age_spread(); ++i)
    age_length_key[i].resize(model_->length_bins().size(),0.0);

  vector<processes::census_data>& census_data = mortality_process_->get_census_data();
  vector<unsigned> census_stratum_ndx;

  bool apply_ageing_error = false;
  if (!ageing_error_) {
    apply_ageing_error = true;
  }

  for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    LOG_FINE() << "About to sort our info for year " << years_[year_ndx];
    for (unsigned stratum_ndx = 0; stratum_ndx < cells_.size(); ++stratum_ndx) {
      LOG_FINE() << "About to sort our info for stratum " << cells_[stratum_ndx];
      census_stratum_ndx.clear();
      for (unsigned i = 0; i < model_->age_spread(); ++i)
        age_length_key[i].resize(model_->length_bins().size(),0.0);

      unsigned census_ndx = 0;
      for (processes::census_data& census : census_data) {
        if ((census.year_ == years_[year_ndx]) && (find(stratum_rows_[cells_[stratum_ndx]].begin(),stratum_rows_[cells_[stratum_ndx]].end(), census.row_) != stratum_rows_[cells_[stratum_ndx]].end()) && (find(stratum_cols_[cells_[stratum_ndx]].begin(),stratum_cols_[cells_[stratum_ndx]].end(), census.col_) != stratum_cols_[cells_[stratum_ndx]].end()))
          census_stratum_ndx.push_back(census_ndx);
        ++census_ndx;
      }
      unsigned samples_to_take = samples_by_year_and_stratum_[years_[year_ndx]][cells_[stratum_ndx]];
      // Some containers to make sure that we are doing sampling WITHOUT replacement
      // This will slow it down, but worth it I believe
      vector<unsigned> census_sampled;
      vector<unsigned> agents_sampled;
      census_ndx = 0;
      unsigned agent_ndx;
      unsigned age_ndx = 0;
      unsigned length_ndx = 0;
      /*
       * Generate an age-length-key for this stratum
      */
      for (unsigned sample_attempt = 0; sample_attempt < samples_to_take;) {
        // Randomly select a cell that a stratum belons to
        census_ndx = rng.chance() * census_stratum_ndx.size();
        processes::census_data& this_census = census_data[census_ndx];
        // Randomly select an agent in that cell
        agent_ndx = this_census.age_ndx_.size() * rng.chance();
        if ((find(census_sampled.begin(), census_sampled.end(),census_ndx) != census_sampled.end()) && (find(agents_sampled.begin(), agents_sampled.end(),agent_ndx) != agents_sampled.end())) {
          // If we have sampled this fish in this cell try again with out counting this as a sample
          continue;
        }
        // else lets remember that we have sampled this fish
        census_sampled.push_back(census_ndx);
        agents_sampled.push_back(agent_ndx);

        // Are we applying ageing error which will be a multinomial process
        age_ndx = this_census.age_ndx_[agent_ndx];
        length_ndx = this_census.length_ndx_[agent_ndx];
        if (apply_ageing_error) {
          vector<vector<float>> &mis_matrix = ageing_error_->mis_matrix();
          vector<float> prob_mis_classification = mis_matrix[age_ndx];
          float temp_prob = 0.0;
          for (unsigned mis_ndx = 0; mis_ndx < prob_mis_classification.size(); ++mis_ndx) {
            temp_prob += prob_mis_classification[mis_ndx];
            if (rng.chance() <= temp_prob)
              age_ndx = mis_ndx;
          }
        }
        age_length_key[age_ndx][length_ndx]++;
        ++sample_attempt;
      }
      /*
       * Calculate the age-frequency by passing the length frequency through the Age-Length key.
       * if Bootstrap=true, do a sample with replacement procedure to calculate C.V for each age bin
      */
    }
  }
}

} /* namespace observations */
} /* namespace niwa */

