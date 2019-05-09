/**
 * @file ProcessRemovalsByLength.cpp
 * @author  C Marsh
 * @version 1.0
 * @date 25/07/18
 * @section LICENSE
 *
 */

// Headers
#include "ProcessRemovalsByLength.h"

#include <algorithm>

#include "Model/Model.h"
#include "Processes/Manager.h"
#include "Likelihoods/Manager.h"
#include "Layers/Manager.h"
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
ProcessRemovalsByLength::ProcessRemovalsByLength(Model* model) : Observation(model) {
  error_values_table_ = new parameters::Table(PARAM_ERROR_VALUES);
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years for which there are observations", "");
  parameters_.Bind<string>(PARAM_PROCESS_LABEL, &process_label_, "Label of of removal process", "", "");
  parameters_.BindTable(PARAM_ERROR_VALUES, error_values_table_, "Table of error values of the observed values (note the units depend on the likelihood)", "", false);
  parameters_.Bind<string>(PARAM_LAYER_OF_CELLS, &layer_label_, "The layer that indicates what area to summarise observations over.", "");
  parameters_.Bind<string>(PARAM_CELLS, &cells_, "The cells we want to generate observations for from the layer of cells supplied", "");
  parameters_.Bind<string>(PARAM_SIMULATION_LIKELIHOOD, &simulation_likelihood_label_, "Simulation likelihood to use", "");

  allowed_likelihood_types_.push_back(PARAM_LOGNORMAL);
  allowed_likelihood_types_.push_back(PARAM_MULTINOMIAL);
}

/**
 * Destructor
 */
ProcessRemovalsByLength::~ProcessRemovalsByLength() {
  delete error_values_table_;
}

/**
 * Validate configuration file parameters
 */
void ProcessRemovalsByLength::DoValidate() {
  for (auto year : years_) {
    LOG_FINEST() << "year : " << year;
  	if((year < model_->start_year()) || (year > model_->final_year()))
  		LOG_ERROR_P(PARAM_YEARS) << "Years can't be less than start_year (" << model_->start_year() << "), or greater than final_year (" << model_->final_year() << "). Please fix this.";
  }
  /**
   * Build our error value map
   */
  vector<vector<string>>& error_values_data = error_values_table_->data();
  if (error_values_data.size() != years_.size()) {
    LOG_ERROR_P(PARAM_ERROR_VALUES) << " has " << error_values_data.size() << " rows defined, but we expected " << years_.size()
        << " to match the number of years provided";
  }

  for (vector<string>& error_values_data_line : error_values_data) {
    if (error_values_data_line.size() != 2 && error_values_data_line.size() != model_->length_bin_mid_points().size()) {
      LOG_FATAL_P(PARAM_ERROR_VALUES) << " has " << error_values_data_line.size() << " values defined, but we expected " << model_->length_bin_mid_points().size()
          << " to match the model length bins";
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
      error_values_by_year_[year].assign(model_->length_bin_mid_points().size(), error_values_by_year_[year][0]);
    }
    LOG_FINEST() << "number of error values in year " << year << " = " << error_values_by_year_[year].size();
    if (error_values_by_year_[year].size() != model_->length_bin_mid_points().size())
      LOG_FATAL_P(PARAM_ERROR_VALUES) << "We counted " << error_values_by_year_[year].size() << " error values by year but expected " << model_->length_bin_mid_points().size() << " based on the obs table";
  }

}

/**
 * Build any runtime relationships we may have and ensure
 * the labels for other objects are valid.
 */
void ProcessRemovalsByLength::DoBuild() {


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

    mortality_process_ = model_->managers().process()->GetMortalityProcess(process_label_);
    if (!mortality_process_)
      LOG_FATAL_P(PARAM_PROCESS_LABEL) << "could not find the process " << process_label_ << ", please make sure it exists";

    // Build and validate layers
    layer_ = model_->managers().layer()->GetCategoricalLayer(layer_label_);
    if (!layer_)
      LOG_FATAL_P(PARAM_LAYER_OF_CELLS) << "could not find layer " << layer_label_ << " does it exist?, if it exists is of type categorical?";

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
}

/**
 * We have to have a pre-execute and execute but all the information is stored on the process so I am just going to do
 * all the calculations in teh Simulate call
 */
void ProcessRemovalsByLength::PreExecute() {

}

/**
 *
 */
void ProcessRemovalsByLength::Execute() {

}

/**
 * This method is called at the end of a model iteration
 * to generate simulated data
 */
void ProcessRemovalsByLength::Simulate() {
  LOG_MEDIUM() << "Simulating data for observation = " << label_;
  /**
   * Simulate or generate results
   * During simulation mode we'll simulate results for this observation
   */
  vector<processes::composition_data>& length_frequency = mortality_process_->get_removals_by_length();
  LOG_FINEST() << "number of years for this observation = " << length_frequency.size();
  vector<float> length_mid_points = model_->length_bin_mid_points();
  // iterate over all the years that we want
  for (unsigned year : years_) {
    for (string cell : cells_) {
      bool cell_found = false;
      vector<float> accumulated_length_frequency(model_->length_bin_mid_points().size(), 0.0);
      for (auto length_comp_data : length_frequency) {
        if ((length_comp_data.year_ == year) && (layer_->get_value(length_comp_data.row_,length_comp_data.col_) == cell)) {
          // Lets accumulate the information for this cell and year
          for(unsigned i = 0; i < length_comp_data.frequency_.size(); ++i) {
            accumulated_length_frequency[i] += length_comp_data.frequency_[i];
          }
          cell_found = true;
        }
      }
      if (not cell_found)
        continue; // to next cell
      /*
       *  Now collapse the number_age into the expected_values for the observation
       */
      for (unsigned k = 0; k < model_->length_bin_mid_points().size(); ++k) {
        SaveComparison(0, length_mid_points[k], cell, accumulated_length_frequency[k], 0.0, error_values_by_year_[year][k], year);
      }
    }
  }
  // Convert to propotions before simulating for each year and cell sum = 1
/*  for (auto& iter : comparisons_) {  // year
    for (auto& second_iter : iter.second) {  // cell
      float total_expec = 0.0;
      for (auto& comparison : second_iter.second)
        total_expec += comparison.expected_;
      for (auto& comparison : second_iter.second)
        comparison.expected_ /= total_expec;
    }
  }*/
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
  }
}

} /* namespace observations */
} /* namespace niwa */
