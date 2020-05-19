/**
 * @file MortalityEventComposition.cpp
 * @author  C.Marsh github.com/Craig44
 * @date 4/11/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "MortalityEventComposition.h"


#include "Processes/Manager.h"
#include "Layers/Manager.h"
#include "AgeingErrors/Manager.h"
#include "Likelihoods/Manager.h"

#include "World/WorldView.h"
#include "World/WorldCell.h"

#include <omp.h>

#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/To.h"
#include "Utilities/Map.h"
#include "Utilities/Math.h"
#include "Utilities/DoubleCompare.h"

// namespaces
namespace niwa {
namespace observations {

namespace utils = niwa::utilities;
namespace math = niwa::utilities::math;


/**
 * Default constructor
 */
MortalityEventComposition::MortalityEventComposition(Model* model) : Observation(model) {
  error_values_table_ = new parameters::Table(PARAM_ERROR_VALUES);
  parameters_.BindTable(PARAM_ERROR_VALUES, error_values_table_, "Table of error values of the observed values (note the units depend on the likelihood)", "", false);
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "The years of the observed values", "");
  parameters_.Bind<string>(PARAM_AGEING_ERROR, &ageing_error_label_, "Label of ageing error to use", "", PARAM_NONE);
  parameters_.Bind<string>(PARAM_PROCESS_LABEL, &process_label_, "Label of of removal process", "", "");
  parameters_.Bind<string>(PARAM_FISHERY_LABEL, &fishery_label_, "Label of of removal process", "");
  parameters_.Bind<string>(PARAM_STRATUM_WEIGHT_METHOD, &stratum_weight_method_, "Method to weight stratum estimates by", "", PARAM_BIOMASS)->set_allowed_values({PARAM_BIOMASS, PARAM_AREA, PARAM_NONE});
  parameters_.Bind<bool>(PARAM_SEXED, &sexed_flag_, "You can ask to 'ignore' sex (only option for unsexed model), or generate composition for a particular sex, either 'male' or 'female", "", false);
  parameters_.Bind<string>(PARAM_COMPOSITION_TYPE, &comp_type_, "Is the composition Age or Length", "", PARAM_AGE)->set_allowed_values({PARAM_AGE, PARAM_LENGTH});
  parameters_.Bind<bool>(PARAM_NORMALISE, &are_obs_props_, "Are the compositions normalised to sum to one", "", true);
  parameters_.Bind<string>(PARAM_SIMULATION_LIKELIHOOD, &simulation_likelihood_label_, "Simulation likelihood to use", "");
  parameters_.Bind<unsigned>(PARAM_MIN_AGE, &min_age_, "Minimum age", "", model_->min_age());
  parameters_.Bind<unsigned>(PARAM_MAX_AGE, &max_age_, "Maximum age", "", model_->max_age());
  parameters_.Bind<bool>(PARAM_PLUS_GROUP, &plus_group_, "max age is a plus group", "", true);

  parameters_.Bind<string>(PARAM_LAYER_OF_STRATUM_DEFINITIONS, &layer_label_, "The layer that indicates what the stratum boundaries are.", "");
  parameters_.Bind<string>(PARAM_STRATUMS_TO_INCLUDE, &cells_, "The cells which represent individual stratum to be included in the analysis, default is all cells are used from the layer", "", true);
  allowed_likelihood_types_.push_back(PARAM_LOGNORMAL);
  allowed_likelihood_types_.push_back(PARAM_MULTINOMIAL);
  allowed_likelihood_types_.push_back(PARAM_DIRICHLET);
  allowed_likelihood_types_.push_back(PARAM_LOGISTIC_NORMAL);
  allowed_likelihood_types_.push_back(PARAM_PSEUDO);
}
/**
 * Destructor
 */
MortalityEventComposition::~MortalityEventComposition() {
  delete error_values_table_;
}
/**
 *
 */
void MortalityEventComposition::DoValidate() {
  LOG_TRACE();
  for (auto year : years_) {
    LOG_FINE() << "year : " << year;
    if ((year < model_->start_year()) || (year > model_->final_year()))
      LOG_ERROR_P(PARAM_YEARS) << "Years can't be less than start_year (" << model_->start_year() << "), or greater than final_year (" << model_->final_year()
          << "). Please fix this.";
  }

  if(!model_->get_sexed()) {
    if (sexed_flag_) {
      LOG_WARNING() << "you asked for a sexed observation but the model isn't sexed so I am ignoring this and giving you unsexed results.";
      sexed_flag_ = false;
    }
  }

  if (obs_type_ == PARAM_NUMBERS)
    are_obs_props_ = false;
  LOG_MEDIUM() << "comp type = " << obs_type_ << " bool = " << are_obs_props_ << " 0 = numbers, 1 = props";

  if (comp_type_ == PARAM_AGE) {
    n_bins_ = model_->age_spread();
    n_unsexed_bins_ = n_bins_;
    is_age_ = true;
  } else {
    n_bins_ = model_->number_of_length_bins();
    n_unsexed_bins_ = n_bins_;

    is_age_ = false;
  }
  if (sexed_flag_)
    n_bins_ *= 2;

  stratum_comp_.resize(n_bins_,0.0);
}

/**
 *
 */
void MortalityEventComposition::DoBuild() {
  LOG_TRACE();
  // Create a pointer to misclassification matrix
  if (is_age_) {
    if (ageing_error_label_ != PARAM_NONE) {
      ageing_error_ = model_->managers().ageing_error()->GetAgeingError(ageing_error_label_);
      if (!ageing_error_)
        LOG_ERROR_P(PARAM_AGEING_ERROR) << "(" << ageing_error_label_ << ") could not be found. Have you defined it?";
       ageing_mis_matrix_ = ageing_error_->mis_matrix();
    }
    if (ageing_error_label_ == PARAM_NONE) {
      LOG_WARNING() << "You are suppling an age based observation with no ageing misclassification error";
    }
  }

  world_ = model_->world_view();
  if (!world_)
    LOG_CODE_ERROR() << "!world_ could not create pointer to world viw model, something is wrong";

  mortality_process_ = model_->managers().process()->GetMortalityEventBiomassProcess(process_label_);
  if (!mortality_process_)
    LOG_FATAL_P(PARAM_PROCESS_LABEL)<< "could not find the process " << process_label_ << ", please make sure it exists and is of type " << PARAM_MORTALITY_EVENT_BIOMASS;

    // Build and validate layers
  layer_ = model_->managers().layer()->GetCategoricalLayer(layer_label_);
  if (!layer_)
    LOG_FATAL_P(PARAM_LAYER_OF_STRATUM_DEFINITIONS)<< "could not find layer " << layer_label_ << " does it exist?, if it exists is of type categorical?";


  // Build Likelihood
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

  LOG_FINE() << "Check stratum are consistent";
  if (parameters_.Get(PARAM_STRATUMS_TO_INCLUDE)->has_been_defined()) {
    // Check all the cells supplied are in the layer
    for (auto cell : cells_) {
      LOG_FINE() << "checking cell " << cell << " exists in layer";
      stratum_area_[cell] = 0.0;
      bool cell_found = false;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          LOG_FINE() << "checking row = " << row << " col = " << col;
          LOG_FINE() << "value = " << layer_->get_value(row, col);
          if (layer_->get_value(row, col) == cell) {
            cell_found = true;
            stratum_rows_[cell].push_back(row);
            stratum_cols_[cell].push_back(col);

            if (stratum_weight_method_ == PARAM_AREA) {
              WorldCell* world_cell = world_->get_base_square(row, col);
              if (world_cell->is_enabled())
                stratum_area_[cell] = world_cell->get_area();
            }
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
        if (stratum_weight_method_ == PARAM_AREA) {
          WorldCell* world_cell = world_->get_base_square(row, col);
          if (world_cell->is_enabled())
            stratum_area_[temp_cell] = world_cell->get_area();
        }
      }
    }
  }

  LOG_FINE() << "Check mortality process is consistent with observation";
  if (not mortality_process_->check_years(years_)) {
    LOG_ERROR_P(PARAM_YEARS)
        << "there was a year that the mortality process doesn't not execute in, can you please check that the years you have supplied for this observation are years that the mortality process occurs in cheers.";
  }

  if (not mortality_process_->check_fishery_exists(fishery_label_)) {
    LOG_FATAL_P(PARAM_FISHERY_LABEL)
        << "could not find the fishery label " << fishery_label_ << " in the mortality process " << process_label_ << " please check it exists and or is spelt correctly";
  }

  fishery_years_ = mortality_process_->get_fishery_years();

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

  /**
   * Build our error value map with dimensions year x stratum
   */
  // Validate stratum errorvalues table
  vector<vector<string>>& error_data = error_values_table_->data();
  if ((error_data.size() - 1) != years_.size()) {
    LOG_ERROR_P(PARAM_ERROR_VALUES) << " has " << (error_data.size() - 1) << " rows defined, but we expected " << years_.size()
        << " to match the number of years provided";
  }
  vector<string>  stratum_order_from_error_table = error_data[0];

  if (stratum_order_from_error_table[0] != PARAM_YEAR)
    LOG_FATAL_P(PARAM_ERROR_VALUES) << "Expected the first column to have the following header '" << PARAM_YEAR << "'";

  for (unsigned col_ndx = 1; col_ndx < stratum_order_from_error_table.size(); ++col_ndx) {
    if (find(cells_.begin(),cells_.end(),stratum_order_from_error_table[col_ndx]) == cells_.end())
      LOG_FATAL_P(PARAM_ERROR_VALUES) << "Could not find the stratum '" << stratum_order_from_error_table[col_ndx] << "' (colum header '" << col_ndx + 1 << "') in the parameter " << PARAM_STRATUMS_TO_INCLUDE << " can you please check that the column headers are consistent with this parameter, chairs";
  }

  LOG_MEDIUM() << "rows in error table = " << error_data.size();
  for (unsigned row_counter = 1; row_counter < error_data.size(); ++row_counter) {
    if (error_data[row_counter].size() != (cells_.size() + 1)) {
      LOG_FATAL_P(PARAM_ERROR_VALUES) << " has " << error_data[row_counter].size() << " values defined, but we expected " << cells_.size() + 1
          << " to match a number for each stratum";
    }
    unsigned year = 0;
    if (!utilities::To<unsigned>(error_data[row_counter][0], year))
      LOG_ERROR_P(PARAM_ERROR_VALUES) << " value " << error_data[row_counter][0] << " could not be converted in to an unsigned integer. It should be the year for this line";
    if (std::find(years_.begin(), years_.end(), year) == years_.end())
      LOG_ERROR_P(PARAM_ERROR_VALUES) << " value " << year << " is not a valid year for this observation";

    LOG_MEDIUM() << "Year = " << year << " in table = " << error_data[row_counter][0];
    for (unsigned i = 1; i < error_data[row_counter].size(); ++i) {
      unsigned value = 0;

      if (!utilities::To<unsigned>(error_data[row_counter][i], value))
        LOG_FATAL_P(PARAM_ERROR_VALUES) << " value (" << error_data[row_counter][i] << " at row = " << row_counter << " and column " << i << " could not be converted to an integer";
      if (value < 0) {
        LOG_ERROR_P(PARAM_ERROR_VALUES) << "at row = " << row_counter << " and column " << i << " the value given = " << value << " this needs to be an integer greater than 0.";
      }

      LOG_MEDIUM() << "Year = " << year << " stratum  = " << stratum_order_from_error_table[i] << " = " << value;

      error_value_by_year_and_stratum_[year][stratum_order_from_error_table[i]] = value;
    }
  }

}

/**
 *
 */
void MortalityEventComposition::PreExecute() {

}

/**
 *
 */
void MortalityEventComposition::Execute() {

}
/**
 *  Reset dynamic containers, in between simulations
 */
void MortalityEventComposition::ResetPreSimulation() {
  LOG_FINE() << "ResetPreSimulation";
/*
  for(unsigned i = 0; i < years_.size(); ++i) {
    for(unsigned j = 0; j < cells_.size(); ++j) {

    }
  }
*/
  ClearComparison(); // Clear comparisons
  LOG_FINE() << "Exit: ResetPreSimulation";
}

/**
 *
 */
void MortalityEventComposition::Simulate() {
  LOG_MEDIUM() << "Simulating data for observation = " << label_;
  ResetPreSimulation();
  vector<vector<processes::census_data>> fishery_age_data = mortality_process_->get_fishery_census_data(fishery_label_);

  if (fishery_age_data.size() == 0) {
    LOG_CODE_ERROR() << "could not return information for fishery " << fishery_label_ << " this error should be dealt with earlier in the code";
  }
  //vector<processes::composition_data>& length_frequency = mortality_process_->get_removals_by_length();
  LOG_FINE() << "length of census data (years * cells) = " << fishery_age_data.size();

  vector<unsigned> census_stratum_ndx;
  if (is_age_) {
    fishery_comp_data_ =  mortality_process_->get_fishery_age_comp(fishery_label_);
  } else {
    fishery_comp_data_ =  mortality_process_->get_fishery_length_comp(fishery_label_);
  }

  /*
   * Year loop
   */
  for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    // find equivalent fishery index
    auto iter = find(fishery_years_.begin(), fishery_years_.end(), years_[year_ndx]);
    unsigned fishery_year_ndx = distance(fishery_years_.begin(), iter);
    LOG_MEDIUM() << "About to sort our info for year " << years_[year_ndx] << " fishery index " << fishery_year_ndx;
    /*
     * Stratum Loop
     */
    for (unsigned stratum_ndx = 0; stratum_ndx < cells_.size(); ++stratum_ndx) {
      stratum_age_frequency_[cells_[stratum_ndx]].resize(model_->age_spread(), 0.0);
      LOG_FINE() << "About to sort our info for stratum " << cells_[stratum_ndx];
      stratum_biomass_[cells_[stratum_ndx]] = 0.0;
      fill(stratum_comp_.begin(),stratum_comp_.end(), 0.0);
      vector<float> biomass_by_cell;

      for (processes::composition_data& comp_data : (*fishery_comp_data_)[fishery_year_ndx]) {
        // Find census elements that are in this stratum
        if ((find(stratum_rows_[cells_[stratum_ndx]].begin(),stratum_rows_[cells_[stratum_ndx]].end(), comp_data.row_) != stratum_rows_[cells_[stratum_ndx]].end()) && (find(stratum_cols_[cells_[stratum_ndx]].begin(),stratum_cols_[cells_[stratum_ndx]].end(), comp_data.col_) != stratum_cols_[cells_[stratum_ndx]].end())) {
          LOG_MEDIUM() << "year = " << comp_data.year_ << " length of freq = " << comp_data.frequency_.size() << " and for female = " << comp_data.female_frequency_.size();
          if (sexed_flag_) {
            for (unsigned i = 0; i < comp_data.frequency_.size(); ++i)
              stratum_comp_[i] += comp_data.frequency_[i];
            for (unsigned i = 0; i < comp_data.female_frequency_.size(); ++i)
              stratum_comp_[i + n_bins_] += comp_data.female_frequency_[i];
          } else {
            for (unsigned i = 0; i < comp_data.frequency_.size(); ++i)
              stratum_comp_[i] += comp_data.frequency_[i];
          }
        }
      }
      /*
       * Now do some age-missclassification if age or store the results.
       *
       */
      if (is_age_ & !ageing_error_) {
        vector<float> vector_for_ageing(stratum_comp_.size(), 0.0); // = rep_vector(0.0, A);

        if (sexed_flag_) {
          for (unsigned i = 0; i < n_unsexed_bins_; ++i) {
            for (unsigned j = 0; j < n_unsexed_bins_; ++j)
              vector_for_ageing[j] += stratum_comp_[i] * ageing_mis_matrix_[i][j];
          }
          for (unsigned i = 0; i < n_unsexed_bins_; ++i)
            stratum_comp_[i] = vector_for_ageing[i];

          fill(vector_for_ageing.begin(), vector_for_ageing.end(), 0.0);
          for (unsigned i = 0; i < n_unsexed_bins_; ++i) {
            for (unsigned j = 0; j < n_unsexed_bins_; ++j)
              vector_for_ageing[j] += stratum_comp_[i + n_bins_] * ageing_mis_matrix_[i][j];
          }
          for (unsigned i = 0; i < n_unsexed_bins_; ++i)
            stratum_comp_[i + n_bins_] = vector_for_ageing[i];
        } else {
          for (unsigned i = 0; i < n_unsexed_bins_; ++i) {
            for (unsigned j = 0; j < n_unsexed_bins_; ++j)
              vector_for_ageing[j] += stratum_comp_[i] * ageing_mis_matrix_[i][j];
          }
          for (unsigned i = 0; i < n_unsexed_bins_; ++i)
            stratum_comp_[i] = vector_for_ageing[i];
        }
      }

      // Save the observationss
      if (sexed_flag_) {
        if (is_age_) {
          for (unsigned i = 0; i < n_unsexed_bins_; ++i)
            SaveComparison(i + model_->min_age(), 0, 0.0, cells_[stratum_ndx], stratum_comp_[i], 0.0, error_value_by_year_and_stratum_[ years_[year_ndx]][cells_[stratum_ndx]], years_[year_ndx] );
          for (unsigned i = 0; i < n_unsexed_bins_; ++i)
            SaveComparison(i + model_->min_age(), 1, 0.0, cells_[stratum_ndx], stratum_comp_[i + n_bins_], 0.0, error_value_by_year_and_stratum_[ years_[year_ndx]][cells_[stratum_ndx]], years_[year_ndx] );
        } else {
          for (unsigned i = 0; i < n_unsexed_bins_; ++i)
            SaveComparison(0, 0, model_->length_bin_mid_points()[i], cells_[stratum_ndx], stratum_comp_[i], 0.0, error_value_by_year_and_stratum_[ years_[year_ndx]][cells_[stratum_ndx]], years_[year_ndx] );
          for (unsigned i = 0; i < n_unsexed_bins_; ++i)
            SaveComparison(0, 1,  model_->length_bin_mid_points()[i], cells_[stratum_ndx], stratum_comp_[i + n_bins_], 0.0, error_value_by_year_and_stratum_[ years_[year_ndx]][cells_[stratum_ndx]], years_[year_ndx] );
        }
      } else {
        if (is_age_) {
          for (unsigned i = 0; i < stratum_comp_.size(); ++i) {
            SaveComparison(i + model_->min_age(), 0, 0.0, cells_[stratum_ndx], stratum_comp_[i], 0.0, error_value_by_year_and_stratum_[ years_[year_ndx]][cells_[stratum_ndx]], years_[year_ndx] );
          }
         } else {
          for (unsigned i = 0; i < stratum_comp_.size(); ++i)
            SaveComparison(0, 0, model_->length_bin_mid_points()[i], cells_[stratum_ndx], stratum_comp_[i], 0.0, error_value_by_year_and_stratum_[ years_[year_ndx]][cells_[stratum_ndx]], years_[year_ndx] );
         }
      }
    } // stratum
  } // year

  /**
   * Simulate or generate results
   * During simulation mode we'll simulate results for this observation
   */
  LOG_MEDIUM() << "Calculating score for observation = " << label_;
  // Convert to propotions before simulating for each year and cell sum = 1
  float total = 0.0;
  vector<float> total_by_cell(cells_.size() * years_.size(), 0.0);
  unsigned counter = 0;
  for (auto& iter : comparisons_) { // cell
    for (auto& second_iter : iter.second) {  // year
      total = 0.0;
      for (auto& comparison : second_iter.second) {
        total += comparison.expected_;
      }
      total_by_cell[counter] = total;
      ++counter;
      for (auto& comparison : second_iter.second)
        comparison.expected_ /= total;
    }
  }

  likelihood_->SimulateObserved(comparisons_);
  // Simualte numbers at age, but we want proportion
  counter = 0;
  for (auto& iter : comparisons_) {
    for (auto& second_iter : iter.second) {  // cell
      total = 0.0;
      for (auto& comparison : second_iter.second)
        total += comparison.simulated_;
      if (are_obs_props_ & ((likelihood_->type() == PARAM_MULTINOMIAL) | (likelihood_->type() == PARAM_DIRICHLET))) {
        for (auto& comparison : second_iter.second)
          comparison.simulated_ /= total;
      }
      if (not are_obs_props_& (likelihood_->type() == PARAM_LOGISTIC_NORMAL)) {
        for (auto& comparison : second_iter.second)
          comparison.simulated_ *=  total_by_cell[counter];
      }
      ++counter;
    }
  }
} // Simulate



void MortalityEventComposition::FillReportCache(ostringstream& cache) {
  LOG_MEDIUM() << "we are here";


}

} /* namespace observations */
} /* namespace niwa */

