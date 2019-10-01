/**
 * @file TagRecaptureByLength.cpp
 * @author  C.Marsh github.com/Craig44
 * @date 4/11/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "TagRecaptureByLength.h"

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
TagRecaptureByLength::TagRecaptureByLength(Model* model) : Observation(model) {

  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "The years of the observed values", "");
  parameters_.Bind<bool>(PARAM_SEXED, &sexed_, "Seperate observation by sex", "", true);
  parameters_.Bind<string>(PARAM_PROCESS_LABEL, &process_label_, "Label of of removal process", "", "");
  parameters_.Bind<string>(PARAM_LAYER_OF_STRATUM_DEFINITIONS, &layer_label_, "The layer that indicates what the stratum boundaries are.", "");
  parameters_.Bind<string>(PARAM_STRATUMS_TO_INCLUDE, &cells_, "The cells which represent individual stratum to be included in the analysis, default is all cells are used from the layer", "", true);
}
/**
 * Destructor
 */
TagRecaptureByLength::~TagRecaptureByLength() {

}
/**
 *
 */
void TagRecaptureByLength::DoValidate() {
  LOG_TRACE();
  for (auto year : years_) {
    LOG_FINE() << "year : " << year;
    if ((year < model_->start_year()) || (year > model_->final_year()))
      LOG_ERROR_P(PARAM_YEARS) << "Years can't be less than start_year (" << model_->start_year() << "), or greater than final_year (" << model_->final_year()
          << "). Please fix this.";
  }

}

/**
 *
 */
void TagRecaptureByLength::DoBuild() {
  LOG_TRACE();
  world_ = model_->world_view();
  if (!world_)
    LOG_CODE_ERROR() << "!world_ could not create pointer to world viw model, something is wrong";

  mortality_process_ = model_->managers().process()->GetMortalityProcess(process_label_);
  if (!mortality_process_)
    LOG_FATAL_P(PARAM_PROCESS_LABEL)<< "could not find the process " << process_label_ << ", please make sure it exists";

  if (mortality_process_->type() == PARAM_MORTALITY_EVENT_BIOMASS) {
    LOG_ERROR_P(PARAM_PROCESS_LABEL) << "We suggest that age observations for the mortality process of type " << PARAM_MORTALITY_EVENT_BIOMASS << " you use the observation type " << PARAM_MORTALITY_EVENT_BIOMASS_SCALED_AGE_FREQUENCY;
  }


    // Build and validate layers
  layer_ = model_->managers().layer()->GetCategoricalLayer(layer_label_);
  if (!layer_)
    LOG_FATAL_P(PARAM_LAYER_OF_STRATUM_DEFINITIONS)<< "could not find layer " << layer_label_ << " does it exist?, if it exists is of type categorical?";

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
}

/**
 *
 */
void TagRecaptureByLength::PreExecute() {

}

/**
 *
 */
void TagRecaptureByLength::Execute() {

}

/**
 *
 */
void TagRecaptureByLength::Simulate() {
  LOG_MEDIUM() << "Simulating data for observation = " << label_;
  ClearComparison(); // Clear comparisons

  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  vector<processes::tag_recapture>& tag_recapture_data = mortality_process_->get_tag_recapture_info();
  //vector<processes::composition_data>& length_frequency = mortality_process_->get_removals_by_length();
  LOG_FINE() << "length of census data = " << tag_recapture_data.size();
  vector<unsigned> census_stratum_ndx;

  for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    LOG_FINE() << "About to sort our info for year " << years_[year_ndx];
    for (unsigned stratum_ndx = 0; stratum_ndx < cells_.size(); ++stratum_ndx) {

    }
  }
} // DoExecute

} /* namespace observations */
} /* namespace niwa */

