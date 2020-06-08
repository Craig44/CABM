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
  parameters_.Bind<unsigned>(PARAM_TAG_RELEASE_YEAR, &tag_release_year_, "The years that the tagged fish were released", "");

  parameters_.Bind<bool>(PARAM_SEXED, &sexed_, "Seperate observation by sex", "", true);
  parameters_.Bind<string>(PARAM_PROCESS_LABEL, &process_label_, "Label of of removal process", "", "");

  parameters_.Bind<string>(PARAM_LAYER_OF_STRATA_DEFINITIONS, &layer_label_, "The layer that indicates what the stratum boundaries are.", "");
  parameters_.Bind<string>(PARAM_RELEASE_STRATUM, &release_stratum_, "Stratum that individuals were released in", "");
  parameters_.Bind<string>(PARAM_RECAPTURE_STRATUM, &recapture_stratum_, "Stratum to record recaptures, that were released in this year and in the release stratum", "", true);
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

  if (mortality_process_->type() != PARAM_MORTALITY_EVENT_BIOMASS) {
    LOG_ERROR_P(PARAM_PROCESS_LABEL) << "We suggest that this observation only be applied for the process of type " << PARAM_MORTALITY_EVENT_BIOMASS;
  }


    // Build and validate layers
  layer_ = model_->managers().layer()->GetCategoricalLayer(layer_label_);
  if (!layer_)
    LOG_FATAL_P(PARAM_LAYER_OF_STRATUM_DEFINITIONS)<< "could not find layer " << layer_label_ << " does it exist?, if it exists is of type categorical?";


  LOG_FINE() << "Check stratum are consistent";
  if (parameters_.Get(PARAM_RECAPTURE_STRATUM)->has_been_defined()) {
    // Check all the cells supplied are in the layer
    LOG_FINE() << "checking cell " << release_stratum_ << " exists in layer";
    bool cell_found = false;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        LOG_FINE() << "checking row = " << row << " col = " << col;
        LOG_FINE() << "value = " << layer_->get_value(row, col);
        if (layer_->get_value(row, col) == release_stratum_) {
          cell_found = true;
          release_stratum_rows_.push_back(row);
          release_stratum_cols_.push_back(col);
        }
      }
    }
    if (not cell_found)
      LOG_ERROR_P(PARAM_STRATUMS_TO_INCLUDE) << "could not find the cell '" << release_stratum_ << "' in the layer " << layer_label_
          << " please make sure that you supply cell labels that are consistent with the layer.";

    // Check all the cells supplied are in the layer
    for (auto cell : recapture_stratum_) {
      LOG_FINE() << "checking cell " << cell << " exists in layer";
      cell_found = false;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          LOG_FINE() << "checking row = " << row << " col = " << col;
          LOG_FINE() << "value = " << layer_->get_value(row, col);
          if (layer_->get_value(row, col) == cell) {
            cell_found = true;
            recapture_stratum_rows_[cell].push_back(row);
            recapture_stratum_cols_[cell].push_back(col);
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
        if (temp_cell == release_stratum_) {
          release_stratum_rows_.push_back(row);
          release_stratum_cols_.push_back(col);
        }
        recapture_stratum_rows_[temp_cell].push_back(row);
        recapture_stratum_cols_[temp_cell].push_back(col);
        recapture_stratum_.push_back(temp_cell);
        if (find(recapture_stratum_.begin(), recapture_stratum_.end(), temp_cell) == recapture_stratum_.end())
          recapture_stratum_.push_back(temp_cell);
      }
    }
  }
  LOG_FINE() << "Check mortality process is consistent with observation";
  if (not mortality_process_->check_years(years_)) {
    LOG_ERROR_P(PARAM_YEARS)
        << "there was a year that the mortality process doesn't not execute in, can you please check that the years you have supplied for this observation are years that the mortality process occurs in cheers.";
  }

  length_freq_.resize(model_->number_of_length_bins(), 0.0);
  scanned_length_freq_.resize(model_->number_of_length_bins(), 0.0);
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

  //utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  tag_recapture_data_ = mortality_process_->get_tag_recapture_info();
  //vector<processes::composition_data>& length_frequency = mortality_process_->get_removals_by_length();
  LOG_FINE() << "length of census data = " << (*tag_recapture_data_).size();
  for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    LOG_MEDIUM() << "About to sort our info for year " << years_[year_ndx];
    for (unsigned recapture_area_ndx = 0; recapture_area_ndx < recapture_stratum_.size(); ++recapture_area_ndx) {
      LOG_MEDIUM() << "About to sort our info for recapture stratum " << recapture_stratum_[recapture_area_ndx];
      // Find out how many agents are available to be used for an ALK in this stratum
      // -- loop over all cells that the fishery occured in (census data)
      // -- if that area and fishery belongs to this stratum then summarise some numbers for use later.
      //
      fill(length_freq_.begin(), length_freq_.end(), 0.0);
      fill(scanned_length_freq_.begin(), scanned_length_freq_.end(), 0.0);

      // Calculate Length frequency for the strata
      unsigned tag_recap_ndx = 0; // links back to the tag-recapture data
      for (processes::tag_recapture &tag_recap : (*tag_recapture_data_)) {
        // Find tag_recap elements that are in this year and stratum
        if ((tag_recap.year_ == years_[year_ndx])
          && (find(recapture_stratum_rows_[recapture_stratum_[recapture_area_ndx]].begin(), recapture_stratum_rows_[recapture_stratum_[recapture_area_ndx]].end(), tag_recap.row_) != recapture_stratum_rows_[recapture_stratum_[recapture_area_ndx]].end())
          && (find(recapture_stratum_cols_[recapture_stratum_[recapture_area_ndx]].begin(), recapture_stratum_cols_[recapture_stratum_[recapture_area_ndx]].end(), tag_recap.col_) != recapture_stratum_cols_[recapture_stratum_[recapture_area_ndx]].end())) {
          if (tag_recap.age_.size() > 0) {
            for (unsigned length_ndx = 0; length_ndx < tag_recap.length_ndx_.size(); ++length_ndx) {
              if (tag_release_year_ == tag_recap.tag_release_year_[length_ndx] &&
                  (find(release_stratum_rows_.begin(),release_stratum_rows_.end(), tag_recap.tag_row_[length_ndx]) != release_stratum_rows_.end()) &&
                  (find(release_stratum_cols_.begin(),release_stratum_cols_.end(), tag_recap.tag_col_[length_ndx]) != release_stratum_cols_.end())) {
                length_freq_[tag_recap.length_ndx_[length_ndx]]++;
              }
            }
          }
          for(unsigned j = 0; j < tag_recap.scanned_length_comp_.size(); ++j)
            scanned_length_freq_[j] += tag_recap.scanned_length_comp_[j];
        }
        ++tag_recap_ndx;
      }
      for (unsigned length_bin_ndx = 0; length_bin_ndx < length_freq_.size(); ++length_bin_ndx)
        SaveComparison(0, model_->length_bin_mid_points()[length_bin_ndx], recapture_stratum_[recapture_area_ndx], length_freq_[length_bin_ndx], 0.0, scanned_length_freq_[length_bin_ndx], years_[year_ndx]);
    }
  }
  for (auto& iter : comparisons_) {
    for (auto& second_iter : iter.second) {  // cell
      for (auto& comparison : second_iter.second)
        comparison.simulated_ += comparison.expected_;
    }
  }
} // DoExecute

void TagRecaptureByLength::FillReportCache(ostringstream& cache) {
  cache << "release_stratum: " << release_stratum_ << "\n";
}

} /* namespace observations */
} /* namespace niwa */

