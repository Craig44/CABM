/**
 * @file AgeLength.cpp
 * @author  C.Marsh
 * @date 01/8/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "AgeLength.h"

//#include "Catchabilities/Manager.h"

#include "TimeSteps/Manager.h"
#include "Layers/Manager.h"
#include "Selectivities/Manager.h"

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
AgeLength::AgeLength(Model* model) : Observation(model) {
  parameters_.Bind<string>(PARAM_TIME_STEP, &time_step_label_, "The label of time-step that the observation occurs in", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "The years of the observed values", "");
  parameters_.Bind<string>(PARAM_SELECTIVITIES, &selectivity_labels_, "Labels of the selectivities", "", true);
  parameters_.Bind<string>(PARAM_LAYER_OF_CELLS, &layer_label_, "The layer that indicates what area to summarise observations over.", "");
  parameters_.Bind<string>(PARAM_CELLS, &cells_, "The cells we want to generate observations for from the layer of cells supplied", "");
  parameters_.Bind<unsigned>(PARAM_NUMBER_OF_SAMPLES, &n_samples_, "The number of samples to collect from each cell", "");

  allowed_likelihood_types_.push_back(PARAM_PSEUDO);
}

/**
 *
 */
void AgeLength::DoValidate() {
  LOG_TRACE();
  // Check there are not duplicate cells
  for (unsigned i = 0; i < cells_.size(); ++i) {
    for (unsigned j = 0; j < cells_.size(); ++j) {
      if (j == i)
        continue;
      if (cells_[i] == cells_[j])
        LOG_ERROR_P(PARAM_CELLS) << "found a non-unique cell '" << cells_[i] << "' at element " << i + 1 << " which matches " << cells_[j] << " at element " << j + 1 << " please provide unique cell labels, otherwise you may create undefined behaviour from me =)";
    }
  }

  if (n_samples_.size() != cells_.size())
    LOG_ERROR_P(PARAM_CELLS) << "number of samples must equal number of cells please sort this out, Chairs";
}

/**
 *
 */
void AgeLength::DoBuild() {
  LOG_TRACE();

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

  world_ = model_->world_view();
  if (!world_)
    LOG_CODE_ERROR() << "!world_ could not create pointer to world viw model, something is wrong";

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

  // Build map of indexes for each cell
  for (auto cell : cells_) {
    vector<unsigned> col_index;
    vector<unsigned> row_index;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        if (layer_->get_value(row,col) == cell) {
          col_index.push_back(col);
          row_index.push_back(row);
        }
      }
    }
    col_index_for_each_cell_[cell] = col_index;
    row_index_for_each_cell_[cell] = row_index;
  }


  // Subscribe this observation to the timestep
  auto time_step = model_->managers().time_step()->GetTimeStep(time_step_label_);
  if (!time_step) {
    LOG_ERROR_P(PARAM_TIME_STEP) << time_step_label_ << " could not be found. Have you defined it?";
  } else {
    for (unsigned year : years_)
      time_step->SubscribeToBlock(this, year);
  }
  for (auto year : years_) {
    if((year < model_->start_year()) || (year > model_->final_year()))
      LOG_ERROR_P(PARAM_YEARS) << "Years can't be less than start_year (" << model_->start_year() << "), or greater than final_year (" << model_->final_year() << "). Please fix this.";
  }
}

/**
 *
 */
void AgeLength::PreExecute() {
  //LOG_FINE();

}

/**
 *
 */
void AgeLength::Execute() {
  LOG_FINE();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  unsigned samples_to_take;
  unsigned n_cells;

  if (selectivity_length_based_) {
    for (unsigned cell_ndx = 0; cell_ndx < cells_.size(); ++cell_ndx) {
      samples_to_take = n_samples_[cell_ndx];
      n_cells = col_index_for_each_cell_[cells_[cell_ndx]].size();
      WorldCell* sample_cell;
      while(samples_to_take > 0) {
        // randomly select a cell in this region
        auto sample_index = rng.chance() * n_cells;
        // get a pointer to the world and select a fish
        sample_cell = world_->get_base_square(row_index_for_each_cell_[cells_[cell_ndx]][sample_index], col_index_for_each_cell_[cells_[cell_ndx]][sample_index]);
        // Randomly select an individual and see if it is selective
        auto agent = sample_cell->agents_.begin();
        advance(agent, sample_cell->agents_.size() * rng.chance());
        if (rng.chance() <= selectivities_[(*agent).get_sex()]->GetResult((*agent).get_length_bin_index()))
          SaveComparison((*agent).get_age(), (*agent).get_length(), sample_cell->get_cell_label(), 0.0, 0.0, 0.0, model_->current_year());
      }
    }
  } else {
    for (unsigned cell_ndx = 0; cell_ndx < cells_.size(); ++cell_ndx) {
      samples_to_take = n_samples_[cell_ndx];
      n_cells = col_index_for_each_cell_[cells_[cell_ndx]].size();
      WorldCell* sample_cell;
      while(samples_to_take > 0) {
        // randomly select a cell in this region
        auto sample_index = rng.chance() * n_cells;
        // get a pointer to the world and select a fish
        sample_cell = world_->get_base_square(row_index_for_each_cell_[cells_[cell_ndx]][sample_index], col_index_for_each_cell_[cells_[cell_ndx]][sample_index]);
        // Randomly select an individual and see if it is selective
        auto agent = sample_cell->agents_.begin();
        advance(agent, sample_cell->agents_.size() * rng.chance());
        if (rng.chance() <= selectivities_[(*agent).get_sex()]->GetResult((*agent).get_age()))
          SaveComparison((*agent).get_age(), (*agent).get_length(), sample_cell->get_cell_label(), 0.0, 0.0, 0.0, model_->current_year());
      }
    }
  }
}

/**
 *
 */
void AgeLength::Simulate() {
  /**
   * Simulate or generate results
   * During simulation mode we'll simulate results for this observation
   */
	LOG_FINEST() << "Calculating score for observation = " << label_;
  //likelihood_->SimulateObserved(comparisons_);  // this is an usual observation and I am not going to simulate from it
}


} /* namespace observations */
} /* namespace niwa */

