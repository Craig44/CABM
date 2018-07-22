/**
 * @file MovementBoxTransfer.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 19/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 *
 */

// headers
#include "MovementBoxTransfer.h"

#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/To.h"

#include "World/WorldCell.h"
#include "World/WorldView.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <omp.h>
// namespaces
namespace niwa {
namespace processes {

/**
 *  constructor
 */
MovementBoxTransfer::MovementBoxTransfer(Model* model) : Process(model) {
  process_type_ = ProcessType::kTransition;
  parameters_.Bind<string>(PARAM_SELECTIVITY, &selectivity_label_, "Selectivity label", "");
  parameters_.Bind<string>(PARAM_ORIGIN_CELL, &origin_cell_, "The origin cell associated with each spatial layer (should have a one to one relationship with specified layers), format follows row-col (1-2)", "");
  parameters_.Bind<string>(PARAM_PROBABILITY_LAYERS, &probability_layer_labels_, "Spatial layers (one layer for each origin cell) describing the probability of moving from an origin cell to all other cells in the spatial domain. ", "");
  // TODO markovian (current implementation) vs natal i.e origin cell refers to where they are (Markovian) or where they were born (Natal Homing)
}

void MovementBoxTransfer::DoValidate() {
  LOG_TRACE();
  if (origin_cell_.size() != probability_layer_labels_.size())
    LOG_ERROR_P(PARAM_PROBABILITY_LAYERS) << "you must specify a layer label for each origin cell. You have supplied '" << origin_cell_.size() << "' origin cells but '" << probability_layer_labels_.size() << "' probability layer labels, please sort this out.";
}
/**
 *  Build relationships with other classes
 */
void MovementBoxTransfer::DoBuild() {
  LOG_TRACE();
  // Get the layers
  for (auto& label : probability_layer_labels_) {
    layers::NumericLayer* temp_layer = nullptr;
    temp_layer = model_->managers().layer()->GetNumericLayer(label);
    if (!temp_layer) {
      LOG_FATAL_P(PARAM_PROBABILITY_LAYERS) << "could not find the layer '" << label << "', please make sure it exists and that it is type 'numeric'";
    }
    probability_layers_.push_back(temp_layer);
  }

  selectivity_ = model_->managers().selectivity()->GetSelectivity(selectivity_label_);
  if (!selectivity_)
    LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << ": selectivity " << selectivity_label_ << " does not exist. Have you defined it?";


  // Split out cell origins for quick look up at execution
  for (auto origin : origin_cell_) {
    unsigned value = 1;
    vector<string> split_cells;
    boost::split(split_cells, origin, boost::is_any_of("-"));

    LOG_FINEST() << "row = " << split_cells[0] << " col = " << split_cells[1];
    if (!utilities::To<unsigned>(split_cells[0], value))
      LOG_ERROR_P(PARAM_ORIGIN_CELL) << " value (" << split_cells[0] << ") could not be converted to a unsigned";
    origin_rows_.push_back(value - 1);
    if (!utilities::To<unsigned>(split_cells[1], value))
      LOG_ERROR_P(PARAM_ORIGIN_CELL) << " value (" << split_cells[1] << ") could not be converted to a unsigned";
    origin_cols_.push_back(value - 1);
  }

  // Check the
  if (world_->get_enabled_cells() != origin_cell_.size())
    LOG_ERROR_P(PARAM_ORIGIN_CELL) << "you haven't supplied an origin cell for all enabled cells in the spatial domain. you have supplied '" << origin_cell_.size() << "' origin cells, but there are '" << world_->get_enabled_cells() << "' enabled cells in the domain, please check this out";

  possible_rows_ = world_->get_enabled_rows();
  possible_cols_ = world_->get_enabled_cols();

}


/**
 *  Execute process
 */
void MovementBoxTransfer::DoExecute() {
  LOG_TRACE();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  // Iterate over origin cells
  #pragma omp parallel for collapse(2)
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      // Find the right probability layer for this combo
      unsigned origin_element = 0;
      for (; origin_element < origin_rows_.size(); ++origin_element) {
        if ((origin_rows_[origin_element] == row) && (origin_cols_[origin_element] == col))
          break;
      }
      LOG_FINEST() << "row = " << row << " col = " << col << " layer pointer vlaue = " << origin_element;
      WorldCell* origin_cell = world_->get_base_square(row, col);
      WorldCell* destination_cell;
      if (origin_cell->is_enabled()) {
        MovementData store_infor(model_->get_height(), model_->get_width(), origin_cell_[origin_element], model_->current_year(), origin_cell->agents_.size());
        for (auto iter = origin_cell->agents_.begin(); iter != origin_cell->agents_.end(); ++iter) {
          // Pick a random destination cell
          unsigned destination_row = possible_rows_[possible_rows_.size() * rng.chance()];
          unsigned destination_col = possible_cols_[possible_cols_.size() * rng.chance()];
          if (rng.chance() <= (probability_layers_[origin_element]->get_value(destination_row, destination_col) * selectivity_->GetResult((*iter).get_age()))) {
            // We are moving splice this agent to the destination cache cell
            store_infor.destination_of_agents_moved_[destination_row][destination_col]++;
            #pragma omp critical // Make a critical section in case multiple threads writing to same destination cell, thread protection // TODO this might slow things down, we could look at putting mutex's on each cell of the grid, and locking and unlocking each destination cell?
            {
              destination_cell = world_->get_cached_square(destination_row, destination_col);
              auto nx = next(iter); // Need to next the iter else we iter changes scope to cached agents, an annoying stl thing
              destination_cell->agents_.splice(destination_cell->agents_.end(), origin_cell->agents_, iter);
              iter = nx;
            }
          }
        }
        if (model_->state() != State::kInitialise) {
          #pragma omp critical
          {
            moved_agents_by_year_.push_back(store_infor);
          }
        }
      }
    }
  }
  // merge destination agents into the actual grid
  world_->MergeCachedGrid();
}

// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MovementBoxTransfer::FillReportCache(ostringstream& cache) {
  for (auto& values : moved_agents_by_year_) {
    cache << "initial_numbers: " << values.initial_numbers_ << "\n";
    cache << values.year_ << "_" << values.origin_cell_ << " " << REPORT_R_MATRIX << "\n";
    for (unsigned i = 0; i < model_->get_height(); ++i) {
      for (unsigned j = 0; j < model_->get_width(); ++j)
        cache << values.destination_of_agents_moved_[i][j] << " ";
      cache << "\n";
    }
  }
}
} /* namespace processes */
} /* namespace niwa */

