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
// namespaces
namespace niwa {
namespace processes {

/**
 *  constructor
 */
MovementBoxTransfer::MovementBoxTransfer(Model* model) : Movement(model) {
  process_type_ = ProcessType::kTransition;
  parameters_.Bind<string>(PARAM_ORIGIN_CELL, &origin_cell_, "The origin cell associated with each spatial layer (should have a one to one relationship with specified layers), format follows row-col (1-2)", "");
  parameters_.Bind<string>(PARAM_PROBABILITY_LAYERS, &probability_layer_labels_, "Spatial layers (one layer for each origin cell) describing the probability of moving from an origin cell to all other cells in the spatial domain. ", "");
  parameters_.Bind<string>(PARAM_MOVEMENT_TYPE, &movement_type_string_, "What type of movement are you applying?", "", PARAM_MARKOVIAN)->set_allowed_values({PARAM_MARKOVIAN, PARAM_NATAL_HOMING});
  parameters_.Bind<string>(PARAM_SELECTIVITY_LABEL, &selectivity_label_, "Label for the selectivity block", "");
}

void MovementBoxTransfer::DoValidate() {
  LOG_TRACE();
  if (origin_cell_.size() != probability_layer_labels_.size())
    LOG_ERROR_P(PARAM_PROBABILITY_LAYERS) << "you must specify a layer label for each origin cell. You have supplied '" << origin_cell_.size() << "' origin cells but '" << probability_layer_labels_.size() << "' probability layer labels, please sort this out.";

  if (movement_type_string_ == PARAM_MARKOVIAN)
    movement_type_ = MovementType::kMarkovian;
  else if (movement_type_string_ == PARAM_NATAL_HOMING)
    movement_type_ = MovementType::kNatal_homing;
}

/**
 *  Build relationships with other classes
 */
void MovementBoxTransfer::DoBuild() {
  LOG_FINE();
  // Get the layers
  for (auto& label : probability_layer_labels_) {
    layers::NumericLayer* temp_layer = nullptr;
    temp_layer = model_->managers().layer()->GetNumericLayer(label);
    if (!temp_layer) {
      LOG_FATAL_P(PARAM_PROBABILITY_LAYERS) << "could not find the layer '" << label << "', please make sure it exists and that it is type 'numeric'";
    }
    probability_layers_.push_back(temp_layer);
  }

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

  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        possible_rows_.push_back(row);
        possible_cols_.push_back(col);
      }
    }
  }

  cell_offset_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i)
    cell_offset_[i].resize(model_->get_width());

  LOG_FINEST() << "selectivities supplied = " << selectivity_label_.size();
  // Build selectivity links
  if (selectivity_label_.size() == 1)
    selectivity_label_.assign(2, selectivity_label_[0]);

  if (selectivity_label_.size() > 2) {
    LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << "You suppled " << selectivity_label_.size()  << " Selectiviites, you can only have one for each sex max = 2";
  }
  LOG_FINEST() << "selectivities supplied = " << selectivity_label_.size();

  bool first = true;
  for (string label : selectivity_label_) {
    Selectivity* temp_selectivity = model_->managers().selectivity()->GetSelectivity(label);
    if (!temp_selectivity)
      LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << ": selectivity " << label << " does not exist. Have you defined it?";

    selectivity_.push_back(temp_selectivity);
    if (first) {
      first = false;
      selectivity_length_based_ = temp_selectivity->is_length_based();
    } else {
      if (selectivity_length_based_ != temp_selectivity->is_length_based()) {
        LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << "The selectivity  " << label << " was not the same type (age or length based) as the previous selectivity label";
      }
    }
  }
}


/**
 *  Execute process
 */
void MovementBoxTransfer::DoExecute() {
  LOG_FINE();

  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  // Pre-calculate agents in the world to set aside our random numbers needed for the operation
  n_agents_ = 0;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        cell_offset_[row][col] = n_agents_;
        n_agents_ += cell->agents_.size();
      }
    }
  }
  // Allocate a single block of memory rather than each thread temporarily allocating their own memory.
  random_numbers_.resize(n_agents_);
  for (unsigned i = 0; i < n_agents_; ++i)
    random_numbers_[i] = rng.chance();


  if (movement_type_ == MovementType::kMarkovian) {
    LOG_FINE() << "applying markovian movement";
    // Iterate over origin cells
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        // Find the right probability layer for this combo
        unsigned origin_element = 0;
        for (; origin_element < origin_rows_.size(); ++origin_element) {
          if ((origin_rows_[origin_element] == row) && (origin_cols_[origin_element] == col))
            break;
        }
        LOG_FINEST() << "row = " << row << " col = " << col << " layer pointer vlaue = " << origin_element << " thread id = " << omp_get_thread_num();
        WorldCell* origin_cell = world_->get_base_square(row, col);
        WorldCell* destination_cell;
        if (origin_cell->is_enabled()) {
          LOG_FINEST() << "number of agents in this cell = " << origin_cell->agents_.size();
          MovementData store_infor(model_->get_height(), model_->get_width(), origin_cell_[origin_element], model_->current_year());
          float temp_sum;
          unsigned counter = 0;
          unsigned counter_jump = 0;

          if (selectivity_length_based_) {
            LOG_FINE() << "selectivity is length based";
            for (auto iter = origin_cell->agents_.begin(); iter != origin_cell->agents_.end(); ++iter) {
              if ((*iter).is_alive()) {
                if (rng.chance() <= selectivity_[(*iter).get_sex()]->GetResult((*iter).get_length_bin_index())) {

                  LOG_FINEST() << "row = " << row << " home row = "<< (*iter).get_home_col() << " col = " << col << " home col = " << (*iter).get_home_col();
                  ++counter;
                  // Iterate over possible cells compare to chance()
                  temp_sum = 0;
                  // Calcualte a multinomial via the following algorithm
                  // compare a random standard uniform value to the cumulative sums of the probabilities and return the first index for which
                  // the cumulative sum is greater than the random uniform
                  for (unsigned potential_destination = 0; potential_destination < possible_rows_.size(); ++potential_destination) {
                    temp_sum += probability_layers_[origin_element]->get_value(possible_rows_[potential_destination], possible_cols_[potential_destination]);
                    if (temp_sum > random_numbers_[cell_offset_[row][col] + counter]) {
                      ++counter_jump;
                      LOG_FINEST() << counter  << " cum prob = " << temp_sum << " random = " << random_numbers_[cell_offset_[row][col] + counter] << " current row = " << row << " home row = " << (*iter).get_home_row() << " current col = " << col << " home row = " << (*iter).get_home_col() << " destination row = " << possible_rows_[potential_destination] << " destination col = " << possible_cols_[potential_destination];
                      //LOG_FINEST() << counter << " removals = " << counter_junp <<   " iter distance = " << distance(origin_cell->agents_.begin(), iter) << " cum prob = " << temp_sum << " random = " << random << " current row = " << row << " current col = " << col << "destination row = " << possible_rows_[potential_destination] << " destination col = " << possible_cols_[potential_destination];
                      store_infor.destination_of_agents_moved_[possible_rows_[potential_destination]][possible_cols_[potential_destination]]++;
                      // if we are moving to this cell lets not move in memory
                      if ((possible_rows_[potential_destination] == row) && (possible_cols_[potential_destination] == col)) {

                        break;
                      }
                      // Make a synchronisation point don't want multiple threads accessing the same pointer simultaneously and splicing to it

                      destination_cell = world_->get_cached_square(possible_rows_[potential_destination], possible_cols_[potential_destination]);
                      destination_cell->agents_.push_back((*iter));
                      // Essentially we have to kill an individual here as it is gone off to be merged
                      (*iter).dies();
                      break;
                    }
                  }
                }
              }
            }
          } else {
            LOG_FINE() << "selectivity is Age based";
            for (auto iter = origin_cell->agents_.begin(); iter != origin_cell->agents_.end(); ++iter) {
              if ((*iter).is_alive()) {
                if (rng.chance() <= selectivity_[(*iter).get_sex()]->GetResult((*iter).get_age_index())) {

                  LOG_FINEST() << "row = " << row << " home row = "<< (*iter).get_home_col() << " col = " << col << " home col = " << (*iter).get_home_col();
                  ++counter;
                  // Iterate over possible cells compare to chance()
                  temp_sum = 0;
                  // Calcualte a multinomial via the following algorithm
                  // compare a random standard uniform value to the cumulative sums of the probabilities and return the first index for which
                  // the cumulative sum is greater than the random uniform
                  for (unsigned potential_destination = 0; potential_destination < possible_rows_.size(); ++potential_destination) {
                    temp_sum += probability_layers_[origin_element]->get_value(possible_rows_[potential_destination], possible_cols_[potential_destination]);
                    if (temp_sum > random_numbers_[cell_offset_[row][col] + counter]) {
                      ++counter_jump;
                      LOG_FINEST() << counter  << " cum prob = " << temp_sum << " random = " << random_numbers_[cell_offset_[row][col] + counter] << " current row = " << row << " home row = " << (*iter).get_home_row() << " current col = " << col << " home row = " << (*iter).get_home_col() << " destination row = " << possible_rows_[potential_destination] << " destination col = " << possible_cols_[potential_destination];
                      //LOG_FINEST() << counter << " removals = " << counter_junp <<   " iter distance = " << distance(origin_cell->agents_.begin(), iter) << " cum prob = " << temp_sum << " random = " << random << " current row = " << row << " current col = " << col << "destination row = " << possible_rows_[potential_destination] << " destination col = " << possible_cols_[potential_destination];
                      store_infor.destination_of_agents_moved_[possible_rows_[potential_destination]][possible_cols_[potential_destination]]++;
                      // if we are moving to this cell lets not move in memory
                      if ((possible_rows_[potential_destination] == row) && (possible_cols_[potential_destination] == col)) {

                        break;
                      }
                      // Make a synchronisation point don't want multiple threads accessing the same pointer simultaneously and splicing to it

                      destination_cell = world_->get_cached_square(possible_rows_[potential_destination], possible_cols_[potential_destination]);
                      destination_cell->agents_.push_back((*iter));
                      // Essentially we have to kill an individual here as it is gone off to be merged
                      (*iter).dies();
                      break;
                    }
                  }
                }
              }
            }
          }
          for (auto iter = origin_cell->agents_.begin(); iter != origin_cell->agents_.end(); ++iter) {
            if ((*iter).is_alive()) {
              LOG_FINEST() << "row = " << row << " home row = "<< (*iter).get_home_col() << " col = " << col << " home col = " << (*iter).get_home_col();
              ++counter;
              // Iterate over possible cells compare to chance()
              temp_sum = 0;
              // Calcualte a multinomial via the following algorithm
              // compare a random standard uniform value to the cumulative sums of the probabilities and return the first index for which
              // the cumulative sum is greater than the random uniform
              for (unsigned potential_destination = 0; potential_destination < possible_rows_.size(); ++potential_destination) {
                temp_sum += probability_layers_[origin_element]->get_value(possible_rows_[potential_destination], possible_cols_[potential_destination]);
                if (temp_sum > random_numbers_[cell_offset_[row][col] + counter]) {
                  ++counter_jump;
                  LOG_FINEST() << counter  << " cum prob = " << temp_sum << " random = " << random_numbers_[cell_offset_[row][col] + counter] << " current row = " << row << " home row = " << (*iter).get_home_row() << " current col = " << col << " home row = " << (*iter).get_home_col() << " destination row = " << possible_rows_[potential_destination] << " destination col = " << possible_cols_[potential_destination];
                  //LOG_FINEST() << counter << " removals = " << counter_junp <<   " iter distance = " << distance(origin_cell->agents_.begin(), iter) << " cum prob = " << temp_sum << " random = " << random << " current row = " << row << " current col = " << col << "destination row = " << possible_rows_[potential_destination] << " destination col = " << possible_cols_[potential_destination];
                  store_infor.destination_of_agents_moved_[possible_rows_[potential_destination]][possible_cols_[potential_destination]]++;
                  // if we are moving to this cell lets not move in memory
                  if ((possible_rows_[potential_destination] == row) && (possible_cols_[potential_destination] == col)) {

                    break;
                  }
                  // Make a synchronisation point don't want multiple threads accessing the same pointer simultaneously and splicing to it

                  destination_cell = world_->get_cached_square(possible_rows_[potential_destination], possible_cols_[potential_destination]);
                  destination_cell->agents_.push_back((*iter));
                  // Essentially we have to kill an individual here as it is gone off to be merged
                  (*iter).dies();
                  break;
                }
              }
            }
          }
          unsigned total = 0;
          for (unsigned i = 0; i < store_infor.destination_of_agents_moved_.size(); ++i) {
            for (unsigned j = 0; j < store_infor.destination_of_agents_moved_[i].size(); ++j )
              total += store_infor.destination_of_agents_moved_[i][j];
          }
          store_infor.initial_numbers_ = counter;
          LOG_FINEST() << "individuals at teh beginning = " << counter << " but we moved " << counter_jump << " total stored = " << total;
          if (model_->state() != State::kInitialise) {
              moved_agents_by_year_.push_back(store_infor);
          }
        }
      }
    }
  } else if (movement_type_ == MovementType::kNatal_homing) {
    LOG_FINE() << "applying Natal Homeing movement";
    // Iterate over origin cells
    // #pragma omp parallel for collapse(2) // I am not 100% confident this can be threaded
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        LOG_FINEST() << "row = " << row << " col = " << col  << " thread id = " << omp_get_thread_num();
        WorldCell* origin_cell = world_->get_base_square(row, col);
        WorldCell* destination_cell;
        if (origin_cell->is_enabled()) {
          unsigned current_cell = 0;
          for (; current_cell < origin_rows_.size(); ++current_cell) {
            if ((origin_rows_[current_cell] == row) && (origin_cols_[current_cell] == col))
              break;
          }
          MovementData store_infor(model_->get_height(), model_->get_width(), origin_cell_[current_cell], model_->current_year());
          LOG_FINEST() << "number of 'possible' agents in this cell = " << origin_cell->agents_.size() << " in current cell " << row + 1 << " col = " << col + 1;
          //
          float temp_sum;
          unsigned counter = 0;
          unsigned counter_jump = 0;
          unsigned origin_element;
          // iterate over partition
          if (selectivity_length_based_) {
            LOG_FINE() << "selectivity is length based";
            for (auto iter = origin_cell->agents_.begin(); iter != origin_cell->agents_.end(); ++iter) {
              if ((*iter).is_alive()) {
                if (rng.chance() <= selectivity_[(*iter).get_sex()]->GetResult((*iter).get_length_bin_index())) {
                  ++counter;
                  origin_element = 0;
                  // Find the home origin cell matrix
                  for (; origin_element < origin_rows_.size(); ++origin_element) {
                    if ((origin_rows_[origin_element] == (*iter).get_home_row()) && (origin_cols_[origin_element] == (*iter).get_home_col()))
                      break;
                  }
                  //LOG_FINEST() << "original element = " << origin_element << " corresponds to cell = " << origin_cell_[origin_element]<<  " row = " << row << " home row = "<< (*iter).get_home_row() << " row to use = " << origin_rows_[origin_element] << " col = " << col << " home col = " << (*iter).get_home_col() << " col to use = " << origin_cols_[origin_element];

                  // Iterate over possible cells compare to chance()
                  temp_sum = 0;
                  // Calcualte a multinomial via the following algorithm
                  // compare a random standard uniform value to the cumulative sums of the probabilities and return the first index for which
                  // the cumulative sum is greater than the random uniform
                  for (unsigned potential_destination = 0; potential_destination < possible_rows_.size(); ++potential_destination) {
                    temp_sum += probability_layers_[origin_element]->get_value(possible_rows_[potential_destination], possible_cols_[potential_destination]);
                    LOG_FINEST() << "temp_sum = " << temp_sum;

                    if (temp_sum > random_numbers_[cell_offset_[row][col] + counter]) {
                      ++counter_jump;
                      LOG_FINEST() << counter  << " cum prob = " << temp_sum << " random = " << random_numbers_[cell_offset_[row][col] + counter] << " current row = " << row << " current col = " << col << " destination row = " << possible_rows_[potential_destination] << " destination col = " << possible_cols_[potential_destination] << " origin cell = " << origin_cell_[origin_element];
                      store_infor.destination_of_agents_moved_[possible_rows_[potential_destination]][possible_cols_[potential_destination]]++;
                      // if we are moving to this cell lets not just move to the next agent
                      if ((possible_rows_[potential_destination] == row) && (possible_cols_[potential_destination] == col)) {
                        break;
                      }
                      destination_cell = world_->get_cached_square(possible_rows_[potential_destination], possible_cols_[potential_destination]);
                      destination_cell->agents_.push_back((*iter));
                      // Essentially we have to kill an individual here as it is gone off to be merged
                      (*iter).dies();
                      break;
                    }
                  }
                }
              }
            }
          } else {
            LOG_FINE() << "selectivity is Age based";
            for (auto iter = origin_cell->agents_.begin(); iter != origin_cell->agents_.end(); ++iter) {
              if ((*iter).is_alive()) {
                if (rng.chance() <= selectivity_[(*iter).get_sex()]->GetResult((*iter).get_age_index())) {
                  ++counter;
                  origin_element = 0;
                  // Find the home origin cell matrix
                  for (; origin_element < origin_rows_.size(); ++origin_element) {
                    if ((origin_rows_[origin_element] == (*iter).get_home_row()) && (origin_cols_[origin_element] == (*iter).get_home_col()))
                      break;
                  }
                  //LOG_FINEST() << "original element = " << origin_element << " corresponds to cell = " << origin_cell_[origin_element]<<  " row = " << row << " home row = "<< (*iter).get_home_row() << " row to use = " << origin_rows_[origin_element] << " col = " << col << " home col = " << (*iter).get_home_col() << " col to use = " << origin_cols_[origin_element];

                  // Iterate over possible cells compare to chance()
                  temp_sum = 0;
                  // Calcualte a multinomial via the following algorithm
                  // compare a random standard uniform value to the cumulative sums of the probabilities and return the first index for which
                  // the cumulative sum is greater than the random uniform
                  for (unsigned potential_destination = 0; potential_destination < possible_rows_.size(); ++potential_destination) {
                    temp_sum += probability_layers_[origin_element]->get_value(possible_rows_[potential_destination], possible_cols_[potential_destination]);
                    LOG_FINEST() << "temp_sum = " << temp_sum;

                    if (temp_sum > random_numbers_[cell_offset_[row][col] + counter]) {
                      ++counter_jump;
                      LOG_FINEST() << counter  << " cum prob = " << temp_sum << " random = " << random_numbers_[cell_offset_[row][col] + counter] << " current row = " << row << " current col = " << col << " destination row = " << possible_rows_[potential_destination] << " destination col = " << possible_cols_[potential_destination] << " origin cell = " << origin_cell_[origin_element];
                      store_infor.destination_of_agents_moved_[possible_rows_[potential_destination]][possible_cols_[potential_destination]]++;
                      // if we are moving to this cell lets not just move to the next agent
                      if ((possible_rows_[potential_destination] == row) && (possible_cols_[potential_destination] == col)) {
                        break;
                      }
                      destination_cell = world_->get_cached_square(possible_rows_[potential_destination], possible_cols_[potential_destination]);
                      destination_cell->agents_.push_back((*iter));
                      // Essentially we have to kill an individual here as it is gone off to be merged
                      (*iter).dies();
                      break;
                    }
                  }
                }
              }
            }
          }
          unsigned total = 0;
          for (unsigned i = 0; i < store_infor.destination_of_agents_moved_.size(); ++i) {
            for (unsigned j = 0; j < store_infor.destination_of_agents_moved_[i].size(); ++j )
              total += store_infor.destination_of_agents_moved_[i][j];
          }
          store_infor.initial_numbers_ = counter;
          LOG_FINEST() << "individuals at teh beginning = " << counter << " but we moved " << counter_jump << " total stored = " << total;
          if (model_->state() != State::kInitialise) {
            moved_agents_by_year_.push_back(store_infor);
          }
        }
      }
    }
  } // if (movement_type_ == MovementType::kNatal_homing)
  // merge destination agents into the actual grid
  world_->MergeCachedGrid(true);
  LOG_FINE();
}


// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MovementBoxTransfer::FillReportCache(ostringstream& cache) {
  LOG_TRACE();
  for (auto& values : moved_agents_by_year_) {
    cache << "year-area-" << values.year_ << "_" << values.origin_cell_ << " " << REPORT_R_LIST << "\n";
    cache << "initial_numbers_in_cell: " << values.initial_numbers_ << "\n";
    cache << "destination_values " << REPORT_R_MATRIX << "\n";
    for (unsigned i = 0; i < values.destination_of_agents_moved_.size(); ++i) {
      for (unsigned j = 0; j < values.destination_of_agents_moved_[i].size(); ++j )
        cache << values.destination_of_agents_moved_[i][j] << " ";
      cache << "\n";
    }
    cache << REPORT_R_LIST_END << "\n";
  }
}
} /* namespace processes */
} /* namespace niwa */

