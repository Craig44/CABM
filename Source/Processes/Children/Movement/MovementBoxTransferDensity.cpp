/**
 * @file MovementBoxTransferDensity.cpp
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
#include "MovementBoxTransferDensity.h"

#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/To.h"
#include "DerivedQuantities/Manager.h"
#include "InitialisationPhases/Manager.h"

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
MovementBoxTransferDensity::MovementBoxTransferDensity(Model* model) : Movement(model) {
  process_type_ = ProcessType::kTransition;
  parameters_.Bind<string>(PARAM_ORIGIN_CELL, &origin_cell_, "The origin cell associated with each spatial layer (should have a one to one relationship with specified layers), format follows row-col (1-2)", "");
  parameters_.Bind<string>(PARAM_PROBABILITY_LAYERS, &probability_layer_labels_, "Spatial layers (one layer for each origin cell) describing the probability of moving from an origin cell to all other cells in the spatial domain. ", "");
  parameters_.Bind<string>(PARAM_SELECTIVITY_LABEL, &selectivity_label_, "Label for the selectivity block", "");
  parameters_.Bind<string>(PARAM_DERIVED_QUANTITY, &derived_quantity_label_, "A label for the derived quantity", "");
  parameters_.Bind<double>(PARAM_THRESHOLD, &threshold_, "The value that triggers movement", "");
  parameters_.Bind<string>(PARAM_DIRECTION, &direction_, "The direction of threshold trigger", "")->set_allowed_values({PARAM_GREATER_THAN, PARAM_LESS_THAN});
  parameters_.Bind<bool>(PARAM_SCALE_BY_INITIAL_VALUE, &scale_dq_by_initial_value_, "Scale by initial value", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "years to apply movement", "", true);


}

void MovementBoxTransferDensity::DoValidate() {
  LOG_TRACE();
  if (origin_cell_.size() != probability_layer_labels_.size())
    LOG_ERROR_P(PARAM_PROBABILITY_LAYERS) << "you must specify a layer label for each origin cell. You have supplied '" << origin_cell_.size() << "' origin cells but '" << probability_layer_labels_.size() << "' probability layer labels, please sort this out.";
  if (!parameters_.Get(PARAM_YEARS)->has_been_defined())
    years_ = model_->years();
}

/**
 *  Build relationships with other classes
 */
void MovementBoxTransferDensity::DoBuild() {
  LOG_FINE();
  
  derived_quantity_ = model_->managers().derived_quantity()->GetDerivedQuantity(derived_quantity_label_);
  if (!derived_quantity_)
    LOG_FATAL_P(PARAM_DERIVED_QUANTITY) << "could not find the @derived_quantity block " << derived_quantity_label_ << ", please make sure it exists";

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

void MovementBoxTransferDensity::ApplyStochasticMovement(vector<Agent>& agents, MovementData& store_infor, bool tagged_partition, unsigned& origin_element, unsigned& row, unsigned& col) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  WorldCell* destination_cell;
  float temp_sum = 0.0;
  float random_chance = 0.0;
  if (selectivity_length_based_) {
  // iterate over all agents
    for (auto iter = agents.begin(); iter != agents.end(); ++iter) {
      if ((*iter).is_alive()) {
          store_infor.initial_numbers_+= (*iter).get_scalar();
        if (rng.chance() <= selectivity_[(*iter).get_sex()]->GetResult((*iter).get_length_bin_index())) {
          // Iterate over possible cells compare to chance()
          temp_sum = 0;
          random_chance = rng.chance();
          // Calcualte a multinomial via the following algorithm
          // compare a random standard uniform value to the cumulative sums of the probabilities and return the first index for which
          // the cumulative sum is greater than the random uniform
          for (unsigned potential_destination = 0; potential_destination < possible_rows_.size(); ++potential_destination) {
            temp_sum += probability_layers_[origin_element]->get_value(possible_rows_[potential_destination], possible_cols_[potential_destination]);
            if (temp_sum > random_chance) {
              //LOG_FINEST() << counter << " removals = " << counter_junp <<   " iter distance = " << distance(origin_cell->agents_.begin(), iter) << " cum prob = " << temp_sum << " random = " << random << " current row = " << row << " current col = " << col << "destination row = " << possible_rows_[potential_destination] << " destination col = " << possible_cols_[potential_destination];
              store_infor.destination_of_agents_moved_[possible_rows_[potential_destination]][possible_cols_[potential_destination]]+= (*iter).get_scalar();
              // if we are moving to this cell lets not move in memory
              if ((possible_rows_[potential_destination] == row) && (possible_cols_[potential_destination] == col)) {
                break;
              }
              // Make a synchronisation point don't want multiple threads accessing the same pointer simultaneously and splicing to it
              destination_cell = world_->get_cached_square(possible_rows_[potential_destination], possible_cols_[potential_destination]);
              if (tagged_partition) {
                destination_cell->tagged_agents_.push_back((*iter));
              } else {
                destination_cell->agents_.push_back((*iter));
              }
              (*iter).dies();
              break;
            }
          }
        }
      }
    }
  } else {
    for (auto iter = agents.begin(); iter != agents.end(); ++iter) {
      if ((*iter).is_alive()) {
        store_infor.initial_numbers_+= (*iter).get_scalar();
        if (rng.chance() <= selectivity_[(*iter).get_sex()]->GetResult((*iter).get_age_index())) {
          // Iterate over possible cells compare to chance()
          temp_sum = 0;
          random_chance = rng.chance();
          // Calcualte a multinomial via the following algorithm
          // compare a random standard uniform value to the cumulative sums of the probabilities and return the first index for which
          // the cumulative sum is greater than the random uniform
          for (unsigned potential_destination = 0; potential_destination < possible_rows_.size(); ++potential_destination) {
            temp_sum += probability_layers_[origin_element]->get_value(possible_rows_[potential_destination], possible_cols_[potential_destination]);
            if (temp_sum > random_chance) {
              //LOG_FINEST() << counter << " removals = " << counter_junp <<   " iter distance = " << distance(origin_cell->agents_.begin(), iter) << " cum prob = " << temp_sum << " random = " << random << " current row = " << row << " current col = " << col << "destination row = " << possible_rows_[potential_destination] << " destination col = " << possible_cols_[potential_destination];
              store_infor.destination_of_agents_moved_[possible_rows_[potential_destination]][possible_cols_[potential_destination]]+= (*iter).get_scalar();
              // if we are moving to this cell lets not move in memory
              if ((possible_rows_[potential_destination] == row) && (possible_cols_[potential_destination] == col)) {
                break;
              }

              destination_cell = world_->get_cached_square(possible_rows_[potential_destination], possible_cols_[potential_destination]);
              if (tagged_partition) {
                destination_cell->tagged_agents_.push_back((*iter));
              } else {
                destination_cell->agents_.push_back((*iter));
              }
              (*iter).dies();
              break;
            }
          }
        }
      }
    }
  }
  
}

/**
 *  Execute process
 */
void MovementBoxTransferDensity::DoExecute() {
  LOG_MEDIUM() << "label: " << label_;
  if (model_->state() != State::kInitialise) {
    if(find(years_.begin(), years_.end(), model_->current_year()) != years_.end()) {
      double dq_value = derived_quantity_->GetValue(model_->current_year());
      initialisationphases::Manager& init_phase_manager = *model_->managers().initialisation_phase();
      double init_dq_value = derived_quantity_->GetLastValueFromInitialisation(init_phase_manager.last_executed_phase());
      if(scale_dq_by_initial_value_)
        dq_value /= init_dq_value;
      bool apply_movement = false;
      if(direction_ == PARAM_GREATER_THAN) {
        if(dq_value > threshold_)
          apply_movement = true;
      } else {
        if(dq_value < threshold_)
          apply_movement = true;
      }
      dq_by_year_[model_->current_year()] = dq_value;
      if(apply_movement) {
        applied_movement_[model_->current_year()] = "yes";
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            // Find the right probability layer for this combo
            WorldCell* origin_cell = world_->get_base_square(row, col);
            if (origin_cell->is_enabled()) {
              unsigned origin_element = 0;
              for (; origin_element < origin_rows_.size(); ++origin_element) {
                if ((origin_rows_[origin_element] == row) && (origin_cols_[origin_element] == col))
                  break;
              }
              MovementData store_infor(model_->get_height(), model_->get_width(), origin_cell_[origin_element], model_->current_year());
              ApplyStochasticMovement(origin_cell->agents_, store_infor, false, origin_element, row, col);
              ApplyStochasticMovement(origin_cell->tagged_agents_, store_infor, true, origin_element, row, col);

              // Save info for reporting
              if (model_->state() != State::kInitialise) {
                  moved_agents_by_year_.push_back(store_infor);
              }
            }
          }
        }
        // merge destination agents into the actual grid
        world_->MergeCachedGrid(true);
        LOG_FINE();
      } else {
        applied_movement_[model_->current_year()] = "no";
      }
    }
  }
}


// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MovementBoxTransferDensity::FillReportCache(ostringstream& cache) {
  LOG_TRACE();
  cache << "year: ";
  for(auto& value : dq_by_year_) {
    cache << value.first  << " ";
  }
  cache << "\n";
  cache << "derived_quantity: ";
  for(auto& value : dq_by_year_) {
    cache << value.second  << " ";
  }
  cache << "\n";
  cache << "applied_movement: ";
  for(auto& value : applied_movement_) {
    cache << value.second  << " ";
  }
  cache << "\n";

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
