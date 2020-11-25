/**
 * @file MovementPreference.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 26/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 *
 */

// headers
#include "MovementPreference.h"

#include "Layers/Manager.h"
#include "PreferenceFunctions/Manager.h"

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
MovementPreference::MovementPreference(Model* model) : Movement(model) {
  process_type_ = ProcessType::kTransition;
  //parameters_.Bind<string>(PARAM_SELECTIVITY, &selectivity_label_, "Selectivity label", "");
  parameters_.Bind<double>(PARAM_D_MAX, &d_max_, "The maxiumu diffusion value", "")->set_lower_bound(0.0, false);
  parameters_.Bind<double>(PARAM_TIME_INTERVAL, &time_interval_, "The time interval that this movement process occurs over (k) in the documentation formula's", "", 1.0)->set_lower_bound(0.0, false);
  parameters_.Bind<double>(PARAM_ZETA, &zeta_, "An arbiituary parameter that controls curvature between preference and diffusion", "", 1.0)->set_lower_bound(0.0, false);
  parameters_.Bind<bool>(PARAM_BROWNIAN_MOTION, &brownian_motion_, "Just apply random movement, undirected movement", "", false);
  parameters_.Bind<string>(PARAM_PREFERENCE_FUNCTIONS, &preference_function_labels_, "The preference functions to apply", "", true);
  parameters_.Bind<string>(PARAM_PREFERENCE_LAYERS, &preference_layer_labels_, "The preference functions to apply", "", true);
  parameters_.Bind<string>(PARAM_SELECTIVITY_LABEL, &selectivity_label_, "Label for the selectivity block", "");
  //parameters_.Bind<unsigned>(PARAM_NUMBER_OF_CELLS_FOR_GRADIENT, &number_of_cells_gradient_, "Number of cells to calculate the (E-W & N-S) gradient over", "", true);

}

/*
 * DoValidate
*/
void MovementPreference::DoValidate() {
  LOG_TRACE();
  // number_of_cells_gradient_ should be an odd number
  if(number_of_cells_gradient_ <= 0)
    LOG_ERROR_P(PARAM_NUMBER_OF_CELLS_FOR_GRADIENT) << "Needs to be bigger than zero";
  if(number_of_cells_gradient_ > model_->get_height())
    LOG_ERROR_P(PARAM_NUMBER_OF_CELLS_FOR_GRADIENT) << "Needs to be less than number of rows in model";
  if(number_of_cells_gradient_ > model_->get_width())
    LOG_ERROR_P(PARAM_NUMBER_OF_CELLS_FOR_GRADIENT) << "Needs to be less than number of cols in model";
  if (number_of_cells_gradient_ % 2 == 0)
    LOG_ERROR_P(PARAM_NUMBER_OF_CELLS_FOR_GRADIENT) << "is even, needs to be odd";

  cell_offset_gradient_ = (number_of_cells_gradient_ - 1) / 2;
  LOG_MEDIUM() << "cell_offset_gradient = " << cell_offset_gradient_ << " number_of_cells_gradient " << number_of_cells_gradient_;
}

/**
 *  Build relationships with other classes
 */
void MovementPreference::DoBuild() {
  LOG_MEDIUM();
  // Get the preference functions
  for (auto& label : preference_function_labels_) {
    PreferenceFunction* temp_func = nullptr;
    temp_func = model_->managers().preference_function()->GetPreferenceFunction(label);
    if (!temp_func) {
      LOG_FATAL_P(PARAM_PREFERENCE_FUNCTIONS) << "could not find the preference function '" << label << "', please make sure it exists'";
    }
    preference_functions_.push_back(temp_func);
  }
  // Get the layers
  unsigned layer_ndx = 0;
  for (auto& label : preference_layer_labels_) {
    layers::NumericLayer* temp_layer = nullptr;
    temp_layer = model_->managers().layer()->GetNumericLayer(label);
    if (!temp_layer) {
      LOG_FATAL_P(PARAM_PREFERENCE_LAYERS) << "could not find the layer '" << label << "', please make sure it exists, and if it does exist make sure it is of type 'numeric''";
    }

    if (not temp_layer->is_static()) {
      non_static_layer_ndx_.push_back(layer_ndx);
      calculate_on_the_fly_ = true;
    }
    LOG_MEDIUM() << "brownian motion = " << brownian_motion_ << " calculate on the fly? " << calculate_on_the_fly_;

    preference_layers_.push_back(temp_layer);
    ++layer_ndx;
  }

  calculate_gradients();

  // check that lat and long layers have been supplied
  if (!model_->lat_and_long_supplied()) {
    LOG_ERROR_P(PARAM_LABEL) << "in order to apply the preference movement function, you need to specify lat and long bounds on the @model block";
  }

  cell_offset_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cell_offset_[i].resize(model_->get_width());
  }


  LOG_MEDIUM() << "selectivities supplied = " << selectivity_label_.size();
  // Build selectivity links
  if (selectivity_label_.size() == 1)
    selectivity_label_.assign(2, selectivity_label_[0]);

  if (selectivity_label_.size() > 2) {
    LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << "You suppled " << selectivity_label_.size()  << " Selectiviites, you can only have one for each sex max = 2";
  }
  LOG_MEDIUM() << "selectivities supplied = " << selectivity_label_.size();

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
void MovementPreference::DoExecute() {
  LOG_MEDIUM() << label_;
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  // Pre-calculate agents in the world to set aside our random numbers needed for the operation
  n_agents_ = 0;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        cell_offset_[row][col] = n_agents_;
        LOG_FINEST() << "row = " << row << " col = " << col << " offset = " << n_agents_ << " n-agents in this cell = " << cell->agents_.size();
        n_agents_ += cell->agents_.size();
      }
    }
  }
  // Allocate a single block of memory rather than each thread temporarily allocating their own memory.
  lon_random_numbers_.resize(n_agents_ + 1);
  lat_random_numbers_.resize(n_agents_ + 1);
  for (unsigned i = 0; i <= n_agents_; ++i) {
    lat_random_numbers_[i] = rng.normal();
    lon_random_numbers_[i] = rng.normal();
  }
  LOG_FINEST() << "random numbers generatored = " << lat_random_numbers_.size();
  unsigned current_year = model_->current_year();

  if (calculate_on_the_fly_ && not brownian_motion_) {
    LOG_MEDIUM() << "Need to do preference function calculation on the fly";
    // need to do some extra calculation before moveing through the main algorithm
    // calculate preference value
    for (unsigned layer_ndx = 0; layer_ndx < non_static_layer_ndx_.size(); ++layer_ndx) {
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          if (model_->state() == State::kInitialise)
            initialisation_preference_value_[row][col] *= preference_functions_[layer_ndx]->get_result(preference_layers_[layer_ndx]->get_value(row, col));  // the 0 in the layers call will return the default layer
          else
            preference_by_year_[current_year][row][col] *= preference_functions_[layer_ndx]->get_result(preference_layers_[layer_ndx]->get_value(row, col));  // the 0 in the layers call will return the default layer
        }
      }
    }
    // take power
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        if (model_->state() == State::kInitialise)
          initialisation_preference_value_[row][col] = pow(initialisation_preference_value_[row][col], float(1.0 / preference_functions_.size()));
        else
          preference_by_year_[current_year][row][col] = pow(preference_by_year_[current_year][row][col], float(1.0 / preference_functions_.size()));
      }
    }
    // calculate gradient
    unsigned offset;
    if (model_->state() == State::kInitialise) {
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          // calcualte meridional gradient
          if (row < cell_offset_gradient_) {
            // at the edge
            offset = cell_offset_gradient_ - (row);
            initialisation_meridonal_gradient_[row][col] = (initialisation_preference_value_[0][col] - initialisation_preference_value_[number_of_cells_gradient_ - offset - 1][col]) / (row + 2);
            //std::cerr << "row = " << row << " col = " << col << " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[0][col] << " x2 = " << initialisation_preference_value_[number_of_cells_gradient_ - offset - 1][col] << " denominator " <<(row + 2) <<"\n";
          } else if (row >= (model_->get_height() - cell_offset_gradient_)) {
            offset = (model_->get_height() - (row));
            initialisation_meridonal_gradient_[row][col] = (initialisation_preference_value_[model_->get_height()  - (2 +  offset)][col] - initialisation_preference_value_[model_->get_height() - 1][col]) / (offset + 1);
            //std::cerr << "row = " << row << " col = " << col  << " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[model_->get_height()  - (2 +  offset)][col]  << " x2 = " << initialisation_preference_value_[model_->get_height() - 1][col] <<"\n";
          } else {
            initialisation_meridonal_gradient_[row][col] = (initialisation_preference_value_[row - cell_offset_gradient_][col] - initialisation_preference_value_[row + cell_offset_gradient_][col]) / (number_of_cells_gradient_ - 1);
            //std::cerr << "row = " << row << " col = " << col <<  " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[row - cell_offset_gradient_][col]  << " x2 = " << initialisation_preference_value_[row + cell_offset_gradient_][col] <<"\n";
          }
          // calcualte zonal gradient
          if (col < cell_offset_gradient_) {
            offset = cell_offset_gradient_ - (col);
            initialisation_zonal_gradient_[row][col] = (initialisation_preference_value_[row][number_of_cells_gradient_ - offset - 1] - initialisation_preference_value_[row][0]) / (col + 2);
          } else if (col >= (model_->get_width() - cell_offset_gradient_)) {
            offset = (model_->get_width() - (col));
            initialisation_zonal_gradient_[row][col] = (initialisation_preference_value_[row][model_->get_width()  - 1] - initialisation_preference_value_[row][model_->get_width()  - (2 +  offset)]) / (offset + 1);
          } else {
            initialisation_zonal_gradient_[row][col] = (initialisation_preference_value_[row][col + cell_offset_gradient_] - initialisation_preference_value_[row][col - cell_offset_gradient_]) / (float)(number_of_cells_gradient_ - 1);
          }
        }
      }
    } else {
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          if (row < cell_offset_gradient_) {
            // at the edge
            offset = cell_offset_gradient_ - (row);
            meridonal_gradient_[current_year][row][col] = (preference_by_year_[current_year][0][col] - preference_by_year_[current_year][number_of_cells_gradient_ - offset - 1][col]) / (row + 2);
            //std::cerr << "row = " << row << " col = " << col << " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[0][col] << " x2 = " << initialisation_preference_value_[number_of_cells_gradient_ - offset - 1][col] << " denominator " <<(row + 2) <<"\n";
          } else if (row >= (model_->get_height() - cell_offset_gradient_)) {
            offset = (model_->get_height() - (row));
            meridonal_gradient_[current_year][row][col] = (preference_by_year_[current_year][model_->get_height()  - (2 +  offset)][col] - preference_by_year_[current_year][model_->get_height() - 1][col]) / (offset + 1);
            //std::cerr << "row = " << row << " col = " << col  << " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[model_->get_height()  - (2 +  offset)][col]  << " x2 = " << initialisation_preference_value_[model_->get_height() - 1][col] <<"\n";
          } else {
            meridonal_gradient_[current_year][row][col] = (preference_by_year_[current_year][row - cell_offset_gradient_][col] - preference_by_year_[current_year][row + cell_offset_gradient_][col]) / (number_of_cells_gradient_ - 1);
            //std::cerr << "row = " << row << " col = " << col <<  " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[row - cell_offset_gradient_][col]  << " x2 = " << initialisation_preference_value_[row + cell_offset_gradient_][col] <<"\n";
          }
          // calcualte zonal gradient
          if (col < cell_offset_gradient_) {
            offset = cell_offset_gradient_ - (col);
            zonal_gradient_[current_year][row][col] = (preference_by_year_[current_year][row][number_of_cells_gradient_ - offset - 1] - preference_by_year_[current_year][row][0]) / (col + 2);
          } else if (col >= (model_->get_width() - cell_offset_gradient_)) {
            offset = (model_->get_width() - (col));
            zonal_gradient_[current_year][row][col] = (preference_by_year_[current_year][row][model_->get_width()  - 1] - preference_by_year_[current_year][row][model_->get_width()  - (2 +  offset)]) / (offset + 1);
          } else {
            zonal_gradient_[current_year][row][col] = (preference_by_year_[current_year][row][col + cell_offset_gradient_] - preference_by_year_[current_year][row][col - cell_offset_gradient_]) / (float)(number_of_cells_gradient_ - 1);
          }
        }
      }
    }
  }
  /*
   * Main movement algorithm
   */
  //Model and world information that we need to build and store before we run and execute in threaded mode
  // - cell pointers in matrix
  // - height and width in matrix
  // - current_year/Initialisaiton in a matrix form (this is the only one that is dynamic godamit)
  // - Max-lat/lon in a matrix

  // Iterate over origin cells
  //#pragma omp parallel for collapse(2)
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      LOG_FINE() << "row = " << row << " col = " << col  << " thread id = " << omp_get_thread_num();
      WorldCell* origin_cell = world_->get_base_square(row, col);
      WorldCell* destination_cell;
      if (origin_cell->is_enabled()) {
        LOG_FINEST() << "number of agents in this cell = " << origin_cell->agents_.size();
        MovementData store_infor(model_->get_height(), model_->get_width(), origin_cell->get_cell_label(), model_->current_year());
        float u, v, lat_distance, lon_distance;
        unsigned destination_row, destination_col;
        // get gradients for current cells
        if (brownian_motion_) {
          u = 0;
          v = 0;
          standard_deviation_ = sqrt(2 * d_max_ * 1);  // TODO if you apply this many times, will need to change the temporal multiplier
        } else if (model_->state() == State::kInitialise) {
          v = initialisation_meridonal_gradient_[row][col];
          u = initialisation_zonal_gradient_[row][col];
          #pragma omp critical
          {
            calculate_diffusion_parameter(initialisation_preference_value_[row][col], standard_deviation_);
          }
        } else  {
          v = meridonal_gradient_[model_->current_year()][row][col];
          u = zonal_gradient_[model_->current_year()][row][col];
          #pragma omp critical
          {
            calculate_diffusion_parameter(preference_by_year_[model_->current_year()][row][col], standard_deviation_);
          }
        }

        LOG_FINE() << "n agents in this cell = " <<  origin_cell->agents_.size();
        store_infor.initial_numbers_ = origin_cell->agents_.size();

        float total_v = 0;
        float total_u = 0;
        float total_jumps = 0;
        float update_lat = 0.0;
        float update_lon = 0.0;

        unsigned counter = 0;
        if (selectivity_length_based_) {
          LOG_FINE() << "Length based preference movement";
          for (auto iter = origin_cell->agents_.begin(); iter != origin_cell->agents_.end(); ++counter, ++iter) {
            // Iterate over possible cells compare to chance()
            if ((*iter).is_alive()) {
              if (rng.chance() <= selectivity_[(*iter).get_sex()]->GetResult((*iter).get_length_bin_index())) {
                lat_distance = v + lat_random_numbers_[cell_offset_[row][col] + counter] * standard_deviation_;
                lon_distance = u + lon_random_numbers_[cell_offset_[row][col] + counter] * standard_deviation_;
                total_v += lat_distance;
                total_u += lon_distance;
                ++total_jumps;
                LOG_FINEST() << counter<< " " << (*iter).get_lat() << " distance = " << lat_distance << " lon = " << (*iter).get_lon() << " distance = " << lon_distance << " Z = " << lon_random_numbers_[cell_offset_[row][col] + counter] << " sigma = " << standard_deviation_;
                // Check bounds and find cell destination
                // this is crude and should be changed into the future, I am thinking some sort of buffer from the edge
                update_lat = (*iter).get_lat() + lat_distance;
                update_lon = (*iter).get_lon() + lon_distance;

                if ((update_lat <= model_->max_lat()) && (update_lat >= model_->min_lat())) {
                  (*iter).set_lat(update_lat);
                } // else they stay as it would be jumping out of bounds

                if ((update_lon <= model_->max_lon()) && (update_lon >= model_->min_lon())) {
                  (*iter).set_lon(update_lon);
                } // else they stay as it would be jumping out of bounds
                (*iter).save_lon_hist(update_lon);
                (*iter).save_lat_hist(update_lat);

                world_->get_cell_element(destination_row, destination_col, (*iter).get_lat(), (*iter).get_lon()); // very difficult to thread this...

                LOG_FINEST() << (*iter).get_lat() << " " << (*iter).get_lon() << " " << destination_row << " " << destination_col << " " << row << " " << col;

                if (destination_row == row && destination_col == col) {
                  store_infor.destination_of_agents_moved_[destination_row][destination_col]++;
                } else {
                  destination_cell = world_->get_cached_square(destination_row, destination_col);
                  if (destination_cell->is_enabled()) {
                    // We are moving 'splice' this agent to the destination cache cell
                    store_infor.destination_of_agents_moved_[destination_row][destination_col]++;
                    //#pragma omp critical
                    {
                      destination_cell->agents_.push_back(*iter);
                    }
                    (*iter).dies();
                  }
                }
              }
            }
          }
        } else {
          LOG_FINE() << "Age based preference movement";
          for (auto iter = origin_cell->agents_.begin(); iter != origin_cell->agents_.end(); ++counter, ++iter) {
            // Iterate over possible cells compare to chance()
            if ((*iter).is_alive()) {
              if (rng.chance() <= selectivity_[(*iter).get_sex()]->GetResult((*iter).get_age_index())) {

                lat_distance = v + lat_random_numbers_[cell_offset_[row][col] + counter] * standard_deviation_;
                lon_distance = u + lon_random_numbers_[cell_offset_[row][col] + counter] * standard_deviation_;
                total_v += lat_distance;
                total_u += lon_distance;
                ++total_jumps;
                update_lat = (*iter).get_lat() + lat_distance;
                update_lon = (*iter).get_lon() + lon_distance;

                if ((update_lat <= model_->max_lat()) && (update_lat >= model_->min_lat())) {
                  (*iter).set_lat(update_lat);
                } // else they stay as it would be jumping out of bounds

                if ((update_lon <= model_->max_lon()) && (update_lon >= model_->min_lon())) {
                  (*iter).set_lon(update_lon);
                } // else they stay as it would be jumping out of bounds
                (*iter).save_lon_hist(update_lon);
                (*iter).save_lat_hist(update_lat);

                world_->get_cell_element(destination_row, destination_col, (*iter).get_lat(), (*iter).get_lon()); // very difficult to thread this...

                LOG_FINEST() << (*iter).get_lat() << " " << (*iter).get_lon() << " " << destination_row << " " << destination_col << " " << row << " " << col;

                if (destination_row == row && destination_col == col) {
                  store_infor.destination_of_agents_moved_[destination_row][destination_col]++;
                } else {
                  destination_cell = world_->get_cached_square(destination_row, destination_col);
                  if (destination_cell->is_enabled()) {
                    destination_cell->add_agent_alive((*iter).get_scalar());
                    origin_cell->remove_agent_alive((*iter).get_scalar());
                    // We are moving 'splice' this agent to the destination cache cell
                    store_infor.destination_of_agents_moved_[destination_row][destination_col]++;
                    //#pragma omp critical
                    {
                      destination_cell->agents_.push_back(*iter);
                    }
                    (*iter).dies();
                  }
                }
              }
            }
          }
        }
        if (model_->state() != State::kInitialise) {
          if (total_jumps != 0) {
            average_zonal_jump_[current_year][row][col] = total_u / total_jumps;
            average_meridonal_jump_[current_year][row][col] = total_v / total_jumps;
          }
        }
        LOG_FINE() << "average latitude jump = " << total_v / total_jumps << " row = " << row + 1 << " col = " << col + 1 << " total jumps = " << total_jumps << " total lat = " << total_v;
        LOG_FINE() << "average longitude jump = " << total_u / total_jumps << " row = " << row + 1 << " col = " << col + 1;

/*        if (model_->state() != State::kInitialise) {
          #pragma omp critical
          {
            moved_agents_by_year_.push_back(store_infor);
          }
        }*/

      } // is enabled
    } // col
  } // row
  // merge destination agents into the actual grid
  LOG_MEDIUM() << "Mergind cached world";
  world_->MergeCachedGrid();
  LOG_MEDIUM();
}

/*
 * Calculate gradient using the difference method, for each year and initialisation
*/
void  MovementPreference::calculate_gradients() {
  LOG_MEDIUM() << "calculate_gradients()";
  if (brownian_motion_)
    return void();
  LOG_FINE();
  // Start with setting the initialisation sets
  initialisation_meridonal_gradient_.resize(model_->get_height());
  initialisation_zonal_gradient_.resize(model_->get_height());
  vector<vector<double>> preference(model_->get_height());
  initialisation_preference_value_.resize(model_->get_height());
  vector<vector<double>> average_zonal(model_->get_height());
  vector<vector<double>> average_meridonal(model_->get_height());
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    initialisation_meridonal_gradient_[row].resize(model_->get_width());
    initialisation_zonal_gradient_[row].resize(model_->get_width(), 1.0);
    preference[row].resize(model_->get_width(), 1.0);
    initialisation_preference_value_[row].resize(model_->get_width(), 1.0);
    average_meridonal[row].resize(model_->get_width(), 0.0);
    average_zonal[row].resize(model_->get_width(),0.0);
  }

  for (auto year : model_->years()) {
    average_zonal_jump_[year] = average_zonal;
    average_meridonal_jump_[year] = average_meridonal;
  }

  // Do iniitalisation
  for (unsigned layer_iter = 0; layer_iter < preference_layer_labels_.size(); ++layer_iter) {
    if (!preference_layers_[layer_iter]->is_static())
      continue; // does this layer need to be called during DoExecute
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        preference[row][col] *= preference_functions_[layer_iter]->get_result(preference_layers_[layer_iter]->get_value(row, col, 0));  // the 0 in the layers call will return the default layer
        LOG_FINE() << "preference layer value = " <<  preference_layers_[layer_iter]->get_value(row, col, 0) << " preference function value = " << preference_functions_[layer_iter]->get_result(preference_layers_[layer_iter]->get_value(row, col, 0));
      }
    }
  }
  if (not calculate_on_the_fly_) {
    // take the preference to the power
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        initialisation_preference_value_[row][col] = pow(preference[row][col], float(1.0 / preference_function_labels_.size()));
      }
    }
    unsigned offset;
    // calculate gradient
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        // calcualte meridional gradient
        if (row < cell_offset_gradient_) {
          // at the edge
          offset = cell_offset_gradient_ - (row);
          initialisation_meridonal_gradient_[row][col] = (initialisation_preference_value_[0][col] - initialisation_preference_value_[number_of_cells_gradient_ - offset - 1][col]) / (row + 2);
          //std::cerr << "row = " << row << " col = " << col << " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[0][col] << " x2 = " << initialisation_preference_value_[number_of_cells_gradient_ - offset - 1][col] << " denominator " <<(row + 2) <<"\n";
        } else if (row >= (model_->get_height() - cell_offset_gradient_)) {
          offset = (model_->get_height() - (row));
          initialisation_meridonal_gradient_[row][col] = (initialisation_preference_value_[model_->get_height()  - (2 +  offset)][col] - initialisation_preference_value_[model_->get_height() - 1][col]) / (offset + 1);
          //std::cerr << "row = " << row << " col = " << col  << " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[model_->get_height()  - (2 +  offset)][col]  << " x2 = " << initialisation_preference_value_[model_->get_height() - 1][col] <<"\n";
        } else {
          initialisation_meridonal_gradient_[row][col] = (initialisation_preference_value_[row - cell_offset_gradient_][col] - initialisation_preference_value_[row + cell_offset_gradient_][col]) / (number_of_cells_gradient_ - 1);
          //std::cerr << "row = " << row << " col = " << col <<  " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[row - cell_offset_gradient_][col]  << " x2 = " << initialisation_preference_value_[row + cell_offset_gradient_][col] <<"\n";
        }
        // calcualte zonal gradient
        if (col < cell_offset_gradient_) {
          offset = cell_offset_gradient_ - (col);
          initialisation_zonal_gradient_[row][col] = (initialisation_preference_value_[row][number_of_cells_gradient_ - offset - 1] - initialisation_preference_value_[row][0]) / (col + 2);
        } else if (col >= (model_->get_width() - cell_offset_gradient_)) {
          offset = (model_->get_width() - (col));
          initialisation_zonal_gradient_[row][col] = (initialisation_preference_value_[row][model_->get_width()  - 1] - initialisation_preference_value_[row][model_->get_width()  - (2 +  offset)]) / (offset + 1);
        } else {
          initialisation_zonal_gradient_[row][col] = (initialisation_preference_value_[row][col + cell_offset_gradient_] - initialisation_preference_value_[row][col - cell_offset_gradient_]) / (float)(number_of_cells_gradient_ - 1);
        }
      }
    }
  } else {
    initialisation_preference_value_ = preference;
  }

  if (calculate_on_the_fly_) {
    // Some prefernce attributes need to be calculated on the fly, usually associated with biomass or density preferences
    // For this we will calculate the preference for all the static layers, to reduce run-time calculations

    // We will allocate memory for all containers used at runtime here
    // Now do it for all years
    for (auto year : model_->years()) {
      LOG_FINE() << "year = " << year;
      // initialise temporary containers
      vector<vector<double>> zonal_gradient(model_->get_height());
      vector<vector<double>> meredional_gradient(model_->get_height());
      vector<vector<double>> preference(model_->get_height());
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        zonal_gradient[row].resize(model_->get_width(), 1.0);
        meredional_gradient[row].resize(model_->get_width(), 1.0);
        preference[row].resize(model_->get_width(),1.0);
      }

      for (unsigned layer_iter = 0; layer_iter < preference_layer_labels_.size(); ++layer_iter) {
        if (!preference_layers_[layer_iter]->is_static())
          continue; // this should not be called if all preference layers are static thus is shares the same code as on the fly calculation

        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            preference[row][col] *= preference_functions_[layer_iter]->get_result(preference_layers_[layer_iter]->get_value(row, col, year));  // the 0 in the layers call will return the default layer

          }
        }
      }
      // don't calculate the power yet do in Execute
      preference_by_year_[year] = preference;

      // store empty containers, but have allocated memory
      meridonal_gradient_[year] = meredional_gradient;
      zonal_gradient_[year] = zonal_gradient;
    }
  } else {
    // We can calculate gradients now as nothing is unknown to us at this point
    for (auto year : model_->years()) {
      LOG_FINE() << "year = " << year << " number of layer labels = " << preference_layer_labels_.size();

      // initialise temporary containers

      vector<vector<double>> zonal_gradient(model_->get_height());
      vector<vector<double>> meredional_gradient(model_->get_height());
      vector<vector<double>> preference(model_->get_height());
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        zonal_gradient[row].resize(model_->get_width(), 1.0);
        meredional_gradient[row].resize(model_->get_width(), 1.0);
        preference[row].resize(model_->get_width(), 1.0);

      }
      for (unsigned layer_iter = 0; layer_iter < preference_layer_labels_.size(); ++layer_iter) {
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            preference[row][col] *= preference_functions_[layer_iter]->get_result(preference_layers_[layer_iter]->get_value(row, col, year));  // the 0 in the layers call will return the default layer

          }
        }
      }
      // take the preference to the power
      LOG_MEDIUM() << "power = " << float(1.0 / preference_function_labels_.size());
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          double val = pow(preference[row][col], float(1.0 / preference_function_labels_.size()));
          preference[row][col] = val;
        }
      }
      preference_by_year_[year] = preference;

      unsigned offset;
      // calculate gradient
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          // calcualte meridional gradient
          if (row < cell_offset_gradient_) {
            // at the edge
            offset = cell_offset_gradient_ - (row);
            meredional_gradient[row][col] = (preference[0][col] - preference[number_of_cells_gradient_ - offset - 1][col]) / (row + 2);
            //std::cerr << "row = " << row << " col = " << col << " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[0][col] << " x2 = " << initialisation_preference_value_[number_of_cells_gradient_ - offset - 1][col] << " denominator " <<(row + 2) <<"\n";
          } else if (row >= (model_->get_height() - cell_offset_gradient_)) {
            offset = (model_->get_height() - (row));
            meredional_gradient[row][col] = (preference[model_->get_height()  - (2 +  offset)][col] - preference[model_->get_height() - 1][col]) / (offset + 1);
            //std::cerr << "row = " << row << " col = " << col  << " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[model_->get_height()  - (2 +  offset)][col]  << " x2 = " << initialisation_preference_value_[model_->get_height() - 1][col] <<"\n";
          } else {
            meredional_gradient[row][col] = (preference[row - cell_offset_gradient_][col] - preference[row + cell_offset_gradient_][col]) / (number_of_cells_gradient_ - 1);
            //std::cerr << "row = " << row << " col = " << col <<  " merid = " << initialisation_meridonal_gradient_[row][col] << " x1 = " << initialisation_preference_value_[row - cell_offset_gradient_][col]  << " x2 = " << initialisation_preference_value_[row + cell_offset_gradient_][col] <<"\n";
          }
          // calcualte zonal gradient
          if (col < cell_offset_gradient_) {
            offset = cell_offset_gradient_ - (col);
            zonal_gradient[row][col] = (preference[row][number_of_cells_gradient_ - offset - 1] - preference[row][0]) / (col + 2);
          } else if (col >= (model_->get_width() - cell_offset_gradient_)) {
            offset = (model_->get_width() - (col));
            zonal_gradient[row][col] = (preference[row][model_->get_width()  - 1] - preference[row][model_->get_width()  - (2 +  offset)]) / (offset + 1);
          } else {
            zonal_gradient[row][col] = (preference[row][col + cell_offset_gradient_] - preference[row][col - cell_offset_gradient_]) / (float)(number_of_cells_gradient_ - 1);
          }
        }
      }
      meridonal_gradient_[year] = meredional_gradient;
      zonal_gradient_[year] = zonal_gradient;
    }
  }
}

/*
 * This function just recalculates the random component of the random walk around the preference space, which depends on habitat quality
*/
void  MovementPreference::calculate_diffusion_parameter(double& preference_value, double& standard_deviation) {
  diffusion_parameter_ = d_max_ * (1 - (preference_value / (zeta_ + preference_value)));
  if (brownian_motion_)
    diffusion_parameter_ = d_max_;
  standard_deviation = sqrt(2 * diffusion_parameter_ * 1);  // TODO if you apply this many times, will need to change the temporal multiplier
}

// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MovementPreference::FillReportCache(ostringstream& cache) {
  LOG_TRACE();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  // Print Preference by year
  cache << "standard_dev: " << standard_deviation_ << "\n";
  // Print gradient by year.
  if (not brownian_motion_) {
    cache << "initialisation_preference " << REPORT_R_MATRIX << "\n";
    for (unsigned row = 0; row < initialisation_preference_value_.size(); ++row) {
      for (unsigned col = 0; col < initialisation_preference_value_[row].size(); ++col) {
        cache << initialisation_preference_value_[row][col] << " ";
      }
      cache << "\n";
    }
    cache << "initialisation_meridional " << REPORT_R_MATRIX << "\n";
    for (unsigned row = 0; row < initialisation_meridonal_gradient_.size(); ++row) {
      for (unsigned col = 0; col < initialisation_meridonal_gradient_[row].size(); ++col) {
        cache << initialisation_meridonal_gradient_[row][col] << " ";
      }
      cache << "\n";
    }

    cache << "initialisation_zonal " << REPORT_R_MATRIX << "\n";
    for (unsigned row = 0; row < initialisation_zonal_gradient_.size(); ++row) {
      for (unsigned col = 0; col < initialisation_zonal_gradient_[row].size(); ++col) {
        cache << initialisation_zonal_gradient_[row][col] << " ";
      }
      cache << "\n";
    }
    // Preference by year
    for (auto& values : preference_by_year_) {
      cache << "preference_" << values.first << " " << REPORT_R_MATRIX << "\n";
      for (unsigned i = 0; i < values.second.size(); ++i) {
        for (unsigned j = 0; j < values.second[i].size(); ++j)
          cache << values.second[i][j] << " ";
        cache << "\n";
      }
    }

    for (auto& values : zonal_gradient_) {
      cache << "zonal_" << values.first << " " << REPORT_R_MATRIX << "\n";
      for (unsigned i = 0; i < values.second.size(); ++i) {
        for (unsigned j = 0; j < values.second[i].size(); ++j)
          cache << values.second[i][j] << " ";
        cache << "\n";
      }
    }

    for (auto& values : meridonal_gradient_) {
      cache << "meridonal_" << values.first << " " << REPORT_R_MATRIX << "\n";
      for (unsigned i = 0; i < values.second.size(); ++i) {
        for (unsigned j = 0; j < values.second[i].size(); ++j)
          cache << values.second[i][j] << " ";
        cache << "\n";
      }
    }

    for (auto& values : average_zonal_jump_) {
      cache << "average_zonal_jump_" << values.first << " " << REPORT_R_MATRIX << "\n";
      for (unsigned i = 0; i < values.second.size(); ++i) {
        for (unsigned j = 0; j < values.second[i].size(); ++j)
          cache << values.second[i][j] << " ";
        cache << "\n";
      }
    }
    for (auto& values : average_meridonal_jump_) {
      cache << "average_meridional_jump_" << values.first << " " << REPORT_R_MATRIX << "\n";
      for (unsigned i = 0; i < values.second.size(); ++i) {
        for (unsigned j = 0; j < values.second[i].size(); ++j)
          cache << values.second[i][j] << " ";
        cache << "\n";
      }
    }
  }

  // Grab a couple of agents and track thier history
  int n_agents = 30;
  int agent_ndx;
  for(unsigned row = 0; row < model_->get_height(); ++row) {
    for(unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        agent_ndx = cell->agents_.size() * rng.chance();
        if( cell->agents_[agent_ndx].is_alive()) {
          cache << n_agents <<"_lat: ";
          for(auto lat : cell->agents_[agent_ndx].lat_history_)
            cache << lat << " ";
          cache << "\n" << n_agents <<"_lon: ";
          for(auto lon : cell->agents_[agent_ndx].lon_history_)
            cache << lon << " ";
          cache << "\n";
          n_agents--;
          if (n_agents < 0)
            break;
        }
      }
    }
  }




/*  for (auto& values : moved_agents_by_year_) {
    cache << "initial_numbers_in_cell: " << values.initial_numbers_ << "\n";
    cache << values.year_ << "_" << values.origin_cell_ << "_destination " << REPORT_R_MATRIX << "\n";
    for (unsigned i = 0; i < values.destination_of_agents_moved_.size(); ++i) {
      for (unsigned j = 0; j < values.destination_of_agents_moved_[i].size(); ++j )
        cache << values.destination_of_agents_moved_[i][j] << " ";
      cache << "\n";
    }
  }*/
}
} /* namespace processes */
} /* namespace niwa */

