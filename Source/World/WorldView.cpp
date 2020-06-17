/**
 * @file WorldView.cpp
 * @author C.Marsh
 * @github
 * @date 11/06/2018
 * @section LICENSE
 *
 */

// headers
#include "WorldView.h"

#include "Model/Model.h"
#include "Layers/Manager.h"
#include "Utilities/To.h"
#include "Utilities/Types.h"
#include "WorldCell.h"
#include "omp.h"
// namespaces
namespace niwa {

// Constructor
WorldView::WorldView(Model* model) :
    model_(model) {
  base_grid_ = 0;
  cached_grid_ = 0;
  init_cached_grid_ = 0;
}

// DeConstructor
WorldView::~WorldView() {
  // Clean Our DifferenceGrid
  if (cached_grid_ != 0) {
    for (unsigned i = 0; i < height_; ++i) {
      delete[] cached_grid_[i];
      cached_grid_[i] = 0;
    }
    delete[] cached_grid_;
  }

  // Clean Our Grid
  if (base_grid_ != 0) {
    for (unsigned i = 0; i < height_; ++i) {
      delete[] base_grid_[i];
      base_grid_[i] = 0;
    }
    delete[] base_grid_;
  }

  // delete initialisation grid.
  if (init_cached_grid_ != 0) {
    for (unsigned i = 0; i < height_; ++i) {
      delete[] init_cached_grid_[i];
      init_cached_grid_[i] = 0;
    }
    delete[] init_cached_grid_;
  }
}

// Build Method
void WorldView::Validate() {
  width_ = model_->get_width();
  height_ = model_->get_height();
}

// Build Method
void WorldView::Reset() {
  LOG_FINE() << "resetting world, clearing all world cells of agents, but keeping all the structure the same.";
  for (unsigned i = 0; i < height_; ++i) {
    for (unsigned j = 0; j < width_; ++j) {
      if (base_grid_[i][j].is_enabled()) {
        base_grid_[i][j].agents_.clear();
      }
    }
  }
}
// Build Method
void WorldView::Build() {
  LOG_TRACE();
  base_layer_ = model_->managers().layer()->GetNumericLayer(model_->get_base_layer_label());
  if (!base_layer_) {
    LOG_FATAL() << "The base layer '" << model_->get_base_layer_label() << "' found in the @model block could not be found, please check there is a @layer defined for this layer and that it is type = 'numeric'";
  }

  // Allocate memory for our cells
  // Allocate Space for our X (Height)
  base_grid_ = new WorldCell*[height_];
  for (unsigned i = 0; i < height_; ++i) {
    base_grid_[i] = new WorldCell[width_];
  }

  // Allocate Our Difference Grid
  cached_grid_ = new WorldCell*[height_];
  for (unsigned i = 0; i < height_; ++i) {
    cached_grid_[i] = new WorldCell[width_];
  }
  // Allocate Our Difference Grid
  init_cached_grid_ = new WorldCell*[height_];
  for (unsigned i = 0; i < height_; ++i) {
    init_cached_grid_[i] = new WorldCell[width_];
  }


  if (model_->lat_and_long_supplied()) {
    lat_bounds_ = model_->get_lat_bounds();
    lon_bounds_ = model_->get_lon_bounds();

    lat_midpoint_by_cell_ = model_->get_lat_mid_points();
    lon_midpoint_by_cell_ = model_->get_lon_mid_points();
    LOG_MEDIUM() << "size of lats = " << lat_midpoint_by_cell_.size() << " size of longs = " << lon_midpoint_by_cell_.size() << " heigght = " << height_ << " width = " << width_;
  }

  float lat, lon;
  for (unsigned i = 0; i < height_; ++i) {
    for (unsigned j = 0; j < width_; ++j) {
      if (model_->lat_and_long_supplied()) {
        lat = lat_midpoint_by_cell_[i];
        lon = lon_midpoint_by_cell_[j];
      } else {
        lat = 0.0;
        lon = 0.0;
      }

      LOG_FINE() << "Building cell " << i << "-" << j << " lat = " << lat << " long = " << lon;
      base_grid_[i][j].Build(i, j, lat, lon, model_);
      cached_grid_[i][j].Build(i, j, lat, lon, model_);
      init_cached_grid_[i][j].Build(i, j, lat, lon, model_);
    }
  }

  // Get Our Base Layer
  // Set Enabled Square Count
  enabled_cells_ = height_ * width_;

  // Flag any as Disabled if they don't match our Base_Layer
  for (unsigned i = 0; i < height_; ++i) {
    for (unsigned j = 0; j < width_; ++j) {
      if (base_layer_->get_value(i, j) <= 0) {
        base_grid_[i][j].set_enabled(false);
        LOG_FINEST() << "disabling cell at row " << i + 1 << " and col " << j + 1;
        enabled_cells_--;
      } else if (base_layer_->get_value(i, j) < 0) {
        LOG_FATAL()<< "found a negative value in the base layer '" << model_->get_base_layer_label() << "', value must be equal to or greater than 0";
      } else {
        base_grid_[i][j].set_area(base_layer_->get_value(i,j));
        enabled_rows_.push_back(i);
        enabled_cols_.push_back(j);
      }
    }
  }

  if (enabled_cells_ <= 0)
    LOG_FATAL()<< "found no valid base squares in the base layer '" << model_->get_base_layer_label() << "', see manual but this should be a matrix of numerics with at least one cell greater than the value '0'. please check";

}

/*
 * This method iterates over all cells and combines the cached grid onto the actual grid, used mainly in movement where we need to temporary hold individuals
 *
*/
void WorldView::MergeCachedGrid() {
  LOG_FINE();
  for (unsigned i = 0; i < height_; ++i) {  // Can't thread this, each cell has a pointer to growth and mortality for update agent, so there is a hidden shared resouce....
    for (unsigned j = 0; j < width_; ++j) {
      if (base_grid_[i][j].is_enabled()) {
        LOG_FINE() << "agents to merge into row " << i << " col = " << j  << " = " << cached_grid_[i][j].agents_.size();
        // Are we updateing agents parameters
        cached_grid_[i][j].update_agent_parameters();
        // Iterate over all source agents and find dead agents to override
        unsigned last_agent_ndx = 0;
        for (unsigned cache_agent_ndx = 0; cache_agent_ndx < cached_grid_[i][j].agents_.size(); ++cache_agent_ndx) {
          base_grid_[i][j].add_agent_alive(cached_grid_[i][j].agents_[cache_agent_ndx].get_scalar());

          while(last_agent_ndx < base_grid_[i][j].agents_.size()) {
            if (not base_grid_[i][j].agents_[last_agent_ndx].is_alive()) {
              base_grid_[i][j].agents_[last_agent_ndx] = cached_grid_[i][j].agents_[cache_agent_ndx];
              break;
            } else {
              last_agent_ndx++;
            }
          }
          //LOG_FINEST() << "last_agent_ndx " << last_agent_ndx << " size = " << base_grid_[i][j].agents_.size();
          if (last_agent_ndx >= base_grid_[i][j].agents_.size()) {
              base_grid_[i][j].agents_.push_back(cached_grid_[i][j].agents_[cache_agent_ndx]);
           }
        }
        cached_grid_[i][j].agents_.clear();
        // Now do tagged agents
        last_agent_ndx = 0;
        for (unsigned cache_agent_ndx = 0; cache_agent_ndx < cached_grid_[i][j].tagged_agents_.size(); ++cache_agent_ndx) {
          base_grid_[i][j].add_agent_alive(cached_grid_[i][j].tagged_agents_[cache_agent_ndx].get_scalar());

          while(last_agent_ndx < base_grid_[i][j].tagged_agents_.size()) {
            if (not base_grid_[i][j].tagged_agents_[last_agent_ndx].is_alive()) {
              base_grid_[i][j].tagged_agents_[last_agent_ndx] = cached_grid_[i][j].tagged_agents_[cache_agent_ndx];
              break;
            } else {
              last_agent_ndx++;
            }
          }
          //LOG_FINEST() << "last_agent_ndx " << last_agent_ndx << " size = " << base_grid_[i][j].agents_.size();
          if (last_agent_ndx >= base_grid_[i][j].tagged_agents_.size()) {
              base_grid_[i][j].tagged_agents_.push_back(cached_grid_[i][j].tagged_agents_[cache_agent_ndx]);
           }
        }
        cached_grid_[i][j].tagged_agents_.clear();
        LOG_FINE() << "cell " << i << "-" << j << " = " << base_grid_[i][j].get_total_individuals_alive();
      }
    }
  }
  LOG_FINE() << "finished merging world";
}

void WorldView::MergeCachedGrid(bool update_lat_long) {
  LOG_FINE();
  for (unsigned i = 0; i < height_; ++i) {  // Can't thread this, each cell has a pointer to growth and mortality for update agent, so there is a hidden shared resouce....
    for (unsigned j = 0; j < width_; ++j) {
      if (base_grid_[i][j].is_enabled()) {
        float lat_mid = base_grid_[i][j].get_lat();
        float lon_mid = base_grid_[i][j].get_lon();
        LOG_FINE() << "agents to merge into row " << i << " col = " << j  << " = " << cached_grid_[i][j].agents_.size();
        // Are we updateing agents parameters
        cached_grid_[i][j].update_agent_parameters();
        // Iterate over all source agents and find dead agents to override
        unsigned last_agent_ndx = 0;
        for (unsigned cache_agent_ndx = 0; cache_agent_ndx < cached_grid_[i][j].agents_.size(); ++cache_agent_ndx) {
          base_grid_[i][j].add_agent_alive(cached_grid_[i][j].agents_[cache_agent_ndx].get_scalar());

          if (update_lat_long) {
            cached_grid_[i][j].agents_[cache_agent_ndx].set_lat(lat_mid);
            cached_grid_[i][j].agents_[cache_agent_ndx].set_lon(lon_mid);
          }

          while(last_agent_ndx < base_grid_[i][j].agents_.size()) {
            if (not base_grid_[i][j].agents_[last_agent_ndx].is_alive()) {
              base_grid_[i][j].agents_[last_agent_ndx] = cached_grid_[i][j].agents_[cache_agent_ndx];
              break;
            } else {
              last_agent_ndx++;
            }
          }
          //LOG_FINEST() << "last_agent_ndx " << last_agent_ndx << " size = " << base_grid_[i][j].agents_.size();
          if (last_agent_ndx >= base_grid_[i][j].agents_.size()) {
              base_grid_[i][j].agents_.push_back(cached_grid_[i][j].agents_[cache_agent_ndx]);
           }
        }
        cached_grid_[i][j].agents_.clear();

        last_agent_ndx = 0;
        for (unsigned cache_agent_ndx = 0; cache_agent_ndx < cached_grid_[i][j].tagged_agents_.size(); ++cache_agent_ndx) {
          base_grid_[i][j].add_agent_alive(cached_grid_[i][j].tagged_agents_[cache_agent_ndx].get_scalar());

          if (update_lat_long) {
            cached_grid_[i][j].tagged_agents_[cache_agent_ndx].set_lat(lat_mid);
            cached_grid_[i][j].tagged_agents_[cache_agent_ndx].set_lon(lon_mid);
          }

          while(last_agent_ndx < base_grid_[i][j].tagged_agents_.size()) {
            if (not base_grid_[i][j].tagged_agents_[last_agent_ndx].is_alive()) {
              base_grid_[i][j].tagged_agents_[last_agent_ndx] = cached_grid_[i][j].tagged_agents_[cache_agent_ndx];
              break;
            } else {
              last_agent_ndx++;
            }
          }
          //LOG_FINEST() << "last_agent_ndx " << last_agent_ndx << " size = " << base_grid_[i][j].agents_.size();
          if (last_agent_ndx >= base_grid_[i][j].tagged_agents_.size()) {
              base_grid_[i][j].tagged_agents_.push_back(cached_grid_[i][j].tagged_agents_[cache_agent_ndx]);
           }
        }
        cached_grid_[i][j].tagged_agents_.clear();

      }
    }
  }
  LOG_FINE() << "finished merging world";
}

/*
 * Temporary cache world grid to save re-running initialisation phase
 *
*/
void WorldView::CachedWorldForInit() {
  LOG_MEDIUM() ;
  double individuals = 0.0;
  for (unsigned i = 0; i < height_; ++i) {  // Can't thread this, each cell has a pointer to growth and mortality for update agent, so there is a hidden shared resouce....
    for (unsigned j = 0; j < width_; ++j) {
      if (base_grid_[i][j].is_enabled()) {
        init_cached_grid_[i][j].agents_ = base_grid_[i][j].agents_;
        init_cached_grid_[i][j].tagged_agents_ = base_grid_[i][j].tagged_agents_;

        individuals = base_grid_[i][j].get_total_individuals_alive();
        init_cached_grid_[i][j].set_total_individuals_alive(individuals);
      }
    }
  }
}

/*
 * fill the row and col parameter with the cell index that contains the lat and lon given
 * Make sure that anything that calls this checks the model that lat and longs have been provided else this will cause
 * a SegFault
 *
*/
void WorldView::MergeWorldForInit() {
  LOG_MEDIUM() ;
  double individuals = 0.0;
  for (unsigned i = 0; i < height_; ++i) {  // Can't thread this, each cell has a pointer to growth and mortality for update agent, so there is a hidden shared resouce....
    for (unsigned j = 0; j < width_; ++j) {
      if (base_grid_[i][j].is_enabled()) {
        base_grid_[i][j].agents_ = init_cached_grid_[i][j].agents_;
        base_grid_[i][j].tagged_agents_ = init_cached_grid_[i][j].tagged_agents_;
        individuals = base_grid_[i][j].get_total_individuals_alive();
        base_grid_[i][j].set_total_individuals_alive(individuals);
      }
    }
  }
}
/*
 * fill the row and col parameter with the cell index that contains the lat and lon given
 * Make sure that anything that calls this checks the model that lat and longs have been provided else this will cause
 * a SegFault
 *
*/
void WorldView::get_cell_element(unsigned& row, unsigned& col, const float lat, const float lon) {
  for (unsigned i = 1; i <= height_; ++i) {
    if (lat < lat_bounds_[i]) {
      row = i - 1;
      break;
    }
  }
  for (unsigned j = 1; j <= width_; ++j) {
    if (lon < lon_bounds_[j]) {
      col = j - 1;
      break;
    }
  }
}


/*
 * This method gets the age frequencey of the world, this is called in intialisation to see if we have meet an equilibrium state
 *
*/
void WorldView::get_world_age_frequency(vector<float>& world_age_freq) {
  LOG_TRACE();
  bool do_age = true;
  world_age_freq.clear();
  world_age_freq.resize(model_->age_spread());
  vector<float> temp;
  for (unsigned i = 0; i < height_; ++i) {
    for (unsigned j = 0; j < width_; ++j) {
      if (base_grid_[i][j].is_enabled()) {
        base_grid_[i][j].get_age_frequency(temp, do_age);
        for(unsigned i = 0; i < temp.size(); ++i)
          world_age_freq[i] += temp[i];
      }
    }
  }
}

// Return a cell of the base world
WorldCell* WorldView::get_base_square(int RowIndex, int ColIndex) {
  //LOG_TRACE();
  return &base_grid_[RowIndex][ColIndex];
}

// Return a cell of the cached world
WorldCell* WorldView::get_cached_square(int RowIndex, int ColIndex) {
  //LOG_TRACE();
  return &cached_grid_[RowIndex][ColIndex];
}

// Called in model.cpp RunBasic(); to update agent parameters.
void WorldView::rebuild_agent_time_varying_params() {
  LOG_FINE();
  if (update_growth_params_ || update_mortality_params_) {
    for (unsigned i = 0; i < height_; ++i) {
      for (unsigned j = 0; j < width_; ++j) {
        if (base_grid_[i][j].is_enabled()) {
          if (update_growth_params_) {
            base_grid_[i][j].apply_growth_time_varying();
          }
          if (update_mortality_params_)
            base_grid_[i][j].apply_mortality_time_varying();
        }
      }
    }
  } else {
    LOG_FINE() << "Don't need to update agent values";
  }
  update_growth_params_ = false;
  update_mortality_params_ = false;
}

} /* namespace niwa */
