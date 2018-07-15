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

// namespaces
namespace niwa {

// Constructor
WorldView::WorldView(Model* model) :
    model_(model) {
  base_grid_ = 0;
  cached_grid_ = 0;
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
}

// Build Method
void WorldView::Validate() {
  width_ = model_->get_width();
  height_ = model_->get_height();
}

// Build Method
void WorldView::Build() {
  LOG_TRACE();
  base_layer_ = model_->managers().layer()->GetNumericLayer(model_->get_base_layer_label());
  if (!base_layer_) {
    LOG_ERROR() << "The base layer '" << model_->get_base_layer_label() << "' found in the @model block could not be found, please check there is a @layer defined for this layer and that it is type = 'numeric'";
  }

  if (model_->get_lat_layer_label() != "") {
    lat_layer_ = model_->managers().layer()->GetNumericLayer(model_->get_lat_layer_label());
    if (!lat_layer_) {
      LOG_ERROR() << "The lattude layer '" << model_->get_lat_layer_label() << "' found in the @model block could not be found, please check there is a @layer defined for this layer and that it is type = 'numeric'";
    }
  }

  if (model_->get_long_layer_label() != "") {
    long_layer_ = model_->managers().layer()->GetNumericLayer(model_->get_long_layer_label());
    if (!long_layer_) {
      LOG_ERROR() << "The longitude layer '" << model_->get_long_layer_label() << "' found in the @model block could not be found, please check there is a @layer defined for this layer and that it is type = 'numeric'";
    }
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

  // Set Variables (Can't do it above. Stupid blah ISO C++)
  double lat = 0.0, lon = 0.0;
  for (unsigned i = 0; i < height_; ++i) {
    for (unsigned j = 0; j < width_; ++j) {
      if (lat_layer_)
        lat = lat_layer_->get_value(i,j);
      if (long_layer_)
        lon = long_layer_->get_value(i,j);
      base_grid_[i][j].Build(i,j,lat,lon, model_->min_age(), model_->max_age());
      cached_grid_[i][j].Build(i,j,lat,lon, model_->min_age(), model_->max_age());
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
        enabled_cells_--;
      } else if (base_layer_->get_value(i, j) < 0) {
        LOG_FATAL()<< "found a negative value in the base layer '" << model_->get_base_layer_label() << "', value must be equal to or greater than 0";
      } else
      base_grid_[i][j].set_area(base_layer_->get_value(i,j));
    }
  }

  if (enabled_cells_ <= 0)
    LOG_FATAL()<< "found no valid base squares in the base layer '" << model_->get_base_layer_label() << "', see manual but this should be a matrix of numerics with at least one cell greater than the value '0'. please check";

}

/*
 * This method gets the age frequencey of the world, this is called in intialisation to see if we have meet an equilibrium state
 *
*/
void WorldView::get_world_age_frequency(vector<unsigned>& world_age_freq) {
  world_age_freq.clear();
  world_age_freq.resize(model_->age_spread());
  vector<unsigned> temp;
  for (unsigned i = 0; i < height_; ++i) {
    for (unsigned j = 0; j < width_; ++j) {
      if (base_grid_[i][j].is_enabled()) {
        base_grid_[i][j].get_age_frequency(temp);
        for(unsigned i = 0; i < temp.size(); ++i)
          world_age_freq[i] += temp[i];
      }
    }
  }
}

// Return a cell of the base world
WorldCell* WorldView::get_base_square(int RowIndex, int ColIndex) {
  LOG_TRACE();
  return &base_grid_[RowIndex][ColIndex];
}

// Return a cell of the cached world
WorldCell* WorldView::get_cached_square(int RowIndex, int ColIndex) {
  LOG_TRACE();
  return &cached_grid_[RowIndex][ColIndex];
}


} /* namespace niwa */
