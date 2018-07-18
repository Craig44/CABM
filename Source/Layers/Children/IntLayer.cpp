/**
 * @file IntLayer.cpp
 * @author  C.Marsh
 * @date 12/07/2018
 * @section LICENSE
 * @description
 *
 */

// headers
#include "IntLayer.h"

namespace niwa {
namespace layers {

// Constructor
IntLayer::IntLayer(Model* model) : Layer(model) {
  int_table_ = new parameters::Table(PARAM_LAYER);
  parameters_.BindTable(PARAM_LAYER, int_table_, "Table of layer attributes", "", false, false);
  // Default Variables
  grid_ = 0;
}

// Deconstructor
IntLayer::~IntLayer() {
  LOG_TRACE();
  // Clean Our Grid
  if (grid_ != 0) {
    for (unsigned i = 0; i < height_; ++i) {
      delete [] grid_[i];
      grid_[i] = 0;
    }
    delete [] grid_;
  }
  // remove input table
  delete int_table_;

}

void IntLayer::DoValidate() {
  LOG_TRACE();
  // Allocate Space for ourlayer on the heap
  grid_ = new unsigned*[height_];
  for (unsigned i = 0; i < height_; ++i)
    grid_[i] = new unsigned[width_];
}



void IntLayer::DoBuild() {
  LOG_TRACE();
  /**
   * Build our layer pointer
   */
  vector<vector<string>>& input_data = int_table_->data();

  LOG_FINEST() << "rows = " << input_data.size() << " columns = " << input_data[0].size();
  unsigned row_iter = 0;
  for (vector<string> row : input_data) {

    if (row.size() != width_)
      LOG_ERROR_P(PARAM_LAYER) << "columns supplied '"<< row.size() << "' ncols on the @model = '" << height_ << "' please resolve this descrepency";

    if ((row_iter + 1) > height_)
      LOG_FATAL_P(PARAM_LAYER) << "you supplied at least '"<< row_iter + 1 << "' rows in the layer but, nrows on the @model = '" << width_ << "' these must be the same";

    unsigned value;
    for (unsigned i = 0; i < row.size(); ++i) {
      value = utilities::ToInline<string, unsigned>(row[i]);
      if (value < 0)
        LOG_ERROR_P(PARAM_LAYER) << "at row '" << row_iter << "' and column '" << i + 1 << "' we found a value less than zero this is not allowed for this layer type, they must be >= 0, please sort out";
      grid_[row_iter][i] = value;
    }
    ++row_iter;
  }
}

/*
 * setValue
*/
unsigned IntLayer::get_value(unsigned RowIndex, unsigned ColIndex) {
#ifndef OPTIMIZE
// TODO do some error catching for debugging purposes
#endif
  return grid_[RowIndex][ColIndex];
}

/*
 * setValue
*/
void IntLayer::set_value(unsigned RowIndex, unsigned ColIndex, unsigned Value) {
#ifndef OPTIMIZE
// TODO do some error catching for debugging purposes
#endif
  grid_[RowIndex][ColIndex] = Value;
}
} /* namespace layers */
} /* namespace niwa */
