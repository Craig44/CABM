/**
 * @file NumericLayer.cpp
 * @author  C.Marsh
 * @date 12/07/2018
 * @section LICENSE
 * @description
 *
 */

// headers
#include "NumericLayer.h"

namespace niwa {
namespace layers {

// Constructor
NumericLayer::NumericLayer(Model* model) : Layer(model) {
  data_table_ = new parameters::Table(PARAM_LAYER);
  parameters_.BindTable(PARAM_LAYER, data_table_, "Table of layer attributes", "", false, true);
  // Default Variables
  grid_ = 0;
}

// Deconstructor
NumericLayer::~NumericLayer() {
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
  delete data_table_;

}

void NumericLayer::DoValidate() {
  LOG_TRACE();
  // Allocate Space for ourlayer on the heap
  grid_ = new float*[height_];
  for (unsigned i = 0; i < height_; ++i)
    grid_[i] = new float[width_];
}



void NumericLayer::DoBuild() {
  LOG_TRACE();
  /**
   * Build our layer pointer
   */
  vector<vector<string>>& input_data = data_table_->data();

  LOG_FINEST() << "rows = " << input_data.size() << " columns = " << input_data[0].size();
  unsigned row_iter = 0;
  for (vector<string> row : input_data) {

    if (row.size() != width_)
      LOG_ERROR_P(PARAM_LAYER) << "columns supplied '"<< row.size() << "' ncols on the @model = '" << height_ << "' please resolve this descrepency";

    if ((row_iter + 1) > height_)
      LOG_FATAL_P(PARAM_LAYER) << "you supplied at least '"<< row_iter + 1 << "' rows in the layer but, nrows on the @model = '" << width_ << "' these must be the same";

    float value;
    for (unsigned i = 0; i < row.size(); ++i) {
      value = utilities::ToInline<string, float>(row[i]);
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
void NumericLayer::set_value(int RowIndex, int ColIndex, float Value) {
#ifndef OPTIMIZE
// TODO do some error catching for debugging purposes
#endif
  grid_[RowIndex][ColIndex] = Value;
}
} /* namespace layers */
} /* namespace niwa */
