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
  int_table_ = new parameters::Table(PARAM_TABLE);
  parameters_.BindTable(PARAM_LAYER, int_table_, "Table of layer attributes", "", false);
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
  const vector<string>& columns = int_table_->columns();

  unsigned row_iter = 1;
  for (vector<string> row : input_data) {
    if (row.size() != height_)
      LOG_ERROR_P(PARAM_LAYER) << "row.size '"<< row.size() << "' != height on the @model = '" << height_ << " please resolve this descrepency" << endl;
    if ((columns.size()) != width_)
      LOG_ERROR_P(PARAM_LAYER) << "columns.size '"<< columns.size() << "' != width on the @model = '" << width_ << " please resolve this descdescrepencyrepecy" << endl;

    unsigned value;
    for (unsigned i = 0; i < row.size(); ++i) {
      value = utilities::ToInline<string, unsigned>(row[i]);
      if (value < 0)
        LOG_ERROR_P(PARAM_LAYER) << "at row '" << row_iter << "' and column '" << i + 1 << "' we found a value less than zero this is not allowed for this layer type, they must be >= 0, please sort out";
      grid_[row_iter - 1][i] = value;
    }
    ++row_iter;
  }
}



/*
 * setValue
*/
void IntLayer::setValue(int RowIndex, int ColIndex, unsigned Value) {
#ifndef OPTIMIZE
// TODO do some error catching for debugging purposes
#endif
  grid_[RowIndex][ColIndex] = Value;
}
} /* namespace layers */
} /* namespace niwa */
