/**
 * @file CategoricalLayer.cpp
 * @author  C.Marsh
 * @date 26/07/2018
 * @section LICENSE
 * @description
 *
 */

// headers
#include "CategoricalLayer.h"

#include "Utilities/DoubleCompare.h"
namespace niwa {
namespace layers {

// Constructor
CategoricalLayer::CategoricalLayer(Model* model) : Layer(model) {
  data_table_ = new parameters::Table(PARAM_LAYER);
  parameters_.BindTable(PARAM_LAYER, data_table_, "Table of layer attributes", "", false, true);
  layer_type_ = LayerType::kCategorical;

  // Default Variables
  grid_ = 0;
}

// Deconstructor
CategoricalLayer::~CategoricalLayer() {
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

void CategoricalLayer::DoValidate() {
  LOG_TRACE();
  // Allocate Space for ourlayer on the heap
  grid_ = new string*[height_];
  for (unsigned i = 0; i < height_; ++i)
    grid_[i] = new string[width_];
}



void CategoricalLayer::DoBuild() {
  LOG_TRACE();
  /**
   * Build our layer pointer
   */
  vector<vector<string>>& input_data = data_table_->data();
  LOG_FINE() << "rows = " << input_data.size() << " columns = " << input_data[0].size();
  unsigned row_iter = 0;
  for (vector<string> row : input_data) {

    if (row.size() != width_)
      LOG_ERROR_P(PARAM_LAYER) << "columns supplied '"<< row.size() << "' ncols on the @model = '" << height_ << "' please resolve this descrepency";

    if ((row_iter + 1) > height_)
      LOG_FATAL_P(PARAM_LAYER) << "you supplied at least '"<< row_iter + 1 << "' rows in the layer but, nrows on the @model = '" << width_ << "' these must be the same";

    for (unsigned i = 0; i < row.size(); ++i) {
      LOG_FINE() << "row " <<  row_iter + 1 << "col = " << i + 1 << " value = " << row[i];
      grid_[row_iter][i] = row[i];
    }
    ++row_iter;
  }
}

/*
 * get_value
*/
string CategoricalLayer::get_value(unsigned RowIndex, unsigned ColIndex) {
  LOG_FINE();
  return grid_[RowIndex][ColIndex];
}

/*
 * get_value
*/
string CategoricalLayer::get_value(unsigned RowIndex, unsigned ColIndex, unsigned year) {
  //LOG_TRACE();
  return grid_[RowIndex][ColIndex];
}


/*
 * setValue
*/
void CategoricalLayer::set_value(unsigned RowIndex, unsigned ColIndex, string Value) {
  grid_[RowIndex][ColIndex] = Value;
}
} /* namespace layers */
} /* namespace niwa */
