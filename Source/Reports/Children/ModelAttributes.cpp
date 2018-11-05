/**
 * @file ModelAttributes.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 18/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "ModelAttributes.h"

#include <boost/algorithm/string/join.hpp>

#include "Model/Managers.h"
#include "Model/Model.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"

// namespaces
namespace niwa {
namespace reports {

/**
 * Default constructor
 *
 * @param model Pointer to the current model context
 */
ModelAttributes::ModelAttributes(Model* model) : Report(model) {
  model_state_ = State::kIterationComplete;;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
}

/**
 * Build our relationships between this object and other objects
 */
void ModelAttributes::DoBuild() {
  LOG_TRACE();
}

/**
 * Execute this report
 */
void ModelAttributes::DoExecute() {
  LOG_FINE() <<" printing report " << label_;

  auto scalars = model_->get_scalars();
  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";
  for (auto scalar : scalars)
    cache_ << scalar.first << ": " << scalar.second << "\n";
  ready_for_writing_ = true;

  cache_ << "length_mid_points: ";
  vector<float> length_mid_point = model_->length_bin_mid_points();
  for (unsigned length_ndx = 0; length_ndx < length_mid_point.size(); ++length_ndx)
    cache_ << length_mid_point[length_ndx] << " ";

}

} /* namespace reports */
} /* namespace niwa */
