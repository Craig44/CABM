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

  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";
  cache_ << "global_scalar: " << model_->get_scalar();
  cache_ << "\n";
  ready_for_writing_ = true;

}

} /* namespace reports */
} /* namespace niwa */
