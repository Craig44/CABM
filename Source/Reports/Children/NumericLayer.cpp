/**
 * @file NumericLayer.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 31/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "NumericLayer.h"


#include "Model/Managers.h"
#include "Model/Model.h"
#include "Layers/Manager.h"
#include "Layers/Layer.h"

// namespaces
namespace niwa {
namespace reports {

/**
 * Default constructor
 *
 * @param model Pointer to the current model context
 */
NumericLayer::NumericLayer(Model* model) : Report(model) {
  model_state_ = State::kExecute;;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
  parameters_.Bind<string>(PARAM_LAYER_LABEL, &layer_label_, "The Numeric Layer label that is reported", "", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years", "", true);
  parameters_.Bind<string>(PARAM_TIME_STEP, &time_step_, "Time Step label", "", "");

}

/**
 * Build our relationships between this object and other objects
 */
void NumericLayer::DoBuild() {
  layer_ = model_->managers().layer()->GetNumericLayer(layer_label_);
  if (!layer_) {
    LOG_ERROR_P(PARAM_LAYER_LABEL) << "layer " << layer_label_ << " could not be found. Have you defined it? If you have make sure it is of type Numeric";
  }
  if (layer_->is_static()) {
    LOG_ERROR_P(PARAM_LAYER_LABEL) << "The point of this report is to print non staic layers, that are not user defined. Are you sure you want this report, go to usermanual for more information";
  }
}

/**
 * Execute this report
 */
void NumericLayer::DoExecute() {
  LOG_FINE() <<" printing report " << label_;
  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";
  cache_ << "year: " << model_->current_year() << "\n";
  cache_ << "time_step: " << time_step_ << "\n";
  cache_ << "values "<< REPORT_R_MATRIX<<"\n";
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col)
      cache_ << layer_->get_value(row, col) << " ";
    cache_ << "\n";
  }
  ready_for_writing_ = true;
}


} /* namespace reports */
} /* namespace niwa */
