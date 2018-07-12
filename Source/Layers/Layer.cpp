/**
 * @file Layer.cpp
 * @author  C.Marsh
 * @date 12/7/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "Layer.h"

#include <limits>


// namespaces
namespace niwa {

/**
 * Default constructor
 *
 * Bind any parameters that are allowed to be loaded from the configuration files.
 * Set bounds on registered parameters
 * Set some initial values
 *
 * Note: The constructor is parsed to generate Latex for the documentation.
 */
Layer::Layer(Model* model)
  : model_(model){
  parameters_.Bind<string>(PARAM_LABEL, &label_, "Label of the Layer", "");
  parameters_.Bind<string>(PARAM_TYPE, &type_, "Type of Layer", "");
}

/**
 * Populate any parameters,
 * Validate values are within expected ranges when we cannot use bind<>() overloads
 *
 * Note: all parameters are populated from configuration files
 */
void Layer::Validate() {
  parameters_.Populate(model_);

  height_ = model_->get_height();
  width_ = model_->get_height();

  DoValidate();
}

/**
 * Build any objects that will need to be utilised by this object.
 * Obtain smart_pointers to any objects that will be used by this object.
 */
void Layer::Build() {
  LOG_TRACE();
/*
  *
   * ensure the time steps we have are valid

  TimeStep* time_step = model_->managers().time_step()->GetTimeStep(time_step_label_);
  if (!time_step)
    LOG_FATAL_P(PARAM_TIME_STEP) << " (" << time_step_label_ << ") could not be found. Have you defined it?";
  time_step->SubscribeToBlock(this);
  time_step->SubscribeToInitialisationBlock(this);*/
  DoBuild();
}

/**
 * Reset our derived quantity values
 */
void Layer::Reset() {

}

} /* namespace niwa */
