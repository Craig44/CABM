/**
 * @file PreferenceFunction.cpp
 * @author  C.Marsh www.github/Craig44
 * @version 1.0
 * @date 26/7/2018
 * @section LICENSE
 *
 *
 */

// Headers
#include "PreferenceFunction.h"


#include "Model/Model.h"

// Namesapces
namespace niwa {

/**
 * Explicit Constructor
 */
PreferenceFunction::PreferenceFunction(Model* model)
: model_(model) {

  parameters_.Bind<string>(PARAM_LABEL, &label_, "The label for this selectivity", "");
  parameters_.Bind<string>(PARAM_TYPE, &type_, "The type of selectivity", "");
}

/**
 *
 */
void PreferenceFunction::Validate() {
  parameters_.Populate(model_);
  DoValidate();
}

/**
 *
 */
void PreferenceFunction::Build() {
  parameters_.Populate(model_);
  DoBuild();
}
/**
 *
 */
void PreferenceFunction::Reset() {

}



} /* namespace niwa */
