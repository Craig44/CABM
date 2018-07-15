/**
 * @file RecruitmentBevertonHolt.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 12/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This file exists to keep documentation generator happy
 */

// headers
#include "RecruitmentBevertonHolt.h"

#include "Layers/Manager.h"
//#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
RecruitmentBevertonHolt::RecruitmentBevertonHolt(Model* model) : Process(model) {
  process_type_ = ProcessType::kRecruitment;

  parameters_.Bind<double>(PARAM_B0, &b0_, "B0", "",false);
  parameters_.Bind<double>(PARAM_STEEPNESS, &steepness_, "Steepness", "", 1.0);
  parameters_.Bind<double>(PARAM_YCS_VALUES, &ycs_values_, "YCS (Year-Class-Strength) Values", "");
  parameters_.Bind<string>(PARAM_RECRUITMENT_LAYER_LABEL, &recruitment_layer_label_, "A label for the recruitment layer", "");
  parameters_.Bind<string>(PARAM_SSB, &ssb_label_, "A label for the SSB layer", "");

}


void RecruitmentBevertonHolt::DoBuild() {

}

} /* namespace processes */
} /* namespace niwa */
