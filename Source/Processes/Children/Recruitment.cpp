/**
 * @file RecruitmentConstant.cpp
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
#include "Recruitment.h"

//#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/DoubleCompare.h"
#include <omp.h>
// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
Recruitment::Recruitment(Model* model) : Process(model) {
  process_type_ = ProcessType::kRecruitment;
  parameters_.Bind<float>(PARAM_B0, &b0_, "B0", "",false);
  parameters_.Bind<string>(PARAM_RECRUITMENT_LAYER_LABEL, &recruitment_layer_label_, "A label for the recruitment layer", "");
  parameters_.Bind<string>(PARAM_SSB, &ssb_label_, "A label for the SSB derived quantity", "");
}

// DoBuild
void Recruitment::DoBuild() {
  LOG_MEDIUM() << "Build Parent this will be overided if the child has this method.";


}

} /* namespace processes */
} /* namespace niwa */
