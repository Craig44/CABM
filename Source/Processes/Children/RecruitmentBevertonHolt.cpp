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

// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
RecruitmentBevertonHolt::RecruitmentBevertonHolt(Model* model) : Process(model) {
  process_type_ = ProcessType::kRecruitment;
}

}
}
