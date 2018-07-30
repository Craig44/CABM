/**
 * @file Movement.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 30/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */

// headers
#include "Movement.h"

#include "World/WorldCell.h"
#include "World/WorldView.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
Movement::Movement(Model* model) : Process(model) {
  process_type_ = ProcessType::kTransition;
}

} /* namespace processes */
} /* namespace niwa */
