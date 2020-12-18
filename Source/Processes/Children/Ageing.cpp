/**
 * @file Ageing.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 15/10/2020
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This file exists to keep documentation generator happy
 */

// headers
#include "Ageing.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "Agents/Agent.h"
// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
Ageing::Ageing(Model* model) : Process(model) {
  process_type_ = ProcessType::kNullProcess;
}


void Ageing::DoExecute() {
  LOG_MEDIUM();
 // #pragma omp parallel for collapse(2)
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        for(auto& agent : cell->agents_) {
          if (agent.is_alive())
            agent.birthday();
        }
        for(auto& agent : cell->tagged_agents_) {
          if (agent.is_alive())
            agent.birthday();
        }
      }
    }
  }
}

}
}
