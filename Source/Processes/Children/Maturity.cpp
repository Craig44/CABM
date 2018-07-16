/**
 * @file Maturity.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 15/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */

// headers
#include "Maturity.h"

#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "Selectivities/Manager.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
Maturity::Maturity(Model* model) : Process(model) {
  process_type_ = ProcessType::kTransition;
  parameters_.Bind<string>(PARAM_SELECTIVITY_LABEL, &selectivity_label_, "Label for the selectivity block note that this selectivity represents probability of becoming mature", "");
}

/**
 * Build class relationships
 */
void Maturity::DoBuild() {
  LOG_TRACE();
  // Get selectivity
  selectivity_ = model_->managers().selectivity()->GetSelectivity(selectivity_label_);
  if (!selectivity_)
    LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << ": selectivity " << selectivity_label_ << " does not exist. Have you defined it?";
}


/**
 * Execute process
 */
void Maturity::DoExecute() {
  LOG_TRACE();
  // Iterate over all cells
  float probability_mature_at_age;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        auto& agents = cell->get_agents();
        LOG_FINEST() << "about to convert " << agents.size() << " through the maturity process";
        unsigned counter = 1;
        for (Agent& agent : agents) {
          probability_mature_at_age = selectivity_->GetResult(agent.age());
          agent.maturity(probability_mature_at_age);
          counter++;
        }
      }
    }
  }
  LOG_TRACE();
}
} /* namespace processes */
} /* namespace niwa */
