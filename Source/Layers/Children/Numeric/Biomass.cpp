/**
 * @file Biomass.cpp
 * @author  C.Marsh
 * @date 31/07/2018
 * @section LICENSE
 * @description
 *
 */

// Local headers
#include "Biomass.h"

#include "World/WorldView.h"
#include "World/WorldCell.h"

namespace niwa {
namespace layers {

/**
 * Usual constructor
 */
Biomass::Biomass(Model* model) : NumericLayer(model) {
  // Default Variables
  parameters_.Bind<bool>(PARAM_MATURE, &mature_, "Mature biomass or total biomass?", "", false);
  layer_type_ = LayerType::kNumeric;
  static_layer_ = false;
}


/**
 * Validate Instance
 */
void Biomass::DoValidate() {
  LOG_TRACE();
}


/**
 * Build pointers class
 */

void Biomass::DoBuild() {
  LOG_TRACE();
  world_ = model_->world_view();
  if(!world_)
    LOG_CODE_ERROR() << "if(!world_) can't get a pointer to the world, what is the point in life?";
}


/**
 * Return value of biomass in this cell on the fly
 */
float Biomass::get_value(unsigned RowIndex, unsigned ColIndex) {
  WorldCell* cell = world_->get_base_square(RowIndex, ColIndex);
  if (cell->is_enabled()) {
    float value = 0.0;
    if (mature_) {
      for (Agent& agent : cell->agents_) {
        if (agent.get_maturity())
          value += agent.get_weight() * agent.get_scalar();
      }
    } else {
      for (Agent& agent : cell->agents_)
        value += agent.get_weight() * agent.get_scalar();
    }
    return value;
  }
  return 1.0;
}

/**
 * Return value of biomass in this cell on the fly
 */
float Biomass::get_value(unsigned RowIndex, unsigned ColIndex, unsigned year) {
  WorldCell* cell = world_->get_base_square(RowIndex, ColIndex);
  if (cell->is_enabled()) {
    float value = 0.0;
    if (mature_) {
      for (Agent& agent : cell->agents_) {
        if (agent.get_maturity())
          value += agent.get_weight() * agent.get_scalar();
      }
    } else {
      for (Agent& agent : cell->agents_)
        value += agent.get_weight() * agent.get_scalar();
    }
    return value;
  }
  return 1.0;
}

} /* namespace layers */
} /* namespace niwa */
