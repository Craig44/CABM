/**
 * @file Abundance.cpp
 * @author  C.Marsh
 * @date 10/09/2018
 * @section LICENSE
 * @description
 *
 */

// Local headers
#include "Abundance.h"

#include "World/WorldView.h"
#include "World/WorldCell.h"

namespace niwa {
namespace layers {

/**
 * Usual constructor
 */
Abundance::Abundance(Model* model) : NumericLayer(model) {
  // Default Variables
  parameters_.Bind<bool>(PARAM_MATURE, &mature_, "Mature biomass or total biomass?", "", false);
  layer_type_ = LayerType::kNumeric;
  static_layer_ = false;
}


/**
 * Validate Instance
 */
void Abundance::DoValidate() {
  LOG_TRACE();
}


/**
 * Build pointers class
 */

void Abundance::DoBuild() {
  LOG_TRACE();
  world_ = model_->world_view();
  if(!world_)
    LOG_CODE_ERROR() << "if(!world_) can't get a pointer to the world, what is the point in life?";
}


/**
 * Return value of Abundance in this cell on the fly
 */
float Abundance::get_value(unsigned RowIndex, unsigned ColIndex) {
  WorldCell* cell = world_->get_base_square(RowIndex, ColIndex);
  float value = 0.0;
  if (cell->is_enabled()) {
    if (mature_) {
      for (Agent& agent : cell->agents_) {
        if (agent.is_alive() and agent.get_maturity())
          value += agent.get_scalar();
      }
    } else {
      for (Agent& agent : cell->agents_) {
        if (agent.is_alive())
          value += agent.get_scalar();
      }
    }
  }
  return value;
}

/**
 * Return value of Abundance in this cell on the fly
 */
float Abundance::get_value(unsigned RowIndex, unsigned ColIndex, unsigned year) {
  WorldCell* cell = world_->get_base_square(RowIndex, ColIndex);
  float value = 0.0;
  if (cell->is_enabled()) {
    if (mature_) {
      for (Agent& agent : cell->agents_) {
        if (agent.is_alive() and agent.get_maturity())
          value += agent.get_scalar();
      }
    } else {
      for (Agent& agent : cell->agents_) {
        if (agent.is_alive())
          value += agent.get_scalar();
      }
    }
  }
  return value;
}

} /* namespace layers */
} /* namespace niwa */
