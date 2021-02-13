/**
 * @file Abundance.cpp
 * @author  C.Marsh
 * @date 31/07/2018
 * @section LICENSE
 * @description
 *
 */

// Local headers
#include "Abundance.h"

#include "World/WorldView.h"
#include "World/WorldCell.h"
#include "Selectivities/Manager.h"
#include "Utilities/RandomNumberGenerator.h"

namespace niwa {
namespace layers {

/**
 * Usual constructor
 */
Abundance::Abundance(Model* model) : NumericLayer(model) {
  // Default Variables
  parameters_.Bind<bool>(PARAM_MATURE, &mature_, "Mature Abundance or total Abundance?", "", true);
  parameters_.Bind<string>(PARAM_SELECTIVITY, &selectivity_label_, "A label for the selectivity", "", true);

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

  if (!parameters_.Get(PARAM_MATURE)->has_been_defined() & !parameters_.Get(PARAM_SELECTIVITY)->has_been_defined())
    LOG_FATAL_P(PARAM_LABEL) << "You need to supply either " << PARAM_MATURE << " or " << PARAM_SELECTIVITY << " inputs, you have not supplied either.";

  if (parameters_.Get(PARAM_MATURE)->has_been_defined() & parameters_.Get(PARAM_SELECTIVITY)->has_been_defined())
    LOG_FATAL_P(PARAM_LABEL) << "You supplied both " << PARAM_MATURE << " and " << PARAM_SELECTIVITY << " inputs, you have not need to supply one or the other.";

  if (parameters_.Get(PARAM_SELECTIVITY)->has_been_defined()) {
    apply_selectivity_ = true;

    // Build selectivity links
    if (selectivity_label_.size() == 1)
      selectivity_label_.assign(2, selectivity_label_[0]);

    if (selectivity_label_.size() > 2) {
      LOG_ERROR_P(PARAM_SELECTIVITY) << "You suppled " << selectivity_label_.size()  << " Selectiviites, you can only have one for each sex max = 2";
    }
    LOG_FINEST() << "selectivities supplied = " << selectivity_label_.size();

    bool first = true;
    for (auto label : selectivity_label_) {
      Selectivity* temp_selectivity = model_->managers().selectivity()->GetSelectivity(label);
      if (!temp_selectivity)
        LOG_ERROR_P(PARAM_SELECTIVITY) << ": selectivity " << label << " does not exist. Have you defined it?";

      selectivity_.push_back(temp_selectivity);
      if (first) {
        first = false;
        length_based_selectivity_ = temp_selectivity->is_length_based();
      } else {
        if (length_based_selectivity_ != temp_selectivity->is_length_based()) {
          LOG_ERROR_P(PARAM_SELECTIVITY) << "The selectivity  " << label << " was not the same type (age or length based) as the previous selectivity label";
        }
      }
    }
    mature_ = false;
  } else {
    apply_selectivity_ = false;
  }


  LOG_MEDIUM() << "apply mature = " << mature_ << " apply selectivity = " << apply_selectivity_ << " length based = " << length_based_selectivity_;
  if (parameters_.Get(PARAM_MATURE)->has_been_defined()){
    if (not mature_)
      LOG_ERROR_P(PARAM_MATURE) << "Only specify the subcommand " << PARAM_MATURE << " if you want it true, otherwise select a selectivity";

  }
}


/**
 * Return value of Abundance in this cell on the fly
 */
float Abundance::get_value(unsigned RowIndex, unsigned ColIndex) {
  LOG_MEDIUM() << "get_value(unsigned RowIndex, unsigned ColIndex)";
  WorldCell* cell = world_->get_base_square(RowIndex, ColIndex);
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  if (cell->is_enabled()) {
    float value = 0.0;
    if (mature_) {
      for (Agent& agent : cell->agents_) {
        if (agent.is_alive() and agent.get_maturity())
          value +=  agent.get_scalar();
      }
    } else {
      if (length_based_selectivity_) {
        for (Agent& agent : cell->agents_) {
          if (agent.is_alive()) {
            if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_length_bin_index()))
              value += agent.get_scalar();
          }
        }
        for (Agent& agent : cell->tagged_agents_) {
          if (agent.is_alive()) {
            if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_length_bin_index()))
              value +=  agent.get_scalar();
          }
        }
      } else  {
        for (Agent& agent : cell->agents_) {
          if (agent.is_alive()) {
            if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_age_index()))
              value +=  agent.get_scalar();
          }
        }
        for (Agent& agent : cell->tagged_agents_) {
          if (agent.is_alive()) {
            if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_age_index()))
              value += agent.get_scalar();
          }
        }
      }
    }
    return value;
  }
  return 1.0;
}

/**
 * Return value of biomass in this cell on the fly
 */
float Abundance::get_value(unsigned RowIndex, unsigned ColIndex, unsigned year) {
  LOG_MEDIUM() << "get_value(unsigned RowIndex, unsigned ColIndex, unsigned year)";

  WorldCell* cell = world_->get_base_square(RowIndex, ColIndex);
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  if (cell->is_enabled()) {
    float value = 0.0;
    if (mature_) {
      for (Agent& agent : cell->agents_) {
        if (agent.is_alive() and agent.get_maturity())
          value += agent.get_scalar();
      }
    } else {
      if (length_based_selectivity_) {
        for (Agent& agent : cell->agents_) {
          if( agent.is_alive()) {
            if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_length_bin_index()))
              value += agent.get_scalar();
          }
        }
        for (Agent& agent : cell->tagged_agents_) {
          if( agent.is_alive()) {
            if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_length_bin_index()))
              value +=  agent.get_scalar();
          }
        }
      } else  {
        for (Agent& agent : cell->agents_) {
          if( agent.is_alive()) {
            if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_age_index()))
              value +=  agent.get_scalar();
          }
        }
        for (Agent& agent : cell->tagged_agents_) {
          if( agent.is_alive()) {
            if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_age_index()))
              value +=  agent.get_scalar();
          }
        }
      }
    }
    return value;
  }
  return 1.0;
}


} /* namespace layers */
} /* namespace niwa */
