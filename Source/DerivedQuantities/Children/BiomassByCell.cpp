/**
 * @file BiomassByCell.cpp
 * @author  C.Marsh
 * @date 30/07/2018
 * @section LICENSE
 * @description
 *
 */
// Local headers
#include "BiomassByCell.h"

#include "Selectivities/Manager.h"
#include "InitialisationPhases/Manager.h"

#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/DoubleCompare.h"
#include "omp.h"

namespace niwa {
namespace derivedquantities {

// Constructor
BiomassByCell::BiomassByCell(Model* model) : DerivedQuantity(model) {
  // Default Variables
  parameters_.Bind<string>(PARAM_SELECTIVITY, &selectivity_label_, "A label for the selectivity", "");
  spatial_ = true;
}


/*
 * DoValidate
*/
void BiomassByCell::DoValidate() {


}

/*
 * DoBuild
*/
void BiomassByCell::DoBuild() {
  LOG_TRACE();

  // create offset contianers to help with threading
  cell_offset_.resize(model_->get_height());
  cache_in_space_.resize(model_->get_height());
  value_in_space_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cache_in_space_[i].resize(model_->get_width());
    value_in_space_[i].resize(model_->get_width());
    cell_offset_[i].resize(model_->get_width());
  }


  LOG_FINEST() << "selectivities supplied = " << selectivity_label_.size();
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
}

/**
 * The main function for calculating biomass
 */
void BiomassByCell::CalcBiomass(vector<Agent>& agents, float& value) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  if (not length_based_selectivity_) {
    for(auto& agent : agents) {
      if (agent.is_alive()) {
        if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_age_index())) {
          value += agent.get_weight() * agent.get_scalar();
        }
      }
    }
  } else {
    for(auto& agent : agents) {
      if (agent.is_alive()) {
        if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_length_bin_index())) {
          value += agent.get_weight() * agent.get_scalar();
        }
      }
    }
  }
}
/*
 * PreExecute
*/
void BiomassByCell::PreExecute() {
  LOG_FINE();
  if (utilities::doublecompare::IsOne(time_step_proportion_))
    return;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {

        cache_in_space_[row][col] = 0.0;
        CalcBiomass(cell->agents_, cache_in_space_[row][col]);
        CalcBiomass(cell->tagged_agents_, cache_in_space_[row][col]);
      }
    }
  }
}
/*
 * PreExecute
*/
void BiomassByCell::Execute() {
  LOG_FINE();
  if (!utilities::doublecompare::IsZero(time_step_proportion_)) {
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          value_in_space_[row][col] = 0.0;
          CalcBiomass(cell->agents_, value_in_space_[row][col]);
          CalcBiomass(cell->tagged_agents_, value_in_space_[row][col]);
        }
      }
    }
  }

  if (model_->state() == State::kInitialise) {
    unsigned initialisation_phase = model_->managers().initialisation_phase()->current_initialisation_phase();
    float value;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        if (initialisation_values_by_space_.size() <= initialisation_phase) {
          initialisation_values_by_space_.resize(initialisation_phase + 1);
          initialisation_values_by_space_[initialisation_phase].resize(model_->get_height());
          for (unsigned row = 0; row < model_->get_height(); ++row)
            initialisation_values_by_space_[initialisation_phase][row].resize(model_->get_width());
        }

        if (utilities::doublecompare::IsZero(time_step_proportion_)) {
          value = cache_in_space_[row][col];
          initialisation_values_by_space_[initialisation_phase][row][col].push_back(value);
        } else if (utilities::doublecompare::IsOne(time_step_proportion_)) {
          value = value_in_space_[row][col];
          initialisation_values_by_space_[initialisation_phase][row][col].push_back(value);
        } else {
          value = cache_in_space_[row][col] + ((value_in_space_[row][col] - cache_in_space_[row][col]) * time_step_proportion_);
          initialisation_values_by_space_[initialisation_phase][row][col].push_back(value);
        }
      }
    }
  } else {
    if (utilities::doublecompare::IsZero(time_step_proportion_))
      values_by_space_[model_->current_year()] = cache_in_space_;
    else if (utilities::doublecompare::IsOne(time_step_proportion_))
      values_by_space_[model_->current_year()] = value_in_space_;
    else {
      vector<vector<float>> temp(model_->get_height());
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        temp[row].resize(model_->get_width());
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          temp[row][col] = cache_in_space_[row][col] + ((value_in_space_[row][col] - cache_in_space_[row][col]) * time_step_proportion_);
        }
      }
      values_by_space_[model_->current_year()] = temp;
    }
  }

}
} /* namespace derivedquantities */
} /* namespace niwa */
