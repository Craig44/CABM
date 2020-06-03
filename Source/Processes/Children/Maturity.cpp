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
#include "Utilities/RandomNumberGenerator.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
Maturity::Maturity(Model* model) : Process(model) {
  process_type_ = ProcessType::kTransition;
}

/**
 * Build class relationships
 */
void Maturity::DoBuild() {
  LOG_TRACE();
  // Get selectivity
  selectivity_label_ = model_->get_maturity_ogive();
  for (auto label : selectivity_label_) {
    Selectivity* temp_selectivity = nullptr;
    temp_selectivity = model_->managers().selectivity()->GetSelectivity(label);
    if (!temp_selectivity)
      LOG_CODE_ERROR()<< "this should have been checked on the ModelDoBuild please check out";
    length_based_selectivity_ = temp_selectivity->is_length_based();
    selectivity_.push_back(temp_selectivity);
  }

  cell_offset_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cell_offset_[i].resize(model_->get_width());
  }

}

/**
 * Execute process
 */
void Maturity::ApplyStochasticMaturity(vector<Agent>& agents) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  if (not length_based_selectivity_) {
    for (auto& agent : agents) {
      if (not agent.get_maturity() & agent.is_alive()) {
        if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_age_index())) {
          agent.set_maturity(true);
          ++mature_conversion_;
        }
      }
    }
  } else {
    for (auto& agent : agents) {
      if (not agent.get_maturity() & agent.is_alive()) {
        if (rng.chance() <= selectivity_[agent.get_sex()]->GetResult(agent.get_length_bin_index())) {
          agent.set_maturity(true);
          ++mature_conversion_;
        }
      }
    }
  }

}
/**
 * Execute process
 */
void Maturity::DoExecute() {
  LOG_MEDIUM() << label_;
  mature_conversion_ = 0;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        LOG_FINEST() << "about to convert " << cell->agents_.size() << " through the maturity process";
        ApplyStochasticMaturity(cell->agents_);
        ApplyStochasticMaturity(cell->tagged_agents_);

      }
    }
  }
  if (model_->state() != State::kInitialise)
    mature_individuals_by_year_[model_->current_year()] = mature_conversion_;
}


// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  Maturity::FillReportCache(ostringstream& cache) {
  cache << "year: ";
  for (auto& year : mature_individuals_by_year_)
    cache << year.first << " ";

  cache << "\nagents_mature: ";
  for (auto& year : mature_individuals_by_year_)
    cache << year.second << " ";
  cache << "\n";

}

} /* namespace processes */
} /* namespace niwa */
