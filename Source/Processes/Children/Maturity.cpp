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
    selectivity_.push_back(temp_selectivity);
  }

  cell_offset_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i)
    cell_offset_[i].resize(model_->get_width());
}


/**
 * Execute process
 */
void Maturity::DoExecute() {
  LOG_TRACE();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  // Pre-calculate agents in the world to set aside our random numbers needed for the operation
  n_agents_ = 0;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        cell_offset_[row][col] = n_agents_;
        n_agents_ += cell->agents_.size();
      }
    }
  }
  // Allocate a single block of memory rather than each thread temporarily allocating their own memory.
  random_numbers_.resize(n_agents_);
  for (unsigned i = 0; i < n_agents_; ++i)
    random_numbers_[i] = rng.chance();

  unsigned mature_conversion = 0;
  // Iterate over all cells
  float probability_mature_at_age;
  #pragma omp parallel for collapse(2)
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        LOG_FINEST() << "about to convert " << cell->agents_.size() << " through the maturity process";
        unsigned counter = 0;
        for (Agent& agent : cell->agents_) {
          if (not agent.get_maturity()) {
            probability_mature_at_age = selectivity_[agent.get_sex()]->GetResult(agent.get_age());
            if (random_numbers_[cell_offset_[row][col] + counter] <= probability_mature_at_age) {
              agent.set_maturity(true);
              ++mature_conversion;
            }
          }
          ++counter;
        }
      }
    }
  }

  if (model_->state() != State::kInitialise)
    mature_individuals_by_year_[model_->current_year()] = mature_conversion;
  LOG_TRACE();
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
