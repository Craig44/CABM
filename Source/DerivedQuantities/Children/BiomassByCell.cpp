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
  selectivity_ = model_->managers().selectivity()->GetSelectivity(selectivity_label_);
  if (!selectivity_) {
    LOG_ERROR_P(PARAM_SELECTIVITY) << "Could not find the selectivity '" << selectivity_label_ << "', please check it exists";
  }

  // create offset contianers to help with threading
  cell_offset_for_selectivity_.resize(model_->get_height());
  cell_offset_.resize(model_->get_height());
  cache_in_space_.resize(model_->get_height());
  value_in_space_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cache_in_space_[i].resize(model_->get_width());
    value_in_space_[i].resize(model_->get_width());
    cell_offset_[i].resize(model_->get_width());
    cell_offset_for_selectivity_[i].resize(model_->get_width());
  }

  if (selectivity_->is_length_based()) {
    length_based_selectivity_ = true;
    for (unsigned i = 0; i < model_->get_height(); ++i) {
      for (unsigned j = 0; j < model_->get_width(); ++j) {
        for (auto len : model_->length_bins())
          cell_offset_for_selectivity_[i][j].push_back(selectivity_->GetResult(len));
      }
    }
  } else {
    for (unsigned i = 0; i < model_->get_height(); ++i) {
      for (unsigned j = 0; j < model_->get_width(); ++j) {
        for (auto age = model_->min_age(); age <= model_->max_age(); ++age)
          cell_offset_for_selectivity_[i][j].push_back(selectivity_->GetResult(age));
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

  if (not length_based_selectivity_) {
    #pragma omp parallel for collapse(2)
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          unsigned counter = 0;
          cache_in_space_[row][col] = 0.0;
          for (Agent& agent : cell->agents_) {
            if (random_numbers_[cell_offset_[row][col] + counter] <= cell_offset_for_selectivity_[row][col][agent.get_age() - model_->min_age()]) {
              cache_in_space_[row][col] += agent.get_weight() * agent.get_scalar();
            }
            ++counter;
          }
        }
      }
    }
  } else {
    #pragma omp parallel for collapse(2)
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          unsigned counter = 0;
          cache_in_space_[row][col] = 0.0;
          for (Agent& agent : cell->agents_) {
            if (random_numbers_[cell_offset_[row][col] + counter] <= cell_offset_for_selectivity_[row][col][agent.get_length_bin_index()]) {
              cache_in_space_[row][col] += agent.get_weight() * agent.get_scalar();
            }
            ++counter;
          }
        }
      }
    }
  }
  LOG_TRACE();
}
/*
 * PreExecute
*/
void BiomassByCell::Execute() {
  LOG_FINE();
  if (!utilities::doublecompare::IsZero(time_step_proportion_)) {
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

    if (not length_based_selectivity_) {
      #pragma omp parallel for collapse(2)
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            unsigned counter = 0;
            value_in_space_[row][col] = 0.0;
            for (Agent& agent : cell->agents_) {
              if (random_numbers_[cell_offset_[row][col] + counter] <= cell_offset_for_selectivity_[row][col][agent.get_age() - model_->min_age()]) {
                value_in_space_[row][col] += agent.get_weight() * agent.get_scalar();
              }
              ++counter;
            }
          }
        }
      }
    } else {
      #pragma omp parallel for collapse(2)
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            unsigned counter = 0;
            value_in_space_[row][col] = 0.0;
            for (Agent& agent : cell->agents_) {
              if (random_numbers_[cell_offset_[row][col] + counter] <= cell_offset_for_selectivity_[row][col][agent.get_length_bin_index()]) {
                value_in_space_[row][col] += agent.get_weight() * agent.get_scalar();
              }
              ++counter;
            }
          }
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
