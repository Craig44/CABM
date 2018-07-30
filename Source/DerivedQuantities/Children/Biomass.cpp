/**
 * @file Biomass.cpp
 * @author  C.Marsh
 * @date 18/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "Biomass.h"

#include "InitialisationPhases/Manager.h"
#include "TimeSteps/Manager.h"
#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/DoubleCompare.h"
#include "omp.h"

// namespaces
namespace niwa {
namespace derivedquantities {

/**
 * Usual constructor
 */
Biomass::Biomass(Model* model) : DerivedQuantity(model) {
  parameters_.Bind<string>(PARAM_SELECTIVITY, &selectivity_label_, "A label for the selectivity", "");
  parameters_.Bind<string>(PARAM_BIOMASS_LAYER_LABEL, &biomass_layer_label_, "A label for the layer that indicates which cells to calculate biomass over", "");
}

/**
 * Validate class
 */
void Biomass::DoValidate() {

}


/**
 * Build pointers class
 */
void Biomass::DoBuild() {
  // Build Layers
  LOG_TRACE();
  biomass_layer_ = model_->managers().layer()->GetIntLayer(biomass_layer_label_);
  if (!biomass_layer_) {
    LOG_ERROR_P(PARAM_BIOMASS_LAYER_LABEL) << "Could not find the layer '" << biomass_layer_label_ << "', please check there is a @layer defined for this layer and that it is type = 'integer'";
  }

  selectivity_ = model_->managers().selectivity()->GetSelectivity(selectivity_label_);
  if (!selectivity_) {
    LOG_ERROR_P(PARAM_SELECTIVITY) << "Could not find the selectivity '" << selectivity_label_ << "', please check it exists";
  }

  // create offset contianers to help with threading
  cell_offset_for_selectivity_.resize(model_->get_height());
  cell_offset_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
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

/**
 * Calculate the cached value to use
 * for any interpolation
 */
void Biomass::PreExecute() {
  LOG_TRACE();
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

  unsigned time_step_index = model_->managers().time_step()->current_time_step();
  LOG_FINE() << "Time step for calculating biomass = " << time_step_index;

  cache_value_ = 0.0;
  if (not length_based_selectivity_) {
    #pragma omp parallel for collapse(2)
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        unsigned val = biomass_layer_->get_value(row, col);
        if (val <= 0)
          continue;
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          unsigned counter = 0;
          for (Agent& agent : cell->agents_) {
            if (random_numbers_[cell_offset_[row][col] + counter] <= cell_offset_for_selectivity_[row][col][agent.get_age() - model_->min_age()]) {
              cache_value_ += agent.get_weight() * agent.get_scalar();
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
        unsigned val = biomass_layer_->get_value(row, col);
        if (val <= 0)
          continue;
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          unsigned counter = 0;
          for (Agent& agent : cell->agents_) {
            if (random_numbers_[cell_offset_[row][col] + counter] <= cell_offset_for_selectivity_[row][col][agent.get_length_bin_index()]) {
              cache_value_ += agent.get_weight() * agent.get_scalar();
            }
            ++counter;
          }
        }
      }
    }
  }
  LOG_TRACE();
}

/**
 * Calculate the derived quantity value for the
 * state of the model.
 *
 * This class will calculate a value that is the sum total
 * of the population in the model filtered by category and
 * multiplied by the selectivities.
 *
 */
void Biomass::Execute() {
  LOG_TRACE();
  float value = 0.0;

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

    unsigned time_step_index = model_->managers().time_step()->current_time_step();
    LOG_FINE() << "Time step for calculating biomass = " << time_step_index;

    if (not length_based_selectivity_) {
      #pragma omp parallel for collapse(2)
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          unsigned val = biomass_layer_->get_value(row, col);
          if (val <= 0)
            continue;
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            unsigned counter = 0;
            for (Agent& agent : cell->agents_) {
              if (random_numbers_[cell_offset_[row][col] + counter] <= cell_offset_for_selectivity_[row][col][agent.get_age() - model_->min_age()]) {
                value += agent.get_weight() * agent.get_scalar();
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
          unsigned val = biomass_layer_->get_value(row, col);
          if (val <= 0)
            continue;
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            unsigned counter = 0;
            for (Agent& agent : cell->agents_) {
              if (random_numbers_[cell_offset_[row][col] + counter] <= cell_offset_for_selectivity_[row][col][agent.get_length_bin_index()]) {
                value += agent.get_weight() * agent.get_scalar();
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
    if (initialisation_values_.size() <= initialisation_phase)
      initialisation_values_.resize(initialisation_phase + 1);

    float b0_value = 0;

    if (utilities::doublecompare::IsZero(time_step_proportion_)) {
      b0_value = cache_value_;
      initialisation_values_[initialisation_phase].push_back(b0_value);
    } else if (utilities::doublecompare::IsOne(time_step_proportion_)) {
      b0_value = value;
      initialisation_values_[initialisation_phase].push_back(b0_value);
    } else if (mean_proportion_method_) {
      b0_value = cache_value_ + ((value - cache_value_) * time_step_proportion_);
      initialisation_values_[initialisation_phase].push_back(b0_value);
    } else {
      b0_value = pow(cache_value_, 1 - time_step_proportion_) * pow(value ,time_step_proportion_);
      initialisation_values_[initialisation_phase].push_back(b0_value);
    }
  } else {
    if (utilities::doublecompare::IsZero(time_step_proportion_))
      values_[model_->current_year()] = cache_value_;
    else if (utilities::doublecompare::IsOne(time_step_proportion_))
      values_[model_->current_year()] = value;
    if (mean_proportion_method_)
      values_[model_->current_year()] = cache_value_ + ((value - cache_value_) * time_step_proportion_);
    else
      values_[model_->current_year()] = pow(cache_value_, 1 - time_step_proportion_) * pow(value ,time_step_proportion_);
    LOG_FINEST() << " Pre Exploitation value " <<  cache_value_ << " Post exploitation " << value << " Final value ";
  }
}

} /* namespace derivedquantities */
} /* namespace niwa */
