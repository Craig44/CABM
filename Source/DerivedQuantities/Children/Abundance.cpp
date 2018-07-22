/**
 * @file Abundance.cpp
 * @author  C.Marsh
 * @date 18/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "Abundance.h"

#include "InitialisationPhases/Manager.h"
#include "TimeSteps/Manager.h"
#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/DoubleCompare.h"

// namespaces
namespace niwa {
namespace derivedquantities {

/**
 * Usual constructor
 */
Abundance::Abundance(Model* model) : DerivedQuantity(model) {
  parameters_.Bind<string>(PARAM_SELECTIVITY, &selectivity_label_, "A label for the selectivity", "");
  parameters_.Bind<string>(PARAM_LAYER_LABEL, &abundance_layer_label_, "A label for the layer that indicates which cells to calculate abundnace in over", "");
}

/**
 * Validate class
 */
void Abundance::DoValidate() {

}


/**
 * Build pointers class
 */
void Abundance::DoBuild() {
  // Build Layers
  LOG_TRACE();
  abundance_layer_ = model_->managers().layer()->GetIntLayer(abundance_layer_label_);
  if (!abundance_layer_) {
    LOG_ERROR_P(PARAM_LAYER_LABEL) << "Could not find the layer '" << abundance_layer_label_ << "', please check there is a @layer defined for this layer and that it is type = 'integer'";
  }

  selectivity_ = model_->managers().selectivity()->GetSelectivity(selectivity_label_);
  if (!selectivity_) {
    LOG_ERROR_P(PARAM_SELECTIVITY) << "Could not find the selectivity '" << selectivity_label_ << "', please check it exists";
  }
}

/**
 * Calculate the cached value to use
 * for any interpolation
 */
void Abundance::PreExecute() {
  LOG_TRACE();
  if (utilities::doublecompare::IsOne(time_step_proportion_))
    return;
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  cache_value_ = 0.0;
  float probability_mature_at_age;
  unsigned time_step_index = model_->managers().time_step()->current_time_step();
  LOG_FINE() << "Time step for calculating biomass = " << time_step_index;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      unsigned val = abundance_layer_->get_value(row, col);
      if (val <= 0)
        continue;
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        for (Agent& agent : cell->agents_) {
          probability_mature_at_age = selectivity_->GetResult(agent.get_age());
          if (rng.chance() <= probability_mature_at_age) {
            ++cache_value_;
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
void Abundance::Execute() {
  LOG_TRACE();

  if (utilities::doublecompare::IsZero(time_step_proportion_))
    return;
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  float value = 0.0;
  float probability_mature_at_age;

  unsigned time_step_index = model_->managers().time_step()->current_time_step();
  LOG_FINE() << "Time step for calculating biomass = " << time_step_index;

  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      unsigned val = abundance_layer_->get_value(row, col);
      if (val <= 0)
        continue;
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        for (Agent& agent : cell->agents_) {
          probability_mature_at_age = selectivity_->GetResult(agent.get_age());
          if (rng.chance() <= probability_mature_at_age) {
            ++value;
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
    LOG_FINEST() << " Pre Exploitation value " <<  cache_value_ << " Post exploitation " << value << " Final value " << values_[model_->current_year()];
  }
}

} /* namespace derivedquantities */
} /* namespace niwa */
