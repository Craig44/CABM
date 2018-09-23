/**
 * @file MatureBiomass.cpp
 * @author  C.Marsh
 * @date 15/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "MatureBiomass.h"

#include "InitialisationPhases/Manager.h"
#include "TimeSteps/Manager.h"
#include "Layers/Manager.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "Utilities/DoubleCompare.h"

// namespaces
namespace niwa {
namespace derivedquantities {

/**
 * Usual constructor
 */
MatureBiomass::MatureBiomass(Model* model) : DerivedQuantity(model) {
  parameters_.Bind<string>(PARAM_BIOMASS_LAYER_LABEL, &biomass_layer_label_, "A label for the layer that indicates which cells to calculate biomass over", "");
}

/**
 * Validate class
 */
void MatureBiomass::DoValidate() {

}


/**
 * Build pointers class
 */
void MatureBiomass::DoBuild() {
  // Build Layers
  LOG_TRACE();
  biomass_layer_ = model_->managers().layer()->GetIntLayer(biomass_layer_label_);
  if (!biomass_layer_) {
    LOG_ERROR_P(PARAM_BIOMASS_LAYER_LABEL) << "Could not find the layer '" << biomass_layer_label_ << "', please check there is a @layer defined for this layer and that it is type = 'integer'";
  }
}

/**
 * Calculate the cached value to use
 * for any interpolation
 */
void MatureBiomass::PreExecute() {
  LOG_FINE();
  if (utilities::doublecompare::IsOne(time_step_proportion_))
    return;
  cache_value_ = 0.0;
  unsigned time_step_index = model_->managers().time_step()->current_time_step();
  LOG_FINE() << "Time step for calculating biomass = " << time_step_index;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      unsigned val = biomass_layer_->get_value(row, col);
      if (val <= 0)
        continue;
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        cache_value_ += cell->get_mature_biomass();
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
void MatureBiomass::Execute() {
  LOG_FINE() << "DQ " << label_;
  float value = 0.0;
  if (!utilities::doublecompare::IsZero(time_step_proportion_)) {
    unsigned time_step_index = model_->managers().time_step()->current_time_step();
    LOG_FINE() << "Time step for calculating biomass = " << time_step_index;

    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        unsigned val = biomass_layer_->get_value(row, col);
        if (val <= 0)
          continue;
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          value += cell->get_mature_biomass();
        }
      }
    }
    LOG_FINEST() << "executing and value = " << value;
  }

  if (model_->state() == State::kInitialise) {
    unsigned initialisation_phase = model_->managers().initialisation_phase()->current_initialisation_phase();
    if (initialisation_values_.size() <= initialisation_phase) {
      LOG_FINE() << "resizing initialisation_values_";
      initialisation_values_.resize(initialisation_phase + 1);
    }

    float b0_value = 0;

    if (utilities::doublecompare::IsZero(time_step_proportion_)) {
      b0_value = cache_value_;
      initialisation_values_[initialisation_phase].push_back(b0_value);
    } else if (utilities::doublecompare::IsOne(time_step_proportion_)) {
      b0_value = value;
      initialisation_values_[initialisation_phase].push_back(b0_value);
    } else {
      b0_value = cache_value_ + ((value - cache_value_) * time_step_proportion_);
      initialisation_values_[initialisation_phase].push_back(b0_value);
    }
    LOG_FINEST() << " Pre Exploitation value " <<  cache_value_ << " Post exploitation " << value << " Final value " << *initialisation_values_[initialisation_phase].rbegin();
  } else {
    if (utilities::doublecompare::IsZero(time_step_proportion_))
      values_[model_->current_year()] = cache_value_;
    else if (utilities::doublecompare::IsOne(time_step_proportion_))
      values_[model_->current_year()] = value;
    else
      values_[model_->current_year()] = cache_value_ + ((value - cache_value_) * time_step_proportion_);

    LOG_FINEST() << " Pre Exploitation value " <<  cache_value_ << " Post exploitation " << value << " Final value " << values_[model_->current_year()];

  }

}

} /* namespace derivedquantities */
} /* namespace niwa */
