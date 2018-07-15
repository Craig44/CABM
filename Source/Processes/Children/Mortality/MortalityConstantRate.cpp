/**
 * @file MortalityConstantRate.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */

// headers
#include "MortalityConstantRate.h"

#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
MortalityConstantRate::MortalityConstantRate(Model* model) : Mortality(model) {
  parameters_.Bind<string>(PARAM_M_LAYER_LABEL, &m_layer_label_, "Label for the numeric layer that describes mean instantaneous mortality rate through space", "");
  parameters_.Bind<string>(PARAM_SELECTIVITY_LABEL, &selectivity_label_, "Label for the selectivity block", "");

}

/**
 * DoBuild
 */
void MortalityConstantRate::DoBuild() {
  // Get the layers
  m_layer_ = model_->managers().layer()->GetNumericLayer(m_layer_label_);
  if (!m_layer_) {
    LOG_ERROR_P(PARAM_M_LAYER_LABEL) << "could not find the layer '" << m_layer_label_ << "', please make sure it exists and that it is type 'numeric'";
  }

  selectivity_ = model_->managers().selectivity()->GetSelectivity(selectivity_label_);
  if (!selectivity_)
    LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << ": selectivity " << selectivity_label_ << " does not exist. Have you defined it?";

}


/**
 * DoExecute
 */
void MortalityConstantRate::DoExecute() {
  // Iterate over all cells
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        auto& agents = cell->get_agents();
        for (auto iter = agents.begin(); iter != agents.end();) {
          (*iter).survival();
          if (not (*iter).is_alive()) {
            iter = agents.erase(iter);
          } else
            ++iter;

        }
      }
    }
  }
}


/*
 * This method is called at when ever an agent is created/seeded or moves. Agents will get a new/updated growth parameters
 * based on the spatial cells of the process. This is called in initialisation/Recruitment and movement processes if needed.
*/
void  MortalityConstantRate::draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<double>& vector) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  double mean_m = m_layer_->get_value(row, col);
  vector.clear();
  vector.resize(number_of_draws);

  for (unsigned i = 0; i < number_of_draws; ++i) {
    double value = rng.lognormal(mean_m, cv_);
    vector[i] = value;
  }

}

} /* namespace processes */
} /* namespace niwa */
