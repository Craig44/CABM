/**
 * @file MortalityCull.cpp
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
#include "MortalityCull.h"

#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "TimeSteps/Manager.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
MortalityCull::MortalityCull(Model* model) : Mortality(model) {
  parameters_.Bind<string>(PARAM_LAYER_LABEL, &layer_label_, "Label for the integer layer which we will kill all agents in", ""); // TODO perhaps as a multiplier, 1.2 * 0.2 = 0.24
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "The years of the observed values", "");

}

/**
 * DoBuild
 */
void MortalityCull::DoBuild() {
  // Get the layers
  layer_ = model_->managers().layer()->GetIntLayer(layer_label_);
  if (!layer_) {
    LOG_ERROR_P(PARAM_LAYER_LABEL) << "could not find the layer '" << layer_label_ << "', please make sure it exists and that it is type 'integer'";
  }
 }

/**
 * DoReset
 */
void MortalityCull::DoReset() {
  LOG_FINE() << "clearing containers";

}

/**
 * The main function of this class. pulled out of DoExcute so that I can apply it to many different vectors.
 */
void MortalityCull::ApplyStochasticMortality(vector<Agent>& agents) {

}
/**
 * DoExecute
 */
void MortalityCull::DoExecute() {
  LOG_MEDIUM() << label_;
  unsigned year = model_->current_year();
  if (find(years_.begin(), years_.end(), year) != years_.end()) {
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        // get the ratio to apply first
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          double agents_removed = 0.0;
          if(layer_->get_value(row, col) > 0) {
            for(auto& agent : cell->agents_) {
              if (agent.is_alive()) {
                agent.dies();
                agents_removed += agent.get_scalar();
              }
            }
            for(auto& agent : cell->tagged_agents_) {
              if (agent.is_alive()) {
                agent.dies();
                agents_removed += agent.get_scalar();
              }
            }
            LOG_MEDIUM() << "year = " << year << " row = "<< row << "  col = " << col << " killed agents = " << agents_removed;
          }
        }
      }
    }
  }
}


/*
 * This method is called at when ever an agent is created/seeded or moves. Agents will get a new/updated growth parameters
 * based on the spatial cells of the process. This is called in initialisation/Recruitment and movement processes if needed.
*/
void  MortalityCull::draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) {

}

// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MortalityCull::FillReportCache(ostringstream& cache) {

}

// A Method for telling the world we need to redistribute Mortality parmaeters
void MortalityCull::RebuildCache() {

}

} /* namespace processes */
} /* namespace niwa */
