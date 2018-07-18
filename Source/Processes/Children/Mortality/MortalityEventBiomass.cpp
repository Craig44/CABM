/**
 * @file MortalityEventBiomass.cpp
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
#include "MortalityEventBiomass.h"

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
MortalityEventBiomass::MortalityEventBiomass(Model* model) : Mortality(model) {
  parameters_.Bind<string>(PARAM_SELECTIVITY, &selectivity_label_, "Selectivity label", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "years to apply the process", "");
  parameters_.Bind<string>(PARAM_CATCH_LAYERS, &catch_layer_label_, "Spatial layer describing catch by cell for each year, there is a one to one link with the year specified, so make sure the order is right", "");
  // TODO add length based selectivity
  // TODO add MLS following the ABM model
  // TODO add discard mortality
  // TODO Tagging events
}

/**
 * Do some initial checks of user supplied parameters.
 */
void MortalityEventBiomass::DoValidate() {
  LOG_TRACE();
  if (years_.size() != catch_layer_label_.size())
    LOG_ERROR_P(PARAM_YEARS) << "you must specify a layer label for each year. You have supplied '" << years_.size() << "' years but '" << catch_layer_label_.size() << "' catch layer labels, please sort this out.";
}

/**
 * DoBuild
 */
void MortalityEventBiomass::DoBuild() {
  LOG_TRACE();
  // Get the layers
  for (auto& label : catch_layer_label_) {
    layers::NumericLayer* temp_layer = nullptr;
    temp_layer = model_->managers().layer()->GetNumericLayer(label);
    if (!temp_layer) {
      LOG_FATAL_P(PARAM_CATCH_LAYERS) << "could not find the layer '" << label << "', please make sure it exists and that it is type 'numeric'";
    }
    catch_layer_.push_back(temp_layer);
  }

  selectivity_ = model_->managers().selectivity()->GetSelectivity(selectivity_label_);
  if (!selectivity_)
    LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << ": selectivity " << selectivity_label_ << " does not exist. Have you defined it?";
}


/**
 * DoExecute
 */
void MortalityEventBiomass::DoExecute() {
  LOG_TRACE();
  auto iter = years_.begin();
  if (find(iter, years_.end(), model_->current_year()) != years_.end()) {
    iter = find(years_.begin(), years_.end(), model_->current_year());
    unsigned catch_ndx = distance(years_.begin(), iter);
    LOG_FINEST() << "applying F in year " << model_->current_year() << " catch index = " << catch_ndx;
    // Get the pointer to the right catch layer
    utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
    float catch_taken;
    float actual_catch_taken = 0;;
    float world_catch_to_take = 0;;

    unsigned catch_attempts = 1;
    unsigned catch_max = 1;
    unsigned random_agent;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          catch_taken = catch_layer_[catch_ndx]->get_value(row, col);
          if (catch_taken > 0) {
            LOG_FINEST() << "We are fishing in cell " << row + 1 << " " << col + 1 << " value = " << catch_taken;
            catch_attempts = 1;
            world_catch_to_take += catch_taken;
            auto& agents = cell->get_agents();
            catch_max = agents.size();
            while (catch_taken > 0) {
              // Random access bullshit for lists
              ++catch_attempts;
              auto iter = agents.begin();
              random_agent = rng.chance() * agents.size();
              advance(iter, random_agent);
              // See if this agent is unlucky
              if (rng.chance() <= selectivity_->GetResult((*iter).get_age())) {
                catch_taken -= (*iter).get_weight() * (*iter).get_scalar();
                actual_catch_taken += (*iter).get_weight() * (*iter).get_scalar();
                agents.erase(iter); // erase agent from memory
              }
              // Make sure we don't end up fishing for infinity
              if (catch_attempts >= catch_max) {
                LOG_FATAL_P(PARAM_LABEL) << "Too many attempts to catch an agent in the process " << label_ << " in year " << model_->current_year() << " in row " << row + 1 << " and column " << col + 1 << " this most likely means you have\n" <<
                   " a model that suggests there should be more agents in this space than than the current agent dynamics are putting in this cell, check the user manual for tips to resolve this situation";
              }
            }
          }
        }
      }
    }
    if (model_->state() != State::kInitialise) {
      removals_by_year_[model_->current_year()] = world_catch_to_take;
      actual_removals_by_year_[model_->current_year()] = actual_catch_taken;
    }
  } // find(years_.begin(), years_.end(), model_->current_year()) != years_.end()
}


// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MortalityEventBiomass::FillReportCache(ostringstream& cache) {
  cache << "year: ";
  for (auto& year : actual_removals_by_year_)
    cache << year.first << " ";

  cache << "\nbiomass_removed: ";
  for (auto& year : actual_removals_by_year_)
    cache << year.second << " ";
  cache << "\ncatch_input_removed: ";
  for (auto& year : removals_by_year_)
    cache << year.second << " ";
  cache << "\n";

}

} /* namespace processes */
} /* namespace niwa */
