/**
 * @file TagShedding.cpp
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
#include "TagShedding.h"

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
TagShedding::TagShedding(Model* model) : Process(model) {
  parameters_.Bind<string>(PARAM_SELECTIVITY_LABEL, &selectivity_label_, "Label for the selectivity block", "");
  parameters_.Bind<float>(PARAM_TIME_STEP_RATIO, &ratios_, "Time step ratios for the Shedding rates to apply in each time step. See manual for how this is applied", "");

}

/**
 * DoBuild
 */
void TagShedding::DoBuild() {

  LOG_FINEST() << "selectivities supplied = " << selectivity_label_.size();
  // Build selectivity links
  if (selectivity_label_.size() == 1)
    selectivity_label_.assign(2, selectivity_label_[0]);

  if (selectivity_label_.size() > 2) {
    LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << "You suppled " << selectivity_label_.size()  << " Selectiviites, you can only have one for each sex max = 2";
  }
  LOG_FINEST() << "selectivities supplied = " << selectivity_label_.size();

  bool first = true;
  for (auto label : selectivity_label_) {
    Selectivity* temp_selectivity = model_->managers().selectivity()->GetSelectivity(label);
    if (!temp_selectivity)
      LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << ": selectivity " << label << " does not exist. Have you defined it?";
    temp_selectivity->SubscribeToRebuildCache(this);

    selectivity_.push_back(temp_selectivity);
    if (first) {
      first = false;
      selectivity_length_based_ = temp_selectivity->is_length_based();
    } else {
      if (selectivity_length_based_ != temp_selectivity->is_length_based()) {
        LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << "The selectivity  " << label << " was not the same type (age or length based) as the previous selectivity label";
      }
    }
  }

  /**
   * Organise our time step ratios. Each time step can
   * apply a different ratio of M so here we want to verify
   * we have enough and re-scale them to 1.0
   */
  vector<TimeStep*> time_steps = model_->managers().time_step()->ordered_time_steps();
  LOG_FINEST() << "time_steps.size(): " << time_steps.size();
  vector<unsigned> active_time_steps;
  for (unsigned i = 0; i < time_steps.size(); ++i) {
    if (time_steps[i]->HasProcess(label_))
      active_time_steps.push_back(i);
  }

  if (ratios_.size() == 0) {
    for (unsigned i : active_time_steps)
      time_step_ratios_[i] = 1.0;
  } else {
    if (ratios_.size() != active_time_steps.size())
      LOG_ERROR_P(PARAM_TIME_STEP_RATIO) << " length (" << ratios_.size()
          << ") does not match the number of time steps this process has been assigned to (" << active_time_steps.size() << ")";

    for (float value : ratios_) {
      if (value < 0.0 || value > 1.0)
        LOG_ERROR_P(PARAM_TIME_STEP_RATIO) << " value (" << value << ") must be between 0.0 (exclusive) and 1.0 (inclusive)";
    }

    for (unsigned i = 0; i < ratios_.size(); ++i)
      time_step_ratios_[active_time_steps[i]] = ratios_[i];
  }


}

/**
 * DoReset
 */
void TagShedding::DoReset() {
  LOG_FINE() << "clearing containers";
}

/**
 * DoExecute
 */
void TagShedding::DoExecute() {
  LOG_MEDIUM();
  /*
  agents_removed_ = 0;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      // get the ratio to apply first
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        // Apply mortality to elements in a cell.
        ApplyStochasticMortality(cell->agents_);
        ApplyStochasticMortality(cell->tagged_agents_);
      }
    }
  }
  if (model_->state() != State::kInitialise)
    removals_by_year_[model_->current_year()] += agents_removed_;*/
}


// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  TagShedding::FillReportCache(ostringstream& cache) {
	/*
  cache << "year: ";
  for (auto& year : removals_by_year_)
    cache << year.first << " ";

  cache << "\nagents_removed: ";
  for (auto& year : removals_by_year_)
    cache << year.second << " ";
  cache << "\n";*/
}

// A Method for telling the world we need to redistribute Mortality parmaeters
void TagShedding::RebuildCache() {
  LOG_FINE();

}

} /* namespace processes */
} /* namespace niwa */
