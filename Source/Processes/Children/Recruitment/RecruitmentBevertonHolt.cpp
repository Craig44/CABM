/**
 * @file RecruitmentBevertonHolt.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 12/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This file exists to keep documentation generator happy
 */

// headers
#include "RecruitmentBevertonHolt.h"

#include "DerivedQuantities/Manager.h"
#include "Layers/Manager.h"
#include "TimeSteps/Manager.h"
#include "DerivedQuantities/Manager.h"
#include "InitialisationPhases/Manager.h"
//#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "Utilities/DoubleCompare.h"
#include <omp.h>
// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
RecruitmentBevertonHolt::RecruitmentBevertonHolt(Model* model) : Recruitment(model) {
  parameters_.Bind<float>(PARAM_STEEPNESS, &steepness_, "Steepness", "", 1.0);
  parameters_.Bind<float>(PARAM_YCS_VALUES, &ycs_values_, "YCS (Year-Class-Strength) Values", "");

  RegisterAsAddressable(PARAM_YCS_VALUES, &ycs_values_by_year_);
}

// DoValidate
void RecruitmentBevertonHolt::DoValidate() {
  vector<unsigned> years = model_->years();
  if (model_->get_sexed()) {
    if (!parameters_.Get(PARAM_PROPORTION_MALE)->has_been_defined())
      LOG_FATAL_P(PARAM_LABEL) << "If you have a sexed model you will want to explicitly specify proportion male";
    if (prop_male_.size() == 1)
      prop_male_.assign(years.size(), prop_male_[0]);
    if (prop_male_.size() != years.size())
      LOG_FATAL_P(PARAM_PROPORTION_MALE) << "Either supply a single value to be applied in all years '" << years.size()<< "', or a value for each year of the model, you have supplied '"<< prop_male_.size() << "', please sort this out";
  } else {
    prop_male_.assign(years.size(), 1.0);
  }

  unsigned prop_ndx = 0;
  for (auto& year : years) {
    prop_male_by_year_[year] = prop_male_[prop_ndx];
    ++prop_ndx;
  }
  // plus set the first value for usage in initialisation, the key = 0;
  prop_male_by_year_[0] = prop_male_[0];

  model_->set_male_proportions(prop_male_by_year_);

  // move to validate because ycs_values_by_year_ is addressable and so will be reset for -i runs
  unsigned iter = 0;
  for (unsigned year = model_->start_year(); year <= model_->final_year(); ++year, ++iter) {
    ycs_values_by_year_[year] = ycs_values_[iter];
  }

}


// DoBuild
void RecruitmentBevertonHolt::DoBuild() {
  LOG_FINE();
  if (model_->years().size() != ycs_values_.size()) {
    LOG_ERROR_P(PARAM_YCS_VALUES) << "number of years '" << model_->years().size() <<  "' doesn't match number of ycs values you supplied '" << ycs_values_.size() << "' please make sure these are the same";
  }


  derived_quantity_ = model_->managers().derived_quantity()->GetDerivedQuantity(ssb_label_);
  if (!derived_quantity_)
    LOG_FATAL_P(PARAM_SSB) << "could not find the @derived_quantity block " << ssb_label_ << ", please make sure it exists";

  recruitment_layer_ = model_->managers().layer()->GetNumericLayer(recruitment_layer_label_);
  if (!recruitment_layer_)
    LOG_FATAL_P(PARAM_RECRUITMENT_LAYER_LABEL) << "could not find the @layer block " << recruitment_layer_label_ << ", please make sure it exists, if it does exits make sure it is of type numeric";

  float props = 0.0;
  // Check that recuitment layer is not in conflict with base layer
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      LOG_FINEST() << "check this all works";
      float value = recruitment_layer_->get_value(row, col);
      LOG_FINEST() << "value = " << value;
      if (!cell->is_enabled() & (value > 0))
        LOG_FATAL_P(PARAM_RECRUITMENT_LAYER_LABEL) << "at row " << row + 1 << " and col " << col + 1 << " the base cell was disabled, but you have proportion of recruits in this cell " << value << ", this is not allowed. please check these layers or see manual";
      props += value;
    }
  }
  if ((props - 1) > 0.0001)
    LOG_FATAL_P(PARAM_RECRUITMENT_LAYER_LABEL) << "the recuitment layer does not sum to 1.0 it was " << props << ", we don't want leakage of indiviuals please sort this out";
  model_->set_b0(label_, b0_);

  /**
   * Check order of sequence make sure Spawning happens before recruitment
   */
  unsigned derived_quantity_time_step_index = model_->managers().time_step()->GetTimeStepIndex(derived_quantity_->time_step());
  unsigned recruitment_index = std::numeric_limits<unsigned>::max();
  const vector<TimeStep*> ordered_time_steps = model_->managers().time_step()->ordered_time_steps();
  unsigned derived_quantity_index = std::numeric_limits<unsigned>::max();
  unsigned time_step_index = 0;
  unsigned process_index = 0;
  bool mortailty_block = false;
  // loop through time steps
  for (auto time_step : ordered_time_steps) {
    if (time_step_index == derived_quantity_time_step_index) {
      for (auto process : time_step->processes()) {
        if (process->process_type() == ProcessType::kMortality) {
          mortailty_block = true;
          derived_quantity_index = process_index;
        }
        process_index++;
      }
      LOG_FINEST() << "process_index = " << process_index;
      if (!mortailty_block) {
        //process_index++;
        derived_quantity_index = process_index;
        process_index++;
      }
    } else {
        process_index = time_step->processes().size();

    }
    time_step_index++;
  }
  recruitment_index = model_->managers().time_step()->GetProcessIndex(label_);

  LOG_FINEST() << "recruitment index = " << recruitment_index << " ssb index = " << derived_quantity_index;
  if ((recruitment_index <= derived_quantity_index) & (model_->min_age() < 1))
    LOG_ERROR_P(PARAM_SSB) << "it seems the derived quantity " << ssb_label_ << " occurs after (or in the same time-step) the recruitment process. For models that start with age = 0, this is not allowed. If this doesn't make much sense look at the usermanual under mortality blocks";



}


// DoExecute
void RecruitmentBevertonHolt::DoExecute() {
  LOG_MEDIUM() << label_;

  if (model_->state() == State::kInitialise) {
    LOG_FINEST() << "applying recruitment in initialisation year " << model_->current_year();
    initialisationphases::Manager& init_phase_manager = *model_->managers().initialisation_phase();

    if (first_enter_execute_) {
      initial_recruits_ = model_->get_r0(label_);
      first_enter_execute_ = false;
      // Don't ask for SSB on the first run,
    } else {
      float SSB = derived_quantity_->GetLastValueFromInitialisation(init_phase_manager.last_executed_phase());
      model_->set_ssb(label_, SSB);
      scalar_ = b0_ / SSB;
    }
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          float value = recruitment_layer_->get_value(row, col);
          unsigned new_agents = (unsigned)(initial_recruits_ * value);
          LOG_FINE() << "row = " << row + 1 << " col = " << col + 1 << " prop = " << value << " initial agents = " << initial_recruits_ << " new agents = " << new_agents;
            cell->birth_agents(new_agents, 1.0);
        }
      }
    }
  } else {

    float SSB = derived_quantity_->GetValue(model_->current_year() - model_->min_age());
    if (model_->current_year() - model_->min_age() < model_->start_year())
      SSB *= model_->get_scalar(label_);

    ssb_by_year_[model_->current_year()] = SSB;
    float ssb_ratio = SSB / b0_;
    float SR = ssb_ratio / (1.0 - ((5.0 * steepness_ - 1.0) / (4.0 * steepness_)) * (1.0 - ssb_ratio));
    float true_ycs = ycs_values_by_year_[model_->current_year()] * SR;
    float amount_per = unsigned(initial_recruits_ * true_ycs);
    recruits_by_year_[model_->current_year()] = amount_per;
    ssb_ratio_[model_->current_year()] = ssb_ratio;
    SR_[model_->current_year()] = SR;
    true_ycs_[model_->current_year()] = true_ycs;
    LOG_FINE() << "applying recruitment in year " << model_->current_year() << " SR = " << SR << " ssb_ratio = " << ssb_ratio << " true ycs = " << true_ycs << " intial recruits = " << initial_recruits_;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          float value = recruitment_layer_->get_value(row, col);
          unsigned new_agents = (unsigned)(amount_per * value);
          LOG_FINEST() << "row = " << row + 1 << " col = " << col + 1 << " prop = " << value << " new agents = " << amount_per << " new agents = " << new_agents;

            cell->birth_agents(new_agents, scalar_);

        }
      }
    }
  }
}


// DoExecute
void RecruitmentBevertonHolt::DoReset() {
  LOG_FINE();
  if(model_->is_initialisation_being_re_run())
    first_enter_execute_ = true;

}
// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void RecruitmentBevertonHolt::FillReportCache(ostringstream& cache) {
  LOG_TRACE();
  cache << "initial_recruits: " << initial_recruits_ << "\n";
  cache << "years: ";
  for (auto& iter : recruits_by_year_)
    cache << iter.first << " ";
  cache << "\nrecruits: ";
  for (auto& iter : recruits_by_year_)
    cache << iter.second << " ";
  cache << "\nssb: ";
  for (auto& iter : ssb_by_year_)
    cache << iter.second << " ";
  cache << "\nssb_ratio: ";
  for (auto& iter : ssb_ratio_)
    cache << iter.second << " ";
  cache << "\nSR: ";
  for (auto& iter : SR_)
    cache << iter.second << " ";
  cache << "\ntrue_ycs: ";
  for (auto& iter : true_ycs_)
    cache << iter.second << " ";
  cache << "\nactual_YCS: ";
  for (auto& iter : ycs_values_by_year_)
    cache << iter.second << " ";

  cache << "\n";
}
// Reset containers when accounting for time-varying object changes
void RecruitmentBevertonHolt::RebuildCache() {
  LOG_FINE();
}

} /* namespace processes */
} /* namespace niwa */
