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
//#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "Utilities/DoubleCompare.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
RecruitmentBevertonHolt::RecruitmentBevertonHolt(Model* model) : Process(model) {
  process_type_ = ProcessType::kRecruitment;

  parameters_.Bind<float>(PARAM_B0, &b0_, "B0", "",false);
  parameters_.Bind<float>(PARAM_STEEPNESS, &steepness_, "Steepness", "", 1.0);
  parameters_.Bind<float>(PARAM_YCS_VALUES, &ycs_values_, "YCS (Year-Class-Strength) Values", "");
  parameters_.Bind<string>(PARAM_RECRUITMENT_LAYER_LABEL, &recruitment_layer_label_, "A label for the recruitment layer", "");
  parameters_.Bind<string>(PARAM_SSB, &ssb_label_, "A label for the SSB derived quantity", "");

}

// DoBuild
void RecruitmentBevertonHolt::DoBuild() {
  LOG_TRACE();
  derived_quantity_ = model_->managers().derived_quantity()->GetDerivedQuantity(ssb_label_);
  if (!derived_quantity_)
    LOG_ERROR_P(PARAM_SSB) << "could not find the @derived_quantity block " << ssb_label_ << ", please make sure it exists";

  recruitment_layer_ = model_->managers().layer()->GetNumericLayer(recruitment_layer_label_);
  if (!recruitment_layer_)
    LOG_ERROR_P(PARAM_RECRUITMENT_LAYER_LABEL) << "could not find the @layer block " << recruitment_layer_label_ << ", please make sure it exists, if it does exits make sure it is of type numeric";

  float props = 0.0;
  // Check that recuitment layer is not in conflict with base layer
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      float value = recruitment_layer_->get_value(row, col);
      if (!cell->is_enabled() & (value > 0))
        LOG_FATAL_P(PARAM_RECRUITMENT_LAYER_LABEL) << "at row " << row + 1 << " and col " << col + 1 << " the base cell was disabled, but you have proportion of recruits in this cell, this is not allowed. please check these layers or see manual";
      props += value;
    }
  }
  if (!utilities::doublecompare::IsOne(props))
    LOG_ERROR_P(PARAM_RECRUITMENT_LAYER_LABEL) << "the recuitment layer does not sum to 1.0 it was " << props << ", we don't want leakage of indiviuals please sort this out";

  model_->set_b0(label_, b0_);

}


// DoExecute
void RecruitmentBevertonHolt::DoExecute() {
  LOG_TRACE();

  if (first_enter_execute_) {
    initial_recruits_ = model_->get_r0(label_);
    first_enter_execute_ = false;
  }
  if (model_->state() == State::kInitialise) {
    LOG_FINEST() << "applying recruitment in initialisation";
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {


        }
      }
    }
  } else {
    LOG_FINEST() << "applying recruitment in year " << model_->current_year();

  }
}

// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void RecruitmentBevertonHolt::FillReportCache(ostringstream& cache) {
  cache << "initial_recruits: " << initial_recruits_ << "\n";
  cache << "years: ";
  for (auto& iter : recruits_by_year_)
    cache << iter.first << " ";
  cache << "\nrecruits: ";
  for (auto& iter : recruits_by_year_)
    cache << iter.second << " ";
  cache << "\n";
}

} /* namespace processes */
} /* namespace niwa */
