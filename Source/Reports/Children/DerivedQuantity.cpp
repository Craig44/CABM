/*
 * DerivedQuantity.cpp
 *
 *  Created on: 4/09/2013
 *      Author: Admin
 */

#include "DerivedQuantity.h"

#include "DerivedQuantities/Manager.h"

namespace niwa {
namespace reports {


/**
 *
 */
DerivedQuantity::DerivedQuantity(Model* model) : Report(model) {
  //run_mode_    = (RunMode::Type)(RunMode::kBasic | RunMode::kSimulation);
  model_state_ = (State::Type)(State::kIterationComplete);
}

/**
 * Validate inputs
 */
void DerivedQuantity::DoValidate() {
  if ((model_->run_mode() == RunMode::kMSE) || (model_->run_mode() == RunMode::kBasic)|| (model_->run_mode() == RunMode::kSimulation))
    run_mode_ = model_->run_mode();

  if ((model_->run_mode() == RunMode::kMSE))
    model_state_ = State::kInputIterationComplete;
}
/**
 *
 */
void DerivedQuantity::DoExecute() {
  LOG_TRACE();
  derivedquantities::Manager& manager = *model_->managers().derived_quantity();
  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";

  auto derived_quantities = manager.objects();
  for (auto dq : derived_quantities) {
    string label =  dq->label();
    cache_ << label << " " << REPORT_R_LIST <<" \n";
    cache_ << "type: " << dq->type() << " \n";
    if (not dq->is_spatial()) {
      vector<vector<double>> init_values = dq->initialisation_values();
      for (unsigned i = 0; i < init_values.size(); ++i) {
        cache_ << "initialisation_phase["<< i + 1 << "]: ";
        cache_ << init_values[i].back() << " ";
        cache_ << "\n";
      }


      const map<unsigned, double> values = dq->values();
      cache_ << "values " << REPORT_R_VECTOR <<"\n";
      for (auto iter = values.begin(); iter != values.end(); ++iter) {
          double weight = iter->second;
          cache_ << iter->first << " " << weight << "\n";
      }
      //cache_ <<"\n";
      cache_ << REPORT_R_LIST_END <<"\n";
    } else {

    }
  }
  ready_for_writing_ = true;
}

} /* namespace reports */
} /* namespace niwa */
