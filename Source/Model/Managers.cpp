/**
 * @file Managers.cpp
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 28/09/2015
 * @section LICENSE
 *
 * Copyright NIWA Science �2015 - www.niwa.co.nz
 *
 */

// headers
#include "Managers.h"

#include "Model/Model.h"
#include "AgeingErrors/Manager.h"
#include "Asserts/Manager.h"
#include "DerivedQuantities/Manager.h"
#include "InitialisationPhases/Manager.h"
#include "Layers/Manager.h"
#include "Likelihoods/Manager.h"
#include "Observations/Manager.h"
#include "PreferenceFunctions/Manager.h"
#include "Processes/Manager.h"
#include "Reports/Manager.h"
#include "Selectivities/Manager.h"
#include "TimeSteps/Manager.h"
#include "TimeVarying/Manager.h"

// namespaces
namespace niwa {

/**
 * Default constructor
 */
Managers::Managers(Model* model) {
  LOG_TRACE();

  model_ = model;

  ageing_error_           = new ageingerrors::Manager();
  assert_                 = new asserts::Manager();
  derived_quantity_       = new derivedquantities::Manager();
  initialisation_phase_   = new initialisationphases::Manager();
  layer_                  = new layers::Manager();
  likelihood_             = new likelihoods::Manager();
  observation_            = new observations::Manager();
  preference_function_    = new preference_functions::Manager();
  process_                = new processes::Manager();
  report_                 = new reports::Manager(model_);
  selectivity_            = new selectivities::Manager();
  time_step_              = new timesteps::Manager();
  time_varying_           = new timevarying::Manager();
}

/**
 * Destructor
 */
Managers::~Managers() {
  delete ageing_error_;
  delete assert_;
  delete derived_quantity_;
  delete initialisation_phase_;
  delete layer_;
  delete likelihood_;
  delete observation_;
  delete preference_function_;
  delete process_;
  delete report_;
  delete selectivity_;
  delete time_step_;
  delete time_varying_;
}

void Managers::Validate() {
  LOG_TRACE();
  ageing_error_->Validate();
  time_step_->Validate(model_);
  initialisation_phase_->Validate();
  assert_->Validate();
  derived_quantity_->Validate();
  layer_->Validate();
  likelihood_->Validate();
  observation_->Validate();
  report_->Validate();
  selectivity_->Validate();
  preference_function_->Validate();
  process_->Validate(model_);
  time_varying_->Validate();
}

void Managers::Build() {
  LOG_TRACE();
  ageing_error_->Build();
  time_step_->Build();
  assert_->Build();
  derived_quantity_->Build();
  likelihood_->Build();
  observation_->Build();
  time_varying_->Build();
  process_->BuildRemainingProcesses();
  preference_function_->Build();
  report_->Build();
  initialisation_phase_->Build(model_);  // This calls report and process() so needs to be built after them
}

// bit of a hack to get around dependencies
void Managers::BuildPreWorldView() {
  LOG_TRACE();
  layer_->Build();
  selectivity_->Build();
  process_->BuildGrowthAndMortalityProcesses();

}


void Managers::Reset() {
  LOG_TRACE();
  ageing_error_->Reset();
  assert_->Reset();
  derived_quantity_->Reset();
  initialisation_phase_->Reset();
  layer_->Reset();
  likelihood_->Reset();
  observation_->Reset();
  preference_function_->Reset();
  process_->Reset();
  report_->Reset();
  selectivity_->Reset();
  time_step_->Reset();
  time_varying_->Reset();
}

} /* namespace niwa */
