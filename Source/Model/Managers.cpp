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
#include "Asserts/Manager.h"
#include "DerivedQuantities/Manager.h"
#include "InitialisationPhases/Manager.h"
#include "Layers/Manager.h"
#include "Likelihoods/Manager.h"
#include "Observations/Manager.h"
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

  assert_                 = new asserts::Manager();
  derived_quantity_       = new derivedquantities::Manager();
  initialisation_phase_   = new initialisationphases::Manager();
  layer_                  = new layers::Manager();
  likelihood_             = new likelihoods::Manager();
  observation_            = new observations::Manager();
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

  delete assert_;
  delete derived_quantity_;
  delete initialisation_phase_;
  delete layer_;
  delete likelihood_;
  delete observation_;
  delete process_;
  delete report_;
  delete selectivity_;
  delete time_step_;
  delete time_varying_;
}

void Managers::Validate() {
  LOG_TRACE();
  time_step_->Validate(model_);
  initialisation_phase_->Validate();
  process_->Validate(model_);
  assert_->Validate();
  derived_quantity_->Validate();
  layer_->Validate();
  likelihood_->Validate();
  observation_->Validate();
  report_->Validate();
  selectivity_->Validate();
  time_varying_->Validate();
}

void Managers::Build() {
  LOG_TRACE();
  time_step_->Build();
  initialisation_phase_->Build(model_);
  process_->Build();
  assert_->Build();
  derived_quantity_->Build();
  layer_->Build();
  likelihood_->Build();
  observation_->Build();
  selectivity_->Build();
  time_varying_->Build();
  report_->Build();
}

void Managers::Reset() {
  LOG_TRACE();
  assert_->Reset();
  derived_quantity_->Reset();
  initialisation_phase_->Reset();
  layer_->Reset();
  likelihood_->Reset();
  observation_->Reset();
  process_->Reset();
  report_->Reset();
  selectivity_->Reset();
  time_step_->Reset();
  time_varying_->Reset();
}

} /* namespace niwa */
