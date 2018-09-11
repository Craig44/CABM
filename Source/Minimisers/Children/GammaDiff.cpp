//============================================================================
// Name        : GammaDiff.cpp
// Author      : S.Rasmussen
// Date        : 8/09/2008
// Copyright   : Copyright NIWA Science ©2008 - www.niwa.co.nz
// Description :
// $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
//============================================================================
// Local Headers
#include "GammaDiff.h"

#include "GammaDiff/Callback.h"
#include "GammaDiff/CallbackBaranov.h"
#include "GammaDiff/Engine.h"

// namespaces
namespace niwa {
namespace minimisers {

/**
 * Default constructor
 */
GammaDiff::GammaDiff(Model* model) : Minimiser(model) {
  parameters_.Bind<int>(PARAM_MAX_ITERATIONS, &max_iterations_, "Maximum number of iterations", "", 1000);
  parameters_.Bind<int>(PARAM_MAX_EVALUATIONS, &max_evaluations_, "Maximum number of evaluations", "", 4000);
  parameters_.Bind<double>(PARAM_TOLERANCE, &gradient_tolerance_, "Tolerance of the gradient for convergence", "", 0.02);
  parameters_.Bind<double>(PARAM_STEP_SIZE, &step_size_, "Minimum Step-size before minimisation fails", "", 1e-7);
}

/**
 * Execute the minimiser to solve the model
 */
void GammaDiff::Execute() {
  LOG_FINE();
  // Variables

  gammadiff::CallBack  call_back(model_);

  vector<double>  lower_bounds = {0.0001,0.0001};
  vector<double>  upper_bounds = {100,100};
  vector<double>  start_values = {5,0.1};

  LOG_FINE() << "Launching minimiser";
  int status = 0;
  gammadiff::Engine clGammaDiff;
  clGammaDiff.optimise_finite_differences(call_back,
      start_values, lower_bounds, upper_bounds,
      status, max_iterations_, max_evaluations_, gradient_tolerance_,
      hessian_,1,step_size_);

}

/**
 * Execute the minimiser to solve the Baranov equation
 */
void GammaDiff::SolveBaranov() {
  LOG_FINE() << "Solve Baranov called";
  gammadiff::CallbackBaranov  call_back_baranov(mortality_process_);
  LOG_FINE() << "check callback";
  vector<double> start_val = {0.003};
  vector<double> lower_bound = {0.00000001};
  vector<double> upper_bound = {1000};

  LOG_FINE() << "Launching minimiser";
  int status = 0;
  gammadiff::Engine clGammaDiff;
  clGammaDiff.optimise_finite_differences(call_back_baranov,
      start_val, lower_bound, upper_bound,
      status, max_iterations_, max_evaluations_, gradient_tolerance_,
      hessian_,1,step_size_);

}
} /* namespace minimisers */
} /* namespace niwa */
