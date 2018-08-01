/**
 * @file DerivedQuantity.cpp
 * @author  C.Marsh
 * @date 1.08/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "DerivedQuantity.h"

#include "Model/Model.h"
#include "DerivedQuantities/Manager.h"

// namespaces
namespace niwa {
namespace asserts {

/**
 * Default constructor
 *
 * Bind any parameters that are allowed to be loaded from the configuration files.
 * Set bounds on registered parameters
 * Register any parameters that can be an estimated or utilised in other run modes (e.g profiling, yields, projections etc)
 * Set some initial values
 *
 * Note: The constructor is parsed to generate Latex for the documentation.
 */
DerivedQuantity::DerivedQuantity(Model* model) : Assert(model) {
  parameters_.Bind<float>(PARAM_VALUE, &value_, "Expected value of the derived quantity for the final year", "");
  parameters_.Bind<string>(PARAM_DERIVED_QUANTITY, &dq_label_, "Label of derived quantity testing", "");
}

/**
 * Build any objects that will need to be utilised by this object.
 * Obtain smart_pointers to any objects that will be used by this object.
 */
void DerivedQuantity::DoBuild() {
  dq_ = model_->managers().derived_quantity()->GetDerivedQuantity(dq_label_);
  if (!dq_)
    LOG_FATAL_P(PARAM_DERIVED_QUANTITY) << "could not find derived quantity " << dq_label_ << " make sure it exists";

  model_->Subscribe(State::kFinalise, this);
}

/**
 * Execute/Run/Process the object.
 */
void DerivedQuantity::Execute() {
  if (abs(value_ - dq_->GetValue(model_->final_year())) > 1e-9)
    LOG_ERROR() << "Assert Failure: derived quantity had actual value " << dq_->GetValue(model_->final_year()) << " when we expected " << value_
        << " with difference: " << abs(value_ - dq_->GetValue(model_->final_year()));
}

} /* namespace asserts */
} /* namespace niwa */
