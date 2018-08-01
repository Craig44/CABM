/**
 * @file DerivedQuantity.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 1/09/2014
 * @section LICENSE
 *
 * Copyright NIWA Science ©2014 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * Used in unit tests, to check whether a run has broken working code.
 */
#ifndef ASSERTS_DERIVED_QUANITY_H_
#define ASSERTS_DERIVED_QUANITY_H_

// headers
#include "Asserts/Assert.h"
#include "DerivedQuantities/DerivedQuantity.h"

// namespaces
namespace niwa {
class DerivedQuantity;
namespace asserts {

// classes
class DerivedQuantity : public niwa::Assert {
public:
  // methods
  DerivedQuantity(Model* model);
  virtual                     ~DerivedQuantity() = default;
  void                        Execute() override final;

protected:
  // methods
  void                        DoValidate() { };
  void                        DoBuild() override final;

private:
  // members
  float                      value_;
  string                     dq_label_;
  niwa::DerivedQuantity* dq_ = nullptr;
};

} /* namespace asserts */
} /* namespace niwa */

#endif /* ASSERTS_DERIVED_QUANITY_H_ */
