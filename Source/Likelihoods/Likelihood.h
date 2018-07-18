/**
 * @file Likelihood.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 22/03/2013
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * The time class represents a moment of time.
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

// Headers
#include "BaseClasses/Object.h"
#include "Observations/Comparison.h"
#include "Utilities/Types.h"

// Namespaces
namespace niwa {
class Model;

/**
 * Class definition
 */
class Likelihood : public niwa::base::Object {
public:
  // Methods
  Likelihood(Model* model);
  virtual                     ~Likelihood() = default;
  void                        Validate();
  void                        Build() { };
  void                        Reset() override final { };
  virtual void                SimulateObserved(map<unsigned, vector<observations::Comparison> >& comparisons) { };
  virtual void                DoValidate() { };

  // accessors
  void                        set_type(const string& type) { type_ = type; }

protected:
  // members
  Model*                      model_ = nullptr;
};
} /* namespace niwa */
#endif /* LIKELIHOOD_H_ */
