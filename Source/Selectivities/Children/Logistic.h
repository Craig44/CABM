/**
 * @file Logistic.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 15/01/2013
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * The time class represents a moment of time.
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef LOGISTIC_H_
#define LOGISTIC_H_

// Headers
#include "Selectivities/Selectivity.h"

// Namespaces
namespace niwa {
namespace selectivities {

/**
 * Class Definition
 */
class Logistic : public niwa::Selectivity {
public:
  // Methods
  explicit Logistic(Model* model);
  virtual                     ~Logistic() = default;
  void                        DoValidate() override final;
  void                        RebuildCache() override final;

protected:
  //Methods

private:
  // Members
  float                      a50_;
  float                      ato95_;
  float                      alpha_;
  unsigned                    low_;
  unsigned                    high_;
};

} /* namespace selectivities */
} /* namespace niwa */
#endif /* LOGISTIC_H_ */
