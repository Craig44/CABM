/**
 * @file Increasing.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 14/01/2013
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
#ifndef INCREASING_H_
#define INCREASING_H_

// Headers
#include "Selectivities/Selectivity.h"

// Namespaces
namespace niwa {
namespace selectivities {

/**
 * Class definition
 */
class Increasing : public niwa::Selectivity {
public:
  // Methods
  explicit Increasing(Model* model);
  virtual                     ~Increasing() = default;
  void                        DoValidate() override final;
  void                        RebuildCache() override final;

protected:
  //Methods

private:
  // Members
  unsigned                    low_;
  unsigned                    high_;
  vector<double>              v_;
  double                      alpha_;
};

} /* namespace selectivities */
} /* namespace niwa */
#endif /* INCREASING_H_ */
