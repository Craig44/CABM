/**
 * @file RandomWalk.h
 * @author Craig Marsh
 * @github https://github.com/Zaita
 * @date 2/02/2016
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This child will add a deviate drawn from a normal distribution to the existing parameter value
 */
#ifndef TIMEVARYING_RANDOM_DRAW_H_
#define TIMEVARYING_RANDOM_DRAW_H_

// headers
#include "TimeVarying/TimeVarying.h"
#include "Utilities/Distribution.h"

// namespaces
namespace niwa {
class Estimate;
namespace timevarying {

/**
 * Class definition
 */
class RandomDraw : public TimeVarying {
public:
  explicit RandomDraw(Model* model);
  virtual                     ~RandomDraw() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final;
  void                        DoUpdate() override final;

private:
  // members
  double                      mu_;
  double                      sigma_;
  string                     distribution_label_;
  double                      lower_bound_;
  double                      upper_bound_;
  Distribution               distribution_ = Distribution::kNormal;

};

} /* namespace timevarying */
} /* namespace niwa */

#endif /* TIMEVARYING_RANDOM_DRAW_H_ */
