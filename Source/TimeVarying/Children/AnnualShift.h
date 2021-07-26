/**
 * @file AnnualShift.h
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 2/02/2015
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * << Add Description >>
 */
#ifndef TIMEVARYING_ANNUALSHIFT_H_
#define TIMEVARYING_ANNUALSHIFT_H_

// headers
#include "TimeVarying/TimeVarying.h"

// namespaces
namespace niwa {
namespace timevarying {

/**
 * Class definition
 */
class AnnualShift : public TimeVarying {
public:
  explicit AnnualShift(Model* model);
  virtual                     ~AnnualShift() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final;
  void                        DoUpdate() override final;

private:
  // members
  vector<double>              values_;
  double                      a_;
  double                      b_;
  double                      c_;
  vector<unsigned>            scaling_years_;
  map<unsigned, double>       values_by_year_;
};

} /* namespace timevarying */
} /* namespace niwa */

#endif /* TIMEVARYING_ANNUALSHIFT_H_ */
