/**
 * @file Likelihood.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 26/4/2019
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This report will take a target likelihood and print all of
 * the values that have been assigned to it.
 */
#ifndef SOURCE_REPORTS_CHILDREN_LIKELIHOOD_H_
#define SOURCE_REPORTS_CHILDREN_LIKELIHOOD_H_

// headers
#include "Reports/Report.h"

// namespaces
namespace niwa {
class Likelihood;
namespace reports {

// classes
class Likelihood : public Report {
public:
  Likelihood(Model* model);
  virtual                     ~Likelihood() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final;
  void                        DoExecute() override final;
  void                        DoReset() override final  { };


private:
  string                      likelihood_label_ = "";
  niwa::Likelihood*           likelihood_ = nullptr;
  bool                        first_run_ = true;

};

} /* namespace reports */
} /* namespace niwa */

#endif /* SOURCE_REPORTS_CHILDREN_LIKELIHOOD_H_ */
