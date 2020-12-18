/*
 * PreferenceFunctions.h
 *
 *  Created on: 25/06/2014
 *      Author: Admin
 */

#ifndef REPORTS_PREFERENCE_FUNCTION_H_
#define REPORTS_PREFERENCE_FUNCTION_H_

#include "Reports/Report.h"

namespace niwa {
class PreferenceFunction;

namespace reports {

class PreferenceFunction : public niwa::Report {
public:
  PreferenceFunction(Model* model);
  virtual ~PreferenceFunction() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoPrepare() override final { };
  void                        DoExecute() override final;
  void                        DoFinalise() override final { };
  void                        DoReset() override final  { };

private:
  string                      pref_fun_label_;
  niwa::PreferenceFunction*   preference_function_ = nullptr;
  bool												first_run_ = true;
  vector<float>               pref_vals_;
};

} /* namespace reports */
} /* namespace niwa */

#endif /* REPORTS_PREFERENCE_FUNCTION_H_ */
