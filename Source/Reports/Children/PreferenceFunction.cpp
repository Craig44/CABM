/**
 * @file PreferenceFunction.cpp
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 12/10/2015
 * @section LICENSE
 *
 * Copyright NIWA Science ©2015 - www.niwa.co.nz
 *
 */
#include "Model/Model.h"
#include "PreferenceFunction.h"
#include "PreferenceFunctions/Manager.h"

namespace niwa {
namespace reports {

PreferenceFunction::PreferenceFunction(Model* model) : Report(model) {
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
  model_state_ = (State::Type)(State::kIterationComplete);

  parameters_.Bind<string>(PARAM_PREFERENCE_FUNCTION_LABEL, &pref_fun_label_, "Preference Function label", "");
  parameters_.Bind<float>(PARAM_PREFERENCE_VALUES, &pref_vals_, "Preference values to report the preference function", "");

}

void PreferenceFunction::DoValidate() {

}

void PreferenceFunction::DoBuild() {
  preference_function_ = model_->managers().preference_function()->GetPreferenceFunction(pref_fun_label_);
  if (!preference_function_)
    LOG_FATAL_P(PARAM_PREFERENCE_FUNCTION_LABEL) << " " << pref_fun_label_ << " does not exist. Have you defined it?";

}


void PreferenceFunction::DoExecute() {
  LOG_FINE();
  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";
  const map<string, Parameter*> parameters = preference_function_->parameters().parameters();

  for (auto iter : parameters) {
    Parameter* x = iter.second;
    cache_  << iter.first << ": ";

    vector<string> values = x->current_values();
    for (string value : values)
      cache_ << value << " ";
    cache_ << "\n";
  }

  cache_ << "Values " << REPORT_R_VECTOR << "\n";
  for (unsigned i = 0; i < pref_vals_.size(); ++i)
    cache_ << pref_vals_[i] << " " << preference_function_->get_result(pref_vals_[i] ) << "\n";
   ready_for_writing_ = true;

}


} /* namespace reports */
} /* namespace niwa */
