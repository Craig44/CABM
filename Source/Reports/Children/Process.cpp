/**
 * @file Process.cpp
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 19/11/2015
 * @section LICENSE
 *
 * Copyright NIWA Science ©2015 - www.niwa.co.nz
 *
 */

// headers
#include "Process.h"

#include <boost/algorithm/string/join.hpp>

#include "Model/Managers.h"
#include "Model/Model.h"
#include "Processes/Manager.h"
#include "Processes/Process.h"

// namespaces
namespace niwa {
namespace reports {

/**
 * Default constructor
 *
 * @param model Pointer to the current model context
 */
Process::Process(Model* model) : Report(model) {
  model_state_ = State::kIterationComplete;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);

  parameters_.Bind<string>(PARAM_PROCESS, &process_label_, "Process label that is reported", "", "");
}

/**
 * Build our relationships between this object and other objects
 */
void Process::DoBuild() {
  process_ = model_->managers().process()->GetProcess(process_label_);
  if (!process_) {
    LOG_ERROR_P(PARAM_PROCESS) << "process " << process_label_ << " could not be found. Have you defined it?";
  }
}

/**
 * Execute this report
 */
void Process::DoExecute() {
  LOG_FINE() <<" printing report " << label_ << " of type " << process_->type();
  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";

  auto parameters = process_->parameters().parameters();
  for (auto parameter : parameters) {
    cache_ << parameter.first << ": ";
    string line = boost::algorithm::join(parameter.second->current_values(), " ");
    cache_ << line << "\n";
  }
  LOG_FINEST() << "about to print cache";
  process_->FillReportCache(cache_);
  LOG_FINEST() << "finished printing cache";

  ready_for_writing_ = true;
}


} /* namespace reports */
} /* namespace niwa */
