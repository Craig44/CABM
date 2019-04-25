/**
 * @file Likelihood.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 26/4/2019
 * @section LICENSE
 *
 * Print likelihood specific information, mainly used for logistic normal to view covariance matrix and other things, but could be in handy later down the track
 * THis basically just calls FillReportCache() and the individual likelihood classes are responsible for reporting what ever they want.
 *
 */

// headers
#include "Likelihood.h"

#include <boost/algorithm/string/join.hpp>

#include "Model/Managers.h"
#include "Model/Model.h"
#include "Likelihoods/Manager.h"
#include "Likelihoods/Likelihood.h"

// namespaces
namespace niwa {
namespace reports {

/**
 * Default constructor
 *
 * @param model Pointer to the current model context
 */
Likelihood::Likelihood(Model* model) : Report(model) {
  model_state_ = State::kIterationComplete;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);

  parameters_.Bind<string>(PARAM_LIKELIHOOD, &likelihood_label_, "Likelihood label that is reported", "", "");
}

/**
 * Build our relationships between this object and other objects
 */
void Likelihood::DoBuild() {
  likelihood_ = model_->managers().likelihood()->GetOrCreateLikelihood(model_, "", likelihood_label_);
  if (!likelihood_) {
    LOG_ERROR_P(PARAM_LIKELIHOOD) << "likelihood " << likelihood_label_ << " could not be found. Have you defined it?";
  }
}

/**
 * Execute this report
 */
void Likelihood::DoExecute() {
  LOG_FINE() <<" printing report " << label_;
  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";
  cache_ << "label: " << likelihood_label_ << "\n";

  LOG_FINEST() << "about to print cache";
  likelihood_->FillReportCache(cache_);
  LOG_FINEST() << "finished printing cache";

  ready_for_writing_ = true;
}


} /* namespace reports */
} /* namespace niwa */
