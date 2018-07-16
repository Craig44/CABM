/**
 * @file WorldAgeFrequency.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 15/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "WorldAgeFrequency.h"

#include <boost/algorithm/string/join.hpp>

#include "Model/Managers.h"
#include "Model/Model.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"

// namespaces
namespace niwa {
namespace reports {

/**
 * Default constructor
 *
 * @param model Pointer to the current model context
 */
WorldAgeFrequency::WorldAgeFrequency(Model* model) : Report(model) {
  model_state_ = State::kExecute;;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years", "", true);
  parameters_.Bind<string>(PARAM_TIME_STEP, &time_step_, "Time Step label", "", "");

}

/**
 * Build our relationships between this object and other objects
 */
void WorldAgeFrequency::DoBuild() {
  LOG_TRACE();
  // Generate a pointer to the world
  world_ = model_->world_view();
  if (!world_)
    LOG_CODE_ERROR() << "could not create the world via reference, something is wrong";

  if (!parameters_.Get(PARAM_YEARS)->has_been_defined()) {
    years_ = model_->years();
  }
}

/**
 * Execute this report
 */
void WorldAgeFrequency::DoExecute() {
  LOG_FINE() <<" printing report " << label_;

  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";
  cache_ << "year: " << model_->current_year() << "\n";
  cache_ << "time_step: " << time_step_ << "\n";
  cache_ << "values "<< REPORT_R_DATAFRAME<<"\n";
  for (unsigned i = model_->min_age(); i <=  model_->max_age(); ++i)
    cache_ << " " << i;
  cache_ << "\n";

  vector<unsigned> age_freq;
  world_->get_world_age_frequency(age_freq);
  LOG_FINEST() << "size of age freq = " << age_freq.size();

  for(auto& age : age_freq)
    cache_ << " " << age;

  cache_ << "\n";
  ready_for_writing_ = true;

}

} /* namespace reports */
} /* namespace niwa */
