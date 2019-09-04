/**
 * @file AgeFrequencyByCell.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 23/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "AgeFrequencyByCell.h"

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
AgeFrequencyByCell::AgeFrequencyByCell(Model* model) : Report(model) {
  model_state_ = State::kExecute;;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years", "", true);
  parameters_.Bind<string>(PARAM_TIME_STEP, &time_step_, "Time Step label", "", "");
}

/**
 * Build our relationships between this object and other objects
 */
void AgeFrequencyByCell::DoBuild() {
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
void AgeFrequencyByCell::DoExecute() {
  LOG_FINE() <<" printing report " << label_;
  cache_ << "*"<< type_ << "[" << label_ << "]\n";
  cache_ << "year: " << model_->current_year() << "\n";
  cache_ << "time_step: " << time_step_ << "\n";
  if (model_->get_sexed()) {
    cache_ << "male-values "<< REPORT_R_DATAFRAME<<"\n";
    cache_ << "row-col";
    for (unsigned i = model_->min_age(); i <=  model_->max_age(); ++i)
      cache_ << " " << i;
    cache_ << "\n";

    vector<float> age_freq_male;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          cache_ << row + 1 << "-" << col + 1;
          cell->get_male_frequency(age_freq_male);
          for(auto& age : age_freq_male)
            cache_ << " " << age;
          cache_ << "\n";
        }
      }
    }
    cache_ << "female-values "<< REPORT_R_DATAFRAME<<"\n";
    cache_ << "row-col";
    for (unsigned i = model_->min_age(); i <=  model_->max_age(); ++i)
      cache_ << " " << i;
    cache_ << "\n";

    vector<float> age_freq_female;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          cache_ << row + 1 << "-" << col + 1;
          cell->get_female_frequency(age_freq_female);
          for(auto& age : age_freq_female)
            cache_ << " " << age;
          cache_ << "\n";
        }
      }
    }
  } else {
    cache_ << "values "<< REPORT_R_DATAFRAME<<"\n";
    cache_ << "row-col";
    for (unsigned i = model_->min_age(); i <=  model_->max_age(); ++i)
      cache_ << " " << i;
    cache_ << "\n";

    vector<float> age_freq;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          cache_ << row + 1 << "-" << col + 1;
          cell->get_age_frequency(age_freq);
          for(auto& age : age_freq)
            cache_ << " " << age;
          cache_ << "\n";
        }
      }
    }
  }
  ready_for_writing_ = true;
  call_number_ = false;
}

} /* namespace reports */
} /* namespace niwa */
