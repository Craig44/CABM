/**
 * @file InitialisationPartition.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 15/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "InitialisationPartition.h"

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
InitialisationPartition::InitialisationPartition(Model* model) : Report(model) {
  parameters_.Bind<bool>(PARAM_DO_LENGTH_FREQUENCY, &do_length_frequency_, "Print the report as length frequency not age.", "", false);
  model_state_ = State::kInitialise;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
}

/**
 * Build our relationships between this object and other objects
 */
void InitialisationPartition::DoBuild() {
  // Generate a pointer to the world
  world_ = model_->world_view();
  if (!world_)
    LOG_CODE_ERROR() << "could not create the world via reference, something is wrong";

}

/**
 * Execute this report
 */
void InitialisationPartition::DoExecute() {
  LOG_FINE() <<" printing report " << label_;
  if (not do_length_frequency_)
    do_age_ = true;
  else
    do_age_ = false;

  LOG_MEDIUM() << "do length =  " << do_length_frequency_ << " do age = " << do_age_;
  if (call_number_) {
    cache_ << "*"<< type_ << "[" << label_ << "_1]" << "\n";
    cache_ << "equilibrium_shortcut: " << model_->current_year() << "\n";
  } else {
    cache_ << "*"<< type_ << "[" << label_ << "_2]" << "\n";
    cache_ << "year: " << model_->current_year() << "\n";
  }
  cache_ << "values "<< REPORT_R_DATAFRAME_ROW_LABELS <<"\n";
  cache_ << "row-col";
  if (do_length_frequency_) {
    for (unsigned i = 0; i <  model_->length_bin_mid_points().size(); ++i)
      cache_ << " " << model_->length_bin_mid_points()[i];
  } else {
    for (unsigned i = model_->min_age(); i <=  model_->max_age(); ++i)
      cache_ << " " << i;
  }
  cache_ << "\n";

  vector<double> age_freq;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        cache_ << row + 1 << "-" << col + 1;
        cell->get_age_frequency(age_freq, do_age_);
        for(auto& age : age_freq)
          cache_ << " " << age;
        cache_ << "\n";
      }
    }
  }
  ready_for_writing_ = true;
  call_number_ = false;
  LOG_FINE() << "finished printing report = " << label_;
}

} /* namespace reports */
} /* namespace niwa */
