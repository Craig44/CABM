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
  if (call_number_) {
    cache_ << "*"<< type_ << "[" << label_ << "_1]" << "\n";
    cache_ << "ewuilibrium_shortcut: " << model_->current_year() << "\n";
  } else {
    cache_ << "*"<< type_ << "[" << label_ << "_2]" << "\n";
    cache_ << "year: " << model_->current_year() << "\n";
  }
  cache_ << "values "<< REPORT_R_DATAFRAME<<"\n";
  cache_ << "row-col";
  for (unsigned i = model_->min_age(); i <=  model_->max_age(); ++i)
    cache_ << " " << i;
  cache_ << "\n";

  vector<unsigned> age_freq;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        cache_ << row + 1 << "-" << col + 1;
        cell->get_age_frequency(age_freq);
        for(auto& age : age_freq)
          cache_ << " " << age;
      }
    }
    cache_ << "\n";
  }
  ready_for_writing_ = true;
  call_number_ = false;
}

} /* namespace reports */
} /* namespace niwa */
