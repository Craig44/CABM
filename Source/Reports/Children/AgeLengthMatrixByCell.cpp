/**
 * @file AgeLengthMatrixByCell.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 23/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "AgeLengthMatrixByCell.h"

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
AgeLengthMatrixByCell::AgeLengthMatrixByCell(Model* model) : Report(model) {
  model_state_ = State::kExecute;;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years", "", true);
  parameters_.Bind<string>(PARAM_TIME_STEP, &time_step_, "Time Step label", "", "");
}

/**
 * Build our relationships between this object and other objects
 */
void AgeLengthMatrixByCell::DoBuild() {
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
void AgeLengthMatrixByCell::DoExecute() {
  LOG_FINE() <<" printing report " << label_;

  cache_ << "*"<< type_ << "[" << label_ << "]\n";
  cache_ << "year: " << model_->current_year() << "\n";
  cache_ << "time_step: " << time_step_ << "\n";
  vector<vector<double>> age_length_mat;
  age_length_mat.resize(model_->age_spread());
  for(unsigned age_ndx = 0; age_ndx < age_length_mat.size(); ++age_ndx)
    age_length_mat[age_ndx].resize(model_->number_of_length_bins());

  if (model_->get_sexed()) {
    LOG_WARNING() << "need to write report for seced model";
  } else {
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          for(unsigned age_ndx = 0; age_ndx < age_length_mat.size(); ++age_ndx)
            fill(age_length_mat[age_ndx].begin(), age_length_mat[age_ndx].end(),0.0);
          for (auto &agent : cell->agents_) {
            if (agent.is_alive())
              age_length_mat[agent.get_age_index()][agent.get_length_bin_index()] += agent.get_scalar();
          }
          cache_ << row + 1 << "-"<< col + 1 << " " << REPORT_R_MATRIX<<"\n";
          for(unsigned age_ndx = 0; age_ndx < age_length_mat.size(); ++age_ndx) {
            for(unsigned len_ndx = 0; len_ndx < age_length_mat[age_ndx].size(); ++len_ndx)
              cache_ << age_length_mat[age_ndx][len_ndx] << " ";
            cache_ << "\n";
          }
        }
      }
    }
  }

  ready_for_writing_ = true;
  call_number_ = false;
}

} /* namespace reports */
} /* namespace niwa */
