/**
 * @file TaggingInfo.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 23/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "TaggingInfo.h"

#include <boost/algorithm/string/join.hpp>

#include "Model/Managers.h"
#include "Model/Model.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
// namespaces
namespace niwa {
namespace reports {

/**
 * Default constructor
 *
 * @param model Pointer to the current model context
 */
TaggingInfo::TaggingInfo(Model* model) : Report(model) {
  model_state_ = State::kExecute;;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years", "", true);
  parameters_.Bind<string>(PARAM_TIME_STEP, &time_step_, "Time Step label", "", "");
  parameters_.Bind<unsigned>(PARAM_RELEASE_YEAR, &release_year_, "Release year", "", "");
  parameters_.Bind<string>(PARAM_RELEASE_CELL, &release_region_, "release cell, for example row 1 and cell 1 = '1-1'", "", "");

}

/**
 * Build our relationships between this object and other objects
 */
void TaggingInfo::DoBuild() {
  // Generate a pointer to the world
  world_ = model_->world_view();
  if (!world_)
    LOG_CODE_ERROR() << "could not create the world via reference, something is wrong";

  if (!parameters_.Get(PARAM_YEARS)->has_been_defined()) {
    years_ = model_->years();
  }

  
  // Split out cell origins for quick look up at execution
  for (auto origin : release_region_) {
    unsigned value = 1;
    vector<string> split_cells;
    boost::split(split_cells, origin, boost::is_any_of("-"));

    LOG_FINE() << "row = " << split_cells[0] << " col = " << split_cells[1];
    if (!utilities::To<unsigned>(split_cells[0], value))
      LOG_ERROR_P(PARAM_RELEASE_CELL) << " value (" << split_cells[0] << ") could not be converted to a unsigned";
    origin_rows_.push_back(value - 1);
    if (!utilities::To<unsigned>(split_cells[1], value))
      LOG_ERROR_P(PARAM_RELEASE_CELL) << " value (" << split_cells[1] << ") could not be converted to a unsigned";
    origin_cols_.push_back(value - 1);
  }
  age_freq_.resize(model_->age_spread(),0.0);
}

/**
 * Execute this report
 */
void TaggingInfo::DoExecute() {
  LOG_FINE() << " printing report " << label_;
  cache_ << "*" << type_ << "[" << label_ << "]\n";
  cache_ << "year: " << model_->current_year() << "\n";
  cache_ << "time_step: " << time_step_ << "\n";
  
  if (model_->get_sexed()) {

  }  else {
    for (auto rel_year : release_year_) {
      for (unsigned reg_counter = 0; reg_counter < release_region_.size(); ++reg_counter) {
        LOG_FINE() << "printing infor for tagged  agents released in year " << rel_year << " and cell = " << origin_rows_[reg_counter] << "-" << origin_cols_[reg_counter] << " (" << release_region_[reg_counter] << ")";
        cache_ << "age_comp-" << rel_year << "-" <<  release_region_[reg_counter] << " " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
        cache_ << "cell";
        for (unsigned i = model_->min_age(); i <=  model_->max_age(); ++i)
          cache_ << " " << i;
        cache_ << "\n";
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            LOG_FINE() << "row = " << row << " col = " << col;
            fill(age_freq_.begin(), age_freq_.end(), 0.0);
            WorldCell *cell = world_->get_base_square(row, col);
            if (cell->is_enabled()) {
              cache_ << row + 1 << "-" << col + 1 << " ";
              unsigned counter = 0;
              LOG_FINE() << "tagged agents in this cell (not all will be alive) " << cell->tagged_agents_.size();
              for (auto iter = cell->tagged_agents_.begin(); iter != cell->tagged_agents_.end(); ++iter) {
                if (((*iter).is_alive()) & ((*iter).get_tag_release_year() == rel_year)) {
                  if(((*iter).get_tag_row() == origin_rows_[reg_counter]) & ((*iter).get_tag_col() == origin_cols_[reg_counter])) {
                    // count this tagged agent.
                    age_freq_[(*iter).get_age_index()]++;
                    counter++;
                  }
                }
              }
              LOG_FINE() << "reported on " << counter << " tagged agents";
              // print out this AF;
              for (auto val : age_freq_)
                cache_ << val << " ";
              cache_ << "\n";
            }
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
