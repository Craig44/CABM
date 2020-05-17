/**
 * @file SummariseAgents.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 17/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "SummariseAgents.h"

#include <boost/algorithm/string/join.hpp>

#include "Model/Managers.h"
#include "Model/Model.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "Utilities/RandomNumberGenerator.h"

// namespaces
namespace niwa {
namespace reports {

/**
 * Default constructor
 *
 * @param model Pointer to the current model context
 */
SummariseAgents::SummariseAgents(Model* model) : Report(model) {
  model_state_ = State::kExecute;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);

  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years", "", true);
  parameters_.Bind<string>(PARAM_TIME_STEP, &time_step_, "Time Step label", "");
  parameters_.Bind<unsigned>(PARAM_NUMBER_OF_AGENTS, &n_agents_, "Number of agents to summarise", "")->set_lower_bound(1,true);
}

/**
 * Build our relationships between this object and other objects
 */
void SummariseAgents::DoBuild() {
  LOG_TRACE();
  world_ = model_->world_view();
  if (!world_)
    LOG_CODE_ERROR() << "!world_, something has gone wrong some where, cannot get a pointer to the world via model";

  for (unsigned i = 0; i < model_->get_height(); ++i)
    rows_.push_back(i);
  for (unsigned i = 0; i < model_->get_width(); ++i)
    cols_.push_back(i);

  if (!parameters_.Get(PARAM_YEARS)->has_been_defined()) {
    years_ = model_->years();
  }
}

/**
 * Execute this report
 */
void SummariseAgents::DoExecute() {
  LOG_FINE() <<" printing report " << label_;
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";
  cache_ << "year: " << model_->current_year() << "\n";
  cache_ << "time_step: " << time_step_ << "\n";
  cache_ << "values "<< REPORT_R_DATAFRAME<<"\n";

  unsigned temp_n_agents = n_agents_;

  cache_ << "row-col agent_ndx age age_index length length_index weight scalar sex mature M L_inf k a b tag alive lat long birth-row-col\n";
  while (temp_n_agents != 0) {
    // Randomly select a cell
    unsigned row_index = rng.chance() * rows_.size();
    unsigned col_index = rng.chance() * cols_.size();
    //LOG_FINEST() << "col index = " << col_index << " row index " << row_index;
    WorldCell* cell = world_->get_base_square(rows_[row_index], cols_[col_index]);
    if (cell->is_enabled()) {
      unsigned agent_ndx = rng.chance() * cell->agents_.size();
      auto& agent = cell->agents_[agent_ndx];
      if (agent.is_alive()) {
        cache_ << row_index + 1 << "-" << col_index + 1 << " " << agent_ndx + 1 <<  " " << agent.get_age() <<  " " << agent.get_age_index() << " " << agent.get_length() << " " << agent.get_length_bin_index() << " " << agent.get_weight() << " " << agent.get_scalar() << " " <<  agent.get_sex() << " " << agent.get_maturity()
              << " " << agent.get_m() << " " << agent.get_first_age_length_par() << " " << agent.get_second_age_length_par() << " " << agent.get_first_length_weight_par() << " " << agent.get_second_length_weight_par()  << " " << agent.get_number_tags() << " " << agent.is_alive() << " " << agent.get_lat()  << " " << agent.get_lon()
              << " " << agent.get_home_row() + 1 << "-" << agent.get_home_col() + 1 << "\n";
        --temp_n_agents;
      }
    }
  }

  ready_for_writing_ = true;
}


} /* namespace reports */
} /* namespace niwa */
