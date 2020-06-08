/**
 * @file TagShedding.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */

// headers
#include "TagShedding.h"

#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "TimeSteps/Manager.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
TagShedding::TagShedding(Model *model) :
    Process(model) {
  parameters_.Bind<string>(PARAM_SELECTIVITY_LABEL, &selectivity_label_, "Label for the selectivity block", "");
  parameters_.Bind<float>(PARAM_TIME_STEP_RATIO, &ratios_, "Time step ratios for the Shedding rates to apply in each time step. See manual for how this is applied", "");
  parameters_.Bind<float>(PARAM_SHEDDING_RATE, &shedding_rate_, "Shedding rate per Tag release event", "");
  parameters_.Bind<string>(PARAM_RELEASE_REGION, &release_region_,
      "The Release region for the corresponding shedding rate, in the format 1-1 for the first row and first column, and '5-2' for the fifth row and secnd column", "");
  parameters_.Bind<unsigned>(PARAM_RELEASE_YEAR, &tag_release_year_, "The Release Year for the corresponding shedding rate", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years to execute the tag shedding process", "");
  parameters_.Bind<unsigned>(PARAM_STOP_TRACKING_IN_YEARS, &stop_tracking_tags_in_year_, "Stop Tracking tagged partition in this year, they will just die. But will stop containimation in observations.", "");

  process_type_ = ProcessType::kTransition;

}

/**
 * DoBuild
 */
void TagShedding::DoValidate() {

  if (shedding_rate_.size() != release_region_.size())
    LOG_ERROR_P(PARAM_SHEDDING_RATE) << "You suppled " << shedding_rate_.size() << " shedding rates, but '" << release_region_.size()
        << "' release regions, there must be a shedding rate for each release regiond";
  if (shedding_rate_.size() != tag_release_year_.size())
    LOG_ERROR_P(PARAM_SHEDDING_RATE) << "You suppled " << shedding_rate_.size() << " shedding rates, but '" << tag_release_year_.size()
        << "' release year, there must be a shedding rate for each release year";

  for(unsigned year_ndx = 0; year_ndx < stop_tracking_tags_in_year_.size(); ++year_ndx) {
    if(find(years_.begin(), years_.end(), stop_tracking_tags_in_year_[year_ndx]) == years_.end())
      LOG_ERROR_P(PARAM_STOP_TRACKING_IN_YEARS) << "year = " << stop_tracking_tags_in_year_[year_ndx] << " needs to be in year " << PARAM_YEARS;
  }
  if(stop_tracking_tags_in_year_.size() != release_region_.size())
    LOG_ERROR_P(PARAM_STOP_TRACKING_IN_YEARS) << "You suppled " << stop_tracking_tags_in_year_.size() << " years to stop tracking, but '" << release_region_.size()
        << "' release regions, there must be a stop year for each release region";

}
/**
 * DoBuild
 */
void TagShedding::DoBuild() {

  LOG_FINEST() << "selectivities supplied = " << selectivity_label_.size();
  // Build selectivity links
  if (selectivity_label_.size() == 1)
    selectivity_label_.assign(2, selectivity_label_[0]);

  if (selectivity_label_.size() > 2) {
    LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << "You suppled " << selectivity_label_.size() << " Selectiviites, you can only have one for each sex max = 2";
  }
  LOG_FINEST() << "selectivities supplied = " << selectivity_label_.size();

  bool first = true;
  for (auto label : selectivity_label_) {
    Selectivity *temp_selectivity = model_->managers().selectivity()->GetSelectivity(label);
    if (!temp_selectivity)
      LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << ": selectivity " << label << " does not exist. Have you defined it?";
    temp_selectivity->SubscribeToRebuildCache(this);

    selectivity_.push_back(temp_selectivity);
    if (first) {
      first = false;
      selectivity_length_based_ = temp_selectivity->is_length_based();
    } else {
      if (selectivity_length_based_ != temp_selectivity->is_length_based()) {
        LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << "The selectivity  " << label << " was not the same type (age or length based) as the previous selectivity label";
      }
    }
  }

  /**
   * Organise our time step ratios. Each time step can
   * apply a different ratio of M so here we want to verify
   * we have enough and re-scale them to 1.0
   */
  vector<TimeStep*> time_steps = model_->managers().time_step()->ordered_time_steps();
  LOG_MEDIUM() << "time_steps.size(): " << time_steps.size();
  vector<unsigned> active_time_steps;
  for (unsigned i = 0; i < time_steps.size(); ++i) {
    if (time_steps[i]->HasProcess(label_))
      active_time_steps.push_back(i);
  }

  if (ratios_.size() == 0) {
    for (unsigned i : active_time_steps)
      time_step_ratios_[i] = 1.0;
  } else {
    if (ratios_.size() != active_time_steps.size())
      LOG_ERROR_P(PARAM_TIME_STEP_RATIO) << " length (" << ratios_.size() << ") does not match the number of time steps this process has been assigned to (" << active_time_steps.size() << ")";

    for (float value : ratios_) {
      if (value < 0.0 || value > 1.0)
        LOG_ERROR_P(PARAM_TIME_STEP_RATIO) << " value (" << value << ") must be between 0.0 (exclusive) and 1.0 (inclusive)";
    }

    for (unsigned i = 0; i < ratios_.size(); ++i)
      time_step_ratios_[active_time_steps[i]] = ratios_[i];
  }
  vector<vector<float>> temp_empty_matrix;
  temp_empty_matrix.resize(model_->get_height());
  for (unsigned row_ndx = 0; row_ndx < model_->get_height(); ++row_ndx) {
    temp_empty_matrix[row_ndx].resize(model_->get_width(), 0.0);
  }
  unsigned row_ndx = 0;
  unsigned col_ndx = 0;

  for (unsigned i = 0; i < tag_release_year_.size(); ++i) {
    shedding_rate_by_year_cell_[tag_release_year_[i]] = temp_empty_matrix;
    vector<string> split_cells;
    boost::split(split_cells, release_region_[i], boost::is_any_of("-"));
    LOG_MEDIUM() << "row = " << split_cells[0] << " col = " << split_cells[1];
    if (!utilities::To<unsigned>(split_cells[0], row_ndx))
      LOG_ERROR_P(PARAM_RELEASE_REGION) << " value (" << split_cells[0] << ") could not be converted to a unsigned";
    if (row_ndx > model_->get_height())
      LOG_ERROR_P(PARAM_RELEASE_REGION) << "The release row " << row_ndx << " is large than the world = '" << model_->get_height() << "'";
    if (row_ndx <= 0)
      LOG_ERROR_P(PARAM_RELEASE_REGION) << "The release row " << row_ndx << " is less than equal to 0, it must be greater than 0";
    release_row_.push_back(row_ndx - 1);

    if (!utilities::To<unsigned>(split_cells[1], col_ndx))
      LOG_ERROR_P(PARAM_RELEASE_REGION) << " value (" << split_cells[1] << ") could not be converted to a unsigned";
    if (col_ndx > model_->get_width())
      LOG_ERROR_P(PARAM_RELEASE_REGION) << "The release col " << col_ndx << " is large than the world = '" << model_->get_width() << "'";
    if (col_ndx <= 0)
      LOG_ERROR_P(PARAM_RELEASE_REGION) << "The release col " << col_ndx << " is less than equal to 0, it must be greater than 0";
    release_col_.push_back(col_ndx - 1);

    shedding_rate_by_year_cell_[tag_release_year_[i]][row_ndx - 1][col_ndx - 1] = shedding_rate_[i];

    LOG_MEDIUM() << "double check shedding rate has been entered correctly for " << i << " shedding rate = " << shedding_rate_by_year_cell_[tag_release_year_[i]][row_ndx - 1][col_ndx - 1]
        << " year = " << tag_release_year_[i] << " row = " << row_ndx - 1 << " col = " << col_ndx - 1;
  }
  tag_shedded_per_release_event_.resize(tag_release_year_.size(), 0);

}

/**
 * DoReset
 */
void TagShedding::DoReset() {
  LOG_MEDIUM() << "clearing containers";

}

/**
 * DoExecute
 */
void TagShedding::DoExecute() {
  LOG_MEDIUM();
  std::pair<bool, int> check = utilities::math::findInVector<unsigned>(years_, model_->current_year());
  LOG_MEDIUM() << "year = " << model_->current_year() << " check first = " << check.first;
  if ((model_->state() != State::kInitialise) & check.first) {
    utilities::RandomNumberGenerator &rng = utilities::RandomNumberGenerator::Instance();
    unsigned time_step = model_->managers().time_step()->current_time_step();
    float ratio = time_step_ratios_[time_step];
    LOG_MEDIUM() << "time_step = " << time_step << " ratio = " << ratio;
    /*
     * Are we removing this cohort of tagged agents
     */
    for(unsigned year_ndx = 0; year_ndx < stop_tracking_tags_in_year_.size(); ++year_ndx) {
      if(model_->current_year() == stop_tracking_tags_in_year_[year_ndx]) {
        LOG_FINE() << "year = " << model_->current_year() << " years to track = " << release_row_[year_ndx] << " year ndx = " << year_ndx;
        LOG_FINE() << "row = " << release_row_[year_ndx] << " col = " << release_col_[year_ndx] << " release year = " << tag_release_year_[year_ndx];
        WorldCell *cell = world_->get_base_square(release_row_[year_ndx], release_col_[year_ndx]);
        LOG_FINE() << "number of tagged agents = " << cell->tagged_agents_.size();

        for (unsigned agent_ndx = 0;  agent_ndx < cell->tagged_agents_.size(); ++agent_ndx) {
          if (cell->tagged_agents_[agent_ndx].is_alive()) {
            if(cell->tagged_agents_[agent_ndx].get_tag_release_year() == tag_release_year_[year_ndx])
              cell->tagged_agents_[agent_ndx].dies();
          }
        }
      }
    }
    /*
     // Agent is tagged and survived the process find a slot to add this tagged agent back in
     */
    unsigned tag_slot = 0;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        // get the ratio to apply first
        WorldCell *cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          if (not selectivity_length_based_) {
            // age based selectivity
            LOG_MEDIUM() << "age based tag shedding , num tags =  " << cell->tagged_agents_.size();
            for (auto &tag_agent : cell->tagged_agents_) {
              if (tag_agent.is_alive()) {
                //LOG_FINE() << "tagged fish shedding rate = " << shedding_rate_by_year_cell_[tag_agent.get_tag_release_year()][tag_agent.get_tag_row()][tag_agent.get_tag_col()] << " release year = " << tag_agent.get_tag_release_year() << " release col = " << tag_agent.get_tag_col() << " release row = " << tag_agent.get_tag_row();
                if (rng.chance()  <= (1 - exp( -ratio * selectivity_[tag_agent.get_sex()]->GetResult(tag_agent.get_age_index())  * shedding_rate_by_year_cell_[tag_agent.get_tag_release_year()][tag_agent.get_tag_row()][tag_agent.get_tag_col()]))) {
                  // Taghas been sheded
                  tag_agent.shed_a_tag();
                  if (tag_agent.get_number_tags() >= 1) {
                    // do nothing, was double tagged and so still in the tagged partition
                    LOG_MEDIUM() << "found tagged agent with  " << tag_agent.get_number_tags() << " tags, this can't be how?";
                  } else {
                    // merge tagged fish backinto untagged partition
                    tag_slot = 0;
                    while (tag_slot < cell->agents_.size()) {
                      if (not cell->agents_[tag_slot].is_alive()) {
                        cell->agents_[tag_slot] = tag_agent;
                        break;
                      } else {
                        tag_slot++;
                      }
                    }
                    if (tag_slot >= cell->agents_.size()) {
                      cell->agents_.push_back(tag_agent);
                    }
                    tag_agent.dies();
                  } // the other conditions should not exist i.e. a tagged agent with < 0 tags
                }
              }
            }
          } else {
            // length based selectivity
            LOG_MEDIUM() << "length based tag shedding , num tags =  " << cell->tagged_agents_.size();
            for (auto &tag_agent : cell->tagged_agents_) {
              if (tag_agent.is_alive()) {
                if (rng.chance()  <= (1  - exp(  -ratio * selectivity_[tag_agent.get_sex()]->GetResult(tag_agent.get_length_bin_index())  * shedding_rate_by_year_cell_[tag_agent.get_tag_release_year()][tag_agent.get_tag_row()][tag_agent.get_tag_col()]))) {
                  // Taghas been sheded
                  tag_agent.shed_a_tag();
                  if (tag_agent.get_number_tags() >= 1) {
                    // do nothing, was double tagged and so still in the tagged partition
                  } else {
                    // merge tagged fish backinto untagged partition
                    tag_slot = 0;
                    while (tag_slot < cell->agents_.size()) {
                      if (not cell->agents_[tag_slot].is_alive()) {
                        cell->agents_[tag_slot] = tag_agent;
                        break;
                      } else {
                        tag_slot++;
                      }
                    }
                    if (tag_slot >= cell->agents_.size()) {
                      cell->agents_.push_back(tag_agent);
                    }
                    tag_agent.dies();
                  } // the other conditions should not exist i.e. a tagged agent with < 0 tags
                }
              }
            }
          }
        } // cell enabled
      } // col enabled
    } // row enabled
  } // check year
} // Doexecute

// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void TagShedding::FillReportCache(ostringstream &cache) {
  /*
   cache << "year: ";
   for (auto& year : removals_by_year_)
   cache << year.first << " ";

   cache << "\nagents_removed: ";
   for (auto& year : removals_by_year_)
   cache << year.second << " ";
   cache << "\n";*/
}

// A Method for telling the world we need to redistribute Mortality parmaeters
void TagShedding::RebuildCache() {
  LOG_FINE();

}

} /* namespace processes */
} /* namespace niwa */
