/**
 * @file Tagging.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 26/7/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This file exists to keep documentation generator happy
 */

// headers
#include "Tagging.h"

#include "Selectivities/Manager.h"
#include "Layers/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "Agents/Agent.h"


#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
Tagging::Tagging(Model* model) : Process(model) {
  proportions_table_ = new parameters::Table(PARAM_PROPORTIONS);

  parameters_.Bind<string>(PARAM_TAGGING_LAYERS, &tag_layer_label_, "Spatial layer describing catch by cell for each year, there is a one to one link with the year specified, so make sure the order is right", "", true);
  parameters_.Bind<string>(PARAM_SELECTIVITIES, &selectivity_labels_, "selectivity used to capture agents", "", "");
  parameters_.Bind<float>(PARAM_HANDLING_MORTALITY, &handling_mortality_, "What is the handling mortality assumed for tagged fish, sometimes called initial mortality", "", 0);
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years to execute the transition in", "");
  parameters_.BindTable(PARAM_PROPORTIONS, proportions_table_, "Table of proportions to move", "" , true, true);

  process_type_ = ProcessType::kTransition;
}

/**
 * Destructor
 */
Tagging::~Tagging() {
  delete proportions_table_;
}

void Tagging::DoValidate() {
  LOG_FINE();
  if (years_.size() != tag_layer_label_.size()) {
    LOG_ERROR_P(PARAM_YEARS) << "there needs to be a layer for all years to apply tagging in. You supplied " << years_.size() << " years but " << tag_layer_label_.size() << " labels, sort this descrepency out please";
  }

  if (parameters_.Get(PARAM_SELECTIVITIES)->has_been_defined() & proportions_table_->has_been_defined()) {
    LOG_ERROR_P(PARAM_SELECTIVITIES) << "You can't provide both table of proportions and selectivity, you need to apply either one or the other.";
  }
  if (!parameters_.Get(PARAM_SELECTIVITIES)->has_been_defined() & !proportions_table_->has_been_defined()) {
    LOG_ERROR_P(PARAM_LABEL) << "You need to provide EITHER a table of proportions OR a selectivity, you haven't provided any";
  }

}

void Tagging::DoBuild() {
  LOG_FINE();
// Get the layers
  for (auto& label : tag_layer_label_) {
    layers::IntLayer* temp_layer = nullptr;
    temp_layer = model_->managers().layer()->GetIntLayer(label);
    if (!temp_layer) {
      LOG_FATAL_P(PARAM_TAGGING_LAYERS) << "could not find the layer '" << label << "', please make sure it exists and that it is type 'integer'";
    }
    tag_layer_.push_back(temp_layer);
  }

  if (model_->get_sexed()) {
     if (selectivity_labels_.size() == 1)
       selectivity_labels_.assign(2, selectivity_labels_[0]);
  }

  // Build selectivities
  bool first = true;
  for (auto label : selectivity_labels_) {
    Selectivity* temp_selectivity = model_->managers().selectivity()->GetSelectivity(label);
    if (!temp_selectivity)
      LOG_ERROR_P(PARAM_SELECTIVITIES) << ": selectivity " << label << " does not exist. Have you defined it?";

    selectivities_.push_back(temp_selectivity);
    if (first) {
      first = false;
      selectivity_length_based_ = temp_selectivity->is_length_based();
    } else {
      if (selectivity_length_based_ != temp_selectivity->is_length_based()) {
        LOG_ERROR_P(PARAM_SELECTIVITIES) << "The selectivity  " << label << " was not the same type (age or length based) as the previous selectivity label";
      }
    }
  }

  // Build table of propotions
  if (proportions_table_->row_count() != 0) {
    /**
     * Load data from proportions table using n parameter
     */
    apply_using_proportions_ = true;
    vector<string> columns = proportions_table_->columns();
    if (columns[0] != PARAM_YEAR)
      LOG_ERROR_P(PARAM_PROPORTIONS) << " first column label (" << columns[0] << ") provided must be 'year'";
    if (columns[1] != PARAM_CELL)
      LOG_ERROR_P(PARAM_PROPORTIONS) << " second column label (" << columns[1] << ") provided must be 'cell'";
    unsigned number_bins = columns.size();
    if (model_->length_plus()) {
      if ((number_bins - 2) != model_->length_bins().size())
        LOG_ERROR_P(PARAM_PROPORTIONS) << "Length bins for this observation are defined in the @model block, there must be a column for each length bin '" << model_->length_bins().size() << "' you supplied '"<< number_bins - 2  << "'. please address this";
    } else {
      if ((number_bins - 2) != (model_->length_bins().size() - 1))
        LOG_ERROR_P(PARAM_PROPORTIONS) << "Length bins for this observation are defined in the @model block, there must be a column for each length bin '" << model_->length_bins().size() - 1 << "' you supplied '"<< number_bins - 2  << "'. please address this";
    }

    // load our table data in to our map
    vector<vector<string>> data = proportions_table_->data();
    unsigned row_counter = 0;
    unsigned year = 0;
    string cell = "";
    float proportion = 0.0;
    for (auto iter : data) {
      ++row_counter;
      if (!utilities::To<unsigned>(iter[0], year))
        LOG_ERROR_P(PARAM_PROPORTIONS) << " value (" << iter[0] << ") is not a valid unsigned value that could be converted to a model year";
      if (find(years_.begin(), years_.end(), year) == years_.end())
        LOG_ERROR_P(PARAM_PROPORTIONS) << "at row '" << row_counter << "', the year " << year << " was not in the years provided in the subcommand " << PARAM_YEARS;
      table_years_.push_back(year);

      unsigned cell_ndx = 0;
      vector<string> split_cells;
      boost::split(split_cells, iter[1], boost::is_any_of("-"));
      LOG_FINEST() << "row = " << split_cells[0] << " col = " << split_cells[1];
      if (!utilities::To<unsigned>(split_cells[0], cell_ndx))
        LOG_ERROR_P(PARAM_PROPORTIONS) << " value (" << split_cells[0] << ") could not be converted to a unsigned";
      if (cell_ndx > model_->get_height())
        LOG_ERROR_P(PARAM_PROPORTIONS) << "The cell row " << cell_ndx << " at row " << row_counter << " is large than the world = '"<< model_->get_height() << "'";
      if (cell_ndx <= 0)
        LOG_ERROR_P(PARAM_PROPORTIONS) << "The cell row " << cell_ndx << " at row " << row_counter << " is less than equal to 0, must be greater than 0";
      table_rows_.push_back(cell_ndx);

      if (!utilities::To<unsigned>(split_cells[1], cell_ndx))
        LOG_ERROR_P(PARAM_PROPORTIONS) << " value (" << split_cells[1] << ") could not be converted to a unsigned";
       table_cols_.push_back(cell_ndx);

      if (cell_ndx > model_->get_width())
        LOG_ERROR_P(PARAM_PROPORTIONS) << "The cell col " << cell_ndx << " at row " << row_counter << " is large than the world = '"<< model_->get_width() << "'";

      if (cell_ndx <= 0)
        LOG_ERROR_P(PARAM_PROPORTIONS) << "The cell col " << cell_ndx << " at row " << row_counter << " is less than equal to 0, must be greater than 0";


      vector<float> proportions;
      float total_proportion = 0.0;
      for (unsigned i = 2; i < iter.size(); ++i) {
        if (!utilities::To<float>(iter[i], proportion))
          LOG_ERROR_P(PARAM_PROPORTIONS) << " value (" << iter[i] << ") could not be converted to a double. Please ensure it's a numeric value";
        proportions.push_back(proportion);
        total_proportion += proportion;
      }

      if (fabs(1.0 - total_proportion) > 0.01)
        LOG_ERROR_P(PARAM_PROPORTIONS) << " total (" << total_proportion << ") is not 1.0 (+- 0.01) for year " << year;
      proportions_data_.push_back(proportions);
    }
  }




  // Build cell specific containers in the hope that we can thread
  // Allocate memory so we are doing this during execute
  cell_offset_for_selectivity_.resize(model_->get_height());
  cell_offset_.resize(model_->get_height());
  model_length_bins_.resize(model_->get_height());
  model_age_bins_.resize(model_->get_height());
  handling_mort_by_space_.resize(model_->get_height());
  current_year_by_space_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cell_offset_[i].resize(model_->get_width());
    cell_offset_for_selectivity_[i].resize(model_->get_width());
    model_length_bins_[i].resize(model_->get_width(), model_->length_bin_mid_points().size());
    model_age_bins_[i].resize(model_->get_width(), model_->age_spread());
    current_year_by_space_[i].resize(model_->get_width());
    handling_mort_by_space_[i].resize(model_->get_width(), handling_mortality_);
  }

  length_distribution_of_tagged_fish_by_year_cell_.resize(years_.size());
  length_observed_tag_of_tagged_fish_by_year_cell_.resize(years_.size());
  age_distribution_of_tagged_fish_by_year_cell_.resize(years_.size());
  age_length_param1_of_tagged_fish_by_year_cell_.resize(years_.size());
  age_length_param2_of_tagged_fish_by_year_cell_.resize(years_.size());

  for(unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    age_distribution_of_tagged_fish_by_year_[years_[year_ndx]].resize(model_->age_spread(),0);
    length_distribution_of_tagged_fish_by_year_[years_[year_ndx]].resize(model_->length_bin_mid_points().size(),0);
    length_distribution_of_tagged_fish_by_year_cell_[year_ndx].resize(model_->get_height());
    length_observed_tag_of_tagged_fish_by_year_cell_[year_ndx].resize(model_->get_height());
    age_length_param1_of_tagged_fish_by_year_cell_[year_ndx].resize(model_->get_height());
    age_length_param2_of_tagged_fish_by_year_cell_[year_ndx].resize(model_->get_height());

    age_distribution_of_tagged_fish_by_year_cell_[year_ndx].resize(model_->get_height());
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      length_distribution_of_tagged_fish_by_year_cell_[year_ndx][row].resize(model_->get_width());
      length_observed_tag_of_tagged_fish_by_year_cell_[year_ndx][row].resize(model_->get_width());
      age_length_param1_of_tagged_fish_by_year_cell_[year_ndx][row].resize(model_->get_width());
      age_length_param2_of_tagged_fish_by_year_cell_[year_ndx][row].resize(model_->get_width());

      age_distribution_of_tagged_fish_by_year_cell_[year_ndx][row].resize(model_->get_width());
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        length_observed_tag_of_tagged_fish_by_year_cell_[year_ndx][row][col].resize(model_->number_of_length_bins(),0.0);
        length_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col].resize(model_->number_of_length_bins(),0.0);
        age_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col].resize(model_->age_spread(),0.0);
      }
    }
  }

  LOG_MEDIUM() << "length based = " << selectivity_length_based_ << " using proportions = " << apply_using_proportions_;
}

/*
 * DoReset()
 */
void Tagging::DoReset() {
  // Reset reporting containers
  for(unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        fill(age_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col].begin(),age_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col].end(),0.0);
        fill(length_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col].begin(),length_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col].end(),0.0);
        fill(length_observed_tag_of_tagged_fish_by_year_cell_[year_ndx][row][col].begin(),length_observed_tag_of_tagged_fish_by_year_cell_[year_ndx][row][col].end(),0.0);


      }
    }
  }
}

void Tagging::DoExecute() {
  LOG_MEDIUM();
  std::pair<bool, int> check = utilities::math::findInVector<unsigned>(years_, model_->current_year());
  if ((model_->state() != State::kInitialise) & check.first) {
    utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
    unsigned year_ndx = check.second;
    // Are we applying Tagging either releases or scanning, or both
    LOG_MEDIUM() << "applying tagging in year " << years_[year_ndx] << " should be = " << model_->current_year() << " ndx = " << year_ndx;
    // Pre-calculate agents in the world to set aside our random numbers needed for the operation
    n_agents_ = 0;
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          cell_offset_[row][col] = n_agents_;
          n_agents_ += cell->agents_.size();
          current_year_by_space_[row][col] = model_->current_year();
        }
      }
    }
    fill(age_distribution_of_tagged_fish_by_year_[model_->current_year()].begin(), age_distribution_of_tagged_fish_by_year_[model_->current_year()].end(), 0.0);
    fill(length_distribution_of_tagged_fish_by_year_[model_->current_year()].begin(), length_distribution_of_tagged_fish_by_year_[model_->current_year()].end(), 0.0);

    LOG_FINE() << "allocate";
    // Allocate a single block of memory rather than each thread temporarily allocating their own memory.
    random_numbers_.resize(n_agents_ + 1);
    selectivity_random_numbers_.resize(n_agents_ + 1);
    handling_mortality_random_numbers_.resize(n_agents_ + 1);
    for (unsigned i = 0; i <= n_agents_; ++i) {
      random_numbers_[i] = rng.chance();
      selectivity_random_numbers_[i] = rng.chance();
      handling_mortality_random_numbers_[i] = rng.chance();
    }
    LOG_MEDIUM() << "about to apply tagging";
    if (not selectivity_length_based_ & not apply_using_proportions_) {
      // Thread out each loop
      //#pragma omp parallel for collapse(2)
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          WorldCell* cell = nullptr;
          unsigned tags_to_release = 0;

          cell = world_->get_base_square(row, col); // Shared resource...
          tags_to_release = tag_layer_[year_ndx]->get_value(row, col);

          if (cell->is_enabled()) {
            unsigned tag_attempts = 1;
            unsigned random_agent;
            unsigned counter = 0;
            unsigned tag_max = cell->agents_.size();
            vector<unsigned>  age_freq(model_age_bins_[row][col],0);
            vector<unsigned>  length_freq(model_length_bins_[row][col],0);
            LOG_FINE() << "row " << row + 1 << " col = " << col + 1 << " tags to release = " << tags_to_release;
            unsigned tag_slot = 0;
            while (tags_to_release > 0) {
              ++tag_attempts;
              // pick a random agent
              random_agent = random_numbers_[cell_offset_[row][col] + counter] * cell->agents_.size();
              if (cell->agents_[random_agent].is_alive()) {
                auto& this_agent = cell->agents_[random_agent];
                // See if it is vulnerable to selectivity
                if (selectivity_random_numbers_[cell_offset_[row][col] + counter] <= selectivities_[this_agent.get_sex()]->GetResult(this_agent.get_age_index())) {
                  if (this_agent.get_number_tags() > 1) // This fish is alread tagged so pretend we didn't catch it
                    continue;

                  // Caught this agent we need to split out a tagged fish and return the two types
                  Agent tagged_agent(this_agent);
                  tagged_agent.apply_tagging_event(1, row, col); // Any tagging attribute should be bundled into this method
                  this_agent.set_scalar(this_agent.get_scalar() - 1.0);
                  tagged_agent.set_scalar(1.0);
                  LOG_MEDIUM() << "original tag " << this_agent.get_number_tags() << " tagged version = " << tagged_agent.get_number_tags();
                  age_freq[this_agent.get_age_index()]++;
                  length_freq[this_agent.get_length_bin_index()]++;
                  length_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col][this_agent.get_length_bin_index()]++;
                  age_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col][this_agent.get_age_index()]++;
                  tags_to_release--;
                  age_length_param1_of_tagged_fish_by_year_cell_[year_ndx][row][col].push_back(tagged_agent.get_first_age_length_par());
                  age_length_param2_of_tagged_fish_by_year_cell_[year_ndx][row][col].push_back(tagged_agent.get_second_age_length_par());

                  // Lets see if it survives handling
                  if (handling_mortality_random_numbers_[cell_offset_[row][col] + counter] <= handling_mort_by_space_[row][col]) {
                    // It died we will never see this or know about this tagged fish so I am just going to skip the rest of the algorithm
                    cell->remove_agent_alive(this_agent.get_scalar());
                    tagged_agent.dies();
                    continue;
                  }
                  // Agent is tagged and survived the process find a slot to add this tagged agent back in
                  while(tag_slot < cell->tagged_agents_.size()) {
                    if (not cell->tagged_agents_[tag_slot].is_alive()) {
                      cell->tagged_agents_[tag_slot] = tagged_agent;
                      break;
                    } else {
                      tag_slot++;
                    }
                  }
                  //LOG_FINEST() << "last_agent_ndx " << last_agent_ndx << " size = " << base_grid_[i][j].agents_.size();
                  if (tag_slot >= cell->tagged_agents_.size()) {
                    cell->tagged_agents_.push_back(tagged_agent);
                  }
                }
                // Make sure we don't end up fishing for infinity
                if (tag_attempts >= tag_max) {
                  LOG_FATAL_P(PARAM_LABEL) << "Too many attempts to catch an agent in the process " << label_ << " in year " << current_year_by_space_[row][col] << " in row " << row + 1 << " and column " << col + 1 << " this most likely means you have" <<
                     " a model that suggests there should be more agents in this space than than the current agent dynamics are putting in this cell, check the user manual for tips to resolve this situation, agents in cell = " << tag_max << " attempts made = " << tag_attempts;
                }
              }
              ++counter;
            }
            // Store global information
            for (unsigned age = 0; age < model_->age_spread(); ++age)
              age_distribution_of_tagged_fish_by_year_[model_->current_year()][age] += age_freq[age];
            for (unsigned length_ndx = 0; length_ndx < model_->length_bin_mid_points().size(); ++length_ndx)
              length_distribution_of_tagged_fish_by_year_[model_->current_year()][length_ndx] += length_freq[length_ndx];


          }
        }
      }
    } else if (selectivity_length_based_ & not apply_using_proportions_) {
      LOG_WARNING() << "in process" << label_ << " length based tagging not yet implemented";



    } else if (apply_using_proportions_) {
      LOG_FINE() << "applying tagging via given proportions";
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          WorldCell* cell = nullptr;
          unsigned tags_to_release = 0;
          // find row index
          unsigned row_ndx = 0;
          for ( ; row_ndx < table_years_.size(); ++row_ndx) {
            if ((table_years_[row_ndx] == model_->current_year()) & (table_rows_[row_ndx] == (row + 1)) & (table_cols_[row_ndx] == (col + 1)))
              break;
          }
          LOG_MEDIUM() << "Year = " << model_->current_year() << " row = " << row + 1<< " col = " << col + 1 << " table ndx = " << row_ndx;
          cell = world_->get_base_square(row, col); // Shared resource...
          tags_to_release = tag_layer_[year_ndx]->get_value(row, col);
          vector<unsigned> tags_by_length_bin(proportions_data_[row_ndx].size(), 0);
          for(unsigned i = 0; i < tags_by_length_bin.size(); ++i) {
            tags_by_length_bin[i] = (unsigned)(proportions_data_[row_ndx][i] * tags_to_release);
            length_observed_tag_of_tagged_fish_by_year_cell_[year_ndx][row][col][i] = tags_by_length_bin[i];
          }
          if (cell->is_enabled()) {
            unsigned tag_slot = 0;
            vector<unsigned>  age_freq(model_age_bins_[row][col],0);
            vector<unsigned>  length_freq(model_length_bins_[row][col],0);
            for(unsigned i = 0; i < tags_by_length_bin.size(); ++i) {
              if (tags_by_length_bin[i] > 0) {
                vector<unsigned> agent_ndx_available_to_sample_;
                unsigned agent_counter = 0;
                for (auto& agent : cell->agents_) {
                  if (agent.is_alive() & (agent.get_length_bin_index() == i))
                    agent_ndx_available_to_sample_.push_back(agent_counter);
                  ++agent_counter;
                }
                LOG_FINE() << "length bin = " << i + 1 << " tags " << tags_by_length_bin[i] << " agents alive that are in this length bin = " << agent_ndx_available_to_sample_.size();
                if(agent_ndx_available_to_sample_.size() <= tags_by_length_bin[i]) {
                  // jump out of this length bin
                  // push tags into next length bin if we aren't in the last bin
                  LOG_WARNING() << "couldn't tag " << tags_by_length_bin[i] << " in length bin " << i + 1;
                  if ((i + 1) < tags_by_length_bin.size()) {
                    tags_by_length_bin[i + 1] += tags_by_length_bin[i];
                  }
                  continue;
                }
                LOG_FINE() << "length bin = " << i + 1 << " tags " << tags_by_length_bin[i] << " agents alive that are in this length bin = " << agent_ndx_available_to_sample_.size();

                agent_counter = tags_by_length_bin[i];
                unsigned agent_ndx = 0;
                while(agent_counter > 0) {
                  //LOG_MEDIUM() << "agent counter = " << agent_counter;
                  agent_ndx = agent_ndx_available_to_sample_[agent_ndx_available_to_sample_.size() * rng.chance()];
                  agent_counter--;
                  Agent tagged_agent(cell->agents_[agent_ndx]);
                  tagged_agent.apply_tagging_event(1, row, col); // Any tagging attribute should be bundled into this method
                  cell->agents_[agent_ndx].set_scalar(cell->agents_[agent_ndx].get_scalar() - 1.0);
                  tagged_agent.set_scalar(1.0);
                  age_freq[tagged_agent.get_age_index()]++;
                  length_freq[tagged_agent.get_length_bin_index()]++;
                  tags_to_release--;
                  //LOG_MEDIUM() << "original tag " << cell->agents_[agent_ndx].get_number_tags() << " tagged version = " << tagged_agent.get_number_tags();

                  length_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col][tagged_agent.get_length_bin_index()]++;
                  age_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col][tagged_agent.get_age_index()]++;
                  age_length_param1_of_tagged_fish_by_year_cell_[year_ndx][row][col].push_back(tagged_agent.get_first_age_length_par());
                  age_length_param2_of_tagged_fish_by_year_cell_[year_ndx][row][col].push_back(tagged_agent.get_second_age_length_par());

                  // Lets see if it survives handling
                  if (rng.chance() <= handling_mortality_) {
                    // It died we will never see this or know about this tagged fish so I am just going to skip the rest of the algorithm
                    tagged_agent.dies();
                    continue;
                  }

                  // Agent is tagged and survived the process find a slot to add this tagged agent back in
                  while(tag_slot < cell->tagged_agents_.size()) {
                    if (not cell->tagged_agents_[tag_slot].is_alive()) {
                      cell->tagged_agents_[tag_slot] = tagged_agent;
                      break;
                    } else {
                      tag_slot++;
                    }
                  }
                  //LOG_FINEST() << "last_agent_ndx " << last_agent_ndx << " size = " << base_grid_[i][j].agents_.size();
                  if (tag_slot >= cell->tagged_agents_.size()) {
                    cell->tagged_agents_.push_back(tagged_agent);
                  }
                }
                if (agent_counter < 0) {
                	LOG_WARNING() << "in Process " << label_ << " couldn't apply all tags in year " << model_->current_year() << " tags that didn't get released = " << tags_to_release;
                }
              }
            }
            // Store global information
            for (unsigned age = 0; age < model_->age_spread(); ++age)
              age_distribution_of_tagged_fish_by_year_[model_->current_year()][age] += age_freq[age];
            for (unsigned length_ndx = 0; length_ndx < model_->length_bin_mid_points().size(); ++length_ndx)
              length_distribution_of_tagged_fish_by_year_[model_->current_year()][length_ndx] += length_freq[length_ndx];

          }
        }
      }
    }
  }
}

// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  Tagging::FillReportCache(ostringstream& cache) {
  cache << "global_tag_release_age_distribution " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
  cache << "year";
  for (unsigned age = 0; age < model_->age_spread(); ++age) {
    cache << " " << age + model_->min_age();
  }
  for (auto& year_value : age_distribution_of_tagged_fish_by_year_) {
    cache << "\n" << year_value.first;
    for (unsigned age_ndx = 0; age_ndx < year_value.second.size(); ++age_ndx)
      cache << " " << year_value.second[age_ndx];
  }
  cache << "\n";
  cache << "global_tag_release_length_distribution " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
  cache << "year";
  for (unsigned length_ndx = 0; length_ndx < model_->length_bin_mid_points().size(); ++length_ndx) {
    cache << " " << model_->length_bin_mid_points()[length_ndx];
  }
  for (auto& year_value : length_distribution_of_tagged_fish_by_year_) {
    cache << "\n" << year_value.first;
    for (unsigned length_ndx = 0; length_ndx < year_value.second.size(); ++length_ndx)
      cache << " " << year_value.second[length_ndx];
  }
  cache << "\n";
  for(unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    cache << "tag_release_by_length_observed-" << years_[year_ndx] << " " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
    cache << "cell";
    for (unsigned length_ndx = 0; length_ndx < model_->length_bin_mid_points().size(); ++length_ndx) {
      cache << " " << model_->length_bin_mid_points()[length_ndx];
    }
    cache << "\n";
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        cache << row + 1 << "-" << col + 1 << " ";
        for (unsigned len_ndx = 0; len_ndx < model_->number_of_length_bins(); ++len_ndx) {
          cache << length_observed_tag_of_tagged_fish_by_year_cell_[year_ndx][row][col][len_ndx] << " ";
        }
        cache << "\n";
      }
    }
  }

  for(unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        cache << "age_length_param1-" << years_[year_ndx] << "-" << row + 1 << "-" << col + 1 << ": ";
        for (unsigned agent_ndx = 0; agent_ndx < age_length_param2_of_tagged_fish_by_year_cell_[year_ndx][row][col].size(); ++agent_ndx) {
          cache << age_length_param2_of_tagged_fish_by_year_cell_[year_ndx][row][col][agent_ndx] << " ";
        }
        cache << "\n";
      }
    }
  }

  for(unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        cache << "age_length_param2-" << years_[year_ndx] << "-" << row + 1 << "-" << col + 1 << ": ";
        for (unsigned agent_ndx = 0; agent_ndx < age_length_param2_of_tagged_fish_by_year_cell_[year_ndx][row][col].size(); ++agent_ndx) {
          cache << age_length_param2_of_tagged_fish_by_year_cell_[year_ndx][row][col][agent_ndx] << " ";
        }
        cache << "\n";
      }
    }
  }

  for(unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    cache << "tag_release_by_length-" << years_[year_ndx] << " " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
    cache << "cell";
    for (unsigned length_ndx = 0; length_ndx < model_->length_bin_mid_points().size(); ++length_ndx) {
      cache << " " << model_->length_bin_mid_points()[length_ndx];
    }
    cache << "\n";
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        cache << row + 1 << "-" << col + 1 << " ";
        for (unsigned len_ndx = 0; len_ndx < model_->number_of_length_bins(); ++len_ndx) {
          cache << length_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col][len_ndx] << " ";
        }
        cache << "\n";
      }
    }
  }


  for(unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    cache << "tag_release_by_age-" << years_[year_ndx] << " " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
    cache << "cell";
    for (unsigned age = 0; age < model_->age_spread(); ++age) {
      cache << " " << age + model_->min_age();
    }
    cache << "\n";
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        cache << row + 1 << "-" << col + 1 << " ";
        for (unsigned age_ndx = 0; age_ndx < model_->age_spread(); ++age_ndx) {
          cache << age_distribution_of_tagged_fish_by_year_cell_[year_ndx][row][col][age_ndx] << " ";
        }
        cache << "\n";
      }
    }
  }
}

}
}
