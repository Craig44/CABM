/**
 * @file MortalityEventBiomass.cpp
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
#include "MortalityEventBiomass.h"

#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
MortalityEventBiomass::MortalityEventBiomass(Model* model) : Mortality(model) {
  catch_table_ = new parameters::Table(PARAM_CATCH);
  parameters_.BindTable(PARAM_CATCH, catch_table_, "Table of catch_layer by year for each fishery", "", true, false);
  scanning_table_ = new parameters::Table(PARAM_SCANNING);
  parameters_.BindTable(PARAM_SCANNING, scanning_table_, "Table of scanning proportions by year for each fishery", "", true, true);
  method_table_ = new parameters::Table(PARAM_METHOD_INFO);
  parameters_.BindTable(PARAM_METHOD_INFO, method_table_, "Table of information for each fishery.", "", true, false);
  parameters_.Bind<bool>(PARAM_PRINT_EXTRA_INFO, &print_extra_info_, "if you have process report for this process you can control the amount of information printed to the file.", "", true);
}

MortalityEventBiomass::~MortalityEventBiomass() {
  delete catch_table_;
  delete method_table_;
  delete scanning_table_;
}

/**
 * Do some initial checks of user supplied parameters.
 */
void MortalityEventBiomass::DoValidate() {
  // check headers
  auto catch_columns = catch_table_->columns();
  LOG_FINE() << "check catch table headers, number of columns = " << catch_columns.size();

  if (catch_columns[0] != PARAM_YEAR)
    LOG_ERROR_P(PARAM_CATCH) << "the first column header needs to be " << PARAM_YEAR << ". Please add it =)";
  for (unsigned j = 1; j < catch_columns.size(); ++j) {
    fishery_index_.push_back(j - 1);
    fishery_label_.push_back(catch_columns[j]);
    LOG_FINE() << "fishery ndx = " << j - 1 << " fishery label " << catch_columns[j];
  }

  LOG_FINE() << "check method info";
  auto meth_columns = method_table_->columns();
  if (std::find(meth_columns.begin(), meth_columns.end(), PARAM_HANDLING_MORTALITY) == meth_columns.end()) {
    LOG_FATAL_P(PARAM_METHOD_INFO) << "The column header " << PARAM_HANDLING_MORTALITY << ", is missing can you please add it";
  }
  if (std::find(meth_columns.begin(), meth_columns.end(), PARAM_MINIMUM_LEGAL_LENGTH) == meth_columns.end()) {
    LOG_FATAL_P(PARAM_METHOD_INFO) << "The column header " << PARAM_MINIMUM_LEGAL_LENGTH << ", is missing can you please add it";
  }

  if (meth_columns[0] != PARAM_METHOD) {
    LOG_FATAL_P(PARAM_METHOD_INFO) << "The first column header must be labelled " << PARAM_METHOD << ", can you please sort this out";
  }

  if (!model_->get_sexed()) {
    // find the column header selectivity
    if (std::find(meth_columns.begin(), meth_columns.end(), PARAM_SELECTIVITY) == meth_columns.end()) {
      LOG_FATAL_P(PARAM_METHOD_INFO) << "because the model is unsexed we expect the column header " << PARAM_SELECTIVITY << ", this is missing can you please add it";
    }
  } else {
    if (std::find(meth_columns.begin(), meth_columns.end(), PARAM_MALE_SELECTIVITY) == meth_columns.end()) {
      LOG_FATAL_P(PARAM_METHOD_INFO) << "because the model is sexed we expect the column header " << PARAM_MALE_SELECTIVITY << ", this is missing can you please add it";
    }
    if (std::find(meth_columns.begin(), meth_columns.end(), PARAM_FEMALE_SELECTIVITY) == meth_columns.end()) {
      LOG_FATAL_P(PARAM_METHOD_INFO) << "because the model is sexed we expect the column header " << PARAM_FEMALE_SELECTIVITY << ", this is missing can you please add it";
    }
  }

  auto meth_data = method_table_->data();
  LOG_FINE() << "rows in meth info table " << meth_data.size() << " length fishery label " << fishery_label_.size();
  for (unsigned i = 0; i < meth_data.size(); ++i) {
    LOG_FINE() << "trying to check this " << meth_data[i][0] << " fishery is in the catch table " << i;
    if (std::find(fishery_label_.begin(), fishery_label_.end(), (string)meth_data[i][0]) == fishery_label_.end()) {
      LOG_ERROR_P(PARAM_METHOD_INFO) << "Could not find the row label " <<  meth_data[i][0] << " these need to be consistent with the column headers with the catch table, can you please double check this";
    } else {
      LOG_FINE() << "found " <<meth_data[i][0] << " in the catch table";
    }
  }

  // Scanning if user has defined it
  if (scanning_table_->has_been_defined()) {
    LOG_FINE() << "check scanning headers";
    auto scan_columns = scanning_table_->columns();
    if (scan_columns[0] != PARAM_YEAR)
      LOG_ERROR_P(PARAM_SCANNING) << "the first column header needs to be " << PARAM_YEAR << ". Please add it =)";
    for (unsigned j = 1; j < scan_columns.size(); ++j) {
      if (std::find(fishery_label_.begin(), fishery_label_.end(), scan_columns[j]) == fishery_label_.end()) {
        LOG_ERROR_P(PARAM_SCANNING) << "could not find the column header " << scan_columns[j] << " in the catch table column headers. The headers of this table must match the column headers of the catch table";
      }
    }
  }


  LOG_FINE() << "finish Validation";
}

/**
 * DoBuild
 */
void MortalityEventBiomass::DoBuild() {
  LOG_FINE() << "Build tables";
  LOG_FINE() << "Building catch table";

  // allocate memory for fishery dimension to related objects
  fishery_catch_layer_labels_.resize(fishery_label_.size());
  fishery_actual_catch_taken_.resize(fishery_label_.size());
  fishery_catch_to_take_.resize(fishery_label_.size());
  fishery_catch_layer_.resize(fishery_label_.size());
  fishery_selectivity_label_.resize(fishery_label_.size());
  fishery_selectivity_.resize(fishery_label_.size());
  fishery_mls_.resize(fishery_label_.size());
  fishery_hand_mort_.resize(fishery_label_.size());
  scanning_proportion_by_fishery_.resize(fishery_label_.size());
  age_comp_by_fishery_.resize(fishery_label_.size());
  length_comp_by_fishery_.resize(fishery_label_.size());
  fishery_census_data_.resize(fishery_label_.size());
  cell_ndx_.resize(fishery_label_.size());


  auto model_years = model_->years();
  // load objects
  vector<vector<string>>& catch_values_data = catch_table_->data();

  for (auto row : catch_values_data) {
    unsigned year = 0;
    if (!utilities::To<string, unsigned>(row[0], year))
       LOG_ERROR_P(PARAM_CATCH) << "year value " << row[0] << " is not numeric.";
     if (std::find(model_years.begin(), model_years.end(), year) == model_years.end())
       LOG_ERROR_P(PARAM_CATCH) << "year " << year << " is not a valid year in this model";
     // Check years are consecutive ascending order.
     // This will mean when I reference catch_ndx later in the code we can have faith it is in order.
     if (catch_year_.size() > 1) {
       if ((year - 1) != catch_year_[catch_year_.size() - 1]) {
         LOG_ERROR_P(PARAM_CATCH) << "years need to be in consecutive ascending order, the year " << catch_year_[catch_year_.size() - 1] << " was followed by " << year << " please sort this out";
       }
     }
     catch_year_.push_back(year);
     for (unsigned i = 1; i < row.size(); ++i) {
       fishery_catch_layer_labels_[i - 1].push_back(row[i]);
       layers::NumericLayer* temp_layer = nullptr;
       temp_layer = model_->managers().layer()->GetNumericLayer(row[i]);
       if (!temp_layer) {
         LOG_FATAL_P(PARAM_CATCH) << "could not find the catch layer '" << row[i] << "', in year " << year << " for fishery " << fishery_label_[i - 1] << " please check it is numeric if it exists";
       }
       fishery_catch_layer_[i - 1].push_back(temp_layer);
     }
  }
  // Allocate memory for some of the catch related maps
  for (unsigned i = 0; i < fishery_catch_layer_.size(); ++i) {
    fishery_actual_catch_taken_[i].resize(fishery_catch_layer_[i].size(), 0.0);
    fishery_catch_to_take_[i].resize(fishery_catch_layer_[i].size(), 0.0);
  }

  // Method info
  auto meth_data = method_table_->data();
  auto meth_cols = method_table_->columns();

  // special case with selectivity in this so do it first
  if (!model_->get_sexed()) {
    // find the column header selectivity
    unsigned selec_index = std::find(meth_cols.begin(), meth_cols.end(), PARAM_SELECTIVITY) - meth_cols.begin();
    for (auto meth_row : meth_data) {

      unsigned meth_index = std::find(fishery_label_.begin(), fishery_label_.end(), meth_row[0]) - fishery_label_.begin();

      fishery_selectivity_label_[fishery_index_[meth_index]].push_back(meth_row[selec_index]);
      fishery_selectivity_label_[fishery_index_[meth_index]].push_back(meth_row[selec_index]);
      Selectivity* temp_selectivity = model_->managers().selectivity()->GetSelectivity(meth_row[selec_index]);
      if (!temp_selectivity)
        LOG_FATAL_P(PARAM_METHOD_INFO) << ": selectivity " << meth_row[selec_index] << " does not exist. Have you defined it? located in column " << selec_index;

      if(temp_selectivity->is_length_based())
        selectivity_length_based_ = true;

      fishery_selectivity_[fishery_index_[meth_index]].push_back(temp_selectivity);
      fishery_selectivity_[fishery_index_[meth_index]].push_back(temp_selectivity);

    }
  }  else {
    unsigned male_selec_index = std::find(meth_cols.begin(), meth_cols.end(), PARAM_MALE_SELECTIVITY) - meth_cols.begin();
    unsigned female_selec_index = std::find(meth_cols.begin(), meth_cols.end(), PARAM_FEMALE_SELECTIVITY) - meth_cols.begin();
    for (auto meth_row : meth_data) {
      unsigned meth_index = std::find(fishery_label_.begin(), fishery_label_.end(), meth_row[0]) - fishery_label_.begin();
      fishery_selectivity_label_[fishery_index_[meth_index]].push_back(meth_row[male_selec_index]);
      fishery_selectivity_label_[fishery_index_[meth_index]].push_back(meth_row[female_selec_index]);
      Selectivity* male_temp_selectivity = model_->managers().selectivity()->GetSelectivity(meth_row[male_selec_index]);
      if (!male_temp_selectivity)
        LOG_FATAL_P(PARAM_METHOD_INFO) << ": male_selectivity " << meth_row[male_selec_index] << " does not exist. Have you defined it?";

      Selectivity* female_temp_selectivity = model_->managers().selectivity()->GetSelectivity(meth_row[female_selec_index]);
      if (!female_temp_selectivity)
        LOG_FATAL_P(PARAM_METHOD_INFO) << ": female_selectivity " << meth_row[female_selec_index] << " does not exist. Have you defined it?";

      if (female_temp_selectivity->is_length_based() != male_temp_selectivity->is_length_based())
        LOG_ERROR_P(PARAM_METHOD_INFO) << "One of the male and one of the female selectivities is age based and the other length based. They both have to be one or the other. Can you please resolve this.";

      if(female_temp_selectivity->is_length_based() && male_temp_selectivity->is_length_based())
        selectivity_length_based_ = true;

      fishery_selectivity_[fishery_index_[meth_index]].push_back(male_temp_selectivity);
      fishery_selectivity_[fishery_index_[meth_index]].push_back(female_temp_selectivity);
    }
  }

  unsigned mls_index = std::find(meth_cols.begin(), meth_cols.end(), PARAM_MINIMUM_LEGAL_LENGTH) - meth_cols.begin();
  unsigned hand_mort_index = std::find(meth_cols.begin(), meth_cols.end(), PARAM_HANDLING_MORTALITY) - meth_cols.begin();
  for (auto meth_row : meth_data) {
    unsigned meth_index = std::find(fishery_label_.begin(), fishery_label_.end(), meth_row[0]) - fishery_label_.begin();
    float mls = 0;
    float hand_mort = 0;
    if (!utilities::To<string, float>(meth_row[mls_index], mls))
       LOG_ERROR_P(PARAM_METHOD_INFO) << PARAM_MINIMUM_LEGAL_LENGTH << " value " << meth_row[mls_index] << " is not numeric, please sort this out.";
    if (!utilities::To<string, float>(meth_row[hand_mort_index], hand_mort))
       LOG_ERROR_P(PARAM_METHOD_INFO) << PARAM_HANDLING_MORTALITY << " value " << meth_row[hand_mort_index] << " is not numeric, please sort this out.";
    fishery_mls_[fishery_index_[meth_index]] = mls;
    fishery_hand_mort_[fishery_index_[meth_index]] = hand_mort;
  }

  if (scanning_table_->has_been_defined()) {
    scanning_ = true;
    auto scan_data = scanning_table_->data();
    auto scan_columns = scanning_table_->columns();

    if (scan_data.size() != catch_values_data.size())
      LOG_FATAL_P(PARAM_SCANNING) << "found " << scan_data.size() << " rows in the scanning table but " << catch_values_data.size() << " rows in the catch table, these have to be the same";
    if (scan_data[0].size() != catch_values_data[0].size())
      LOG_FATAL_P(PARAM_SCANNING) << "found " << scan_data[0].size() << " columns in the scanning table but " << catch_values_data[0].size() << " columns in the catch table, these have to be the same";
    vector<unsigned> fish_index_map;
    for(unsigned i = 1; i < scan_columns.size(); ++i) {
      fish_index_map.push_back(std::find(fishery_label_.begin(), fishery_label_.end(), scan_columns[i]) - fishery_label_.begin());
    }

    for (auto scan_row : scan_data) {
      unsigned year = 0;
      if (!utilities::To<string, unsigned>(scan_row[0], year))
        LOG_FATAL_P(PARAM_SCANNING) << "year value " << scan_row[0] << " is not numeric.";
      if (std::find(catch_year_.begin(), catch_year_.end(), year) == catch_year_.end())
        LOG_FATAL_P(PARAM_SCANNING) << "year " << year << " is not found in the catch table, there needs to be a year for every catch year.";
      scanning_years_.push_back(year);

      float prop = 0;
      for (unsigned i = 1; i < scan_row.size(); ++i) {
        if (!utilities::To<string, float>(scan_row[i], prop))
          LOG_FATAL_P(PARAM_SCANNING) << "proportion value " << scan_row[i] << " is not numeric.";
        if (prop < 0 || prop > 1) {
          LOG_FATAL_P(PARAM_SCANNING) << "found a proportion less then zero or greater the one, please sort this out";
        }
        scanning_proportion_by_fishery_[fish_index_map[i - 1]].push_back(prop);
      }
    }
  }

  for (unsigned i = 0; i < age_comp_by_fishery_.size(); ++i) {
    age_comp_by_fishery_[i].resize( catch_year_.size());
    length_comp_by_fishery_[i].resize( catch_year_.size());
    fishery_census_data_[i].resize( catch_year_.size());
  }

  scanning_this_year_.resize(fishery_index_.size(), false);
  LOG_FINE() << "finished building Tables";

  catch_to_take_by_fishery_.resize(fishery_index_.size(),0.0);
}


/**
 * DoReset
 */
void MortalityEventBiomass::DoReset() {
  LOG_FINE() << "clearing containers";
  removals_by_age_and_area_.clear();
  removals_by_length_and_area_.clear();
  removals_census_.clear();
  removals_tag_recapture_.clear();
}

/**
 * DoExecute
 */
void MortalityEventBiomass::DoExecute() {
  LOG_MEDIUM();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance(); // shared resource
  vector<unsigned> global_age_freq(model_->age_spread(), 0);
  unsigned time_step = model_->get_time_step_counter();

  auto iter = catch_year_.begin();
  if (model_->state() != State::kInitialise) {
    if (find(iter, catch_year_.end(), model_->current_year()) != catch_year_.end()) {
      iter = find(catch_year_.begin(), catch_year_.end(), model_->current_year());
      unsigned catch_ndx = distance(catch_year_.begin(), iter);

      if (scanning_years_.size() > 0) {
        LOG_FINE() << "checking scanning";
        LOG_FINE() << "scanning_proportion_by_fishery_.size() " << scanning_proportion_by_fishery_.size() << " scan years " << scanning_years_.size();
        auto scan_iter = scanning_years_.begin();
        scan_iter = find(scanning_years_.begin(), scanning_years_.end(), model_->current_year());
        unsigned scan_ndx = distance(scanning_years_.begin(), scan_iter);
        LOG_FINE() << "scan ndx = " << scan_ndx;

        for(unsigned i = 0; i < fishery_index_.size(); ++i) {
          LOG_FINE() << "i = " << scanning_proportion_by_fishery_[i][scan_ndx];
          if (scanning_proportion_by_fishery_[i][scan_ndx] > 0) {
            scanning_this_year_[i] = true;
            LOG_FINE() << "scanning for fishery " << fishery_label_[i] << " in year " << model_->current_year();
          } else {
            scanning_this_year_[i] = false;
          }
        }
      }

      LOG_FINE() << "about to kick into the gear, are we scanning (1 = yes) " << scanning_;
      float world_catch_to_take = 0;
      if (not selectivity_length_based_) {
        LOG_FINE() << "Age based F";
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            unsigned catch_attempts = 1;
            unsigned catch_max = 1;
            WorldCell* cell = nullptr;
            float catch_taken = 0;
            cell = world_->get_base_square(row, col); // Shared resource...
            // Which fisheries are we taken catches for, not all fisheries are taking catch every year.
            vector<unsigned> fisheries_to_sample_from;
            if (cell->is_enabled()) {
              fill(catch_to_take_by_fishery_.begin(),catch_to_take_by_fishery_.end(),0.0);
              float catch_for_fishery = 0.0;
              LOG_FINE() << "cell: " << row << "-" << col << " what we are storing = " << catch_to_take_by_fishery_.size() << " looping over " << fishery_catch_layer_.size();
              for (unsigned i = 0; i < fishery_catch_layer_.size(); ++i) {
                catch_for_fishery = fishery_catch_layer_[i][catch_ndx]->get_value(row, col);
                fishery_catch_to_take_[i][catch_ndx] += catch_for_fishery;
                LOG_FINE() << "in year " << model_->current_year() << " fishery " << fishery_label_[i] << "need to remove " << catch_for_fishery;
                catch_taken += catch_for_fishery;
                composition_data age_freq_for_fishery(PARAM_AGE, model_->current_year(), row, col, model_->age_spread());
                age_comp_by_fishery_[i][catch_ndx].push_back(age_freq_for_fishery);
                composition_data length_freq_for_fishery(PARAM_LENGTH, model_->current_year(), row, col, model_->number_of_length_bins());
                length_comp_by_fishery_[i][catch_ndx].push_back(length_freq_for_fishery);
                census_data census_fishery_specific(model_->current_year(), row, col);
                fishery_census_data_[i][catch_ndx].push_back(census_fishery_specific);
                catch_to_take_by_fishery_[i] = catch_for_fishery;
                if (catch_for_fishery > 0.0)
                  fisheries_to_sample_from.push_back(i);

                cell_ndx_[i] = age_comp_by_fishery_[i][catch_ndx].size() - 1;
              }

              vector<float> prop_catch_by_fishery(fisheries_to_sample_from.size(),0.0);
              for (unsigned i = 0; i < prop_catch_by_fishery.size(); ++i) {
                prop_catch_by_fishery[i] = catch_to_take_by_fishery_[fisheries_to_sample_from[i]] / catch_taken;
                LOG_FINE() << "fishery ndx = " << fisheries_to_sample_from[i] << " proportion of catch = " << prop_catch_by_fishery[i];
              }

              world_catch_to_take += catch_taken;
              LOG_FINE() << "need to remove " << catch_taken << " weight of fish across all fisheries";
              unsigned counter = 0;
              if (catch_taken > 0) {
                LOG_FINEST() << "We are fishing in cell " << row + 1 << " " << col + 1 << " value = " << catch_taken;
                census_data census_fishery(model_->current_year(), row, col);
                tag_recapture tag_recapture_info(model_->current_year(), row, col, time_step);

                composition_data age_freq(PARAM_AGE,  model_->current_year(), row, col, model_->age_spread());
                composition_data length_freq(PARAM_LENGTH,  model_->current_year(), row, col, model_->number_of_length_bins());

                catch_attempts = 1;
                catch_max = cell->agents_.size() * 2;
                LOG_FINEST() << "individuals = " << catch_max;
                unsigned fishery_ndx = 0;
                float random_fish, temp_sum;
                double individuals_in_cell = cell->get_total_individuals_alive();
                double double_check = 0;
                for(auto& agent : cell->agents_) {
                  if (agent.is_alive())
                    double_check += agent.get_scalar();
                }
                LOG_FINE() <<"individuals_in_cell = " << individuals_in_cell << " but we count = " << double_check;
                if (individuals_in_cell > 0)
                  LOG_FINE() << "do this";
                //unsigned individual_ndx = 0;
                /*
                 *  Main loop
                */
                while (catch_taken > 0) {
                  LOG_FINEST() << "catch taken = " << catch_taken << " attempts = " << catch_attempts;
                  ++catch_attempts;
                  // Account for unequal probability
                  // randomly find agent
                  auto& this_agent = cell->agents_[rng.chance() * cell->agents_.size()];
                  if (this_agent.is_alive()) {
                    // calculate fishery based on proportion of catch
                    temp_sum = 0.0;
                    random_fish = rng.chance();
                    for (unsigned i = 0; i < fisheries_to_sample_from.size(); ++i) {
                      temp_sum += prop_catch_by_fishery[i];
                      if (temp_sum >= random_fish) {
                        fishery_ndx = fisheries_to_sample_from[i];
                        break;
                      }
                    }
                    // Do we need to take catch from this fishery
                    if (catch_to_take_by_fishery_[fishery_ndx] > 0) {
                      if (rng.chance() <= fishery_selectivity_[fishery_ndx][this_agent.get_sex()]->GetResult(this_agent.get_age_index())) {
                        // check under MLS and apply some handling mortality
                        if (this_agent.get_length() < fishery_mls_[fishery_ndx]) {
                          if (rng.chance() <= fishery_hand_mort_[fishery_ndx]) {
                            cell->remove_agent_alive(this_agent.get_scalar());
                            this_agent.dies();
                          }
                        } else {
                          // record information
                          //LOG_FINE() << catch_taken << " " <<  this_agent.get_sex() << " " << this_agent.get_weight() << " " << this_agent.get_scalar();
                          catch_taken -= this_agent.get_weight() * this_agent.get_scalar();
                          catch_to_take_by_fishery_[fishery_ndx] -= this_agent.get_weight() * this_agent.get_scalar();
                          fishery_actual_catch_taken_[fishery_ndx][catch_ndx] += this_agent.get_weight() * this_agent.get_scalar();
                          age_freq.frequency_[this_agent.get_age_index()] += this_agent.get_scalar(); // This actually represents many individuals.
                          length_freq.frequency_[this_agent.get_length_bin_index()] += this_agent.get_scalar();

                          if (this_agent.get_sex() == 0) {
                            age_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].frequency_[this_agent.get_age_index()] += this_agent.get_scalar();
                            length_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].frequency_[this_agent.get_length_bin_index()] += this_agent.get_scalar();
                            fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].biomass_+= this_agent.get_weight() * this_agent.get_scalar();
                          } else {
                            age_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_frequency_[this_agent.get_age_index()] += this_agent.get_scalar();
                            length_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_frequency_[this_agent.get_length_bin_index()] += this_agent.get_scalar();
                            fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_biomass_+= this_agent.get_weight() * this_agent.get_scalar();
                          }
                          global_age_freq[this_agent.get_age_index()] += this_agent.get_scalar();
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].age_ndx_.push_back(this_agent.get_age_index());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].sex_.push_back(this_agent.get_sex());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].fishery_ndx_.push_back(fishery_ndx);
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].length_ndx_.push_back(this_agent.get_length_bin_index());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].scalar_.push_back(this_agent.get_scalar());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].weight_.push_back(this_agent.get_weight());

                          if (scanning_ & scanning_this_year_[fishery_ndx]) {
                            // Probability of scanning agent
                            if (rng.chance() <= scanning_proportion_by_fishery_[fishery_ndx][catch_ndx]) {
                              // We scanned this agent
                              tag_recapture_info.scanned_fish_+= this_agent.get_scalar();
                              if (this_agent.get_number_tags() > 0) {
                                //fish has a tag record it
                                tag_recapture_info.age_.push_back(this_agent.get_age());
                                tag_recapture_info.sex_.push_back(this_agent.get_sex());
                                tag_recapture_info.length_.push_back(this_agent.get_length());
                                tag_recapture_info.fishery_ndx_.push_back(fishery_ndx);
                                tag_recapture_info.time_at_liberty_.push_back(this_agent.get_time_at_liberty(time_step));
                                tag_recapture_info.length_increment_.push_back(this_agent.get_length_increment_since_tag());
                                tag_recapture_info.tag_row_.push_back(this_agent.get_tag_row());
                                tag_recapture_info.tag_col_.push_back(this_agent.get_tag_col());
                              }
                            }
                          }
                        }
                        cell->remove_agent_alive(this_agent.get_scalar());
                        this_agent.dies();
                      }
                    }
                    // Make sure we don't end up fishing for infinity
                    if (catch_attempts >= catch_max) {
                     LOG_WARNING() << "Too many attempts to catch an agent in the process " << label_ << " in year " << model_->current_year() << " in row " << row + 1 << " and column " << col + 1 << ", remaining catch to take = " << catch_taken << " this most likely means you have" <<
                         " a model that suggests there should be more agents in this space than than the current agent dynamics are putting in this cell, check the user manual for tips to resolve this situation. attempts = " << catch_attempts << " max attempts allowed = " << catch_max;
                     // Kick out of this cell
                     break;
                    }
                    ++counter;
                  }
                } // while (catch_taken > 0

                for (auto& catch_ : catch_to_take_by_fishery_) {
                  LOG_FINE() << catch_;
                }
                removals_by_length_and_area_.push_back(length_freq);
                removals_by_age_and_area_.push_back(age_freq);
                removals_census_.push_back(census_fishery);
                if (scanning_ & (tag_recapture_info.age_.size() > 0)) {
                  LOG_FINE() << "saving tag-recapture";
                  removals_tag_recapture_.push_back(tag_recapture_info);
                }
              } //if (catch_taken > 0) {
            } // cell ->is_enabled()
          } // col
        } // row
        } else {
        LOG_FINE() << "Applying length based F";
        float world_catch_to_take = 0;
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            unsigned catch_attempts = 1;
            unsigned catch_max = 1;
            WorldCell* cell = nullptr;
            float catch_taken = 0;
            cell = world_->get_base_square(row, col); // Shared resource...
            // Which fisheries are we taken catches for, not all fisheries are taking catch every year.
            vector<unsigned> fisheries_to_sample_from;
            if (cell->is_enabled()) {
              float catch_for_fishery = 0.0;
              for (unsigned i = 0; i < fishery_catch_layer_.size(); ++i) {
                LOG_FINE() << "in year " << model_->current_year() << " fishery " << fishery_label_[i];
                catch_for_fishery = fishery_catch_layer_[i][catch_ndx]->get_value(row, col);
                fishery_catch_to_take_[i][catch_ndx] += catch_for_fishery;
                catch_taken += catch_for_fishery;
                LOG_FINE() << "need to remove " << catch_for_fishery;
                composition_data age_freq_for_fishery(PARAM_AGE, model_->current_year(), row, col, model_->age_spread());
                age_comp_by_fishery_[i][catch_ndx].push_back(age_freq_for_fishery);
                composition_data length_freq_for_fishery(PARAM_LENGTH, model_->current_year(), row, col, model_->number_of_length_bins());
                length_comp_by_fishery_[i][catch_ndx].push_back(length_freq_for_fishery);
                census_data census_fishery_specific(model_->current_year(), row, col);
                fishery_census_data_[i][catch_ndx].push_back(census_fishery_specific);
                catch_to_take_by_fishery_[i] = catch_for_fishery;
                if (catch_for_fishery > 0.0)
                  fisheries_to_sample_from.push_back(i);

                cell_ndx_[i] = age_comp_by_fishery_[i][catch_ndx].size() - 1;
              }

              vector<float> prop_catch_by_fishery(fisheries_to_sample_from.size(),0.0);
              for (unsigned i = 0; i < prop_catch_by_fishery.size(); ++i) {
                prop_catch_by_fishery[i] = catch_to_take_by_fishery_[fisheries_to_sample_from[i]] / catch_taken;
                LOG_FINE() << "fishery ndx = " << fisheries_to_sample_from[i] << " proportion of catch = " << prop_catch_by_fishery[i];
              }

              world_catch_to_take += catch_taken;
              LOG_FINE() << "need to remove " << catch_taken << " weight of fish across all fisheries";
              unsigned counter = 0;
              if (catch_taken > 0) {
                LOG_FINE() << "We are fishing in cell " << row + 1 << " " << col + 1 << " value = " << catch_taken;
                census_data census_fishery(model_->current_year(), row, col);
                tag_recapture tag_recapture_info(model_->current_year(), row, col, time_step);

                composition_data age_freq(PARAM_AGE,  model_->current_year(), row, col, model_->age_spread());
                composition_data length_freq(PARAM_LENGTH,  model_->current_year(), row, col, model_->number_of_length_bins());

                catch_attempts = 1;
                catch_max = cell->agents_.size() * 50;
                LOG_FINEST() << "individuals = " << catch_max;
                unsigned fishery_ndx = 0;
                float random_fish, temp_sum;

                /*
                 *  Main loop
                */
                while (catch_taken > 0) {
                  LOG_FINEST() << "catch taken = " << catch_taken << " attempts = " << catch_attempts;
                  ++catch_attempts;
                  // randomly find agent
                  auto& this_agent = cell->agents_[rng.chance() * cell->agents_.size()];
                  if (this_agent.is_alive()) {
                    // calculate fishery based on proportion of catch
                    temp_sum = 0.0;
                    random_fish = rng.chance();
                    for (unsigned i = 0; i < fisheries_to_sample_from.size(); ++i) {
                      temp_sum += prop_catch_by_fishery[i];
                      if (temp_sum >= random_fish) {
                        fishery_ndx = fisheries_to_sample_from[i];
                        break;
                      }
                    }
                    // Do we need to take catch from this fishery
                    if (catch_to_take_by_fishery_[fishery_ndx] > 0) {
                      if (rng.chance() <= fishery_selectivity_[fishery_ndx][this_agent.get_sex()]->GetResult(this_agent.get_length_bin_index())) {
                        // check under MLS and apply some handling mortality
                        if (this_agent.get_length() < fishery_mls_[fishery_ndx]) {
                          if (rng.chance() <= fishery_hand_mort_[fishery_ndx]) {
                            cell->remove_agent_alive(this_agent.get_scalar());
                            this_agent.dies();
                          }
                        } else {
                          // record information
                          //LOG_FINE() << catch_taken << " " <<  this_agent.get_sex() << " " << this_agent.get_weight() << " " << this_agent.get_scalar();
                          catch_taken -= this_agent.get_weight() * this_agent.get_scalar();
                          catch_to_take_by_fishery_[fishery_ndx] -= this_agent.get_weight() * this_agent.get_scalar();
                          fishery_actual_catch_taken_[fishery_ndx][catch_ndx] += this_agent.get_weight() * this_agent.get_scalar();
                          age_freq.frequency_[this_agent.get_age_index()] += this_agent.get_scalar(); // This actually represents many individuals.
                          length_freq.frequency_[this_agent.get_length_bin_index()] += this_agent.get_scalar();

                          if (this_agent.get_sex() == 0) {
                            age_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].frequency_[this_agent.get_age_index()] += this_agent.get_scalar();
                            length_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].frequency_[this_agent.get_length_bin_index()] += this_agent.get_scalar();
                            fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].biomass_+= this_agent.get_weight() * this_agent.get_scalar();
                          } else {
                            age_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_frequency_[this_agent.get_age_index()] += this_agent.get_scalar();
                            length_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_frequency_[this_agent.get_length_bin_index()] += this_agent.get_scalar();
                            fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_biomass_+= this_agent.get_weight() * this_agent.get_scalar();
                          }
                          global_age_freq[this_agent.get_age_index()] += this_agent.get_scalar();
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].age_ndx_.push_back(this_agent.get_age_index());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].sex_.push_back(this_agent.get_sex());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].fishery_ndx_.push_back(fishery_ndx);
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].length_ndx_.push_back(this_agent.get_length_bin_index());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].scalar_.push_back(this_agent.get_scalar());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].weight_.push_back(this_agent.get_weight());

                          if (scanning_ & scanning_this_year_[fishery_ndx]) {
                            // Probability of scanning agent
                            if (rng.chance() <= scanning_proportion_by_fishery_[fishery_ndx][catch_ndx]) {
                              // We scanned this agent
                              tag_recapture_info.scanned_fish_+= this_agent.get_scalar();
                              if (this_agent.get_number_tags() > 0) {
                                //fish has a tag record it
                                tag_recapture_info.age_.push_back(this_agent.get_age());
                                tag_recapture_info.sex_.push_back(this_agent.get_sex());
                                tag_recapture_info.length_.push_back(this_agent.get_length());
                                tag_recapture_info.fishery_ndx_.push_back(fishery_ndx);
                                tag_recapture_info.time_at_liberty_.push_back(this_agent.get_time_at_liberty(time_step));
                                tag_recapture_info.length_increment_.push_back(this_agent.get_length_increment_since_tag());
                                tag_recapture_info.tag_row_.push_back(this_agent.get_tag_row());
                                tag_recapture_info.tag_col_.push_back(this_agent.get_tag_col());
                              }
                            }
                          }
                        }
                        cell->remove_agent_alive(this_agent.get_scalar());
                        this_agent.dies();
                      }
                    }
                    // Make sure we don't end up fishing for infinity
                    if (catch_attempts >= catch_max) {
                     LOG_WARNING() << "Too many attempts to catch an agent in the process " << label_ << " in year " << model_->current_year() << " in row " << row + 1 << " and column " << col + 1 << ", remaining catch to take = " << catch_taken << " this most likely means you have" <<
                         " a model that suggests there should be more agents in this space than than the current agent dynamics are putting in this cell, check the user manual for tips to resolve this situation. attempts = " << catch_attempts << " max attempts allowed = " << catch_max;
                     // Kick out of this cell
                     break;
                    }
                    ++counter;
                  }
                } // while (catch_taken > 0

                for (auto& catch_ : catch_to_take_by_fishery_)
                  LOG_FINE() << catch_;
                removals_by_length_and_area_.push_back(length_freq);
                removals_by_age_and_area_.push_back(age_freq);
                removals_census_.push_back(fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]]);
                if (scanning_ & (tag_recapture_info.age_.size() > 0)) {
                  removals_tag_recapture_.push_back(tag_recapture_info);
                }
              } //if (catch_taken > 0) {
            } // cell ->is_enabled()
          } // col
        } // row
      } // length based
    } // find(years_.begin(), years_.end(), model_->current_year()) != years_.end()
  }  //model_->state() != State::kInitialise
  LOG_FINE() << "finished Biomass Mort process.";
}


// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MortalityEventBiomass::FillReportCache(ostringstream& cache) {
  LOG_FINE() << "";
  // Fishery specific info
  // actual catches
  if (fishery_actual_catch_taken_[0].size() > 0) {
    cache << "actual_catches " << REPORT_R_DATAFRAME << "\n";
    cache << "year ";
    for (auto& fishery : fishery_label_)
      cache << fishery << " ";
    cache << "\n";
    for (unsigned i = 0; i < fishery_actual_catch_taken_[0].size(); ++i){
      cache << catch_year_[i] << " ";
      for (auto& fishery_ndx : fishery_index_)
        cache << fishery_actual_catch_taken_[fishery_ndx][i] << " ";
      cache << "\n";
    }
  }

  if (fishery_catch_to_take_[0].size() > 0) {
    cache << "catches " << REPORT_R_DATAFRAME << "\n";
    cache << "year ";
    for (auto& fishery : fishery_label_)
      cache << fishery << " ";
    cache << "\n";
    for (unsigned i = 0; i < fishery_catch_to_take_[0].size(); ++i){
      cache << catch_year_[i] << " ";
      for (auto& fishery_ndx : fishery_index_)
        cache << fishery_catch_to_take_[fishery_ndx][i] << " ";
      cache << "\n";
    }
  }

  // Get rid of this later, but check census data is the same as age-freq which has been validated
  if (fishery_census_data_.size() > 0) {
    if (model_->get_sexed()) {
      for (auto& fishery : fishery_index_) {
        cache << "age_freq-census-male-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year ";
        for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age)
          cache << age << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < fishery_census_data_[fishery].size(); ++year_ndx) {
          vector<unsigned> temp_age_freq(model_->age_spread(), 0.0);
          cache << catch_year_[year_ndx] << " ";
          for (unsigned cell_ndx = 0; cell_ndx < fishery_census_data_[fishery][year_ndx].size(); ++cell_ndx) {
            if (fishery_census_data_[fishery][year_ndx][cell_ndx].year_ == catch_year_[year_ndx]) {
              for (unsigned age_ndx = 0; age_ndx < fishery_census_data_[fishery][year_ndx][cell_ndx].age_ndx_.size(); ++age_ndx) {
                if (fishery_census_data_[fishery][year_ndx][cell_ndx].sex_[age_ndx] == 0)
                  temp_age_freq[fishery_census_data_[fishery][year_ndx][cell_ndx].age_ndx_[age_ndx]]++;
              }
            }
          }
          for (auto& age_freq : temp_age_freq)
            cache << age_freq << " ";
          cache << "\n";
        }
      }
      for (auto& fishery : fishery_index_) {
        cache << "age_freq-census-female-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year ";
        for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age)
          cache << age << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < fishery_census_data_[fishery].size(); ++year_ndx) {
          vector<unsigned> temp_age_freq(model_->age_spread(), 0.0);
          cache << catch_year_[year_ndx] << " ";
          for (unsigned cell_ndx = 0; cell_ndx < fishery_census_data_[fishery][year_ndx].size(); ++cell_ndx) {
            if (fishery_census_data_[fishery][year_ndx][cell_ndx].year_ == catch_year_[year_ndx]) {
              for (unsigned age_ndx = 0; age_ndx < fishery_census_data_[fishery][year_ndx][cell_ndx].age_ndx_.size(); ++age_ndx) {
                if (fishery_census_data_[fishery][year_ndx][cell_ndx].sex_[age_ndx] == 1)
                  temp_age_freq[fishery_census_data_[fishery][year_ndx][cell_ndx].age_ndx_[age_ndx]]++;
              }
            }
          }
          for (auto& age_freq : temp_age_freq)
            cache << age_freq << " ";
          cache << "\n";
        }
      }
    }
  }

  // age frequency by sex fishery and year
  if (age_comp_by_fishery_.size() > 0) {

    if (model_->get_sexed()) {
      for (auto& fishery : fishery_index_) {
        cache << "age_freq-male-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year ";
        for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age)
          cache << age << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < age_comp_by_fishery_[fishery].size(); ++year_ndx) {
          vector<unsigned> temp_age_freq(model_->age_spread(), 0.0);
          cache << catch_year_[year_ndx] << " ";
          // Accumulate over all cells
          for (unsigned cell_ndx = 0; cell_ndx < age_comp_by_fishery_[fishery][year_ndx].size(); ++cell_ndx) {
            if (age_comp_by_fishery_[fishery][year_ndx][cell_ndx].year_ == catch_year_[year_ndx]) {
              for (unsigned age_ndx = 0; age_ndx < age_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_.size(); ++age_ndx)
                temp_age_freq[age_ndx] += age_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_[age_ndx];
            }
          }
          for (auto& age_freq : temp_age_freq)
            cache << age_freq << " ";
          cache << "\n";
        }
      }

      for (auto& fishery : fishery_index_) {
        cache << "age_freq-female-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year ";
        for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age)
          cache << age << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < age_comp_by_fishery_[fishery].size(); ++year_ndx) {
          vector<unsigned> temp_age_freq(model_->age_spread(), 0.0);
          cache << catch_year_[year_ndx] << " ";
          // Accumulate over all cells
          for (unsigned cell_ndx = 0; cell_ndx < age_comp_by_fishery_[fishery][year_ndx].size(); ++cell_ndx) {
            if (age_comp_by_fishery_[fishery][year_ndx][cell_ndx].year_ == catch_year_[year_ndx]) {
              for (unsigned age_ndx = 0; age_ndx < age_comp_by_fishery_[fishery][year_ndx][cell_ndx].female_frequency_.size(); ++age_ndx)
                temp_age_freq[age_ndx] += age_comp_by_fishery_[fishery][year_ndx][cell_ndx].female_frequency_[age_ndx];
            }
          }
          for (auto& age_freq : temp_age_freq)
            cache << age_freq << " ";
          cache << "\n";
        }
      }

    } else {
      for (auto& fishery : fishery_index_) {
        cache << "age_freq-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year ";
        for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age)
          cache << age << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < age_comp_by_fishery_[fishery].size(); ++year_ndx) {
          vector<unsigned> temp_age_freq(model_->age_spread(), 0.0);
          cache << catch_year_[year_ndx] << " ";
          // Accumulate over all cells
          for (unsigned cell_ndx = 0; cell_ndx < age_comp_by_fishery_[fishery][year_ndx].size(); ++cell_ndx) {
            if (age_comp_by_fishery_[fishery][year_ndx][cell_ndx].year_ == catch_year_[year_ndx]) {
              for (unsigned age_ndx = 0; age_ndx < age_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_.size(); ++age_ndx)
                temp_age_freq[age_ndx] += age_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_[age_ndx];
            }
          }
          for (auto& age_freq : temp_age_freq)
            cache << age_freq << " ";
          cache << "\n";
        }
      }
    }
  }

  if (removals_by_age_.size() > 0) {
    cache << "age_frequency " << REPORT_R_DATAFRAME << "\n";
    cache << "year ";
    for (unsigned i = model_->min_age(); i <= model_->max_age(); ++i)
      cache << i << " ";
    cache << "\n";
    for (auto& age_freq : removals_by_age_) {
      cache << age_freq.first << " ";
      for (auto age_value : age_freq.second)
        cache << age_value << " ";
      cache << "\n";
    }
  }

  if (print_extra_info_) {
    if (removals_tag_recapture_.size() > 0) {
      for (auto& tag_recapture : removals_tag_recapture_) {
        cache << "tag_recapture_info-" << tag_recapture.year_ << "-" << tag_recapture.row_ << "-" << tag_recapture.col_ << " " << REPORT_R_LIST << "\n";
        cache << "scanned_fish: " << tag_recapture.scanned_fish_ << "\n";
        cache << "values " << REPORT_R_MATRIX << "\n";

        LOG_MEDIUM() << "ages = " << tag_recapture.age_.size();

        //cache << "age length length-increment time_at_liberty\n";
        for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
          cache << tag_recapture.age_[ndx] << " ";
        cache << "\n";
        for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
          cache << tag_recapture.length_[ndx] << " ";
        cache << "\n";
        for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
          cache << tag_recapture.time_at_liberty_[ndx] << " ";
        cache << "\n";
        for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
          cache << tag_recapture.length_increment_[ndx] << " ";
        cache << "\n";
        for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
          cache << tag_recapture.tag_row_[ndx] << " ";
        cache << "\n";
        for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
          cache << tag_recapture.tag_col_[ndx] << " ";
        cache << "\n" << REPORT_R_LIST_END << "\n";
      }
    }

    LOG_FINE() << fishery_census_data_.size() << " census vals";
    // Print census information
    for(unsigned fish_ndx = 0; fish_ndx < fishery_census_data_.size(); ++fish_ndx) {
      for(unsigned year_ndx = 0; year_ndx < fishery_census_data_[fish_ndx].size(); ++year_ndx) {
        for(unsigned cell_ndx = 0; cell_ndx < fishery_census_data_[fish_ndx][year_ndx].size(); ++cell_ndx) {
          cache << "census_info-" << fishery_label_[fishery_index_[fish_ndx]] << "-" << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].year_ << "-" << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].row_ << "-" << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].col_ << " " << REPORT_R_MATRIX << "\n";
          for (unsigned age_ndx = 0; age_ndx < fishery_census_data_[fish_ndx][year_ndx][cell_ndx].age_ndx_.size(); ++age_ndx)
            cache << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].age_ndx_[age_ndx] + model_->min_age() << " ";
          cache << "\n";
          for (unsigned length_ndx = 0; length_ndx < fishery_census_data_[fish_ndx][year_ndx][cell_ndx].length_ndx_.size(); ++length_ndx)
            cache << model_->length_bin_mid_points()[fishery_census_data_[fish_ndx][year_ndx][cell_ndx].length_ndx_[length_ndx]] << " ";
          cache << "\n";
        }
      }
    }
  }

}

// true means all years are found, false means there is a mismatch in years
bool MortalityEventBiomass::check_years(vector<unsigned> years_to_check_) {
  LOG_FINE();
  for (unsigned year_ndx = 0; year_ndx < years_to_check_.size(); ++year_ndx) {
    if (find(catch_year_.begin(), catch_year_.end(), years_to_check_[year_ndx]) == catch_year_.end())
      return false;
  }

  return true;

}
// Does this fishery exist?
// yes = true
// no = false
bool MortalityEventBiomass::check_fishery_exists(string fishery_label) {
  LOG_FINE();
  if (find(fishery_label_.begin(), fishery_label_.end(), fishery_label) != fishery_label_.end())
    return true;

  return false;

}

// find and return age comp data for this fishery.
vector<vector<census_data>>  MortalityEventBiomass::get_fishery_census_data(string fishery_label) {
  LOG_FINE();
  vector<vector<census_data>> fish_comp;
  for (unsigned i = 0; i < fishery_label_.size(); ++i) {
    if (fishery_label_[i] == fishery_label)
      fish_comp = fishery_census_data_[fishery_index_[i]];
  }
  return fish_comp;
}

} /* namespace processes */
} /* namespace niwa */
