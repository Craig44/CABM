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
MortalityEventBiomass::MortalityEventBiomass(Model *model) :
    Mortality(model) {
  catch_table_ = new parameters::Table(PARAM_CATCH);
  parameters_.BindTable(PARAM_CATCH, catch_table_, "Table of catch_layer by year for each fishery", "", true, false);
  scanning_table_ = new parameters::Table(PARAM_SCANNING);
  parameters_.BindTable(PARAM_SCANNING, scanning_table_, "Table of scanning proportions by year for each fishery", "", true, true);
  method_table_ = new parameters::Table(PARAM_METHOD_INFO);
  parameters_.BindTable(PARAM_METHOD_INFO, method_table_, "Table of information for each fishery.", "", true, false);
  parameters_.Bind<bool>(PARAM_PRINT_CENSUS_INFO, &print_census_info_, "if you have process report for this process you can control the amount of information printed to the file.", "", true);
  parameters_.Bind<bool>(PARAM_PRINT_TAG_RECAP_INFO, &print_tag_recap_info_, "if you have process report for this process you can control the amount of information printed to the file.", "", true);

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
    LOG_FATAL_P(PARAM_METHOD_INFO)
    << "The column header " << PARAM_HANDLING_MORTALITY << ", is missing can you please add it";
  }
  if (std::find(meth_columns.begin(), meth_columns.end(), PARAM_MINIMUM_LEGAL_LENGTH) == meth_columns.end()) {
    LOG_FATAL_P(PARAM_METHOD_INFO)
    << "The column header " << PARAM_MINIMUM_LEGAL_LENGTH << ", is missing can you please add it";
  }

  if (meth_columns[0] != PARAM_METHOD) {
    LOG_FATAL_P(PARAM_METHOD_INFO)
    << "The first column header must be labelled " << PARAM_METHOD << ", can you please sort this out";
  }

  if (!model_->get_sexed()) {
    // find the column header selectivity
    if (std::find(meth_columns.begin(), meth_columns.end(), PARAM_SELECTIVITY) == meth_columns.end()) {
      LOG_FATAL_P(PARAM_METHOD_INFO)
      << "because the model is unsexed we expect the column header " << PARAM_SELECTIVITY << ", this is missing can you please add it";
    }
  } else {
    if (std::find(meth_columns.begin(), meth_columns.end(), PARAM_MALE_SELECTIVITY) == meth_columns.end()) {
      LOG_FATAL_P(PARAM_METHOD_INFO)
      << "because the model is sexed we expect the column header " << PARAM_MALE_SELECTIVITY << ", this is missing can you please add it";
    }
    if (std::find(meth_columns.begin(), meth_columns.end(), PARAM_FEMALE_SELECTIVITY) == meth_columns.end()) {
      LOG_FATAL_P(PARAM_METHOD_INFO)
      << "because the model is sexed we expect the column header " << PARAM_FEMALE_SELECTIVITY << ", this is missing can you please add it";
    }
  }

  auto meth_data = method_table_->data();
  LOG_FINE() << "rows in meth info table " << meth_data.size() << " length fishery label " << fishery_label_.size();
  for (unsigned i = 0; i < meth_data.size(); ++i) {
    LOG_FINE() << "trying to check this " << meth_data[i][0] << " fishery is in the catch table " << i;
    if (std::find(fishery_label_.begin(), fishery_label_.end(), (string) meth_data[i][0]) == fishery_label_.end()) {
      LOG_ERROR_P(PARAM_METHOD_INFO) << "Could not find the row label " << meth_data[i][0]
          << " these need to be consistent with the column headers with the catch table, can you please double check this";
    } else {
      LOG_FINE() << "found " << meth_data[i][0] << " in the catch table";
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
        LOG_ERROR_P(PARAM_SCANNING) << "could not find the column header " << scan_columns[j]
            << " in the catch table column headers. The headers of this table must match the column headers of the catch table";
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
  only_mature_partition_.resize(fishery_label_.size());
  auto model_years = model_->years();

  // load objects
  vector<vector<string>> &catch_values_data = catch_table_->data();

  for (auto row : catch_values_data) {
    unsigned year = 0;
    if (!utilities::To<string, unsigned>(row[0], year))
      LOG_ERROR_P(PARAM_CATCH) << "year value " << row[0] << " is not numeric.";
    if (std::find(model_years.begin(), model_years.end(), year) == model_years.end())
      LOG_ERROR_P(PARAM_CATCH) << "year " << year << " is not a valid year in this model";
    // Check years are consecutive ascending order.
    // This will mean when I reference catch_ndx later in the code we can have faith it is in order.
    if (years_.size() > 1) {
      if ((year - 1) != years_[years_.size() - 1]) {
        LOG_ERROR_P(PARAM_CATCH) << "years need to be in consecutive ascending order, the year " << years_[years_.size() - 1] << " was followed by " << year << " please sort this out";
      }
    }
    years_.push_back(year);
    LOG_FINE() << "year = " << year;
    for (unsigned i = 1; i < row.size(); ++i) {
      fishery_catch_layer_labels_[i - 1].push_back(row[i]);
      layers::NumericLayer *temp_layer = nullptr;
      temp_layer = model_->managers().layer()->GetNumericLayer(row[i]);
      if (!temp_layer) {
        LOG_FATAL_P(PARAM_CATCH)
        << "could not find the catch layer '" << row[i] << "', in year " << year << " for fishery " << fishery_label_[i - 1] << " please check it is numeric if it exists";
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

      if (std::find(meth_cols.begin(), meth_cols.end(), PARAM_ONLY_MATURE_PARTITION) == meth_cols.end()) {
        only_mature_partition_[meth_index] = false;
      } else {
        unsigned mature_ndx = std::find(meth_cols.begin(), meth_cols.end(), PARAM_ONLY_MATURE_PARTITION) - meth_cols.begin();
        bool mature_result;

        if (!utilities::To<string, bool>(meth_row[mature_ndx], mature_result))
          LOG_ERROR_P(PARAM_METHOD_INFO) << "couldn't convert column " << PARAM_ONLY_MATURE_PARTITION << " for fishery = " << meth_row[0] << " from string to boolean, check spelling";
        only_mature_partition_[meth_index] = mature_result;
        LOG_MEDIUM() << "mature_ndx = " << mature_ndx << " fishing method = " << meth_row[0] << " result = " << mature_result;

      }
      fishery_selectivity_label_[fishery_index_[meth_index]].push_back(meth_row[selec_index]);
      fishery_selectivity_label_[fishery_index_[meth_index]].push_back(meth_row[selec_index]);
      Selectivity *temp_selectivity = model_->managers().selectivity()->GetSelectivity(meth_row[selec_index]);
      if (!temp_selectivity)
        LOG_FATAL_P(PARAM_METHOD_INFO)
        << ": selectivity " << meth_row[selec_index] << " does not exist. Have you defined it? located in column " << selec_index;

      if (temp_selectivity->is_length_based())
        selectivity_length_based_ = true;

      fishery_selectivity_[fishery_index_[meth_index]].push_back(temp_selectivity);
      fishery_selectivity_[fishery_index_[meth_index]].push_back(temp_selectivity);

    }
  } else {
    unsigned male_selec_index = std::find(meth_cols.begin(), meth_cols.end(), PARAM_MALE_SELECTIVITY) - meth_cols.begin();
    unsigned female_selec_index = std::find(meth_cols.begin(), meth_cols.end(), PARAM_FEMALE_SELECTIVITY) - meth_cols.begin();
    for (auto meth_row : meth_data) {
      unsigned meth_index = std::find(fishery_label_.begin(), fishery_label_.end(), meth_row[0]) - fishery_label_.begin();
      fishery_selectivity_label_[fishery_index_[meth_index]].push_back(meth_row[male_selec_index]);
      fishery_selectivity_label_[fishery_index_[meth_index]].push_back(meth_row[female_selec_index]);
      Selectivity *male_temp_selectivity = model_->managers().selectivity()->GetSelectivity(meth_row[male_selec_index]);
      if (!male_temp_selectivity)
        LOG_FATAL_P(PARAM_METHOD_INFO)
        << ": male_selectivity " << meth_row[male_selec_index] << " does not exist. Have you defined it?";

      if (std::find(meth_cols.begin(), meth_cols.end(), PARAM_ONLY_MATURE_PARTITION) == meth_cols.end()) {
        only_mature_partition_[meth_index] = false;
      } else {
        unsigned mature_ndx = std::find(meth_cols.begin(), meth_cols.end(), PARAM_ONLY_MATURE_PARTITION) - meth_cols.begin();
        bool mature_result;
        if (!utilities::To<string, bool>(meth_row[mature_ndx], mature_result))
          LOG_ERROR_P(PARAM_METHOD_INFO) << "couldn't convert column " << PARAM_ONLY_MATURE_PARTITION << " for fishery = " << meth_row[0] << " from string to boolean, check spelling";
        only_mature_partition_[meth_index] = mature_result;
        LOG_MEDIUM() << "mature_ndx = " << mature_ndx << " fishing method = " << meth_row[0] << " result = " << mature_result;

      }
      Selectivity *female_temp_selectivity = model_->managers().selectivity()->GetSelectivity(meth_row[female_selec_index]);
      if (!female_temp_selectivity)
        LOG_FATAL_P(PARAM_METHOD_INFO)
        << ": female_selectivity " << meth_row[female_selec_index] << " does not exist. Have you defined it?";

      if (female_temp_selectivity->is_length_based() != male_temp_selectivity->is_length_based())
        LOG_ERROR_P(PARAM_METHOD_INFO)
            << "One of the male and one of the female selectivities is age based and the other length based. They both have to be one or the other. Can you please resolve this.";

      if (female_temp_selectivity->is_length_based() && male_temp_selectivity->is_length_based())
        selectivity_length_based_ = true;

      fishery_selectivity_[fishery_index_[meth_index]].push_back(male_temp_selectivity);
      fishery_selectivity_[fishery_index_[meth_index]].push_back(female_temp_selectivity);
    }
  }

  unsigned mls_index = std::find(meth_cols.begin(), meth_cols.end(), PARAM_MINIMUM_LEGAL_LENGTH) - meth_cols.begin();
  unsigned hand_mort_index = std::find(meth_cols.begin(), meth_cols.end(), PARAM_HANDLING_MORTALITY) - meth_cols.begin();
  for (auto meth_row : meth_data) {
    unsigned meth_index = std::find(fishery_label_.begin(), fishery_label_.end(), meth_row[0]) - fishery_label_.begin();
    double mls = 0;
    double hand_mort = 0;
    if (!utilities::To<string, double>(meth_row[mls_index], mls))
      LOG_ERROR_P(PARAM_METHOD_INFO) << PARAM_MINIMUM_LEGAL_LENGTH << " value " << meth_row[mls_index] << " is not numeric, please sort this out.";
    if (!utilities::To<string, double>(meth_row[hand_mort_index], hand_mort))
      LOG_ERROR_P(PARAM_METHOD_INFO) << PARAM_HANDLING_MORTALITY << " value " << meth_row[hand_mort_index] << " is not numeric, please sort this out.";
    fishery_mls_[fishery_index_[meth_index]] = mls;
    fishery_hand_mort_[fishery_index_[meth_index]] = hand_mort;
  }

  if (scanning_table_->has_been_defined()) {
    scanning_ = true;
    auto scan_data = scanning_table_->data();
    auto scan_columns = scanning_table_->columns();

    if (scan_data.size() != catch_values_data.size())
      LOG_FATAL_P(PARAM_SCANNING)
      << "found " << scan_data.size() << " rows in the scanning table but " << catch_values_data.size() << " rows in the catch table, these have to be the same";
    if (scan_data[0].size() != catch_values_data[0].size())
      LOG_FATAL_P(PARAM_SCANNING)
      << "found " << scan_data[0].size() << " columns in the scanning table but " << catch_values_data[0].size() << " columns in the catch table, these have to be the same";
    vector<unsigned> fish_index_map;
    for (unsigned i = 1; i < scan_columns.size(); ++i) {
      fish_index_map.push_back(std::find(fishery_label_.begin(), fishery_label_.end(), scan_columns[i]) - fishery_label_.begin());
    }

    for (auto scan_row : scan_data) {
      unsigned year = 0;
      if (!utilities::To<string, unsigned>(scan_row[0], year))
        LOG_FATAL_P(PARAM_SCANNING)
        << "year value " << scan_row[0] << " is not numeric.";
      if (std::find(years_.begin(), years_.end(), year) == years_.end())
        LOG_FATAL_P(PARAM_SCANNING)
        << "year " << year << " is not found in the catch table, there needs to be a year for every catch year.";
      scanning_years_.push_back(year);

      double prop = 0;
      for (unsigned i = 1; i < scan_row.size(); ++i) {
        if (!utilities::To<string, double>(scan_row[i], prop))
          LOG_FATAL_P(PARAM_SCANNING)
          << "proportion value " << scan_row[i] << " is not numeric.";
        if (prop < 0 || prop > 1) {
          LOG_FATAL_P(PARAM_SCANNING)
          << "found a proportion less then zero or greater the one, please sort this out";
        }
        scanning_proportion_by_fishery_[fish_index_map[i - 1]].push_back(prop);
      }
    }
  }

  for (unsigned i = 0; i < age_comp_by_fishery_.size(); ++i) {
    age_comp_by_fishery_[i].resize(years_.size());
    length_comp_by_fishery_[i].resize(years_.size());
    fishery_census_data_[i].resize(years_.size());
  }

  scanning_this_year_.resize(fishery_index_.size(), false);
  LOG_FINE() << "finished building Tables";

  catch_to_take_by_fishery_.resize(fishery_index_.size(), 0.0);
  actual_catch_by_area_.resize(years_.size());
  for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    actual_catch_by_area_[year_ndx].resize(fishery_label_.size());
    for (unsigned fishery_ndx = 0; fishery_ndx < fishery_label_.size(); ++fishery_ndx) {
      actual_catch_by_area_[year_ndx][fishery_ndx].resize(model_->get_height());
      for (unsigned row = 0; row < model_->get_height(); ++row)
        actual_catch_by_area_[year_ndx][fishery_ndx][row].resize(model_->get_width(), 0.0);
    }
  }
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
  for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    
    for (unsigned fishery_ndx = 0; fishery_ndx < fishery_label_.size(); ++fishery_ndx) {
      fishery_census_data_[fishery_ndx][year_ndx].clear();
      age_comp_by_fishery_[fishery_ndx][year_ndx].clear();
      length_comp_by_fishery_[fishery_ndx][year_ndx].clear();
      for (unsigned row = 0; row < model_->get_height(); ++row)
        fill(actual_catch_by_area_[year_ndx][fishery_ndx][row].begin(), actual_catch_by_area_[year_ndx][fishery_ndx][row].end(), 0.0);
    }
  }
}

/**
 * DoExecute
 */
void MortalityEventBiomass::DoExecute() {
  utilities::RandomNumberGenerator &rng = utilities::RandomNumberGenerator::Instance(); // shared resource
  vector<unsigned> global_age_freq(model_->age_spread(), 0);
  unsigned time_step = model_->get_time_step_counter();
  account_for_tag_population_ = false;
  auto iter = years_.begin();
  scanning_ = false;
  if (model_->state() != State::kInitialise) {
    if (std::find(iter, years_.end(), model_->current_year()) != years_.end()) {
      iter = find(years_.begin(), years_.end(), model_->current_year());
      unsigned catch_ndx = distance(years_.begin(), iter);

      if (scanning_years_.size() > 0) {
        LOG_MEDIUM() << "checking scanning";
        LOG_MEDIUM() << "scanning_proportion_by_fishery_.size() " << scanning_proportion_by_fishery_.size() << " scan years " << scanning_years_.size();
        auto scan_iter = scanning_years_.begin();
        scan_iter = find(scanning_years_.begin(), scanning_years_.end(), model_->current_year());
        unsigned scan_ndx = distance(scanning_years_.begin(), scan_iter);
        LOG_FINE() << "scan ndx = " << scan_ndx;

        for (unsigned i = 0; i < fishery_index_.size(); ++i) {
          LOG_FINE() << "i = " << scanning_proportion_by_fishery_[i][scan_ndx];
          if (scanning_proportion_by_fishery_[i][scan_ndx] > 0) {
            scanning_this_year_[i] = true;
            account_for_tag_population_ = true;
            scanning_ = true;
            LOG_MEDIUM() << "scanning for fishery " << fishery_label_[i] << " in year " << model_->current_year();
          } else {
            scanning_this_year_[i] = false;
          }
        }
      }

      LOG_FINE() << "about to kick into the gear, are we scanning (1 = yes) " << scanning_;
      double world_catch_to_take = 0;
      if (not selectivity_length_based_) {
        LOG_FINE() << "Age based F";
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            unsigned catch_attempts = 1;
            unsigned catch_max = 1;
            WorldCell *cell = nullptr;
            double catch_taken = 0;
            cell = world_->get_base_square(row, col); // Shared resource...
            // Which fisheries are we taken catches for, not all fisheries are taking catch every year.
            vector<unsigned> fisheries_to_sample_from;
            if (cell->is_enabled()) {
              fill(catch_to_take_by_fishery_.begin(), catch_to_take_by_fishery_.end(), 0.0);
              double catch_for_fishery = 0.0;
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
                LOG_MEDIUM() << "ndx = " << cell_ndx_[i] << " for fishery = " << fishery_label_[i];
              }

              vector<double> prop_catch_by_fishery(fisheries_to_sample_from.size(), 0.0);
              for (unsigned i = 0; i < prop_catch_by_fishery.size(); ++i) {
                prop_catch_by_fishery[i] = catch_to_take_by_fishery_[fisheries_to_sample_from[i]] / catch_taken;
                LOG_FINE() << "fishery ndx = " << fisheries_to_sample_from[i] << " proportion of catch = " << prop_catch_by_fishery[i];
              }

              world_catch_to_take += catch_taken;
              LOG_FINE() << "need to remove " << catch_taken << " weight of fish across all fisheries";
              //unsigned agent_counter = 0;
              if (catch_taken > 0) {
                LOG_MEDIUM() << "We are fishing in cell " << row + 1 << " " << col + 1 << " value = " << catch_taken;
                census_data census_fishery(model_->current_year(), row, col);
                tag_recapture tag_recapture_info(model_->current_year(), row, col, time_step);
                tag_recapture_info.scanned_age_comp_.resize(model_->age_spread(), 0.0);
                tag_recapture_info.scanned_length_comp_.resize(model_->number_of_length_bins(), 0.0);
                composition_data age_freq(PARAM_AGE, model_->current_year(), row, col, model_->age_spread());
                composition_data length_freq(PARAM_LENGTH, model_->current_year(), row, col, model_->number_of_length_bins());

                catch_attempts = 1;
                catch_max = cell->agents_.size() * 2;
                LOG_FINEST() << "individuals = " << catch_max;
                unsigned fishery_ndx = 0;
                double random_fish, temp_sum;
                double vulnerable_total_in_cell = 0.0; // an approximation. assumes catches proportion for entire population, it's also stochastic with selectivity
                double tagged_fish = 0;
                double tagged_fish_vulnerable = 0;
                double average_counter = 0.0;
                double prob_of_untagged_fish = 0.0;
                double prob_of_tagged_fish = 0.0;
                double vulnerable_agent_in_cell = 0.0;
                double tagged_agent_vulnerable = 0.0;

                // vulnerable biomass an approximation
                // we want to find the proportion of vulnerable tagged fish (can argue they are all vulnerable) to vulnerable untagged fish
                // TODO: probably should be fishery specific, because of different selectivities, for now it is a mixed value
                if (account_for_tag_population_) {
                  tag_recapture_info.tagged_fish_available_.resize(model_->age_spread(), 0);
                  tag_recapture_info.all_fish_available_.resize(model_->age_spread(), 0);
                  LOG_MEDIUM() << "number of elements in nontagged partition = " << cell->agents_.size();
                  for (auto &agent : cell->agents_) {
                    if (agent.is_alive()) {
                      temp_sum = 0.0;
                      random_fish = rng.chance();
                      for (unsigned i = 0; i < fisheries_to_sample_from.size(); ++i) {
                        temp_sum += prop_catch_by_fishery[i];
                        if (temp_sum >= random_fish) {
                          fishery_ndx = fisheries_to_sample_from[i];
                          break;
                        }
                      }
                      if (rng.chance() <= fishery_selectivity_[fishery_ndx][agent.get_sex()]->GetResult(agent.get_age_index())) {
                        vulnerable_total_in_cell += agent.get_scalar();                // * agent.get_weight();
                        tag_recapture_info.all_fish_available_[agent.get_age_index()] += agent.get_scalar();                // * agent.get_weight();
                        average_counter++;
                        vulnerable_agent_in_cell++;
                      }
                    }
                  }
                  // Tagged partition
                  LOG_MEDIUM() << "number of elements in tagged partition = " << cell->tagged_agents_.size();
                  for (auto &agent : cell->tagged_agents_) {
                    if (agent.is_alive()) {
                      tagged_fish += agent.get_scalar(); // if all tagged fish are vulnerable to fishing
                      temp_sum = 0.0;
                      random_fish = rng.chance();
                      for (unsigned i = 0; i < fisheries_to_sample_from.size(); ++i) {
                        temp_sum += prop_catch_by_fishery[i];
                        if (temp_sum >= random_fish) {
                          fishery_ndx = fisheries_to_sample_from[i];
                          break;
                        }
                      }
                      // Assuming law of large numbers allows this to go for
                      if (rng.chance() <= fishery_selectivity_[fishery_ndx][agent.get_sex()]->GetResult(agent.get_age_index())) {
                        tagged_fish_vulnerable += agent.get_scalar(); // * agent.get_weight();
                        tag_recapture_info.tagged_fish_available_[agent.get_age_index()] += agent.get_scalar(); // * agent.get_weight();
                        tag_recapture_info.all_fish_available_[agent.get_age_index()] += agent.get_scalar(); // * agent.get_weight();
                        tagged_agent_vulnerable++;
                      }
                    }
                  }

                  tag_recapture_info.expected_scanned_ = vulnerable_total_in_cell + tagged_fish_vulnerable;
                  tag_recapture_info.scanned_agents_ = vulnerable_agent_in_cell + tagged_agent_vulnerable;

                  prob_of_untagged_fish = vulnerable_total_in_cell / (vulnerable_total_in_cell + tagged_fish_vulnerable);
                  prob_of_tagged_fish = 1 - prob_of_untagged_fish;

                  LOG_FINE() << "prob untagged = " << prob_of_untagged_fish << " prob_of_tagged_fish = " << prob_of_tagged_fish << " avg scalar = " << (vulnerable_total_in_cell / average_counter);
                  LOG_FINE() << "individuals_in_cell = " << vulnerable_total_in_cell << " tagged fish = " << tagged_fish << " prob of untagged fish = " << prob_of_untagged_fish
                      << " tagged fish that are vulnerable = " << tagged_fish_vulnerable << " prob untagged = " << prob_of_untagged_fish << " prob tagged = " << prob_of_tagged_fish;
                  // Proportion of individuals tagged
                  tag_recapture_info.proportion_inidividuals_tagged_ = prob_of_tagged_fish;
                  // Account for sampling of agents NOT individuals
                  // see line 551: ValidationCode/Agent based tag simulation.R for confirmation
                  prob_of_untagged_fish /= (vulnerable_total_in_cell / average_counter);
                  prob_of_untagged_fish = prob_of_untagged_fish / (prob_of_untagged_fish + prob_of_tagged_fish);
                  prob_of_tagged_fish = 1 - prob_of_untagged_fish;
                  tag_recapture_info.prob_sample_tagged_agents_ = prob_of_tagged_fish;
                }

                unsigned agent_has_tag = 0;
                unsigned agent_counter = 0;

                bool already_checked_selec = false;
                unsigned sample_agent_ndx = 0;
                // initialise agent that we will be assigning values to
                Agent *this_agent = nullptr;

                /*
                 *  Main loop
                 *  - Randomly select Fishery, weighted by catch (multinomial)
                 *  - Randomly select an agent
                 *  - if scanning for tags, do a binomial try
                 *  -- check if we are only fishing mature fish
                 *  -- check tagging info
                 *  -- chance() kill entier agent, and record information
                 */
                LOG_MEDIUM() << "begin main loop, catch to take = " << catch_taken << " max attempts allowed = " << catch_max;
                while (catch_taken > 0) {
                  /*
                   if (catch_attempts % 1000 == 0) {
                   cerr << catch_attempts << " ";
                   }
                   */
                  agent_counter = 0;
                  //LOG_MEDIUM() << "catch taken = " << catch_taken << " attempts = " << catch_attempts;
                  ++catch_attempts;
                  // Make sure we don't end up fishing for infinity if there isn't enough to catch
                  if (catch_attempts >= catch_max) {
                    LOG_WARNING() << "Too many attempts to catch an agent in the process " << label_ << " in year " << model_->current_year() << " in row " << row + 1 << " and column " << col + 1
                        << ", remaining catch to take = " << catch_taken << " this most likely means you have"
                        << " a model that suggests there should be more agents in this space than than the current agent dynamics are putting in this cell, check the user manual for tips to resolve this situation. attempts = "
                        << catch_attempts << " max attempts allowed = " << catch_max;
                    // Kick out of this cell
                    break;
                  }
                  already_checked_selec = false;   // this is only used in tagging.
                  // randomly select a fishery based on proportion of catch to take
                  temp_sum = 0.0;
                  random_fish = rng.chance();
                  for (unsigned i = 0; i < fisheries_to_sample_from.size(); ++i) {
                    temp_sum += prop_catch_by_fishery[i];
                    if (temp_sum >= random_fish) {
                      fishery_ndx = fisheries_to_sample_from[i];
                      break;
                    }
                  }
                  if (catch_to_take_by_fishery_[fishery_ndx] > 0) {
                    // randomly find agent
                    sample_agent_ndx = rng.chance() * cell->agents_.size();
                    this_agent = &cell->agents_[sample_agent_ndx];

                    // if scanning for tag-recoveries check that
                    if (scanning_this_year_[fishery_ndx]) {
                      agent_has_tag = (unsigned) rng.binomial(prob_of_tagged_fish, 1);
                      LOG_FINEST() << "find value " << agent_counter << " " << cell->agents_.size() << " fishery ndx = " << fishery_ndx << " agents has tag = " << agent_has_tag;
                      tag_recapture_info.tag_draws_ += agent_has_tag;
                      // find an untagged fish that is vulnerable and alive
                      if (agent_has_tag == 0) {
                        already_checked_selec = false;
                        while (not already_checked_selec) {
                          ++agent_counter;
                          if (agent_counter > (cell->agents_.size() * 5)) {
                            LOG_WARNING() << "couldn't find agent, that is vulnerable and alive in year " << model_->current_year() << " for process = " << label_ << " fishery = "
                                << fishery_label_[fishery_ndx];
                            break;
                          }
                          this_agent = &cell->agents_[rng.chance() * cell->agents_.size()];
                          if ((*this_agent).is_alive() & (rng.chance() <= (fishery_selectivity_[fishery_ndx][(*this_agent).get_sex()]->GetResult((*this_agent).get_age_index()))))
                            already_checked_selec = true;
                        }
                      } else if (agent_has_tag == 1) {
                        // rabdomly find an tagged agent that is vulnerable and alive
                        already_checked_selec = false;
                        while (not already_checked_selec) {
                          ++agent_counter;
                          if (agent_counter > (cell->tagged_agents_.size() * 10)) {
                            LOG_WARNING() << "couldn't find agent, that is vulnerable and alive in year " << model_->current_year() << " for process = " << label_ << " taged (1 = yes) "
                                << agent_has_tag << " fishery = " << fishery_label_[fishery_ndx] << " cell: " << row << "-" << col;
                            break;
                          }
                          this_agent = &cell->tagged_agents_[rng.chance() * cell->tagged_agents_.size()];
                          if ((*this_agent).is_alive() & (rng.chance() <= (fishery_selectivity_[fishery_ndx][(*this_agent).get_sex()]->GetResult((*this_agent).get_age_index()))))
                            already_checked_selec = true;
                        }
                      }
                    }
                    // Do the Moratlity check and remove from partition aka die
                    if ((*this_agent).is_alive()) {
                      if ((rng.chance() <= fishery_selectivity_[fishery_ndx][(*this_agent).get_sex()]->GetResult((*this_agent).get_age_index())) | already_checked_selec) {

                        if (only_mature_partition_[fishery_ndx]) { // for most models this will be false
                          if (not (*this_agent).get_maturity()) {
                            continue;
                          }
                        }
/*

                        if (scanning_this_year_[fishery_ndx]) {
                          LOG_FINE() << "agent age = " << (*this_agent).get_age() << " tag = " << agent_has_tag << " (" << (*this_agent).get_number_tags() << ") contributing "
                              << (*this_agent).get_scalar();
                        }
*/

                        // check under MLS and apply some handling mortality: NOTE no checking for tags here, TO ADD
                        if ((*this_agent).get_length() < fishery_mls_[fishery_ndx]) {
                          if (rng.chance() <= fishery_hand_mort_[fishery_ndx]) {
                            cell->remove_agent_alive((*this_agent).get_scalar());
                            (*this_agent).dies();
                          }
                        } else {
                          // record information
                          //LOG_FINE() << catch_taken << " " <<  (*this_agent).get_sex() << " " << (*this_agent).get_weight() << " " << (*this_agent).get_scalar();
                          catch_taken -= (*this_agent).get_weight() * (*this_agent).get_scalar();
                          catch_to_take_by_fishery_[fishery_ndx] -= (*this_agent).get_weight() * (*this_agent).get_scalar();
                          fishery_actual_catch_taken_[fishery_ndx][catch_ndx] += (*this_agent).get_weight() * (*this_agent).get_scalar();
                          age_freq.frequency_[(*this_agent).get_age_index()] += (*this_agent).get_scalar(); // This actually represents many individuals.
                          length_freq.frequency_[(*this_agent).get_length_bin_index()] += (*this_agent).get_scalar();

                          if ((*this_agent).get_sex() == 0) {
                            age_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].frequency_[(*this_agent).get_age_index()] += (*this_agent).get_scalar();
                            length_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].frequency_[(*this_agent).get_length_bin_index()] += (*this_agent).get_scalar();
                            fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].biomass_ += (*this_agent).get_weight() * (*this_agent).get_scalar();
                          } else {
                            age_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_frequency_[(*this_agent).get_age_index()] += (*this_agent).get_scalar();
                            length_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_frequency_[(*this_agent).get_length_bin_index()] += (*this_agent).get_scalar();
                            fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_biomass_ += (*this_agent).get_weight() * (*this_agent).get_scalar();
                          }
                          global_age_freq[(*this_agent).get_age_index()] += (*this_agent).get_scalar();
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].age_ndx_.push_back((*this_agent).get_age_index());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].sex_.push_back((*this_agent).get_sex());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].fishery_ndx_.push_back(fishery_ndx);
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].length_ndx_.push_back((*this_agent).get_length_bin_index());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].scalar_.push_back((*this_agent).get_scalar());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].weight_.push_back((*this_agent).get_weight());
                          actual_catch_by_area_[catch_ndx][fishery_ndx][row][col] += (*this_agent).get_scalar() * (*this_agent).get_weight();
                          tag_recapture_info.agents_caught_++;
                          tag_recapture_info.individuals_caught_ += (*this_agent).get_scalar();

                          if (scanning_this_year_[fishery_ndx]) {
                            // Probability of scanning agent
                            if (rng.chance() <= scanning_proportion_by_fishery_[fishery_ndx][catch_ndx]) {
                              // We scanned this agent
                              tag_recapture_info.scanned_age_comp_[(*this_agent).get_age_index()] += (*this_agent).get_scalar();
                              tag_recapture_info.scanned_length_comp_[(*this_agent).get_length_bin_index()] += (*this_agent).get_scalar();
                              tag_recapture_info.scanned_fish_ += (*this_agent).get_scalar() * (*this_agent).get_weight();
                              tag_recapture_info.agents_sampled_++;
                              if ((agent_has_tag == 1) & ((*this_agent).get_number_tags() <= 0)) {
                                LOG_MEDIUM() << "found a tagged agent, that has no tags??, tag release " << (*this_agent).get_tag_release_year() << " release region = " <<  (*this_agent).get_tag_row() + 1 << "-" <<  (*this_agent).get_tag_col() << " is it alive ? " << (*this_agent).is_alive();

                              }
                              if (agent_has_tag == 1) {
                                //fish has a tag record it
                                tag_recapture_info.age_.push_back((*this_agent).get_age());
                                tag_recapture_info.sex_.push_back((*this_agent).get_sex());
                                tag_recapture_info.length_ndx_.push_back((*this_agent).get_length_bin_index());
                                tag_recapture_info.length_.push_back((*this_agent).get_length());
                                tag_recapture_info.fishery_ndx_.push_back(fishery_ndx);
                                tag_recapture_info.time_at_liberty_.push_back((*this_agent).get_time_at_liberty(time_step));
                                tag_recapture_info.length_increment_.push_back((*this_agent).get_length_increment_since_tag());
                                tag_recapture_info.tag_row_.push_back((*this_agent).get_tag_row());
                                tag_recapture_info.tag_col_.push_back((*this_agent).get_tag_col());
                                tag_recapture_info.tag_release_year_.push_back((*this_agent).get_tag_release_year());
                              }
                            }
                          }
                        }
                        (*this_agent).dies();
                      }
                    }
                  }
                } // while (catch_taken > 0

                for (auto &catch_ : catch_to_take_by_fishery_) {
                  LOG_FINE() << catch_;
                }
                removals_by_length_and_area_.push_back(length_freq);
                removals_by_age_and_area_.push_back(age_freq);
                removals_census_.push_back(census_fishery);
                if (scanning_ & (tag_recapture_info.age_.size() > 0)) {
                  LOG_FINE() << "saving tag-recapture";
                  removals_tag_recapture_.push_back(tag_recapture_info);
                }
                LOG_MEDIUM() << "row = " << row + 1 << " col = " << col + 1 << " tags drawn = " << tag_recapture_info.tag_draws_ << " tags saved " << tag_recapture_info.length_ndx_.size();
                LOG_MEDIUM() << "catch attempts for agents counter =  " << catch_attempts;

              } //if (catch_taken > 0) {
            } // cell ->is_enabled()
          } // col
        } // row
      } else {
        LOG_FINE() << "Applying length based F";
        double world_catch_to_take = 0;
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            unsigned catch_attempts = 1;
            unsigned catch_max = 1;
            WorldCell *cell = nullptr;
            double catch_taken = 0;
            cell = world_->get_base_square(row, col); // Shared resource...
            // Which fisheries are we taken catches for, not all fisheries are taking catch every year.
            vector<unsigned> fisheries_to_sample_from;
            if (cell->is_enabled()) {
              double catch_for_fishery = 0.0;
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

              vector<double> prop_catch_by_fishery(fisheries_to_sample_from.size(), 0.0);
              for (unsigned i = 0; i < prop_catch_by_fishery.size(); ++i) {
                prop_catch_by_fishery[i] = catch_to_take_by_fishery_[fisheries_to_sample_from[i]] / catch_taken;
                LOG_FINE() << "fishery ndx = " << fisheries_to_sample_from[i] << " proportion of catch = " << prop_catch_by_fishery[i];
              }

              world_catch_to_take += catch_taken;
              LOG_FINE() << "need to remove " << catch_taken << " weight of fish across all fisheries";
              if (catch_taken > 0) {
                LOG_FINE() << "We are fishing in cell " << row + 1 << " " << col + 1 << " value = " << catch_taken;
                census_data census_fishery(model_->current_year(), row, col);
                tag_recapture tag_recapture_info(model_->current_year(), row, col, time_step);
                tag_recapture_info.scanned_age_comp_.resize(model_->age_spread(), 0.0);
                tag_recapture_info.scanned_length_comp_.resize(model_->number_of_length_bins(), 0.0);
                composition_data age_freq(PARAM_AGE, model_->current_year(), row, col, model_->age_spread());
                composition_data length_freq(PARAM_LENGTH, model_->current_year(), row, col, model_->number_of_length_bins());

                catch_attempts = 1;
                catch_max = cell->agents_.size() * 50;
                LOG_FINEST() << "individuals = " << catch_max;
                unsigned fishery_ndx = 0;

                double random_fish, temp_sum;
                double vulnerable_total_in_cell = 0.0; // an approximation. assumes catches proportion for entire population, it's also stochastic with selectivity
                double tagged_fish = 0;
                double tagged_fish_vulnerable = 0;
                double vulnerable_agent_in_cell = 0;
                double tagged_agent_vulnerable = 0;
                double average_counter = 0.0;
                double prob_of_untagged_fish = 0.0;
                double prob_of_tagged_fish = 0.0;

                // vulnerable biomass an approximation
                // we want to find the proportion of vulnerable tagged fish (can argue they are all vulnerable) to vulnerable untagged fish
                // TODO: probably should be fishery specific, because of different selectivities, for now it is a mixed value
                if (account_for_tag_population_) {
                  tag_recapture_info.tagged_fish_available_.resize(model_->age_spread(), 0);
                  tag_recapture_info.all_fish_available_.resize(model_->age_spread(), 0);
                  LOG_MEDIUM() << "number of elements in nontagged partition = " << cell->agents_.size();
                  for (auto &agent : cell->agents_) {
                    if (agent.is_alive()) {
                      temp_sum = 0.0;
                      random_fish = rng.chance();
                      for (unsigned i = 0; i < fisheries_to_sample_from.size(); ++i) {
                        temp_sum += prop_catch_by_fishery[i];
                        if (temp_sum >= random_fish) {
                          fishery_ndx = fisheries_to_sample_from[i];
                          break;
                        }
                      }
                      if (rng.chance() <= fishery_selectivity_[fishery_ndx][agent.get_sex()]->GetResult(agent.get_length_bin_index())) {
                        vulnerable_total_in_cell += agent.get_scalar();                // * agent.get_weight();
                        tag_recapture_info.all_fish_available_[agent.get_age_index()] += agent.get_scalar();                // * agent.get_weight();
                        average_counter++;
                        vulnerable_agent_in_cell++;
                      }
                    }
                  }
                  // Tagged partition
                  LOG_MEDIUM() << "number of elements in tagged partition = " << cell->tagged_agents_.size();
                  for (auto &agent : cell->tagged_agents_) {
                    if (agent.is_alive()) {
                      tagged_fish += agent.get_scalar(); // if all tagged fish are vulnerable to fishing
                      temp_sum = 0.0;
                      random_fish = rng.chance();
                      for (unsigned i = 0; i < fisheries_to_sample_from.size(); ++i) {
                        temp_sum += prop_catch_by_fishery[i];
                        if (temp_sum >= random_fish) {
                          fishery_ndx = fisheries_to_sample_from[i];
                          break;
                        }
                      }
                      // Assuming law of large numbers allows this to go for
                      if (rng.chance() <= fishery_selectivity_[fishery_ndx][agent.get_sex()]->GetResult(agent.get_length_bin_index())) {
                        tagged_fish_vulnerable += agent.get_scalar(); // * agent.get_weight();
                        tag_recapture_info.tagged_fish_available_[agent.get_age_index()] += agent.get_scalar(); // * agent.get_weight();
                        tag_recapture_info.all_fish_available_[agent.get_age_index()] += agent.get_scalar(); // * agent.get_weight();
                        tagged_agent_vulnerable++;
                      }
                    }
                  }

                  tag_recapture_info.expected_scanned_ = vulnerable_total_in_cell + tagged_fish_vulnerable;
                  tag_recapture_info.scanned_agents_ = vulnerable_agent_in_cell + tagged_agent_vulnerable;

                  prob_of_untagged_fish = vulnerable_total_in_cell / (vulnerable_total_in_cell + tagged_fish_vulnerable);
                  prob_of_tagged_fish = 1 - prob_of_untagged_fish;

                  LOG_FINE() << "prob untagged = " << prob_of_untagged_fish << " prob_of_tagged_fish = " << prob_of_tagged_fish << " avg scalar = " << (vulnerable_total_in_cell / average_counter);
                  LOG_FINE() << "individuals_in_cell = " << vulnerable_total_in_cell << " tagged fish = " << tagged_fish << " prob of untagged fish = " << prob_of_untagged_fish
                      << " tagged fish that are vulnerable = " << tagged_fish_vulnerable << " prob untagged = " << prob_of_untagged_fish << " prob tagged = " << prob_of_tagged_fish;
                  // Proportion of individuals tagged
                  tag_recapture_info.proportion_inidividuals_tagged_ = prob_of_tagged_fish;
                  // Account for sampling of entire agents NOT individuals
                  // see line 551: ValidationCode/Agent based tag simulatio.R for confirmation
                  prob_of_untagged_fish /= (vulnerable_total_in_cell / average_counter);
                  prob_of_untagged_fish = prob_of_untagged_fish / (prob_of_untagged_fish + prob_of_tagged_fish);
                  prob_of_tagged_fish = 1 - prob_of_untagged_fish;
                  tag_recapture_info.prob_sample_tagged_agents_ = prob_of_tagged_fish;
                }

                unsigned agent_has_tag = 0;
                unsigned agent_counter = 0;

                bool already_checked_selec = false;

                // initialise agent that we will be assigning values to
                Agent *this_agent = nullptr;

                /*
                 *  Main loop
                 */
                while (catch_taken > 0) {
                  agent_counter = 0;
                  LOG_FINEST() << "catch taken = " << catch_taken << " attempts = " << catch_attempts;
                  ++catch_attempts;
                  already_checked_selec = false;
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
                  if (catch_to_take_by_fishery_[fishery_ndx] > 0) {
                    if (scanning_this_year_[fishery_ndx])
                      LOG_FINE() << "about to calculate agent";
                    // randomly find agent
                    this_agent = &cell->agents_[rng.chance() * cell->agents_.size()];

                    if (scanning_this_year_[fishery_ndx]) {
                      agent_has_tag = (unsigned) rng.binomial(prob_of_tagged_fish, 1);
                      LOG_FINEST() << "find value " << agent_counter << " " << cell->agents_.size() << " fishery ndx = " << fishery_ndx << " agents has tag = " << agent_has_tag;
                      tag_recapture_info.tag_draws_ += agent_has_tag;
                      // find an untagged fish that is vulnerable and alive
                      if (agent_has_tag == 0) {
                        already_checked_selec = false;
                        while (not already_checked_selec) {
                          ++agent_counter;
                          if (agent_counter > (cell->agents_.size() * 5)) {
                            LOG_WARNING() << "couldn't find agent, that is vulnerable and alive in year " << model_->current_year() << " for process = " << label_ << " fishery = "
                                << fishery_label_[fishery_ndx];
                            break;
                          }
                          this_agent = &cell->agents_[rng.chance() * cell->agents_.size()];
                          if ((*this_agent).is_alive() & (rng.chance() <= (fishery_selectivity_[fishery_ndx][(*this_agent).get_sex()]->GetResult((*this_agent).get_length_bin_index()))))
                            already_checked_selec = true;
                        }
                      } else if (agent_has_tag == 1) {
                        // find an tagged fish that is vulnerable and alive
                        already_checked_selec = false;
                        while (not already_checked_selec) {
                          ++agent_counter;
                          if (agent_counter > (cell->tagged_agents_.size() * 10)) {
                            LOG_WARNING() << "couldn't find agent, that is vulnerable and alive in year " << model_->current_year() << " for process = " << label_ << " taged (1 = yes) "
                                << agent_has_tag << " fishery = " << fishery_label_[fishery_ndx] << " cell: " << row << "-" << col;
                            break;
                          }
                          this_agent = &cell->tagged_agents_[rng.chance() * cell->tagged_agents_.size()];
                          if ((*this_agent).is_alive() & (rng.chance() <= (fishery_selectivity_[fishery_ndx][(*this_agent).get_sex()]->GetResult((*this_agent).get_length_bin_index()))))
                            already_checked_selec = true;
                        }
                      }
                    }
                    if ((*this_agent).is_alive()) {
                      if ((rng.chance() <= fishery_selectivity_[fishery_ndx][(*this_agent).get_sex()]->GetResult((*this_agent).get_length_bin_index())) | already_checked_selec) {
/*

                        if (scanning_this_year_[fishery_ndx]) {
                          LOG_FINE() << "agent age = " << (*this_agent).get_age() << " tag = " << agent_has_tag << " (" << (*this_agent).get_number_tags() << ") contributing "
                              << (*this_agent).get_scalar();
                        }
*/

                        // check under MLS and apply some handling mortality
                        if ((*this_agent).get_length() < fishery_mls_[fishery_ndx]) {
                          if (rng.chance() <= fishery_hand_mort_[fishery_ndx]) {
                            cell->remove_agent_alive((*this_agent).get_scalar());
                            (*this_agent).dies();
                          }
                        } else {
                          // record information
                          //LOG_FINE() << catch_taken << " " <<  (*this_agent).get_sex() << " " << (*this_agent).get_weight() << " " << (*this_agent).get_scalar();
                          catch_taken -= (*this_agent).get_weight() * (*this_agent).get_scalar();
                          catch_to_take_by_fishery_[fishery_ndx] -= (*this_agent).get_weight() * (*this_agent).get_scalar();
                          fishery_actual_catch_taken_[fishery_ndx][catch_ndx] += (*this_agent).get_weight() * (*this_agent).get_scalar();
                          age_freq.frequency_[(*this_agent).get_age_index()] += (*this_agent).get_scalar(); // This actually represents many individuals.
                          length_freq.frequency_[(*this_agent).get_length_bin_index()] += (*this_agent).get_scalar();

                          if ((*this_agent).get_sex() == 0) {
                            age_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].frequency_[(*this_agent).get_age_index()] += (*this_agent).get_scalar();
                            length_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].frequency_[(*this_agent).get_length_bin_index()] += (*this_agent).get_scalar();
                            fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].biomass_ += (*this_agent).get_weight() * (*this_agent).get_scalar();
                          } else {
                            age_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_frequency_[(*this_agent).get_age_index()] += (*this_agent).get_scalar();
                            length_comp_by_fishery_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_frequency_[(*this_agent).get_length_bin_index()] += (*this_agent).get_scalar();
                            fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].female_biomass_ += (*this_agent).get_weight() * (*this_agent).get_scalar();
                          }
                          global_age_freq[(*this_agent).get_age_index()] += (*this_agent).get_scalar();
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].age_ndx_.push_back((*this_agent).get_age_index());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].sex_.push_back((*this_agent).get_sex());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].fishery_ndx_.push_back(fishery_ndx);
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].length_ndx_.push_back((*this_agent).get_length_bin_index());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].scalar_.push_back((*this_agent).get_scalar());
                          fishery_census_data_[fishery_ndx][catch_ndx][cell_ndx_[fishery_ndx]].weight_.push_back((*this_agent).get_weight());
                          actual_catch_by_area_[catch_ndx][fishery_ndx][row][col] += (*this_agent).get_scalar() * (*this_agent).get_weight();
                          tag_recapture_info.agents_caught_++;
                          tag_recapture_info.individuals_caught_ += (*this_agent).get_scalar();
                          if (scanning_ & scanning_this_year_[fishery_ndx]) {
                            // Probability of scanning agent
                            if (rng.chance() <= scanning_proportion_by_fishery_[fishery_ndx][catch_ndx]) {
                              // We scanned this agent
                              tag_recapture_info.scanned_age_comp_[(*this_agent).get_age_index()] += (*this_agent).get_scalar();
                              tag_recapture_info.scanned_length_comp_[(*this_agent).get_length_bin_index()] += (*this_agent).get_scalar();
                              tag_recapture_info.scanned_fish_ += (*this_agent).get_scalar() * (*this_agent).get_weight();
                              tag_recapture_info.agents_sampled_++;
                              if ((*this_agent).get_number_tags() > 0) {
                                //fish has a tag record it
                                tag_recapture_info.age_.push_back((*this_agent).get_age());
                                tag_recapture_info.sex_.push_back((*this_agent).get_sex());
                                tag_recapture_info.length_ndx_.push_back((*this_agent).get_length_bin_index());
                                tag_recapture_info.length_.push_back((*this_agent).get_length());
                                tag_recapture_info.fishery_ndx_.push_back(fishery_ndx);
                                tag_recapture_info.time_at_liberty_.push_back((*this_agent).get_time_at_liberty(time_step));
                                tag_recapture_info.length_increment_.push_back((*this_agent).get_length_increment_since_tag());
                                tag_recapture_info.tag_row_.push_back((*this_agent).get_tag_row());
                                tag_recapture_info.tag_col_.push_back((*this_agent).get_tag_col());
                                tag_recapture_info.tag_release_year_.push_back((*this_agent).get_tag_release_year());
                              }
                            }
                          }
                        }
                        (*this_agent).dies();
                      }
                    }
                    // Make sure we don't end up fishing for infinity
                    if (catch_attempts >= catch_max) {
                      LOG_WARNING() << "Too many attempts to catch an agent in the process " << label_ << " in year " << model_->current_year() << " in row " << row + 1 << " and column " << col + 1
                          << ", remaining catch to take = " << catch_taken << " this most likely means you have"
                          << " a model that suggests there should be more agents in this space than than the current agent dynamics are putting in this cell, check the user manual for tips to resolve this situation. attempts = "
                          << catch_attempts << " max attempts allowed = " << catch_max;
                      // Kick out of this cell
                      break;
                    }
                  }
                } // while (catch_taken > 0

                for (auto &catch_ : catch_to_take_by_fishery_) {
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
      } // length based
    } // find(years_.begin(), years_.end(), model_->current_year()) != years_.end()
  }  //model_->state() != State::kInitialise
  LOG_FINE() << "finished Biomass Mort process. NUmber of tag-recapture obs " << removals_tag_recapture_.size();
}

// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void MortalityEventBiomass::FillReportCache(ostringstream &cache) {
  LOG_FINE() << "";
  // Fishery specific info
  // age frequency by sex fishery and year
  if (age_comp_by_fishery_.size() > 0) {
    vector<unsigned> temp_age_freq(model_->age_spread(), 0.0);
    if (model_->get_sexed()) {
      for (auto &fishery : fishery_index_) {
        cache << "age_freq-male-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year cell ";
        for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age)
          cache << age << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < age_comp_by_fishery_[fishery].size(); ++year_ndx) {
          fill(temp_age_freq.begin(), temp_age_freq.end(), 0.0);
          // Print by cell
          LOG_FINE() << "cells = " << age_comp_by_fishery_[fishery][year_ndx].size();
          for (unsigned cell_ndx = 0; cell_ndx < age_comp_by_fishery_[fishery][year_ndx].size(); ++cell_ndx) {
            cache << years_[year_ndx] << " " << age_comp_by_fishery_[fishery][year_ndx][cell_ndx].row_ + 1 << "-" << age_comp_by_fishery_[fishery][year_ndx][cell_ndx].col_ + 1 << " ";
            if (age_comp_by_fishery_[fishery][year_ndx][cell_ndx].year_ == years_[year_ndx]) {
              for (unsigned age_ndx = 0; age_ndx < age_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_.size(); ++age_ndx) {
                cache << age_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_[age_ndx] << " ";
                temp_age_freq[age_ndx] += age_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_[age_ndx];
              }
              cache << "\n";
            }
          }
          cache << years_[year_ndx] << " Total ";
          for (auto &age_freq : temp_age_freq)
            cache << age_freq << " ";
          cache << "\n";
        }
      }

      for (auto &fishery : fishery_index_) {
        cache << "age_freq-female-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year cell ";
        for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age)
          cache << age << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < age_comp_by_fishery_[fishery].size(); ++year_ndx) {
          fill(temp_age_freq.begin(), temp_age_freq.end(), 0.0);
          // By cell
          for (unsigned cell_ndx = 0; cell_ndx < age_comp_by_fishery_[fishery][year_ndx].size(); ++cell_ndx) {
            if (age_comp_by_fishery_[fishery][year_ndx][cell_ndx].year_ == years_[year_ndx]) {
              cache << years_[year_ndx] << " " << age_comp_by_fishery_[fishery][year_ndx][cell_ndx].row_ + 1 << "-" << age_comp_by_fishery_[fishery][year_ndx][cell_ndx].col_ + 1 << " ";

              for (unsigned age_ndx = 0; age_ndx < age_comp_by_fishery_[fishery][year_ndx][cell_ndx].female_frequency_.size(); ++age_ndx) {
                cache << age_comp_by_fishery_[fishery][year_ndx][cell_ndx].female_frequency_[age_ndx] << " ";
                temp_age_freq[age_ndx] += age_comp_by_fishery_[fishery][year_ndx][cell_ndx].female_frequency_[age_ndx];

              }
              cache << "\n";
            }
          }
          cache << years_[year_ndx] << " Total ";
          for (auto &age_freq : temp_age_freq)
            cache << age_freq << " ";
          cache << "\n";
        }
      }
    } else {
      for (auto &fishery : fishery_index_) {
        cache << "age_freq-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year cell ";
        for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age)
          cache << age << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < age_comp_by_fishery_[fishery].size(); ++year_ndx) {
          fill(temp_age_freq.begin(), temp_age_freq.end(), 0.0);
          // By cell
          for (unsigned cell_ndx = 0; cell_ndx < age_comp_by_fishery_[fishery][year_ndx].size(); ++cell_ndx) {
            if (age_comp_by_fishery_[fishery][year_ndx][cell_ndx].year_ == years_[year_ndx]) {
              cache << years_[year_ndx] << " " << age_comp_by_fishery_[fishery][year_ndx][cell_ndx].row_ + 1 << "-" << age_comp_by_fishery_[fishery][year_ndx][cell_ndx].col_ + 1 << " ";
              for (unsigned age_ndx = 0; age_ndx < age_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_.size(); ++age_ndx) {
                cache << age_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_[age_ndx] << " ";
                temp_age_freq[age_ndx] += age_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_[age_ndx];
              }
              cache << "\n";
            }
          }
          cache << years_[year_ndx] << " Total ";
          for (auto &age_freq : temp_age_freq)
            cache << age_freq << " ";
          cache << "\n";
        }
      }
    }
  }
  // Length comp by fishery
  if (length_comp_by_fishery_.size() > 0) {
    vector<unsigned> temp_len_freq(model_->length_bin_mid_points().size(), 0.0);
    if (model_->get_sexed()) {
      for (auto &fishery : fishery_index_) {
        cache << "length_freq-male-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year cell ";
        for (auto len_bin : model_->length_bin_mid_points())
          cache << len_bin << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < length_comp_by_fishery_[fishery].size(); ++year_ndx) {
          fill(temp_len_freq.begin(), temp_len_freq.end(), 0.0);
          // By cell
          for (unsigned cell_ndx = 0; cell_ndx < length_comp_by_fishery_[fishery][year_ndx].size(); ++cell_ndx) {
            if (length_comp_by_fishery_[fishery][year_ndx][cell_ndx].year_ == years_[year_ndx]) {
              cache << years_[year_ndx] << " " << length_comp_by_fishery_[fishery][year_ndx][cell_ndx].row_ + 1 << "-" << length_comp_by_fishery_[fishery][year_ndx][cell_ndx].col_ + 1 << " ";
              for (unsigned len_ndx = 0; len_ndx < length_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_.size(); ++len_ndx) {
                cache << length_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_[len_ndx] << " ";
                temp_len_freq[len_ndx] += length_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_[len_ndx];
              }
              cache << "\n";
            }
          }
          cache << years_[year_ndx] << " Total ";
          for (auto &len_freq : temp_len_freq)
            cache << len_freq << " ";
          cache << "\n";
        }
      }

      for (auto &fishery : fishery_index_) {
        cache << "length_freq-female-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year cell ";
        for (auto len_bin : model_->length_bin_mid_points())
          cache << len_bin << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < length_comp_by_fishery_[fishery].size(); ++year_ndx) {
          fill(temp_len_freq.begin(), temp_len_freq.end(), 0.0);
          // By cell
          for (unsigned cell_ndx = 0; cell_ndx < length_comp_by_fishery_[fishery][year_ndx].size(); ++cell_ndx) {
            if (length_comp_by_fishery_[fishery][year_ndx][cell_ndx].year_ == years_[year_ndx]) {
              cache << years_[year_ndx] << " " << length_comp_by_fishery_[fishery][year_ndx][cell_ndx].row_ + 1 << "-" << length_comp_by_fishery_[fishery][year_ndx][cell_ndx].col_ + 1 << " ";
              for (unsigned len_ndx = 0; len_ndx < length_comp_by_fishery_[fishery][year_ndx][cell_ndx].female_frequency_.size(); ++len_ndx) {
                cache << length_comp_by_fishery_[fishery][year_ndx][cell_ndx].female_frequency_[len_ndx] << " ";
                temp_len_freq[len_ndx] += length_comp_by_fishery_[fishery][year_ndx][cell_ndx].female_frequency_[len_ndx];
              }
              cache << "\n";
            }
          }
          cache << years_[year_ndx] << " Total ";
          for (auto &len_freq : temp_len_freq)
            cache << len_freq << " ";
          cache << "\n";
        }
      }
    } else {
      for (auto &fishery : fishery_index_) {
        cache << "length_freq-" << fishery_label_[fishery] << " " << REPORT_R_DATAFRAME << "\n";
        cache << "year cell ";
        for (auto len_bin : model_->length_bin_mid_points())
          cache << len_bin << " ";
        cache << "\n";
        for (unsigned year_ndx = 0; year_ndx < length_comp_by_fishery_[fishery].size(); ++year_ndx) {
          fill(temp_len_freq.begin(), temp_len_freq.end(), 0.0);
          // By cell
          for (unsigned cell_ndx = 0; cell_ndx < length_comp_by_fishery_[fishery][year_ndx].size(); ++cell_ndx) {
            if (length_comp_by_fishery_[fishery][year_ndx][cell_ndx].year_ == years_[year_ndx]) {
              cache << years_[year_ndx] << " " << length_comp_by_fishery_[fishery][year_ndx][cell_ndx].row_ + 1 << "-" << length_comp_by_fishery_[fishery][year_ndx][cell_ndx].col_ + 1 << " ";
              for (unsigned len_ndx = 0; len_ndx < length_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_.size(); ++len_ndx) {
                cache << length_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_[len_ndx] << " ";
                temp_len_freq[len_ndx] += length_comp_by_fishery_[fishery][year_ndx][cell_ndx].frequency_[len_ndx];
              }
              cache << "\n";
            }
          }
          cache << years_[year_ndx] << " Total ";
          for (auto &len_freq : temp_len_freq)
            cache << len_freq << " ";
          cache << "\n";
        }
      }
    }
  }
  // Print actual catches by year x fishery x area
  for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    cache << "actual_catches-" << years_[year_ndx] << " " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
    cache << "cell ";
    for (unsigned fishery_ndx = 0; fishery_ndx < fishery_label_.size(); ++fishery_ndx)
      cache << fishery_label_[fishery_ndx] << " ";
    cache << "\n";
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        cache << row + 1 << "-" << col + 1 << " ";
        for (unsigned fishery_ndx = 0; fishery_ndx < fishery_label_.size(); ++fishery_ndx) {
          cache << actual_catch_by_area_[year_ndx][fishery_ndx][row][col] << " ";
        }
        cache << "\n";
      }
    }
  }
  // Print Recapture information
  if (removals_tag_recapture_.size() > 0) {
    for (auto &tag_recapture : removals_tag_recapture_) {
      if (tag_recapture.age_.size() > 0) {
        cache << "tag_recapture_info-" << tag_recapture.year_ << "-" << tag_recapture.row_ + 1 << "-" << tag_recapture.col_ + 1 << " " << REPORT_R_LIST << "\n";
        cache << "scanned_weight: " << tag_recapture.scanned_fish_ << "\n"; // over all fisheries
        // before fishing
        cache << "approximate_individuals_vulnerable: " << tag_recapture.expected_scanned_ << "\n";  // both tagged and untagged
        cache << "approximate_agents_vulnerable: " << tag_recapture.scanned_agents_ << "\n";// both tagged and untagged
        cache << "approximate_proportion_tagged_individuals_vulnerable: " << tag_recapture.proportion_inidividuals_tagged_ << "\n";
        cache << "prob_sample_tagged_agent: " << tag_recapture.prob_sample_tagged_agents_ << "\n"; // accounted for agents vs indivudals of tagged
        cache << "agents_scanned: " << tag_recapture.agents_sampled_ << "\n";
        cache << "binomial_samples: " << tag_recapture.tag_draws_ << "\n";
        cache << "agents_caught: " << tag_recapture.agents_caught_ << "\n";
        cache << "individuals_caught: " << tag_recapture.individuals_caught_ << "\n";
        cache << "tagged_fish_available: ";
        for (auto &tag_age : tag_recapture.tagged_fish_available_)
          cache << tag_age << " ";
        cache << "\n";
        cache << "all_fish_available: ";
        for (auto &tag_age : tag_recapture.all_fish_available_)
          cache << tag_age << " ";
        cache << "\n";
        if (print_tag_recap_info_) {
          cache << "values " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
          LOG_MEDIUM() << "ages = " << tag_recapture.age_.size();
          cache << "rowlabs ";
          for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
            cache << ndx + 1 << " ";
          cache << "\n";
          cache << "release_year ";
          for (unsigned ndx = 0; ndx < tag_recapture.tag_release_year_.size(); ++ndx)
            cache << tag_recapture.tag_release_year_[ndx] << " ";
          cache << "\n";
          cache << "age ";
          //cache << "age length length-increment time_at_liberty\n";
          for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
            cache << tag_recapture.age_[ndx] << " ";
          cache << "\n";
          cache << "length_bin ";
          for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
            cache << tag_recapture.length_ndx_[ndx] << " ";
          cache << "\n";
          cache << "length ";
          for (unsigned ndx = 0; ndx < tag_recapture.length_.size(); ++ndx)
            cache << tag_recapture.length_[ndx] << " ";
          cache << "\n";
          cache << "at_liberty ";
          for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
            cache << tag_recapture.time_at_liberty_[ndx] << " ";
          cache << "\n";
          cache << "length_increment ";
          for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
            cache << tag_recapture.length_increment_[ndx] << " ";
          cache << "\n";
          cache << "tag_row ";
          for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
            cache << tag_recapture.tag_row_[ndx] + 1 << " ";
          cache << "\n";
          cache << "fishery ";
          for (unsigned ndx = 0; ndx < tag_recapture.fishery_ndx_.size(); ++ndx)
            cache << fishery_label_[tag_recapture.fishery_ndx_[ndx]] << " ";
          cache << "\n";
          cache << "tag_col ";
          for (unsigned ndx = 0; ndx < tag_recapture.age_.size(); ++ndx)
            cache << tag_recapture.tag_col_[ndx] + 1 << " ";
          cache << "\n" << REPORT_R_LIST_END << "\n";
        } else {
          cache << REPORT_R_LIST_END << "\n";
        }
      }
    }
  }
  if (print_census_info_) {
    if (removals_by_age_.size() > 0) {
      cache << "age_frequency " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
      cache << "year ";
      for (unsigned i = model_->min_age(); i <= model_->max_age(); ++i)
        cache << i << " ";
      cache << "\n";
      for (auto &age_freq : removals_by_age_) {
        cache << age_freq.first << " ";
        for (auto age_value : age_freq.second)
          cache << age_value << " ";
        cache << "\n";
      }
    }
     
  

    // Print Census information (warning this will print alot of information

    if (fishery_census_data_.size() > 0) {
      LOG_FINE() << fishery_census_data_.size() << " census vals";
      // Print census information
      for (unsigned fish_ndx = 0; fish_ndx < fishery_census_data_.size(); ++fish_ndx) {
        for (unsigned year_ndx = 0; year_ndx < fishery_census_data_[fish_ndx].size(); ++year_ndx) {
          for (unsigned cell_ndx = 0; cell_ndx < fishery_census_data_[fish_ndx][year_ndx].size(); ++cell_ndx) {
            if (fishery_census_data_[fish_ndx][year_ndx][cell_ndx].age_ndx_.size() > 0) {
              cache << "census_info-" << fishery_label_[fishery_index_[fish_ndx]] << "-" << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].year_ << "-"
                  << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].row_ + 1 << "-" << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].col_ + 1 << " " << REPORT_R_MATRIX << "\n";
              for (unsigned age_ndx = 0; age_ndx < fishery_census_data_[fish_ndx][year_ndx][cell_ndx].age_ndx_.size(); ++age_ndx)
                cache << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].age_ndx_[age_ndx] + model_->min_age() << " ";
              cache << "\n";
              for (unsigned length_ndx = 0; length_ndx < fishery_census_data_[fish_ndx][year_ndx][cell_ndx].length_ndx_.size(); ++length_ndx)
                cache << model_->length_bin_mid_points()[fishery_census_data_[fish_ndx][year_ndx][cell_ndx].length_ndx_[length_ndx]] << " ";
              cache << "\n";
              for (unsigned sex_ndx = 0; sex_ndx < fishery_census_data_[fish_ndx][year_ndx][cell_ndx].sex_.size(); ++sex_ndx)
                cache << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].sex_[sex_ndx] << " ";
              cache << "\n";
			  for (unsigned scalar_ndx = 0; scalar_ndx < fishery_census_data_[fish_ndx][year_ndx][cell_ndx].scalar_.size(); ++scalar_ndx)
                cache << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].scalar_[scalar_ndx] << " ";
              cache << "\n";
			  for (unsigned weight_ndx = 0; weight_ndx < fishery_census_data_[fish_ndx][year_ndx][cell_ndx].scalar_.size(); ++weight_ndx)
                cache << fishery_census_data_[fish_ndx][year_ndx][cell_ndx].weight_[weight_ndx] << " ";
              cache << "\n";
            }
          }
        }
      }
    }
  }
}

/*
 * Set Catch based on R code
 */
void MortalityEventBiomass::set_HCR(map<unsigned, map<string, double>> future_catches) {
  // find
  for (auto year_map : future_catches) {
    if(find(years_.begin(), years_.end(), year_map.first) == years_.end())
      LOG_FATAL() << "could not find year " << year_map.first << " for setting HCR rule";
    harvest_control_years_.push_back(year_map.first);
    for (auto fish_map : year_map.second) {
      if (find(fishery_label_.begin(), fishery_label_.end(), fish_map.first) == fishery_label_.end())
        LOG_FATAL() << "could not find fishery " << fish_map.first << " for setting HCR rule";
      // set save the values
      harvest_control_Fs_[year_map.first][fish_map.first] = fish_map.second;
    }
  }

}

} /* namespace processes */
} /* namespace niwa */
