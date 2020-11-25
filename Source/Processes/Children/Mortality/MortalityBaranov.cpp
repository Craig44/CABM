/**
 * @file MortalityBaranov.cpp
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
#include "MortalityBaranov.h"

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
MortalityBaranov::MortalityBaranov(Model *model) :
    Mortality(model) {
  f_table_ = new parameters::Table(PARAM_F_TABLE);
  parameters_.BindTable(PARAM_F_TABLE, f_table_, "Table of f_layer by year for each fishery", "", true, false);
  method_table_ = new parameters::Table(PARAM_METHOD_INFO);
  parameters_.BindTable(PARAM_METHOD_INFO, method_table_, "Table of information for each fishery.", "", true, false);
  parameters_.Bind<bool>(PARAM_PRINT_CENSUS_INFO, &print_census_info_, "if you have process report for this process you can control the amount of information printed to the file.", "", true);

}

MortalityBaranov::~MortalityBaranov() {
  delete f_table_;
  delete method_table_;
}

/**
 * Do some initial checks of user supplied parameters.
 */
void MortalityBaranov::DoValidate() {
  // check headers
  auto f_columns = f_table_->columns();
  LOG_FINE() << "check f table headers, number of columns = " << f_columns.size();

  if (f_columns[0] != PARAM_YEAR)
    LOG_ERROR_P(PARAM_F_TABLE) << "the first column header needs to be " << PARAM_YEAR << ". Please add it =)";
  for (unsigned j = 1; j < f_columns.size(); ++j) {
    fishery_index_.push_back(j - 1);
    fishery_label_.push_back(f_columns[j]);
    LOG_FINE() << "fishery ndx = " << j - 1 << " fishery label " << f_columns[j];
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
    LOG_FINE() << "trying to check this " << meth_data[i][0] << " fishery is in the f table " << i;
    if (std::find(fishery_label_.begin(), fishery_label_.end(), (string) meth_data[i][0]) == fishery_label_.end()) {
      LOG_ERROR_P(PARAM_METHOD_INFO) << "Could not find the row label " << meth_data[i][0]
          << " these need to be consistent with the column headers with the f table, can you please double check this";
    } else {
      LOG_FINE() << "found " << meth_data[i][0] << " in the f table";
    }
  }
  LOG_FINE() << "finish Validation";
}

/**
 * DoBuild
 */
void MortalityBaranov::DoBuild() {
  LOG_FINE() << "Build tables";
  LOG_FINE() << "Building f table";
  // allocate memory for fishery dimension to related objects
  fishery_f_layer_labels_.resize(fishery_label_.size());
  fishery_actual_catch_taken_.resize(fishery_label_.size());
  fishery_f_to_take_.resize(fishery_label_.size());
  fishery_f_layer_.resize(fishery_label_.size());
  fishery_selectivity_label_.resize(fishery_label_.size());
  fishery_selectivity_.resize(fishery_label_.size());
  fishery_mls_.resize(fishery_label_.size());
  fishery_hand_mort_.resize(fishery_label_.size());
  age_comp_by_fishery_.resize(fishery_label_.size());
  length_comp_by_fishery_.resize(fishery_label_.size());
  fishery_census_data_.resize(fishery_label_.size());
  cell_ndx_.resize(fishery_label_.size());
  auto model_years = model_->years();
  // load objects
  vector<vector<string>> &f_values_data = f_table_->data();

  for (auto row : f_values_data) {
    unsigned year = 0;
    if (!utilities::To<string, unsigned>(row[0], year))
      LOG_ERROR_P(PARAM_F_TABLE) << "year value " << row[0] << " is not numeric.";
    if (std::find(model_years.begin(), model_years.end(), year) == model_years.end())
      LOG_ERROR_P(PARAM_F_TABLE) << "year " << year << " is not a valid year in this model";
    // Check years are consecutive ascending order.
    // This will mean when I reference year_ndx later in the code we can have faith it is in order.
    if (years_.size() > 1) {
      if ((year - 1) != years_[years_.size() - 1]) {
        LOG_ERROR_P(PARAM_F_TABLE) << "years need to be in consecutive ascending order, the year " << years_[years_.size() - 1] << " was followed by " << year << " please sort this out";
      }
    }
    years_.push_back(year);
    for (unsigned i = 1; i < row.size(); ++i) {
      fishery_f_layer_labels_[i - 1].push_back(row[i]);
      layers::NumericLayer *temp_layer = nullptr;
      temp_layer = model_->managers().layer()->GetNumericLayer(row[i]);
      if (!temp_layer) {
        LOG_FATAL_P(PARAM_F_TABLE)
        << "could not find the f layer '" << row[i] << "', in year " << year << " for fishery " << fishery_label_[i - 1] << " please check it is numeric if it exists";
      }
      fishery_f_layer_[i - 1].push_back(temp_layer);
    }
  }
  // Allocate memory for some of the f related maps
  for (unsigned i = 0; i < fishery_f_layer_.size(); ++i) {
    fishery_actual_catch_taken_[i].resize(fishery_f_layer_[i].size(), 0.0);
    fishery_f_to_take_[i].resize(fishery_f_layer_[i].size(), 0.0);
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
    float mls = 0;
    float hand_mort = 0;
    if (!utilities::To<string, float>(meth_row[mls_index], mls))
      LOG_ERROR_P(PARAM_METHOD_INFO) << PARAM_MINIMUM_LEGAL_LENGTH << " value " << meth_row[mls_index] << " is not numeric, please sort this out.";
    if (!utilities::To<string, float>(meth_row[hand_mort_index], hand_mort))
      LOG_ERROR_P(PARAM_METHOD_INFO) << PARAM_HANDLING_MORTALITY << " value " << meth_row[hand_mort_index] << " is not numeric, please sort this out.";
    fishery_mls_[fishery_index_[meth_index]] = mls;
    fishery_hand_mort_[fishery_index_[meth_index]] = hand_mort;
  }


  for (unsigned i = 0; i < age_comp_by_fishery_.size(); ++i) {
    age_comp_by_fishery_[i].resize(years_.size());
    length_comp_by_fishery_[i].resize(years_.size());
    fishery_census_data_[i].resize(years_.size());
  }

  LOG_FINE() << "finished building Tables";


  prop_F_fishery_and_bin_.resize(fishery_label_.size());
  f_to_take_by_fishery_.resize(fishery_label_.size());
  for (unsigned fishery_ndx = 0; fishery_ndx < fishery_label_.size(); ++fishery_ndx) {
    prop_F_fishery_and_bin_[fishery_ndx].resize(model_->get_height());
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      prop_F_fishery_and_bin_[fishery_ndx][row].resize(model_->get_width());
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        if (model_->get_sexed()) {
          prop_F_fishery_and_bin_[fishery_ndx][row][col].resize(2);
          if (selectivity_length_based_) {
            prop_F_fishery_and_bin_[fishery_ndx][row][col][0].resize(model_->number_of_length_bins(),0.0);
            prop_F_fishery_and_bin_[fishery_ndx][row][col][1].resize(model_->number_of_length_bins(),0.0);
           } else {
             prop_F_fishery_and_bin_[fishery_ndx][row][col][0].resize(model_->age_spread(),0.0);
             prop_F_fishery_and_bin_[fishery_ndx][row][col][1].resize(model_->age_spread(),0.0);
           }
        } else {
          prop_F_fishery_and_bin_[fishery_ndx][row][col].resize(1);
          if (selectivity_length_based_) {
            prop_F_fishery_and_bin_[fishery_ndx][row][col][0].resize(model_->number_of_length_bins(),0.0);
           } else {
             prop_F_fishery_and_bin_[fishery_ndx][row][col][0].resize(model_->age_spread(),0.0);
           }
        }
      }
    }
  }
  F_by_year_bin_.resize(years_.size());
  actual_catch_by_area_.resize(years_.size());

  for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
     F_by_year_bin_[year_ndx].resize(model_->get_height());
     for (unsigned row = 0; row < model_->get_height(); ++row) {
       F_by_year_bin_[year_ndx][row].resize(model_->get_width());
       for (unsigned col = 0; col < model_->get_width(); ++col) {
         if (model_->get_sexed()) {
           F_by_year_bin_[year_ndx][row][col].resize(2);
           if (selectivity_length_based_) {
             F_by_year_bin_[year_ndx][row][col][0].resize(model_->number_of_length_bins(),0.0);
             F_by_year_bin_[year_ndx][row][col][1].resize(model_->number_of_length_bins(),0.0);
            } else {
              F_by_year_bin_[year_ndx][row][col][0].resize(model_->age_spread(),0.0);
              F_by_year_bin_[year_ndx][row][col][1].resize(model_->age_spread(),0.0);
            }
         } else {
           F_by_year_bin_[year_ndx][row][col].resize(1);
           if (selectivity_length_based_) {
             F_by_year_bin_[year_ndx][row][col][0].resize(model_->number_of_length_bins(),0.0);
            } else {
              F_by_year_bin_[year_ndx][row][col][0].resize(model_->age_spread(),0.0);
            }
         }
       }
     }

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
void MortalityBaranov::DoReset() {
  LOG_FINE() << "clearing containers";
  removals_by_age_and_area_.clear();
  removals_by_length_and_area_.clear();
  removals_census_.clear();
  removals_tag_recapture_.clear();
  for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for(unsigned sex_ndx = 0; sex_ndx < F_by_year_bin_[year_ndx][col][row].size(); ++sex_ndx)
          fill(F_by_year_bin_[year_ndx][col][row][sex_ndx].begin(), F_by_year_bin_[year_ndx][col][row][sex_ndx].end(), 0.0);
      }
    }
    for (unsigned fishery_ndx = 0; fishery_ndx < fishery_label_.size(); ++fishery_ndx) {
      fishery_f_to_take_[fishery_ndx][year_ndx] = 0.0;
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
void MortalityBaranov::DoExecute() {
  utilities::RandomNumberGenerator &rng = utilities::RandomNumberGenerator::Instance(); // shared resource
  vector<unsigned> global_age_freq(model_->age_spread(), 0);
  auto iter = years_.begin();
  if (model_->state() != State::kInitialise) {
    if (std::find(iter, years_.end(), model_->current_year()) != years_.end()) {
      iter = find(years_.begin(), years_.end(), model_->current_year());
      unsigned year_ndx = distance(years_.begin(), iter);

      if (not selectivity_length_based_) {
        LOG_FINE() << "Age based F";
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            WorldCell *cell = nullptr;
            float f_taken = 0;
            cell = world_->get_base_square(row, col); // Shared resource...
            // Which fisheries are we taken fes for, not all fisheries are taking f every year.
            vector<unsigned> fisheries_to_sample_from;
            if (cell->is_enabled()) {
              unsigned age_iter = 0;
              fill(f_to_take_by_fishery_.begin(), f_to_take_by_fishery_.end(), 0.0);
              float f_for_fishery = 0.0;
              LOG_FINE() << "cell: " << row << "-" << col << " what we are storing = " << f_to_take_by_fishery_.size() << " looping over " << fishery_f_layer_.size() << " fishery label = " <<  fishery_label_.size();
              LOG_FINE() << " prop_F_fishery_and_bin_.size() " << prop_F_fishery_and_bin_.size();
              for (unsigned i = 0; i < fishery_label_.size(); ++i) {
                f_for_fishery = fishery_f_layer_[i][year_ndx]->get_value(row, col);
                fishery_f_to_take_[i][year_ndx] += f_for_fishery;
                f_taken += f_for_fishery;
                LOG_FINE() << "in year " << model_->current_year() << " fishery " << fishery_label_[i] << "need to remove " << f_for_fishery << " cell_ndx_.size() " << cell_ndx_.size();
                if (model_->get_sexed()) {
                  age_iter = 0;
                  for (auto age = model_->min_age(); age <= model_->max_age(); ++age, ++age_iter) {
                    F_by_year_bin_[year_ndx][row][col][0][age_iter] += f_for_fishery * fishery_selectivity_[i][0]->GetResult(age_iter);
                    F_by_year_bin_[year_ndx][row][col][1][age_iter] += f_for_fishery * fishery_selectivity_[i][1]->GetResult(age_iter);
                    prop_F_fishery_and_bin_[i][row][col][0][age_iter] = f_for_fishery * fishery_selectivity_[i][0]->GetResult(age_iter);
                    prop_F_fishery_and_bin_[i][row][col][1][age_iter] = f_for_fishery * fishery_selectivity_[i][1]->GetResult(age_iter);
                  }
                } else {
                  age_iter = 0;
                  for (auto age = model_->min_age(); age <= model_->max_age(); ++age, ++age_iter) {
                    LOG_FINE() << "age = " << age << " age iter = " << age_iter;
                    F_by_year_bin_[year_ndx][row][col][0][age_iter] += f_for_fishery * fishery_selectivity_[i][0]->GetResult(age_iter);
                    prop_F_fishery_and_bin_[i][row][col][0][age_iter] = f_for_fishery * fishery_selectivity_[i][0]->GetResult(age_iter);
                  }
                }
                LOG_FINE() << "appropriated F";

                composition_data age_freq_for_fishery(PARAM_AGE, model_->current_year(), row, col, model_->age_spread());
                age_comp_by_fishery_[i][year_ndx].push_back(age_freq_for_fishery);
                composition_data length_freq_for_fishery(PARAM_LENGTH, model_->current_year(), row, col, model_->number_of_length_bins());
                length_comp_by_fishery_[i][year_ndx].push_back(length_freq_for_fishery);
                census_data census_fishery_specific(model_->current_year(), row, col);
                fishery_census_data_[i][year_ndx].push_back(census_fishery_specific);
                f_to_take_by_fishery_[i] = f_for_fishery;
                LOG_MEDIUM() << "ndx = " << cell_ndx_[i];
                cell_ndx_[i] = age_comp_by_fishery_[i][year_ndx].size() - 1;
                LOG_MEDIUM() << "ndx = " << cell_ndx_[i] << " for fishery = " << fishery_label_[i];
              }
              // Convert prop_F_fishery_and_bin_ to a proportion of F among fisheries.
              for (unsigned i = 0; i < fishery_label_.size(); ++i) {
                LOG_FINE() << "fishery " << fishery_label_[i];
                age_iter = 0;
                if (model_->get_sexed()) {
                  for (auto age = model_->min_age(); age <= model_->max_age(); ++age, ++age_iter) {
                    prop_F_fishery_and_bin_[i][row][col][0][age_iter] /= F_by_year_bin_[year_ndx][row][col][0][age_iter];
                    prop_F_fishery_and_bin_[i][row][col][1][age_iter] /= F_by_year_bin_[year_ndx][row][col][1][age_iter];
                  }
                } else {
                  for (auto age = model_->min_age(); age <= model_->max_age(); ++age, ++age_iter) {
                    prop_F_fishery_and_bin_[i][row][col][0][age_iter] /= F_by_year_bin_[year_ndx][row][col][0][age_iter];
                    LOG_FINE() << "age = " << age << " prop F " << prop_F_fishery_and_bin_[i][row][col][0][age_iter];
                  }
                }
              }
              LOG_FINE() << "need to remove " << f_taken << " weight of fish across all fisheries";
              //unsigned agent_counter = 0;
              if (f_taken > 0) {
                LOG_MEDIUM() << "We are fishing in cell " << row + 1 << " " << col + 1 << " value = " << f_taken;
                census_data census_fishery(model_->current_year(), row, col);
                composition_data age_freq(PARAM_AGE, model_->current_year(), row, col, model_->age_spread());
                composition_data length_freq(PARAM_LENGTH, model_->current_year(), row, col, model_->number_of_length_bins());

                /*
                 *  Main loop
                 *  - loop over all agents
                 *  - apply total F over all agents chance() < exp(-F_by_year_bin_[year_ndx][row][col][0][age_iter])
                 *  - If caught then do a multinomial draw, to see which fishery caught this agent.
                 */
                float temp_sum = 0.0;
                float random_chance = 0.0;
                for(auto& agent : cell->agents_) {
                  if (agent.is_alive()) {
                    //LOG_FINEST() << "selectivity = " << selectivity_at_age << " m = " << (*iter).get_m();
                    if (rng.chance() <= (1 - std::exp(- F_by_year_bin_[year_ndx][row][col][agent.get_sex()][agent.get_age_index()]))) {
                      // caught this fish
                      // Which fishery
                      temp_sum = 0;
                      random_chance = rng.chance();
                      for (unsigned fish_ndx = 0; fish_ndx < fishery_label_.size(); ++fish_ndx) {
                        temp_sum += prop_F_fishery_and_bin_[fish_ndx][row][col][agent.get_sex()][agent.get_age_index()];
                        if (temp_sum > random_chance) {
                          // caught by this fishery
                          // check under MLS and apply some handling mortality: NOTE no checking for tags here, TO ADD
                          if (agent.get_length() < fishery_mls_[fish_ndx]) {
                            // released
                            if (rng.chance() <= fishery_hand_mort_[fish_ndx])
                              agent.dies();
                          } else {
                            // save obs
                            fishery_actual_catch_taken_[fish_ndx][year_ndx] += agent.get_weight() * agent.get_scalar();
                            age_freq.frequency_[agent.get_age_index()] += agent.get_scalar(); // This actually represents many individuals.
                            length_freq.frequency_[agent.get_length_bin_index()] += agent.get_scalar();

                            if (agent.get_sex() == 0) {
                              age_comp_by_fishery_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].frequency_[agent.get_age_index()] += agent.get_scalar();
                              length_comp_by_fishery_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].frequency_[agent.get_length_bin_index()] += agent.get_scalar();
                              fishery_census_data_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].biomass_ +=agent.get_weight() * agent.get_scalar();
                            } else {
                              age_comp_by_fishery_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].female_frequency_[agent.get_age_index()] += agent.get_scalar();
                              length_comp_by_fishery_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].female_frequency_[agent.get_length_bin_index()] += agent.get_scalar();
                              fishery_census_data_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].female_biomass_ += agent.get_weight() * agent.get_scalar();
                            }
                            global_age_freq[agent.get_age_index()] += agent.get_scalar();
                            fishery_census_data_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].age_ndx_.push_back(agent.get_age_index());
                            fishery_census_data_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].sex_.push_back(agent.get_sex());
                            fishery_census_data_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].fishery_ndx_.push_back(fish_ndx);
                            fishery_census_data_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].length_ndx_.push_back(agent.get_length_bin_index());
                            fishery_census_data_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].scalar_.push_back(agent.get_scalar());
                            fishery_census_data_[fish_ndx][year_ndx][cell_ndx_[fish_ndx]].weight_.push_back(agent.get_weight());
                            actual_catch_by_area_[year_ndx][fish_ndx][row][col] += agent.get_scalar() * agent.get_weight();
                            agent.dies();
                          }
                          break;
                        }
                      }
                    }
                  }
                }
                removals_by_length_and_area_.push_back(length_freq);
                removals_by_age_and_area_.push_back(age_freq);
                removals_census_.push_back(census_fishery);
              } //if (f_taken > 0) {
            } // cell ->is_enabled()
          } // col
        } // row
      } else {
        LOG_FINE() << "Applying length based F";
      }
    } // find(years_.begin(), years_.end(), model_->current_year()) != years_.end()
  }  //model_->state() != State::kInitialise
  LOG_FINE() << "finished Biomass Mort process.";
}

// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void MortalityBaranov::FillReportCache(ostringstream &cache) {
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
  // Print actual fes by year x fishery x area
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

  // Print actual fes by year x fishery x area
  if( not model_->get_sexed()) {
    for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
      cache << "F_by_bin-" << years_[year_ndx] << " " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
      cache << "cell ";
      if(selectivity_length_based_) {
        for (auto len_bin : model_->length_bin_mid_points())
          cache << len_bin << " ";
      } else {
        for (auto age_bin = model_->min_age(); age_bin <= model_->max_age(); ++age_bin)
          cache << age_bin << " ";
      }
      cache << "\n";
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          cache << row + 1 << "-" << col + 1 << " ";
          for (unsigned bin_ndx = 0; bin_ndx < F_by_year_bin_[year_ndx][row][col][0].size(); ++bin_ndx) {
            cache <<  F_by_year_bin_[year_ndx][row][col][0][bin_ndx] << " ";
          }
          cache << "\n";
        }
      }
    }
  } else {
    LOG_WARNING() << "need to complete fishing report for sexed models.";
  }

  // Print actual fes by year x fishery x area
  WorldCell *cell = nullptr;
  for (unsigned year_ndx = 0; year_ndx < years_.size(); ++year_ndx) {
    cache << "F-" << years_[year_ndx] << " " << REPORT_R_DATAFRAME_ROW_LABELS << "\n";
    cache << "cell ";
    for (unsigned fishery_ndx = 0; fishery_ndx < fishery_label_.size(); ++fishery_ndx)
      cache << fishery_label_[fishery_ndx] << " ";
    cache << "\n";
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        cache << row + 1 << "-" << col + 1 << " ";
        cell = world_->get_base_square(row, col); // Shared resource...
        if (cell->is_enabled()) {
          for (unsigned i = 0; i < fishery_label_.size(); ++i) {
            cache << fishery_f_layer_[i][year_ndx]->get_value(row, col) << " ";
          }
          cache << "\n";
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
            }
          }
        }
      }
    }
  }
}

} /* namespace processes */
} /* namespace niwa */
