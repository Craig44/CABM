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
  parameters_.Bind<string>(PARAM_SELECTIVITY, &selectivity_label_, "Selectivity label", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "years to apply the process", "");
  parameters_.Bind<string>(PARAM_CATCH_LAYERS, &catch_layer_label_, "Spatial layer describing catch by cell for each year, there is a one to one link with the year specified, so make sure the order is right", "");
  parameters_.Bind<float>(PARAM_MINIMUM_LEGAL_LENGTH, &mls_, "The minimum legal length for this fishery, any individual less than this will be returned using some discard mortality", "", 0);
  parameters_.Bind<float>(PARAM_HANDLING_MORTALITY, &discard_mortality_, "if discarded due to being under the minimum legal length, what is the probability the individual will die when released", "", 1.0);
  parameters_.Bind<bool>(PARAM_PRINT_EXTRA_INFO, &print_extra_info_, "if you have process report for this process you can control the amount of information printed to the file.", "", true);
  // Tag-recapture inputs
  parameters_.Bind<unsigned>(PARAM_SCANNING_YEARS, &scanning_years_, "Years to scan for tags in the fishery", "", true);
  parameters_.Bind<float>(PARAM_SCANNING_PROPROTION, &scanning_proportion_, "The proportion of catch scanned in each year", "", true);


}

/**
 * Do some initial checks of user supplied parameters.
 */
void MortalityEventBiomass::DoValidate() {
  LOG_TRACE();
  if (years_.size() != catch_layer_label_.size())
    LOG_ERROR_P(PARAM_YEARS) << "you must specify a layer label for each year. You have supplied '" << years_.size() << "' years but '" << catch_layer_label_.size() << "' catch layer labels, please sort this out.";

  if (scanning_proportion_.size() == 1)
    scanning_proportion_.resize(scanning_years_.size(), scanning_proportion_[0]);
  if (scanning_years_.size() != scanning_proportion_.size()) {
    LOG_ERROR_P(PARAM_SCANNING_YEARS) << "there needs to be a proportion for all years to apply tagging in. You supplied " << scanning_years_.size() << " years but " << scanning_proportion_.size() << " proportions, sort this descrepency out please";
  }

  for (auto& prop : scanning_proportion_) {
    if (prop < 0.0)
      LOG_ERROR_P(PARAM_SCANNING_PROPROTION) << "you cannot supply a proportion less than 0";
  }

}

/**
 * DoBuild
 */
void MortalityEventBiomass::DoBuild() {
  LOG_FINE();
  // Get the layers
  for (auto& label : catch_layer_label_) {
    layers::NumericLayer* temp_layer = nullptr;
    temp_layer = model_->managers().layer()->GetNumericLayer(label);
    if (!temp_layer) {
      LOG_FATAL_P(PARAM_CATCH_LAYERS) << "could not find the layer '" << label << "', please make sure it exists and that it is type 'numeric'";
    }
    catch_layer_.push_back(temp_layer);
  }

  // Users don't have to apply tagging
  if (parameters_.Get(PARAM_SCANNING_YEARS)->has_been_defined()) {
    for (auto year : scanning_years_) {
      if (find(years_.begin(), years_.end(), year) == years_.end())
        LOG_ERROR_P(PARAM_SCANNING_YEARS) << "could not find the year " << year << " in the process years, please sort this out";
    }
  }

  LOG_FINEST() << "selectivities supplied = " << selectivity_label_.size();
  // Build selectivity links
  if (selectivity_label_.size() == 1)
    selectivity_label_.assign(2, selectivity_label_[0]);

  if (selectivity_label_.size() > 2) {
    LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << "You suppled " << selectivity_label_.size()  << " Selectiviites, you can only have one for each sex max = 2";
  }
  LOG_FINEST() << "selectivities supplied = " << selectivity_label_.size();

  bool first = true;
  for (auto label : selectivity_label_) {
    Selectivity* temp_selectivity = model_->managers().selectivity()->GetSelectivity(label);
    if (!temp_selectivity)
      LOG_ERROR_P(PARAM_SELECTIVITY_LABEL) << ": selectivity " << label << " does not exist. Have you defined it?";

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

  cell_offset_.resize(model_->get_height());
  model_length_bins_.resize(model_->get_height());
  model_age_bins_.resize(model_->get_height());
  mls_by_space_.resize(model_->get_height());
  current_year_by_space_.resize(model_->get_height());
  scanning_prop_year_by_space_.resize(model_->get_height());
  discard_by_space_.resize(model_->get_height());
  current_time_step_by_space_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cell_offset_[i].resize(model_->get_width());
    model_length_bins_[i].resize(model_->get_width(), model_->length_bins().size());
    model_age_bins_[i].resize(model_->get_width(), model_->age_spread());
    mls_by_space_[i].resize(model_->get_width(), mls_);
    current_year_by_space_[i].resize(model_->get_width());
    discard_by_space_[i].resize(model_->get_width(), discard_mortality_);
    scanning_prop_year_by_space_[i].resize(model_->get_width(), 0.0);
    current_time_step_by_space_[i].resize(model_->get_width(), 1);
  }
}

/**
 * DoExecute
 */
void MortalityEventBiomass::DoExecute() {
  LOG_MEDIUM();
  vector<unsigned> global_age_freq(model_->age_spread(), 0);
  auto iter = years_.begin();
  if (model_->state() != State::kInitialise) {
    if (find(iter, years_.end(), model_->current_year()) != years_.end()) {
      iter = find(years_.begin(), years_.end(), model_->current_year());
      unsigned catch_ndx = distance(years_.begin(), iter);

      // Are we applying Tagging either releases or scanning, or both
      auto scan_iter = scanning_years_.begin();
      bool scanning = false; // A bool used to see if we are scanning
      unsigned scanning_ndx = 0;
      if (find(scan_iter , scanning_years_.end(), model_->current_year()) != scanning_years_.end()) {
        scanning = true;
        scanning_ndx = distance(scanning_years_.begin(),scan_iter);
      }
      LOG_FINE() << "scanning index = " << scanning_ndx;
      LOG_FINE() << "applying F in year " << model_->current_year() << " catch index = " << catch_ndx;
      // Pre-calculate agents in the world to set aside our random numbers needed for the operation
      n_agents_ = 0;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            LOG_FINEST() << "Setting spatial areas, row = " << row + 1 << " col = " << col + 1;
            cell_offset_[row][col] = n_agents_;
            LOG_FINEST() << "cell_offset done";
            n_agents_ += cell->agents_.size();
            current_year_by_space_[row][col] = model_->current_year();
            LOG_FINEST() << "year by size done";
            if (scanning) {
              scanning_prop_year_by_space_[row][col] = scanning_proportion_[scanning_ndx];
              LOG_FINEST() << "scanning by size done";
            }
            current_time_step_by_space_[row][col] = model_->get_time_step_counter();
            LOG_FINEST() << "time-step by size done";

          }
        }
      }




      // Allocate a single block of memory rather than each thread temporarily allocating their own memory.
      random_numbers_.resize(n_agents_ + 1);
      discard_random_numbers_.resize(n_agents_ + 1);
      selectivity_random_numbers_.resize(n_agents_ + 1);
      scanning_random_numbers_.resize(n_agents_ + 1);
      utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
      for (unsigned i = 0; i <= n_agents_; ++i) {
        random_numbers_[i] = rng.chance();
        discard_random_numbers_[i] = rng.chance();
        selectivity_random_numbers_[i] = rng.chance();
        scanning_random_numbers_[i] = rng.chance();
      }

      LOG_FINE() << "about to kick into the gear";
      float actual_catch_taken = 0;
      float world_catch_to_take = 0;
      if (not selectivity_length_based_) {
        // Thread out each loop
        // #pragma omp parallel for collapse(2)
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            unsigned catch_attempts = 1;
            unsigned catch_max = 1;
            WorldCell* cell = nullptr;
            float catch_taken = 0;
            float actual_catch_this_cell = 0;
            cell = world_->get_base_square(row, col); // Shared resource...
            catch_taken = catch_layer_[catch_ndx]->get_value(row, col);
            world_catch_to_take += catch_taken;
            if (cell->is_enabled()) {
              unsigned counter = 0;
              if (catch_taken > 0) {
                LOG_FINEST() << "We are fishing in cell " << row + 1 << " " << col + 1 << " value = " << catch_taken;
                composition_data age_freq(PARAM_AGE, current_year_by_space_[row][col], row, col, model_age_bins_[row][col]);
                composition_data length_freq(PARAM_LENGTH, current_year_by_space_[row][col], row, col, model_length_bins_[row][col]);
                census_data census_fishery(current_year_by_space_[row][col], row, col);
                tag_recapture tag_recapture_info(current_year_by_space_[row][col], row, col, current_time_step_by_space_[row][col]);

                catch_attempts = 1;
                catch_max = cell->agents_.size();
                LOG_FINEST() << "individuals = " << catch_max;
                while (catch_taken > 0) {
                  ++catch_attempts;
                  auto& this_agent = cell->agents_[random_numbers_[cell_offset_[row][col] + counter] * cell->agents_.size()];
                  if (this_agent.is_alive()) {
                    //LOG_FINEST() << "attempt "<< catch_attempts << " catch_taken = " << catch_taken << " age = " << this_agent.get_age_index() << " random number = " << selectivity_random_numbers_[cell_offset_[row][col] + counter]  << " selectivity = " << cell_offset_for_selectivity_[row][col][this_agent.get_sex() * model_->age_spread() + this_agent.get_age_index()];
                    // See if this agent is unlucky
                    if (selectivity_random_numbers_[cell_offset_[row][col] + counter] <= selectivity_[this_agent.get_sex()]->GetResult(this_agent.get_age_index())) {
                      //LOG_FINEST() << "vulnerable to gear catch_taken = " << catch_taken << " age = " << this_agent.get_age_index() << " individual weight = " << this_agent.get_weight() << " scalar = " <<  this_agent.get_scalar();
                      if (this_agent.get_length() < mls_by_space_[row][col]) {
                        if (discard_random_numbers_[cell_offset_[row][col] + counter] <= discard_by_space_[row][col]) {
                          this_agent.dies();
                        }
                      } else {
                        // record information
                        catch_taken -= this_agent.get_weight() * this_agent.get_scalar();
                        actual_catch_this_cell += this_agent.get_weight() * this_agent.get_scalar();
                        age_freq.frequency_[this_agent.get_age_index()] += this_agent.get_scalar(); // This catch actually represents many individuals.
                        length_freq.frequency_[this_agent.get_length_bin_index()] += this_agent.get_scalar();
                        census_fishery.age_.push_back(this_agent.get_age());
                        census_fishery.length_.push_back(this_agent.get_length());
                        census_fishery.scalar_.push_back(this_agent.get_scalar());
                        if (scanning) {
                          // Probability of scanning agent
                          if (scanning_random_numbers_[cell_offset_[row][col] + counter] <= scanning_prop_year_by_space_[row][col]) {
                            // We scanned this agent
                            tag_recapture_info.scanned_fish_++;
                            if (this_agent.get_number_tags() > 0) {
                              // fish has a tag record it
                              tag_recapture_info.age_.push_back(this_agent.get_age());
                              tag_recapture_info.length_.push_back(this_agent.get_length());
                              tag_recapture_info.time_at_liberty_.push_back(this_agent.get_time_at_liberty(current_time_step_by_space_[row][col]));
                              tag_recapture_info.length_increment_.push_back(this_agent.get_length_increment_since_tag());
                              tag_recapture_info.tag_row_.push_back(this_agent.get_tag_row());
                              tag_recapture_info.tag_col_.push_back(this_agent.get_tag_col());

                            }
                          }
                        }

                        this_agent.dies();
                      }
                    }
                    // Make sure we don't end up fishing for infinity
                    if (catch_attempts >= catch_max) {
                      LOG_FATAL_P(PARAM_TYPE) << "Too many attempts to catch an agent in the process " << label_ << " in year " << current_year_by_space_[row][col] << " in row " << row + 1 << " and column " << col + 1 << ", remaining catch to take = " << catch_taken << " this most likely means you have" <<
                         " a model that suggests there should be more agents in this space than than the current agent dynamics are putting in this cell, check the user manual for tips to resolve this situation";
                    }
                  }
                  ++counter;
                }
                #pragma omp critical
                {
                  for (unsigned i = 0; i < model_age_bins_[row][col]; ++i)
                    global_age_freq[i] += age_freq.frequency_[i];
                  removals_by_length_and_area_.push_back(length_freq);
                  removals_by_age_and_area_.push_back(age_freq);
                  actual_catch_taken += actual_catch_this_cell;
                  removals_census_.push_back(census_fishery);
                  if (scanning) {
                    removals_tag_recapture_.push_back(tag_recapture_info);
                  }
                }
              }
              LOG_FINEST() << "individuals = " << cell->agents_.size();
            }
          }
        }
      } else { // Selectivity = length based.
        iter = find(years_.begin(), years_.end(), model_->current_year());
        unsigned catch_ndx = distance(years_.begin(), iter);
        // Get the pointer to the right catch layer
        // Thread out each loop
        #pragma omp parallel for collapse(2)
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            unsigned catch_attempts = 1;
            unsigned catch_max = 1;
            float actual_catch_this_cell = 0;
            WorldCell* cell = nullptr;
            float catch_taken = 0;
            #pragma omp critical
            {
              cell = world_->get_base_square(row, col); // Shared resource...
              catch_taken = catch_layer_[catch_ndx]->get_value(row, col);
              world_catch_to_take += catch_taken;
            }
            if (cell->is_enabled()) {
              unsigned counter = 0;
              if (catch_taken > 0) {
                LOG_FINEST() << "We are fishing in cell " << row + 1 << " " << col + 1 << " value = " << catch_taken;
                // create reporting class
                composition_data age_freq(PARAM_AGE, current_year_by_space_[row][col], row, col, model_age_bins_[row][col]);
                composition_data length_freq(PARAM_LENGTH, current_year_by_space_[row][col], row, col, model_length_bins_[row][col]);
                census_data census_fishery(current_year_by_space_[row][col], row, col);
                tag_recapture tag_recapture_info(current_year_by_space_[row][col], row, col,current_time_step_by_space_[row][col]);

                catch_attempts = 1;
                catch_max = cell->agents_.size();
                while (catch_taken > 0) {
                  // Random access bullshit for lists
                  ++catch_attempts;
                  auto& this_agent = cell->agents_[random_numbers_[cell_offset_[row][col] + counter] * cell->agents_.size()];
                  if (this_agent.is_alive()) {
                      // See if this agent is unlucky
                    if (selectivity_random_numbers_[cell_offset_[row][col] + counter] <= selectivity_[this_agent.get_sex()]->GetResult(this_agent.get_length_bin_index())) {
                      if (this_agent.get_length() < mls_by_space_[row][col]) {
                        if (discard_random_numbers_[cell_offset_[row][col] + counter] <= discard_by_space_[row][col]) {
                          this_agent.dies();
                        }
                      } else {
                        catch_taken -= this_agent.get_weight() * this_agent.get_scalar();
                        actual_catch_this_cell += this_agent.get_weight() * this_agent.get_scalar();
                        age_freq.frequency_[this_agent.get_age_index()]+= this_agent.get_scalar(); // This catch actually represents many individuals.
                        length_freq.frequency_[this_agent.get_length_bin_index()]+= this_agent.get_scalar();
                        census_fishery.age_.push_back(this_agent.get_age());
                        census_fishery.length_.push_back(this_agent.get_length());
                        census_fishery.scalar_.push_back(this_agent.get_scalar());
                        if (scanning) {
                          // Probability of scanning agent
                          if (scanning_random_numbers_[cell_offset_[row][col] + counter] <= scanning_prop_year_by_space_[row][col]) {
                            // We scanned this agent
                            tag_recapture_info.scanned_fish_++;
                            if (this_agent.get_number_tags() > 0) {
                              // fish has a tag record it
                              tag_recapture_info.age_.push_back(this_agent.get_age());
                              tag_recapture_info.length_.push_back(this_agent.get_length());
                              tag_recapture_info.time_at_liberty_.push_back(this_agent.get_time_at_liberty(current_time_step_by_space_[row][col]));
                              tag_recapture_info.length_increment_.push_back(this_agent.get_length_increment_since_tag());
                              tag_recapture_info.tag_row_.push_back(this_agent.get_tag_row());
                              tag_recapture_info.tag_col_.push_back(this_agent.get_tag_col());
                            }
                          }
                        }
                        this_agent.dies();
                      }
                    }
                    // Make sure we don't end up fishing for infinity if there are not enough fish here
                    if (catch_attempts >= catch_max) {
                      LOG_FATAL_P(PARAM_TYPE) << "Too many attempts to catch an agent in the process " << label_ << " in year " << current_year_by_space_[row][col] << " in row " << row + 1 << " and column " << col + 1 << ", remaining catch to take = " << catch_taken << " this most likely means you have" <<
                         " a model that suggests there should be more agents in this space than than the current agent dynamics are putting in this cell, check the user manual for tips to resolve this situation";
                    }
                  }
                  ++counter;
                }
                #pragma omp critical
                {
                  for (unsigned i = 0; i < model_age_bins_[row][col]; ++i)
                    global_age_freq[i] += age_freq.frequency_[i];
                  removals_by_length_and_area_.push_back(length_freq);
                  removals_by_age_and_area_.push_back(age_freq);
                  actual_catch_taken += actual_catch_this_cell;
                  removals_census_.push_back(census_fishery);
                  if (scanning) {
                    removals_tag_recapture_.push_back(tag_recapture_info);
                  }
                }
              } // if catch > 0
            } // is enabled
          } //col
        } // row
      } // length_based_selectivity
      LOG_FINE() << "world catch taken = " << world_catch_to_take << " actual catch taken = " << actual_catch_taken;
      removals_by_year_[model_->current_year()] = world_catch_to_take;
      actual_removals_by_year_[model_->current_year()] = actual_catch_taken;
      removals_by_age_[model_->current_year()] = global_age_freq;
      LOG_MEDIUM();
    } // find(years_.begin(), years_.end(), model_->current_year()) != years_.end()
  }  //model_->state() != State::kInitialise
}


// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MortalityEventBiomass::FillReportCache(ostringstream& cache) {
  cache << "biomass_removed: ";
  for (auto& year : actual_removals_by_year_)
    cache << year.second << " ";
  cache << "\ncatch_input_removed: ";
  for (auto& year : removals_by_year_)
    cache << year.second << " ";
  cache << "\n";

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

    // Print census information
    for (auto& census : removals_census_) {
      cache << "census_info-" << census.year_ << "-" << census.row_ << "-" << census.col_ << " " << REPORT_R_MATRIX << "\n";
      //cache << "age length scalar\n";
      for (unsigned ndx = 0; ndx < census.age_.size(); ++ndx)
        cache << census.age_[ndx] << " ";
      cache << "\n";
      for (unsigned ndx = 0; ndx < census.age_.size(); ++ndx)
        cache << census.length_[ndx] << " ";
      cache << "\n";
      for (unsigned ndx = 0; ndx < census.age_.size(); ++ndx)
        cache << census.scalar_[ndx] << " ";
      cache << "\n";
    }
  }
}

// true means all years are found, false means there is a mismatch in years
bool MortalityEventBiomass::check_years(vector<unsigned> years_to_check_) {
  LOG_FINE();
  for (unsigned year_ndx = 0; year_ndx < years_to_check_.size(); ++year_ndx) {
    if (find(years_.begin(), years_.end(), years_to_check_[year_ndx]) == years_.end())
      return false;
  }
  return true;
}

} /* namespace processes */
} /* namespace niwa */
