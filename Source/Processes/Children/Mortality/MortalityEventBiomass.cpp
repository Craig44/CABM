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
  parameters_.Bind<float>(PARAM_MINIMUM_LEGAL_LENGTH, &mls_, "The minimum legal length for this fishery, any individual less than this will be returned using some discard mortality", "");
  parameters_.Bind<float>(PARAM_HANDLING_MORTALITY, &discard_mortality_, "if discarded due to being under the minimum legal length, what is the probability the individual will die when released", "");
  // TODO Tagging events
}

/**
 * Do some initial checks of user supplied parameters.
 */
void MortalityEventBiomass::DoValidate() {
  LOG_TRACE();
  if (years_.size() != catch_layer_label_.size())
    LOG_ERROR_P(PARAM_YEARS) << "you must specify a layer label for each year. You have supplied '" << years_.size() << "' years but '" << catch_layer_label_.size() << "' catch layer labels, please sort this out.";
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

  cell_offset_for_selectivity_.resize(model_->get_height());
  cell_offset_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cell_offset_[i].resize(model_->get_width());
    cell_offset_for_selectivity_[i].resize(model_->get_width());
  }

  if (selectivity_length_based_) {
    for (unsigned i = 0; i < model_->get_height(); ++i) {
      for (unsigned j = 0; j < model_->get_width(); ++j) {
        for (unsigned ogive = 0; ogive < selectivity_label_.size(); ++ogive) {
          for (auto len : model_->length_bins())
            cell_offset_for_selectivity_[i][j].push_back(selectivity_[ogive]->GetResult(len));
        }
      }
    }
  } else {
    for (unsigned i = 0; i < model_->get_height(); ++i) {
      for (unsigned j = 0; j < model_->get_width(); ++j) {
        for (unsigned ogive = 0; ogive < selectivity_label_.size(); ++ogive) {
          for (auto age = model_->min_age(); age <= model_->max_age(); ++age)
          cell_offset_for_selectivity_[i][j].push_back(selectivity_[ogive]->GetResult(age));
        }
      }
    }
  }
}


/**
 * DoExecute
 */
void MortalityEventBiomass::DoExecute() {
  LOG_FINE();
  vector<unsigned> global_age_freq(model_->age_spread());
  auto iter = years_.begin();
  if (model_->state() != State::kInitialise) {
    if (find(iter, years_.end(), model_->current_year()) != years_.end()) {
      iter = find(years_.begin(), years_.end(), model_->current_year());
      unsigned catch_ndx = distance(years_.begin(), iter);
      LOG_FINEST() << "applying F in year " << model_->current_year() << " catch index = " << catch_ndx;
      // Pre-calculate agents in the world to set aside our random numbers needed for the operation
      n_agents_ = 0;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            cell_offset_[row][col] = n_agents_;
            n_agents_ += cell->agents_.size();
          }
        }
      }

      // Allocate a single block of memory rather than each thread temporarily allocating their own memory.
      random_numbers_.resize(n_agents_ + 1);
      discard_random_numbers_.resize(n_agents_ + 1);
      selectivity_random_numbers_.resize(n_agents_ + 1);
      utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
      for (unsigned i = 0; i <= n_agents_; ++i) {
        random_numbers_[i] = rng.chance();
        discard_random_numbers_[i] = rng.chance();
        selectivity_random_numbers_[i] = rng.chance();
      }


      float actual_catch_taken = 0;;
      float world_catch_to_take = 0;;
      if (selectivity_length_based_) {
        // Thread out each loop
        #pragma omp parallel for collapse(2)
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            unsigned catch_attempts = 1;
            unsigned catch_max = 1;
            unsigned random_agent;
            WorldCell* cell = world_->get_base_square(row, col);
            if (cell->is_enabled()) {
              float catch_taken = catch_layer_[catch_ndx]->get_value(row, col);
              unsigned counter = 0;
              if (catch_taken > 0) {
                LOG_FINEST() << "We are fishing in cell " << row + 1 << " " << col + 1 << " value = " << catch_taken;
                composition_data age_freq(PARAM_AGE, model_->current_year(), row, col, model_->age_spread());
                composition_data length_freq(PARAM_LENGTH, model_->current_year(), row, col, model_->length_bins().size());
                catch_attempts = 1;
                world_catch_to_take += catch_taken;
                catch_max = cell->agents_.size();
                while (catch_taken > 0) {
                  ++catch_attempts;
                  auto iter = cell->agents_.begin();
                  random_agent = random_numbers_[cell_offset_[row][col] + counter] * cell->agents_.size();
                  advance(iter, random_agent);
                  // See if this agent is unlucky
                  if (selectivity_random_numbers_[cell_offset_[row][col] + counter] <= cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()]) {
                    if ((*iter).get_length() < mls_) {
                      if (discard_random_numbers_[cell_offset_[row][col] + counter] <= discard_mortality_) {
                        LOG_MEDIUM();
                        cell->agents_.erase(iter); // erase agent from discard mortality
                      }
                    } else {
                      catch_taken -= (*iter).get_weight() * (*iter).get_scalar();
                      actual_catch_taken += (*iter).get_weight() * (*iter).get_scalar();
                      age_freq.frequency_[(*iter).get_age() - model_->min_age()]++;
                      length_freq.frequency_[(*iter).get_length_bin_index()]++;
                      global_age_freq[(*iter).get_age() - model_->min_age()]++;
                      cell->agents_.erase(iter); // erase agent from memory
                    }
                  }
                  // Make sure we don't end up fishing for infinity
                  if (catch_attempts >= catch_max) {
                    LOG_FATAL_P(PARAM_LABEL) << "Too many attempts to catch an agent in the process " << label_ << " in year " << model_->current_year() << " in row " << row + 1 << " and column " << col + 1 << " this most likely means you have" <<
                       " a model that suggests there should be more agents in this space than than the current agent dynamics are putting in this cell, check the user manual for tips to resolve this situation, agents in cell = " << catch_max << " attempts made = " << catch_attempts;
                  }
                  ++counter;
                }
                removals_by_length_and_area_.push_back(length_freq);
                removals_by_age_and_area_.push_back(age_freq);
              }

            }
          }
        }
      } else {
        iter = find(years_.begin(), years_.end(), model_->current_year());
        unsigned catch_ndx = distance(years_.begin(), iter);
        LOG_FINEST() << "applying F in year " << model_->current_year() << " catch index = " << catch_ndx;
        // Get the pointer to the right catch layer
        float actual_catch_taken = 0;;
        float world_catch_to_take = 0;;
        // Thread out each loop
        #pragma omp parallel for collapse(2)
        for (unsigned row = 0; row < model_->get_height(); ++row) {
          for (unsigned col = 0; col < model_->get_width(); ++col) {
            unsigned catch_attempts = 1;
            unsigned catch_max = 1;
            unsigned random_agent;
            WorldCell* cell = world_->get_base_square(row, col);
            if (cell->is_enabled()) {

              float catch_taken = catch_layer_[catch_ndx]->get_value(row, col);
              unsigned counter = 0;

              if (catch_taken > 0) {
                LOG_FINEST() << "We are fishing in cell " << row + 1 << " " << col + 1 << " value = " << catch_taken;
                // create reporting class
                composition_data age_freq(PARAM_AGE, model_->current_year(), row, col, model_->age_spread());
                composition_data length_freq(PARAM_LENGTH, model_->current_year(), row, col, model_->length_bins().size());

                catch_attempts = 1;
                world_catch_to_take += catch_taken;
                catch_max = cell->agents_.size();
                while (catch_taken > 0) {
                  // Random access bullshit for lists
                  ++catch_attempts;
                  auto iter = cell->agents_.begin();
                  random_agent = random_numbers_[cell_offset_[row][col] + counter] * cell->agents_.size();
                  advance(iter, random_agent);
                  // See if this agent is unlucky
                  if (selectivity_random_numbers_[cell_offset_[row][col] + counter] <= cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->length_bins().size() + (*iter).get_length_bin_index()]) {
                    if ((*iter).get_length() < mls_) {
                      if (discard_random_numbers_[cell_offset_[row][col] + counter] <= discard_mortality_) {
                        cell->agents_.erase(iter); // erase agent from discard mortality
                      }
                    } else {
                      catch_taken -= (*iter).get_weight() * (*iter).get_scalar();
                      actual_catch_taken += (*iter).get_weight() * (*iter).get_scalar();
                      age_freq.frequency_[(*iter).get_age() - model_->min_age()]++;
                      length_freq.frequency_[(*iter).get_length_bin_index()]++;
                      global_age_freq[(*iter).get_age() - model_->min_age()]++;
                      cell->agents_.erase(iter); // erase agent from memory
                    }
                  }
                  // Make sure we don't end up fishing for infinity if there are not enough fish here
                  if (catch_attempts >= catch_max) {
                    LOG_FATAL_P(PARAM_LABEL) << "Too many attempts to catch an agent in the process " << label_ << " in year " << model_->current_year() << " in row " << row + 1 << " and column " << col + 1 << " this most likely means you have\n" <<
                       " a model that suggests there should be more agents in this space than than the current agent dynamics are putting in this cell, check the user manual for tips to resolve this situation";
                  }
                  ++counter;
                }
                removals_by_length_and_area_.push_back(length_freq);
                removals_by_age_and_area_.push_back(age_freq);
              }

            }
          }
        }
      }
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
  cache << "year: ";
  for (auto& year : actual_removals_by_year_)
    cache << year.first << " ";

  cache << "\nbiomass_removed: ";
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
}

} /* namespace processes */
} /* namespace niwa */
