/**
 * @file MortalityEffortBasedWithCovar.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 10/10/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */

// headers
#include "MortalityEffortBasedWithCovar.h"

#include "Layers/Manager.h"
#include "PreferenceFunctions/Manager.h"
#include "Selectivities/Manager.h"
#include "Minimisers/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/Math.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"



// namespaces
namespace niwa {
namespace processes {

namespace math = niwa::utilities::math;

/**
 * constructor
 */
MortalityEffortBasedWithCovar::MortalityEffortBasedWithCovar(Model* model) : Mortality(model) {
  parameters_.Bind<string>(PARAM_SELECTIVITY, &selectivity_label_, "Selectivity label", "");
  parameters_.Bind<string>(PARAM_MINIMIZER, &minimiser_label_, "Label of the minimser to solve the problem", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "years to apply the process", "");
  parameters_.Bind<float>(PARAM_CATCHES, &catches_, "Total catch by year","");
  parameters_.Bind<float>(PARAM_EFFORT_VALUES, &effort_input_, "A vector of effort values, one for each enabled cell of the model, these should represent the variability of effort of the fishery mimicking","");
 // parameters_.Bind<string>(PARAM_EFFORT_LAYER_LABEL, &effort_layer_label_, "A layer label that is a numeric layer label that contains effort values for each year.","","");
  parameters_.Bind<double>(PARAM_STARTING_VALUE_FOR_LAMBDA, &start_value_for_lambda_, "Total catch by year","", true);
  parameters_.Bind<string>(PARAM_PREFERENCE_FUNCTIONS, &preference_function_labels_, "The preference functions to apply to each layer", "");
  parameters_.Bind<string>(PARAM_PREFERENCE_LAYERS, &preference_layer_labels_, "The preference functions to apply", "");
  parameters_.Bind<float>(PARAM_PREFERENCE_WEIGHTS, &preference_weights_, "The weight to each preference when calculating mean spatial effort distribution", "");
  //parameters_.Bind<float>(PARAM_CATCHABILITY, &catchability_, "An arbiturary scalar to get the effort value","");
}

/**
 * Do some initial checks of user supplied parameters.
 */
void MortalityEffortBasedWithCovar::DoValidate() {
  LOG_TRACE();
  if (years_.size() != catches_.size())
    LOG_ERROR_P(PARAM_YEARS) << "you must specify a layer label for each year. You have supplied '" << years_.size() << "' years but '" << catches_.size() << "' catches, please sort this out.";
}

/**
 * DoBuild
 */
void MortalityEffortBasedWithCovar::DoBuild() {
  LOG_FINE();

  // Get the preference functions
  for (auto& label : preference_function_labels_) {
    PreferenceFunction* temp_func = nullptr;
    temp_func = model_->managers().preference_function()->GetPreferenceFunction(label);
    if (!temp_func) {
      LOG_FATAL_P(PARAM_PREFERENCE_FUNCTIONS) << "could not find the preference function '" << label << "', please make sure it exists'";
    }
    preference_functions_.push_back(temp_func);
  }
  LOG_FINEST();

  base_preference_.resize(model_->get_height());
  vary_preference_.resize(model_->get_height());
  temp_preference_.resize(model_->get_height());

  for (unsigned row = 0; row < model_->get_height(); ++row) {
    base_preference_[row].resize(model_->get_width(), 0.0);
    vary_preference_[row].resize(model_->get_width(), 0.0);
    temp_preference_[row].resize(model_->get_width(), 0.0);
  }

  if(preference_function_labels_.size() != preference_layer_labels_.size())
    LOG_ERROR_P(PARAM_PREFERENCE_LAYERS)  << "You need to supply the same number of preference functions as layers, reciegved " << preference_function_labels_.size() << " functions but " << preference_layer_labels_.size() << " layers";
  LOG_FINEST() << " number of preference functions = " << preference_function_labels_.size();


  // Get the layers
  unsigned layer_ndx = 0;
  for (auto& label : preference_layer_labels_) {
    LOG_FINEST() << "trying to get layer = " << label << " layer_ndx = " << layer_ndx;
    layers::NumericLayer* temp_layer = nullptr;
    temp_layer = model_->managers().layer()->GetNumericLayer(label);
    if (!temp_layer) {
      LOG_FATAL_P(PARAM_PREFERENCE_LAYERS) << "could not find the layer '" << label << "', please make sure it exists, and if it does exist make sure it is of type 'numeric''";
    }

    if (not temp_layer->is_static()) {
      LOG_FINEST() << "layer is not static";
      non_static_layer_ndx_.push_back(layer_ndx);
      calculate_on_the_fly_ = true;
    } else {
      LOG_FINEST() << "layer is static";
      // Can calculate preference values here and not worry about them later on as they wont change during model run time
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          LOG_FINEST() << "row = " << row + 1 << " col = " << col + 1 << " init value = " << base_preference_[row][col];

          // rest these containers
          base_preference_[row][col] += preference_weights_[layer_ndx] *  preference_functions_[layer_ndx]->get_result(temp_layer->get_value(row, col));  // the 0 in the layers call will return the default layer

        }
      }
    }
    //LOG_MEDIUM() << "brownian motion = " << brownian_motion_ << " calculate on the fly? " << calculate_on_the_fly_;
    LOG_FINEST();

    preference_layers_.push_back(temp_layer);
    ++layer_ndx;
  }


  LOG_MEDIUM() << "enabled cells = " << world_->get_enabled_cells();
  // Order the effort vector so it what?
  sort(effort_input_.begin(), effort_input_.end());  // low -> high

  // allocate memory for runtime containers.
  cell_offset_.resize(model_->get_height());
  effort_by_cell_.resize(model_->get_height());
  vulnerable_by_cell_.resize(model_->get_height());
  F_by_cell_.resize(model_->get_height());
  removals_by_cell_.resize(model_->get_height());

  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cell_offset_[i].resize(model_->get_width());
    effort_by_cell_[i].resize(model_->get_width(),0.0);
    F_by_cell_[i].resize(model_->get_width(),0.0);
    vulnerable_by_cell_[i].resize(model_->get_width(),0.0);
    removals_by_cell_[i].resize(model_->get_width(),0.0);
  }

  // Get the layers
  unsigned catch_ndx = 0;
  for (auto& year : years_) {
    catches_by_year_[year] = catches_[catch_ndx];
    ++catch_ndx;
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

  // Build GammaDiff
  minimiser_ =  model_->managers().minimiser()->get_minimiser(minimiser_label_);
  if (!minimiser_) {
    LOG_FATAL_P(PARAM_MINIMIZER) << "Could not find the @minimiser block labelled " << minimiser_label_ << " please make sure that it exists.";
  }
  LOG_FINE() << "about to call execute on the minimsier";
  minimiser_->PassBaranovProcess(label_);

  if (parameters_.Get(PARAM_STARTING_VALUE_FOR_LAMBDA)->has_been_defined()) {
    if (start_value_for_lambda_.size() != 1)
      LOG_FATAL_P(PARAM_STARTING_VALUE_FOR_LAMBDA) << "Can only supply a single value for this parameter, you have specified " << start_value_for_lambda_.size();
    minimiser_->PassStartValueBaranov(start_value_for_lambda_);
  }

  age_freq_census_.resize(model_->age_spread());


}

/**
 * DoReset
 */
void MortalityEffortBasedWithCovar::DoReset() {
  LOG_FINE() << "clearing containers";
  removals_by_age_and_area_.clear();
  removals_by_length_and_area_.clear();
  removals_census_.clear();
  removals_tag_recapture_.clear();
}

/**
 * DoExecute
 */
void MortalityEffortBasedWithCovar::DoExecute() {
  LOG_MEDIUM() << "DoExecute";
  if (first_execute_) {
    // Bit of an annoying hack, but cannot check at build because of Building mortality classes before World, blah TODO move this process into the PreWorldProcessBuild in Process Manager.
    if (world_->get_enabled_cells() != effort_input_.size()) {
      LOG_FATAL_P(PARAM_EFFORT_VALUES) << "You must specify an effort value for each enabled (Base layer > 0) cell, there are " << world_->get_enabled_cells() << " Base cells, but you have supplied " << effort_input_.size() << ". Can you please sort this descrepency out";
    }

    vulnerable_biomass_vector_format_.resize(world_->get_enabled_cells(), 0.0);
    effort_organised_vector_format_.resize(world_->get_enabled_cells(), 0.0);
    pref_organised_vector_format_.resize(world_->get_enabled_cells(), 0.0);
    first_execute_ = false;
  }
  auto iter = years_.begin();
  if (model_->state() != State::kInitialise) {
    if (find(iter, years_.end(), model_->current_year()) != years_.end()) {
      iter = find(years_.begin(), years_.end(), model_->current_year());
      unsigned catch_ndx = distance(years_.begin(), iter);
      LOG_MEDIUM() << "applying F in year " << model_->current_year() << " catch index = " << catch_ndx;
      // Pre-calculate agents in the world to set aside our random numbers needed for the operation
      n_agents_ = 0;
      unsigned cell_counter = 0;

      // Calculate preference for time-varying layers
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        fill(vary_preference_[row].begin(), vary_preference_[row].end(), 0.0);
        fill(removals_by_cell_[row].begin(), removals_by_cell_[row].end(), 0.0);
        fill(vulnerable_by_cell_[row].begin(), vulnerable_by_cell_[row].end(), 0.0);
        fill(effort_by_cell_[row].begin(), effort_by_cell_[row].end(), 0.0);
        fill(F_by_cell_[row].begin(), F_by_cell_[row].end(), 0.0);


        for (unsigned col = 0; col < model_->get_width(); ++col) {
          effort_by_cell_[row][col] = 0.0;
          removals_by_cell_[row][col] = 0.0;
          vulnerable_by_cell_[row][col] = 0.0;
          WorldCell *cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            cell_offset_[row][col] = n_agents_;
            n_agents_ += cell->agents_.size();
            for (unsigned layer_ndx = 0; layer_ndx < non_static_layer_ndx_.size(); ++layer_ndx) {
              vary_preference_[row][col] += preference_weights_[non_static_layer_ndx_[layer_ndx]] * preference_functions_[non_static_layer_ndx_[layer_ndx]]->get_result(preference_layers_[non_static_layer_ndx_[layer_ndx]]->get_value(row, col, model_->current_year()));
            }
          }
        }
      }


      fill(age_freq_census_.begin(), age_freq_census_.end(), 0.0);

      vulnerable_biomass_by_year_[model_->current_year()] = 0;
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
      cell_counter = 0;
      max_vulnerable_ = 0.0;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {

          WorldCell* cell = world_->get_base_square(row, col);

          if (cell->is_enabled()) {
            LOG_MEDIUM() << "checking cell in row " << row + 1 << " col = " << col + 1;
            // iterate through and calcualte vulnerable biomass in each cell exactly, nothing random here
            unsigned counter = 0;
            if (selectivity_length_based_) {
              for (auto agent_iter = cell->agents_.begin(); agent_iter != cell->agents_.end(); ++counter,++agent_iter) {
                //LOG_FINE() << "counter = " << counter;
                if ((*agent_iter).is_alive()) {
                  if (random_numbers_[cell_offset_[row][col] + counter] <= selectivity_[(*agent_iter).get_sex()]->GetResult((*agent_iter).get_length_bin_index()))
                    vulnerable_by_cell_[row][col] += (*agent_iter).get_weight() * (*agent_iter).get_scalar();
                }
              }
              for (auto agent_iter = cell->tagged_agents_.begin(); agent_iter != cell->tagged_agents_.end(); ++counter,++agent_iter) {
                if ((*agent_iter).is_alive()) {
                  if (random_numbers_[cell_offset_[row][col] + counter] <= selectivity_[(*agent_iter).get_sex()]->GetResult((*agent_iter).get_length_bin_index()))
                    vulnerable_by_cell_[row][col] += (*agent_iter).get_weight() * (*agent_iter).get_scalar();
                }
              }
            } else {
              for (auto agent_iter = cell->agents_.begin(); agent_iter != cell->agents_.end(); ++counter,++agent_iter) {
                //LOG_FINE() << "counter = " << counter;
                if ((*agent_iter).is_alive()) {
                  if (random_numbers_[cell_offset_[row][col] + counter] <= selectivity_[(*agent_iter).get_sex()]->GetResult((*agent_iter).get_age_index()))
                    vulnerable_by_cell_[row][col] += (*agent_iter).get_weight() * (*agent_iter).get_scalar();
                }
              }
              for (auto agent_iter = cell->tagged_agents_.begin(); agent_iter != cell->tagged_agents_.end(); ++counter,++agent_iter) {
                //LOG_FINE() << "counter = " << counter;
                if ((*agent_iter).is_alive()) {
                  if (random_numbers_[cell_offset_[row][col] + counter] <= selectivity_[(*agent_iter).get_sex()]->GetResult((*agent_iter).get_age_index()))
                    vulnerable_by_cell_[row][col] += (*agent_iter).get_weight() * (*agent_iter).get_scalar();
                }
              }
            }
            LOG_MEDIUM() << "Biomass = " << vulnerable_by_cell_[row][col] << " cell counter = " << cell_counter << " size = " << vulnerable_biomass_vector_format_.size();

            vulnerable_biomass_by_year_[model_->current_year()] += vulnerable_by_cell_[row][col];

            vulnerable_biomass_vector_format_[cell_counter] = vulnerable_by_cell_[row][col];
            ++cell_counter;

            //effort_by_cell_[row][col] = vulnerable_by_cell_[row][col];// * catchability_;
            if(vulnerable_by_cell_[row][col] > max_vulnerable_)
              max_vulnerable_ = vulnerable_by_cell_[row][col];
          }
        }
      }
      LOG_MEDIUM() << "year = " << model_->current_year() << " max vulnerable biomasss = " << max_vulnerable_;
      cell_counter = 0;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          temp_preference_[row][col] = base_preference_[row][col] + vary_preference_[row][col] + vulnerable_by_cell_[row][col] / max_vulnerable_;
          temp_preference_[row][col] /= (preference_layer_labels_.size() + 1);
          pref_organised_vector_format_[cell_counter] = temp_preference_[row][col];
          cell_counter++;
        }
      }
      preference_by_year_[model_->current_year()] = temp_preference_;
      //
      // Final attribute effort based on an ideal free distribution The highest vulnerable biomass gets the most effort.
      effort_index_ = math::sort_indexes(pref_organised_vector_format_);
      for(unsigned i = 0; i < vulnerable_biomass_vector_format_.size(); ++i) {
        LOG_MEDIUM() << "i  = " << i << " biomass ndx = " << effort_index_[i]<< " effort value = " << effort_input_[i] << " biomass = " << pref_organised_vector_format_[effort_index_[i]];
        effort_organised_vector_format_[effort_index_[i]] = effort_input_[i];
      }

      cell_counter = 0;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            effort_by_cell_[row][col] = effort_organised_vector_format_[cell_counter];
            cell_counter++;
          }
        }
      }

      LOG_MEDIUM() << "about to check baranov";
      time(&start_time_);
      minimiser_->SolveBaranov();
      time_by_year_[model_->current_year()] = (time(NULL) - start_time_); // seconds
      message_by_year_[model_->current_year()] = minimiser_->get_message();
      LOG_MEDIUM() << "Finished minimsation lambda value = " << lambda_  << " time now = " << time(NULL) << " start time = " << start_time_;
      catch_based_on_baranov_by_year_[model_->current_year()] = catch_based_on_baranov_;
      lambda_by_year_[model_->current_year()] = lambda_;

      // Now apply Mortality with probability with F as we do with M in mortality constant
      actual_catch_ = 0.0;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            LOG_FINE() << "checking cell in row " << row + 1 << " col = " << col + 1;
            composition_data age_freq(PARAM_AGE, model_->current_year(), row, col,  model_->age_spread());
            composition_data length_freq(PARAM_LENGTH, model_->current_year(), row, col,  model_->length_bin_mid_points().size());
            census_data census_fishery(model_->current_year(), row, col);

            // iterate through and calcualte vulnerable biomass in each cell exactly, nothing random here
            unsigned counter = 0;
            double F_this_cell = lambda_ * effort_by_cell_[row][col];
            F_by_cell_[row][col] = F_this_cell;
            LOG_FINE() << "About to kill agents.";
            if (selectivity_length_based_) {
              for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++counter, ++iter) {
                //LOG_MEDIUM() << "rand number = " << random_numbers_[cell_offset_[row][col] + counter] << " selectivity = " << cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()] << " is alove = " << (*iter).is_alive() << " counter = " << counter;
                if ((*iter).is_alive()) {
                  if (random_numbers_[cell_offset_[row][col] + counter] <= (1.0 - std::exp(-F_this_cell * selectivity_[(*iter).get_sex()]->GetResult((*iter).get_length_bin_index())))) {
                    actual_catch_ += (*iter).get_weight() * (*iter).get_scalar();
                    removals_by_cell_[row][col] += (*iter).get_weight() * (*iter).get_scalar();
                    age_freq.frequency_[(*iter).get_age_index()] += (*iter).get_scalar(); // This catch actually represents many individuals.
                    age_freq.biomass_ += (*iter).get_weight() * (*iter).get_scalar();
                    length_freq.frequency_[(*iter).get_length_bin_index()] += (*iter).get_scalar();
                    length_freq.biomass_ += (*iter).get_weight() * (*iter).get_scalar();
                    census_fishery.age_ndx_.push_back((*iter).get_age_index());
                    census_fishery.length_ndx_.push_back((*iter).get_length_bin_index());
                    census_fishery.scalar_.push_back((*iter).get_scalar());
                    census_fishery.biomass_+= (*iter).get_weight() * (*iter).get_scalar();
                    census_fishery.weight_.push_back((*iter).get_weight());
                    census_fishery.sex_.push_back((*iter).get_sex());
                    age_freq_census_[(*iter).get_age_index()] += (*iter).get_scalar();
                    cell->remove_agent_alive((*iter).get_scalar());
                    (*iter).dies();
                  }
                }
              }
            } else {
              for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++counter, ++iter) {
                //LOG_MEDIUM() << "rand number = " << random_numbers_[cell_offset_[row][col] + counter] << " selectivity = " << cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()] << " is alove = " << (*iter).is_alive() << " counter = " << counter;
                if ((*iter).is_alive()) {
                  if (random_numbers_[cell_offset_[row][col] + counter]
                      <= (1.0 - std::exp(-F_this_cell * selectivity_[(*iter).get_sex()]->GetResult((*iter).get_age_index())))) {
                    actual_catch_ += (*iter).get_weight() * (*iter).get_scalar();
                    removals_by_cell_[row][col] += (*iter).get_weight() * (*iter).get_scalar();
                    age_freq.frequency_[(*iter).get_age_index()] += (*iter).get_scalar(); // This catch actually represents many individuals.
                    age_freq.biomass_ += (*iter).get_weight() * (*iter).get_scalar();
                    length_freq.frequency_[(*iter).get_length_bin_index()] += (*iter).get_scalar();
                    length_freq.biomass_ += (*iter).get_weight() * (*iter).get_scalar();
                    census_fishery.age_ndx_.push_back((*iter).get_age_index());
                    census_fishery.length_ndx_.push_back((*iter).get_length_bin_index());
                    census_fishery.scalar_.push_back((*iter).get_scalar());
                    census_fishery.biomass_+= (*iter).get_weight() * (*iter).get_scalar();
                    census_fishery.weight_.push_back((*iter).get_weight());
                    census_fishery.sex_.push_back((*iter).get_sex());
                    age_freq_census_[(*iter).get_age_index()] += (*iter).get_scalar();
                    cell->remove_agent_alive((*iter).get_scalar());
                    (*iter).dies();
                  }
                }
              }
            }
            LOG_FINE() << "Applied mortality";
            removals_by_length_and_area_.push_back(length_freq);
            removals_by_age_and_area_.push_back(age_freq);
            removals_census_.push_back(census_fishery);
          }
        }
      }
      actual_catch_by_year_[model_->current_year()] = actual_catch_;
      actual_removals_by_year_and_cell_[model_->current_year()] = removals_by_cell_;
      F_by_year_and_cell_[model_->current_year()] = F_by_cell_;
      vulnerable_by_year_and_cell_[model_->current_year()] = vulnerable_by_cell_;
      effort_by_year_and_cell_[model_->current_year()] = effort_by_cell_;
      age_freq_census_by_year_[model_->current_year()] = age_freq_census_;

    } // year
  } // initialisation
} // DoExecute

/**
 * Evaluate SSE for catch using the baranov catch equation
 * There is stochastic fishing in the sense we generate random numbers once to evaluate the vulnerability of agents to fishing to generate vulnerable biomass
 * After that, we solve for (observed_catch - est_catch)^2  where est_catch =
 *
 */
double MortalityEffortBasedWithCovar::SolveBaranov() {
  LOG_FINE();
  catch_based_on_baranov_ = 0;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        unsigned counter = 0;
        double F_this_cell = lambda_ * effort_by_cell_[row][col];
        if (selectivity_length_based_) {
          for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++iter, ++counter) {
            //LOG_FINEST() << "rand number = " << random_numbers_[cell_offset_[row][col] + counter] << " selectivity = " << cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()] << " survivorship = " << (1 - std::exp(-(*iter).get_m() *  cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()])) << " M = " << (*iter).get_m();
            if ((*iter).is_alive()) {
              if (random_numbers_[cell_offset_[row][col] + counter] <= (1.0 - std::exp(-F_this_cell * selectivity_[(*iter).get_sex()]->GetResult((*iter).get_length_bin_index())))) {
                catch_based_on_baranov_ += (*iter).get_weight() * (*iter).get_scalar();
              }
            }
          }
        } else {
          for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++iter, ++counter) {
            //LOG_FINEST() << "rand number = " << random_numbers_[cell_offset_[row][col] + counter] << " selectivity = " << cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()] << " survivorship = " << (1 - std::exp(-(*iter).get_m() *  cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()])) << " M = " << (*iter).get_m();
            if ((*iter).is_alive()) {
              if (random_numbers_[cell_offset_[row][col] + counter] <= (1.0 - std::exp(-F_this_cell * selectivity_[(*iter).get_sex()]->GetResult((*iter).get_age_index())))) {
                catch_based_on_baranov_ += (*iter).get_weight() * (*iter).get_scalar();
              }
            }
          }
        }
      }
    }
  }
  //
  LOG_MEDIUM() << "proposed catch = " << catch_based_on_baranov_ << " catch = " << catches_by_year_[model_->current_year()] <<  " SSE = " << pow(catches_by_year_[model_->current_year()] - catch_based_on_baranov_,2) << " year = " << model_->current_year() << " lambda = " << lambda_;

  return pow(catches_by_year_[model_->current_year()] - catch_based_on_baranov_,2);
}


// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MortalityEffortBasedWithCovar::FillReportCache(ostringstream& cache) {
  cache << "actual_catch: ";
  for(auto& val : actual_catch_by_year_)
    cache << val.second << " ";
  cache << "\ncatch_based_on_baranov: ";
  for(auto& val : catch_based_on_baranov_by_year_)
    cache << val.second << " ";
  cache << "\nvulnerable_biomass: ";
  for(auto& val : vulnerable_biomass_by_year_)
    cache << val.second << " ";
  cache << "\nseconds_to_minimise: ";
  for(auto& val : time_by_year_)
    cache << val.second << " ";
  cache << "\nlambda: ";
  for(auto& val : lambda_by_year_)
    cache << val.second << " ";
  cache << "\n";

  for(auto& val : message_by_year_) {
    cache << val.first  << "_message " << REPORT_R_STRING_VECTOR << "\n";
    cache << val.second << "\n";
  }

  for (auto& values : F_by_year_and_cell_) {
    cache << "F_by_cell_" << values.first << " " << REPORT_R_MATRIX << "\n";
    for (unsigned i = 0; i < values.second.size(); ++i) {
      for (unsigned j = 0; j < values.second[i].size(); ++j)
        cache << values.second[i][j] << " ";
      cache << "\n";
    }
  }


  cache << "AgeFrequency" << " " << REPORT_R_DATAFRAME_ROW_LABELS << "\nyear ";
  for (unsigned age = model_->min_age(); age <= model_->max_age(); ++age)
    cache << age << " ";
  cache << "\n";
  for (auto& values : age_freq_census_by_year_) {
    cache << values.first << " ";
    for (auto& age_val : values.second)
      cache << age_val << " ";
    cache << "\n";
  }

  cache << "base_preference " << REPORT_R_MATRIX << "\n";
  for (unsigned i = 0; i < base_preference_.size(); ++i) {
    for (unsigned j = 0; j < base_preference_[i].size(); ++j)
      cache << base_preference_[i][j] << " ";
    cache << "\n";
  }



  for (auto& values : preference_by_year_) {
    cache << "preference_by_cell_" << values.first << " " << REPORT_R_MATRIX << "\n";
    for (unsigned i = 0; i < values.second.size(); ++i) {
      for (unsigned j = 0; j < values.second[i].size(); ++j)
        cache << values.second[i][j] << " ";
      cache << "\n";
    }
  }


  for (auto& values : effort_by_year_and_cell_) {
    cache << "effort_by_cell_" << values.first << " " << REPORT_R_MATRIX << "\n";
    for (unsigned i = 0; i < values.second.size(); ++i) {
      for (unsigned j = 0; j < values.second[i].size(); ++j)
        cache << values.second[i][j] << " ";
      cache << "\n";
    }
  }

  for (auto& values : vulnerable_by_year_and_cell_) {
    cache << "vulnerable_by_cell_" << values.first << " " << REPORT_R_MATRIX << "\n";
    for (unsigned i = 0; i < values.second.size(); ++i) {
      for (unsigned j = 0; j < values.second[i].size(); ++j)
        cache << values.second[i][j] << " ";
      cache << "\n";
    }
  }

  for (auto& values : actual_removals_by_year_and_cell_) {
    cache << "removals_" << values.first << " " << REPORT_R_MATRIX << "\n";
    for (unsigned i = 0; i < values.second.size(); ++i) {
      for (unsigned j = 0; j < values.second[i].size(); ++j)
        cache << values.second[i][j] << " ";
      cache << "\n";
    }
  }

}



} /* namespace processes */
} /* namespace niwa */
