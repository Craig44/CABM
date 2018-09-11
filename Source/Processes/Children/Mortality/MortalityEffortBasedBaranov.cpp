/**
 * @file MortalityEffortBasedBaranov.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 10/09/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */

// headers
#include "MortalityEffortBasedBaranov.h"

#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "Minimisers/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
MortalityEffortBasedBaranov::MortalityEffortBasedBaranov(Model* model) : Mortality(model) {
  parameters_.Bind<string>(PARAM_SELECTIVITY, &selectivity_label_, "Selectivity label", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "years to apply the process", "");
  parameters_.Bind<float>(PARAM_CATCHES, &catches_, "Total catch by year","");
  parameters_.Bind<float>(PARAM_CATCHABILITY, &catchability_, "An arbiturary scalar to get the effort value","");

}

/**
 * Do some initial checks of user supplied parameters.
 */
void MortalityEffortBasedBaranov::DoValidate() {
  LOG_TRACE();
  if (years_.size() != catches_.size())
    LOG_ERROR_P(PARAM_YEARS) << "you must specify a layer label for each year. You have supplied '" << years_.size() << "' years but '" << catches_.size() << "' catches, please sort this out.";
}

/**
 * DoBuild
 */
void MortalityEffortBasedBaranov::DoBuild() {
  LOG_FINE();

  // allocate memory for runtime containers.
  cell_offset_for_selectivity_.resize(model_->get_height());
  cell_offset_.resize(model_->get_height());
  effort_by_cell_.resize(model_->get_height());
  vulnerable_by_cell_.resize(model_->get_height());
  removals_by_cell_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cell_offset_[i].resize(model_->get_width());
    cell_offset_for_selectivity_[i].resize(model_->get_width());
    effort_by_cell_[i].resize(model_->get_width(),0.0);
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



  if (selectivity_length_based_) {
    LOG_FATAL_P(PARAM_SELECTIVITY_LABEL) << "haven't coded for length based selectivity in this process.";
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
  // Build GammaDiff
  minimiser_ =  model_->managers().minimiser()->active_minimiser();
  if (!minimiser_) {
    LOG_FATAL_P(PARAM_LABEL) << "you need to declare a @minimiser block to use this process, and it needs to be of type gamma_diff. See the manual for more information";
  }
  LOG_FINE() << "about to call execute on the minimsier";
  minimiser_->PassBaranovProcess(label_);

  //minimiser_->Execute();
  // Prepare the model to solve the baranov catch equation.
  //model_->prepare_for_effort_based_baranov(label_);

}

/**
 * DoExecute
 */
void MortalityEffortBasedBaranov::DoExecute() {
  LOG_FINE() << "DoExecute";

  auto iter = years_.begin();
  if (model_->state() != State::kInitialise) {
    if (find(iter, years_.end(), model_->current_year()) != years_.end()) {
      iter = find(years_.begin(), years_.end(), model_->current_year());
      unsigned catch_ndx = distance(years_.begin(), iter);
      LOG_FINE() << "applying F in year " << model_->current_year() << " catch index = " << catch_ndx;
      // Pre-calculate agents in the world to set aside our random numbers needed for the operation
      n_agents_ = 0;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          WorldCell* cell = world_->get_base_square(row, col);
          // reset these containers
          effort_by_cell_[row][col] = 0.0;
          removals_by_cell_[row][col] = 0.0;
          vulnerable_by_cell_[row][col] = 0.0;
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
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {

          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            LOG_FINE() << "checking cell in row " << row + 1 << " col = " << col + 1;
            // iterate through and calcualte vulnerable biomass in each cell exactly, nothing random here
            unsigned counter = 0;
            for (auto agent_iter = cell->agents_.begin(); agent_iter != cell->agents_.end(); ++counter,++agent_iter) {
              //LOG_FINE() << "counter = " << counter;
              vulnerable_by_cell_[row][col] += (*agent_iter).get_weight() * (*agent_iter).get_scalar() * cell_offset_for_selectivity_[row][col][(*agent_iter).get_sex() * model_->age_spread() + (*agent_iter).get_age_index()];
            }
            effort_by_cell_[row][col] = vulnerable_by_cell_[row][col] * catchability_;
          }
        }
      }
      LOG_FINE() << "about to ckech baranov";
      minimiser_->SolveBaranov();
      LOG_FINE() << "Finished temp value";
      catch_based_on_baranov__by_year_[model_->current_year()] = catch_based_on_baranov_;
      lambda_by_year_[model_->current_year()] = lambda_;

      // Now apply Mortality with probability with F as we do with M in mortality constant
      actual_catch_ = 0.0;
      for (unsigned row = 0; row < model_->get_height(); ++row) {
        for (unsigned col = 0; col < model_->get_width(); ++col) {
          WorldCell* cell = world_->get_base_square(row, col);
          if (cell->is_enabled()) {
            LOG_FINE() << "checking cell in row " << row + 1 << " col = " << col + 1;
            // iterate through and calcualte vulnerable biomass in each cell exactly, nothing random here
            unsigned counter = 0;
            double F_this_cell = lambda_ * effort_by_cell_[row][col];
            vulnerable_by_cell_[row][col] = F_this_cell;
            unsigned initial_size = cell->agents_.size();
            LOG_FINEST() << initial_size << " initial agents";
            for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++counter) {
              //LOG_FINEST() << "rand number = " << random_numbers_[cell_offset_[row][col] + counter] << " selectivity = " << cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()] << " survivorship = " << (1 - std::exp(-(*iter).get_m() *  cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()])) << " M = " << (*iter).get_m();
              if (random_numbers_[cell_offset_[row][col] + counter]
                  <= (1.0 - std::exp(-F_this_cell * cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()]))) {
                actual_catch_ += (*iter).get_weight() * (*iter).get_scalar();
                removals_by_cell_[row][col] += (*iter).get_weight() * (*iter).get_scalar();
                iter = cell->agents_.erase(iter);
                initial_size--;
              } else
                ++iter;
            }
            LOG_FINEST() << initial_size << " after mortality";
          }
        }
      }
      actual_catch_by_year_[model_->current_year()] = actual_catch_;
      actual_removals_by_year_and_cell_[model_->current_year()] = removals_by_cell_;
      F_by_year_and_cell_[model_->current_year()] = vulnerable_by_cell_;
    } // year
  } // initialisation
} // DoExecute

/**
 * Evaluate SSE for catch using the baranov catch equation
 */
double MortalityEffortBasedBaranov::SolveBaranov() {
  LOG_FINE();
  catch_based_on_baranov_ = 0;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        unsigned counter = 0;
        double F_this_cell = lambda_ * effort_by_cell_[row][col];
        for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++iter, ++counter) {
          //LOG_FINEST() << "rand number = " << random_numbers_[cell_offset_[row][col] + counter] << " selectivity = " << cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()] << " survivorship = " << (1 - std::exp(-(*iter).get_m() *  cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()])) << " M = " << (*iter).get_m();
          if (random_numbers_[cell_offset_[row][col] + counter]
              <= (1.0 - std::exp(-F_this_cell * cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()]))) {
            catch_based_on_baranov_ += (*iter).get_weight() * (*iter).get_scalar();
          }
        }
      }
    }
  }
  //
  return pow(catches_by_year_[model_->current_year()] - catch_based_on_baranov_,2);
}


// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MortalityEffortBasedBaranov::FillReportCache(ostringstream& cache) {
  cache << "actual_catch: ";
  for(auto& val : actual_catch_by_year_)
    cache << val.second << " ";
  cache << "\ncatch_based_on_baranov: ";
  for(auto& val : catch_based_on_baranov__by_year_)
    cache << val.second << " ";
  cache << "\nlambda: ";
  for(auto& val : lambda_by_year_)
    cache << val.second << " ";
  cache << "\n";

  for (auto& values : F_by_year_and_cell_) {
    cache << "F_by_cell_" << values.first << " " << REPORT_R_MATRIX << "\n";
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
