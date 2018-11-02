/**
 * @file MortalityConstantRate.cpp
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
#include "MortalityConstantRate.h"

#include "Layers/Manager.h"
#include "Selectivities/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "TimeSteps/Manager.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
MortalityConstantRate::MortalityConstantRate(Model* model) : Mortality(model) {
  parameters_.Bind<string>(PARAM_DISTRIBUTION, &distribution_, "the distribution to allocate the parameters to the agents","", PARAM_NORMAL)->set_allowed_values({PARAM_NORMAL, PARAM_LOGNORMAL});
  parameters_.Bind<float>(PARAM_CV, &cv_, "The cv of the distribution", "");
  parameters_.Bind<string>(PARAM_M_MULTIPLIER_LAYER_LABEL, &m_layer_label_, "Label for the numeric layer that describes a multiplier of M through space", "", ""); // TODO perhaps as a multiplier, 1.2 * 0.2 = 0.24
  parameters_.Bind<string>(PARAM_SELECTIVITY_LABEL, &selectivity_label_, "Label for the selectivity block", "");
  parameters_.Bind<float>(PARAM_M, &m_, "Natural mortality for the model", "");
  parameters_.Bind<bool>(PARAM_UPDATE_MORTALITY_PARAMETERS, &update_natural_mortality_parameters_, "If an agent/individual moves do you want it to take on the natural mortality parameters of the new spatial cell", "");
  parameters_.Bind<float>(PARAM_TIME_STEP_RATIO, &ratios_, "Time step ratios for the mortality rates to apply in each time step. See manual for how this is applied", "");

  RegisterAsAddressable(PARAM_M, &m_);

}

/**
 * DoBuild
 */
void MortalityConstantRate::DoBuild() {
  // Get the layers
  if (m_layer_label_ != "") {
    m_layer_ = model_->managers().layer()->GetNumericLayer(m_layer_label_);
    if (!m_layer_) {
      LOG_ERROR_P(PARAM_M_MULTIPLIER_LAYER_LABEL) << "could not find the layer '" << m_layer_label_ << "', please make sure it exists and that it is type 'numeric'";
    }
    // Do the multiplication
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        float multiplier = m_layer_->get_value(row, col);
        LOG_FINEST() << "multiplier = " << multiplier << " m value = " << m_ <<  " check we set the right value = " << m_layer_->get_value(row, col, multiplier * m_);
      }
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

  /**
   * Organise our time step ratios. Each time step can
   * apply a different ratio of M so here we want to verify
   * we have enough and re-scale them to 1.0
   */
  vector<TimeStep*> time_steps = model_->managers().time_step()->ordered_time_steps();
  LOG_FINEST() << "time_steps.size(): " << time_steps.size();
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
      LOG_ERROR_P(PARAM_TIME_STEP_RATIO) << " length (" << ratios_.size()
          << ") does not match the number of time steps this process has been assigned to (" << active_time_steps.size() << ")";

    for (float value : ratios_) {
      if (value < 0.0 || value > 1.0)
        LOG_ERROR_P(PARAM_TIME_STEP_RATIO) << " value (" << value << ") must be between 0.0 (exclusive) and 1.0 (inclusive)";
    }

    for (unsigned i = 0; i < ratios_.size(); ++i)
      time_step_ratios_[active_time_steps[i]] = ratios_[i];
  }


  model_->set_m(m_);

  cell_offset_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cell_offset_[i].resize(model_->get_width());
  }

}

/**
 * DoExecute
 */
void MortalityConstantRate::DoExecute() {
  LOG_MEDIUM();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
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
  random_numbers_.resize(n_agents_);
  for (unsigned i = 0; i < n_agents_; ++i)
    random_numbers_[i] = rng.chance();

  LOG_FINE() << "number of agents = " << n_agents_;
  unsigned agents_removed = 0;
  if (selectivity_length_based_) {
    //#pragma omp parallel for collapse(2)
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        // get the ratio to apply first
        unsigned time_step = model_->managers().time_step()->current_time_step();
        float ratio = time_step_ratios_[time_step];
        WorldCell* cell = world_->get_base_square(row, col);

        if (cell->is_enabled()) {
          unsigned counter = 0;
          unsigned initial_size = cell->agents_.size();
          LOG_FINEST() << initial_size << " initial agents";
          for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++counter,++iter) {
            if ((*iter).is_alive()) {
              //LOG_FINEST() << "selectivity = " << selectivity_at_age << " m = " << (*iter).get_m();
              if (random_numbers_[cell_offset_[row][col] + counter] <= (1 - std::exp(- ratio * (*iter).get_m() * selectivity_[(*iter).get_sex()]->GetResult((*iter).get_length_bin_index())))) {
                (*iter).dies();
                initial_size--;
                agents_removed++;
              }
            }
          }
          LOG_FINEST() << initial_size << " after mortality";
        }
      }
    }
  } else {
    //#pragma omp parallel for collapse(2)
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          // need thread safe access to rng
          unsigned time_step = model_->managers().time_step()->current_time_step();
          float ratio = time_step_ratios_[time_step];
          WorldCell* cell = world_->get_base_square(row, col);
          unsigned counter = 0;
          unsigned initial_size = cell->agents_.size();
          LOG_FINEST() << initial_size << " initial agents, ratio = " << ratio;
          for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++counter,++iter) {
            if ((*iter).is_alive()) {
              //LOG_FINEST() << "rand number = " << random_numbers_[cell_offset_[row][col] + counter] << " selectivity = " << cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()] << " survivorship = " << (1 - std::exp(-(*iter).get_m() *  cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()])) << " M = " << (*iter).get_m();
              if (random_numbers_[cell_offset_[row][col] + counter] <= (1.0 - std::exp(- ratio * (*iter).get_m() * selectivity_[(*iter).get_sex()]->GetResult((*iter).get_age_index())))) {
                (*iter).dies();
                initial_size--;
                agents_removed++;
              }
            }
          }
          LOG_FINEST() << initial_size << " after mortality";
        }
      }
    }
  }
  if (model_->state() != State::kInitialise)
    removals_by_year_[model_->current_year()] += agents_removed;
}


/*
 * This method is called at when ever an agent is created/seeded or moves. Agents will get a new/updated growth parameters
 * based on the spatial cells of the process. This is called in initialisation/Recruitment and movement processes if needed.
*/
void  MortalityConstantRate::draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  float mean_m;
  if (m_layer_)
    mean_m = m_layer_->get_value(row, col) * m_;
  else
    mean_m = m_;

  LOG_FINE() << "mean M = " << mean_m;
  vector.clear();
  vector.resize(number_of_draws);

  //LOG_FINEST() << "mean M = " << mean_m;
  if (distribution_ == PARAM_NORMAL) {
    for (unsigned i = 0; i < number_of_draws; ++i) {
      float value = rng.normal(mean_m, cv_ * mean_m);
      //LOG_FINEST() << "value of M = " << value;
      vector[i] = value;
    }
  } else {
    for (unsigned i = 0; i < number_of_draws; ++i) {
      float value = rng.lognormal(mean_m, cv_);
      //LOG_FINEST() << "value of M = " << value;
      vector[i] = value;
    }
  }
}

// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void  MortalityConstantRate::FillReportCache(ostringstream& cache) {
  cache << "year: ";
  for (auto& year : removals_by_year_)
    cache << year.first << " ";

  cache << "\nagents_removed: ";
  for (auto& year : removals_by_year_)
    cache << year.second << " ";
  cache << "\n";
}

// A Method for telling the world we need to redistribute Mortality parmaeters
void MortalityConstantRate::RebuildCache() {
  LOG_FINE();
  world_->rebuild_mort_params();

}

} /* namespace processes */
} /* namespace niwa */
