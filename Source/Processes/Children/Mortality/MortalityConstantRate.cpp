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
  parameters_.Bind<double>(PARAM_CV, &cv_, "The cv of the distribution", "");
  parameters_.Bind<string>(PARAM_M_MULTIPLIER_LAYER_LABEL, &m_layer_label_, "Label for the numeric layer that describes a multiplier of M through space", "", ""); // TODO perhaps as a multiplier, 1.2 * 0.2 = 0.24
  parameters_.Bind<string>(PARAM_SELECTIVITY_LABEL, &selectivity_label_, "Label for the selectivity block", "");
  parameters_.Bind<double>(PARAM_M, &m_, "Natural mortality for the model", "");
  parameters_.Bind<bool>(PARAM_UPDATE_MORTALITY_PARAMETERS, &update_natural_mortality_parameters_, "If an agent/individual moves do you want it to take on the natural mortality parameters of the new spatial cell", "");
  parameters_.Bind<double>(PARAM_TIME_STEP_RATIO, &ratios_, "Time step ratios for the mortality rates to apply in each time step. See manual for how this is applied", "");

  RegisterAsAddressable(PARAM_M, &m_, addressable::kAll, addressable::kno);

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
        double multiplier = m_layer_->get_value(row, col);
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
  LOG_FINEST() << "time_steps.size(): " << time_steps.size();
  vector<unsigned> active_time_steps;
  for (unsigned i = 0; i < time_steps.size(); ++i) {
    if (time_steps[i]->HasProcess(label_))
      active_time_steps.push_back(i);
  }

  LOG_FINE() << "active time-steps = " << active_time_steps.size() << " number of ratios = " << ratios_.size();
  if (ratios_.size() == 0) {
    for (unsigned i : active_time_steps)
      time_step_ratios_[i] = 1.0;
  } else {
    if (ratios_.size() != active_time_steps.size())
      LOG_ERROR_P(PARAM_TIME_STEP_RATIO) << " length (" << ratios_.size()
          << ") does not match the number of time steps this process has been assigned to (" << active_time_steps.size() << ")";

    for (double value : ratios_) {
      if (value < 0.0 || value > 1.0)
        LOG_ERROR_P(PARAM_TIME_STEP_RATIO) << " value (" << value << ") must be between 0.0 (exclusive) and 1.0 (inclusive)";
    }

    for (unsigned i = 0; i < ratios_.size(); ++i) {
      LOG_FINE() << "timestep = " << i << " value = " << ratios_[i] << " timestep ndx = " << active_time_steps[i];
      time_step_ratios_[active_time_steps[i]] = ratios_[i];
    }
  }

  LOG_FINE() << "setting value of m in the model = " << m_;
  model_->set_m(m_);

  cell_offset_.resize(model_->get_height());
  for (unsigned i = 0; i < model_->get_height(); ++i) {
    cell_offset_[i].resize(model_->get_width());
  }

}

/**
 * DoReset
 */
void MortalityConstantRate::DoReset() {
  LOG_FINE() << "clearing containers";
  removals_by_age_and_area_.clear();
  removals_by_length_and_area_.clear();
  removals_census_.clear();
  removals_tag_recapture_.clear();
}

/**
 * The main function of this class. pulled out of DoExcute so that I can apply it to many different vectors.
 */
void MortalityConstantRate::ApplyStochasticMortality(vector<Agent>& agents) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  unsigned time_step = model_->managers().time_step()->current_time_step();
  double ratio = time_step_ratios_[time_step];

  LOG_FINE() << "year = " << model_->current_year() << " time step = " << time_step << " ratio = " << ratio;
  if (selectivity_length_based_) {
    for(auto& agent : agents) {
      if (agent.is_alive()) {
        //LOG_FINEST() << "selectivity = " << selectivity_at_age << " m = " << (*iter).get_m();
        if (rng.chance() <= (1 - std::exp(- ratio * agent.get_m() * selectivity_[agent.get_sex()]->GetResult(agent.get_length_bin_index())))) {
          agent.dies();
          agents_removed_++;
        }
      }
    }
  } else {
    for(auto& agent : agents) {
      if (agent.is_alive()) {
        if (rng.chance() <= (1 - std::exp(- ratio * agent.get_m() * selectivity_[agent.get_sex()]->GetResult(agent.get_age_index())))) {
          agent.dies();
          agents_removed_++;
        }
      }
    }
  }
}
/**
 * DoExecute
 */
void MortalityConstantRate::DoExecute() {
  LOG_MEDIUM() << label_;
  agents_removed_ = 0;
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      // get the ratio to apply first
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        // Apply mortality to elements in a cell.
        ApplyStochasticMortality(cell->agents_);
        ApplyStochasticMortality(cell->tagged_agents_);
      }
    }
  }
  if (model_->state() != State::kInitialise)
    removals_by_year_[model_->current_year()] += agents_removed_;
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

  //LOG_FINE() << "mean M = " << mean_m;
  vector.clear();
  vector.resize(number_of_draws);

  //LOG_FINEST() << "mean M = " << mean_m;
  float value = 0.0;
  if (distribution_ == PARAM_NORMAL) {
    for (unsigned i = 0; i < number_of_draws; ++i) {
      value = rng.normal(mean_m, cv_ * mean_m);
      //LOG_FINEST() << "value of M = " << value;
      vector[i] = value;
    }
  } else {
    for (unsigned i = 0; i < number_of_draws; ++i) {
      value = rng.lognormal(mean_m, cv_);
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
