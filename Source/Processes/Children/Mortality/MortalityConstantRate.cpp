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
// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
MortalityConstantRate::MortalityConstantRate(Model* model) : Mortality(model) {
  parameters_.Bind<string>(PARAM_DISTRIBUTION, &distribution_, "the distribution to allocate the parameters to the agents", "");
  parameters_.Bind<float>(PARAM_CV, &cv_, "The cv of the distribution", "");
  parameters_.Bind<string>(PARAM_M_MULTIPLIER_LAYER_LABEL, &m_layer_label_, "Label for the numeric layer that describes a multiplier of M through space", "", ""); // TODO perhaps as a multiplier, 1.2 * 0.2 = 0.24
  parameters_.Bind<string>(PARAM_SELECTIVITY_LABEL, &selectivity_label_, "Label for the selectivity block", "");
  parameters_.Bind<float>(PARAM_M, &m_, "Natural mortality for the model", "");
  parameters_.Bind<bool>(PARAM_UPDATE_MORTALITY_PARAMETERS, &update_natural_mortality_parameters_, "If an agent/individual moves do you want it to take on the natural mortality parameters of the new spatial cell", "");

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
        LOG_FINEST() << "multiplier = " << multiplier << " m value = " << m_;
        m_layer_->set_value(row, col, multiplier * m_);
        LOG_FINEST() << "check we set the right value = " << m_layer_->get_value(row, col);
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
  model_->set_m(m_);

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
void MortalityConstantRate::DoExecute() {
  LOG_TRACE();
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

  unsigned agents_removed = 0;
  if (selectivity_length_based_) {
    #pragma omp parallel for collapse(2)
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          unsigned counter = 0;
          unsigned initial_size = cell->agents_.size();
          LOG_FINEST() << initial_size << " initial agents";
          for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++counter) {
            //LOG_FINEST() << "selectivity = " << selectivity_at_age << " m = " << (*iter).get_m();
            if (random_numbers_[cell_offset_[row][col] + counter] <= (1 - std::exp(-(*iter).get_m() * cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->length_bins().size() + (*iter).get_length_bin_index()]))) {
              iter = cell->agents_.erase(iter);
              initial_size--;
              agents_removed++;
            } else
              ++iter;
          }
          LOG_FINEST() << initial_size << " after mortality";
        }
      }
    }
  } else {
    #pragma omp parallel for collapse(2)
    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        WorldCell* cell = world_->get_base_square(row, col);
        if (cell->is_enabled()) {
          // need thread safe access to rng
          unsigned counter = 0;
          unsigned initial_size = cell->agents_.size();
          LOG_FINEST() << initial_size << " initial agents";
          for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++counter) {
            //LOG_FINEST() << "selectivity = " << selectivity_at_age << " m = " << (*iter).get_m();
            if (random_numbers_[cell_offset_[row][col] + counter] <= (1 - std::exp(-(*iter).get_m() *  cell_offset_for_selectivity_[row][col][(*iter).get_sex() * model_->age_spread() + (*iter).get_age_index()]))) {
              iter = cell->agents_.erase(iter);
              initial_size--;
              agents_removed++;
            } else
              ++iter;
          }
          LOG_FINEST() << initial_size << " after mortality";
        }
      }
    }
  }
  if (model_->state() != State::kInitialise)
    removals_by_year_[model_->current_year()] = agents_removed;
}


/*
 * This method is called at when ever an agent is created/seeded or moves. Agents will get a new/updated growth parameters
 * based on the spatial cells of the process. This is called in initialisation/Recruitment and movement processes if needed.
*/
void  MortalityConstantRate::draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  float mean_m;
  if (m_layer_)
    mean_m = m_layer_->get_value(row, col);
  else
    mean_m = m_;

  vector.clear();
  vector.resize(number_of_draws);

  //LOG_FINEST() << "mean M = " << mean_m;
  for (unsigned i = 0; i < number_of_draws; ++i) {
    float value = rng.lognormal(mean_m, cv_);
    //LOG_FINEST() << "value of M = " << value;
    vector[i] = value;
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

} /* namespace processes */
} /* namespace niwa */
