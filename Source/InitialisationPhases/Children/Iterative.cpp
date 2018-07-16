/**
 * @file Iterative.cpp
 * @author  C.Marsh
 * @date 12/07/2018
 * @section LICENSE
 *
 *
 */
// headers
#include "Iterative.h"

#include <algorithm>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>


#include "Processes/Manager.h"
#include "TimeSteps/Factory.h"
#include "TimeSteps/Manager.h"
#include "Layers/Manager.h"
#include "Processes/Manager.h"
#include "Reports/Manager.h"
#include "Model/Model.h"
#include "World/WorldCell.h"
// namespaces
namespace niwa {
namespace initialisationphases {

/**
 * Default constructor
 */
Iterative::Iterative(Model* model)
  : InitialisationPhase(model) {
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "The number of iterations (years) over which to execute this initialisation phase", "");
  parameters_.Bind<unsigned>(PARAM_CONVERGENCE_YEARS, &convergence_years_, "The iteration (year) when the test for converegence (lambda) is evaluated", "", true);
  parameters_.Bind<float>(PARAM_LAMBDA, &lambda_, "The maximum value of the absolute sum of differences (lambda) between the partition at year-1 and year that indicates successfull convergence", "", float(0.0));
  parameters_.Bind<unsigned>(PARAM_NUMBER_OF_AGENTS, &number_agents_, "The number of agents to initially seed in the partition", "");
  parameters_.Bind<string>(PARAM_LAYER_LABEL, &intial_layer_label_, "The label of a layer that you want to seed a distribution by.", "", "");
}

/**
 *
 */
void Iterative::DoValidate() {
  LOG_TRACE();
}

/**
 *
 */
void Iterative::DoBuild() {
  LOG_TRACE();

  // Did the user supply a layer
  if (intial_layer_label_ != "") {
    initial_layer_ = model_->managers().layer()->GetNumericLayer(intial_layer_label_);
    if (!initial_layer_)
      LOG_ERROR_P(PARAM_LAYER_LABEL) << "could not find layer does it exist? if it does exist, is of type numeric? because it needs to be";
  }

  // Generate a pointer to the world
  world_ = model_->world_view();
  if (!world_)
    LOG_CODE_ERROR() << "could not create the world via reference, something is wrong";

  time_steps_ = model_->managers().time_step()->ordered_time_steps();


  // Check convergence years are in strictly increasing order
  if (convergence_years_.size() != 0) {
    std::sort(convergence_years_.begin(), convergence_years_.end());
    if ((*convergence_years_.rbegin()) != years_)
      convergence_years_.push_back(years_);
  }

  model_->set_n_agents(number_agents_);
}

/**
 * Execute our iterative initialisation phases.
 */
void Iterative::Execute() {
  LOG_TRACE();

  // Check we get a sensible estimate of M
  if (model_->get_m() <= 0.0)
    LOG_FATAL() << "Could not get a sensible M value from the mortality process defined this is unusual you want to check everything is okay in the mortality process or contact a developer";

  // Calculate the R0 values for the recruitment processes this is a bit crude but will do for now
  unsigned large_age = model_->max_age() * 4;
  float number = 0;
  float m = model_->get_m();
  for (unsigned age = 0; age < large_age; ++age)
    number += exp(- m * age);

  vector<string> recruitment_labels;
  float total_b0 = 0.0;
  for (auto iter : model_->get_b0s()) {
    recruitment_labels.push_back(iter.first);
    total_b0 += iter.second;
  }

  for (auto iter : model_->get_b0s()) {
    unsigned value = number_agents_ / (unsigned)number * (unsigned)iter.second / (unsigned)total_b0;
    model_->set_r0(iter.first, value);
  }


  // Move on and seed our n_agents
  unsigned cells = world_->get_enabled_cells();
  unsigned agents_per_cell = (int)number_agents_ / (int)cells;
  LOG_FINEST() << "number of cells = " << cells << " number of agents per cell = " << agents_per_cell;

  float seed_z = model_->get_initial_seed_z();

  LOG_FINEST() << "check seed = " << seed_z;
  // Seed some individuals in the world
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        cell->seed_agents(agents_per_cell, seed_z);
        LOG_FINEST() << "row " << row + 1 << " col = " << col + 1 << " seeded " << cell->get_agents().size();
      }
    }
  }

  // Ask to print the intialisation phase here;
  model_->managers().report()->Execute(State::kInitialise);



  if (convergence_years_.size() == 0) {
    timesteps::Manager& time_step_manager = *model_->managers().time_step();
    time_step_manager.ExecuteInitialisation(label_, years_);
  } else {
    unsigned total_years = 0;
    for (unsigned years : convergence_years_) {
      timesteps::Manager& time_step_manager = *model_->managers().time_step();
      time_step_manager.ExecuteInitialisation(label_, years - (total_years + 1));
      total_years += years - (total_years + 1);
      if ((total_years + 1) >= years_) {
        time_step_manager.ExecuteInitialisation(label_, 1);
        break;
      }
      vector<unsigned> cached_partition;
      world_->get_world_age_frequency(cached_partition);
      time_step_manager.ExecuteInitialisation(label_, 1);
      ++total_years;

      vector<unsigned> current_partition;
      world_->get_world_age_frequency(current_partition);
      if (CheckConvergence(cached_partition, current_partition)) {
        LOG_FINEST() << " year Convergence was reached = " << years;
        break;
      }
      LOG_FINEST() << "Initial year = " << years;
    }
  }
}

/**
 * Check for convergence on our partition and return true if it exceeds the
 * lambda threshold so we can quit early and save time.
 * @param pre_quantity a vector of partition attributes before we check
 * @param post_quantity a vector of partition attributes at the point of the check
 *
 * @return True if convergence, False otherwise
 */
template<typename T>
bool Iterative::CheckConvergence(vector<T> pre_quantity, vector<T> post_quantity) {
  LOG_TRACE();
  float variance = 0.0;

  auto pre_iter = pre_quantity.begin();
  auto post_iter = post_quantity.begin();

  T sum = 0;
  for (auto iter = post_quantity.begin(); iter != post_quantity.end(); ++iter)
    sum += *iter;

  LOG_FINEST() << "sum = " << sum;
  for (; pre_iter != pre_quantity.end(); ++pre_iter, ++post_iter) {
    variance += fabs((*pre_iter) - (*post_iter)) / sum;
  }

  LOG_FINEST() << "variance = " << variance << " lambda = " << lambda_;

  if (variance < lambda_)
    return true;


  return false;
}


} /* namespace initialisationphases */
} /* namespace niwa */
