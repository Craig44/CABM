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
  parameters_.Bind<double>(PARAM_LAMBDA, &lambda_, "The maximum value of the absolute sum of differences (lambda) between the partition at year-1 and year that indicates successfull convergence", "", Double(0.0));
  parameters_.Bind<unsigned>(PARAM_NUMBER_OF_AGENTS, &number_agents_, "The number of agents to initially seed in the partition", "");
  parameters_.Bind<string>(PARAM_LAYER_LABEL, &intial_layer_label_, "The label of a layer that you want to seed a distribution by.", "", "");
  parameters_.Bind<string>(PARAM_GROWTH_PROCESS_LABEL, &growth_process_label_, "Label for the growth process in the annual cycle", "");
  parameters_.Bind<string>(PARAM_NATURAL_MORTALITY_PROCESS_LABEL, &natural_mortality_label_, "Label for the natural mortality process in the annual cycle", "");
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

  // Link to growth and mortality processes
  processes::Growth* growth = model_->managers().process()->GetGrowthProcess(growth_process_label_);
  if (!growth) {
    LOG_FATAL_P(PARAM_GROWTH_PROCESS_LABEL) << "could not find this process does it exist?";
  }

  processes::Mortality* mortality = model_->managers().process()->GetMortalityProcess(natural_mortality_label_);
  if (!mortality) {
    LOG_FATAL_P(PARAM_NATURAL_MORTALITY_PROCESS_LABEL) << "could not find this process does it exist? if it does exist, is of type mortality? because it needs to be";
  }
  unsigned cells = world_->get_enabled_cells();
  unsigned agents_per_cell = (int)number_agents_ / (int)cells;
  LOG_FINEST() << "number of cells = " << cells << " number of agents per cell = " << agents_per_cell;
  vector<double> Ms;
  vector<vector<double>> growth_pars;
  double seed_z = model_->get_initial_seed_z();

  LOG_FINEST() << "check seed = " << seed_z;
  // Seed some individuals in the world
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      mortality->draw_rate_param(row, col, agents_per_cell, Ms);
      growth->draw_growth_param(row, col, agents_per_cell, growth_pars);
      LOG_FINEST() << "cell = " << row + 1 << " col = " << col + 1 << " size of Ms = " << Ms.size() << " number of seeds = " << agents_per_cell;

      if (cell->is_enabled()) {
        cell->seed_agents(agents_per_cell, Ms, growth_pars, seed_z);
        LOG_FINEST() << "row " << row + 1 << " col = " << col + 1 << " seeded " << cell->get_agents().size();
      }
    }
  }

  // Ask to print the intialisation phase here;
  model_->managers().report()->Execute(State::kInitialise);


  time_steps_ = model_->managers().time_step()->ordered_time_steps();


  // Check convergence years are in strictly increasing order
  if (convergence_years_.size() != 0) {
    std::sort(convergence_years_.begin(), convergence_years_.end());
    if ((*convergence_years_.rbegin()) != years_)
      convergence_years_.push_back(years_);
  }
}

/**
 * Execute our iterative initialisation phases.
 */
void Iterative::Execute() {
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
  Double variance = 0.0;

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
