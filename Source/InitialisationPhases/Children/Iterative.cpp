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

  LOG_FINEST() << "total b0 = " << total_b0 << " number = " << number << " number of agents = " << number_agents_;

  for (auto iter : model_->get_b0s()) {
    unsigned value = (unsigned)(number_agents_ / number * iter.second / total_b0);
    LOG_FINEST() << "setting R0 for " << iter.first << " = " << value;
    model_->set_r0(iter.first, value);
  }


  // Move on and seed our n_agents
  unsigned cells = world_->get_enabled_cells();
  unsigned agents_per_cell = (int)number_agents_ / (int)cells;
  LOG_FINEST() << "number of cells = " << cells << " number of agents per cell = " << agents_per_cell;

  float seed_z = model_->get_m();
  // This is important to set, otherwise age distribution will get all weird
  unsigned init_year = model_->start_year() - years_ - 1;
  model_->set_current_year_in_initialisation(init_year);

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

  // Run out the initialisation phase
  timesteps::Manager& time_step_manager = *model_->managers().time_step();
  model_->set_current_year_in_initialisation(model_->start_year() - years_);

  time_step_manager.ExecuteInitialisation(label_, years_);
  // I removed the convergence check because we seed agents by birth_year so we really need to have an incremental inital years to have a sensible age distribution
  // if we cut at 50 years when we thought we might be running for 100 years then this would cause an improper age distribution, So I leave it for the user to define
  // The correct burin in time to reach equilibrium

  float ssb = 0.0, b0 = 0.0, scalar;
  for (auto iter : model_->get_b0s()) {
    b0 += iter.second;
    ssb += model_->get_ssb(iter.first);
  }
  scalar = b0 / ssb;
  model_->set_scalar(scalar);
  LOG_FINEST() << "scalar = " << scalar;
  // Set scalar before continuing on
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        auto& agents = cell->get_agents();
        for (auto iter = agents.begin(); iter != agents.end(); ++iter) {
          (*iter).set_scalar(scalar);
        }
      }
    }
  }
}


} /* namespace initialisationphases */
} /* namespace niwa */
