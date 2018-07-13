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
  parameters_.Bind<string>(PARAM_LAYER_LABEL, &intial_layer_label_, "The label of a layer that you want to seed a distribution by.", "")->is_optional();
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

  unsigned cells = world_->get_enabled_cells();
  unsigned agents_per_cell = (int)number_agents_ / (int)cells;
  LOG_FINEST() << "number of cells = " << cells << " number of agents per cell = " << agents_per_cell;
  // Seed some individuals in the world
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        cell->seed_agents(agents_per_cell);
        LOG_FINEST() << "row " << row + 1 << " col = " << col + 1 << " seeded " << cell->get_agents().size();
      }
    }
  }

  time_steps_ = model_->managers().time_step()->ordered_time_steps();
/*

  // handle any new processes we want to insert
  for (string insert : insert_processes_) {
    vector<string> pieces;
    boost::split(pieces, insert, boost::is_any_of("()="), boost::token_compress_on);

    string target_process   = pieces.size() == 3 ? pieces[1] : "";
    string new_process      = pieces.size() == 3 ? pieces[2] : pieces[1];

    auto time_step = model_->managers().time_step()->GetTimeStep(pieces[0]);
    vector<string> process_labels = time_step->initialisation_process_labels(label_);

    if (target_process == "") {
      process_labels.push_back(new_process);
    } else {
      vector<string>::iterator iter = std::find(process_labels.begin(), process_labels.end(), target_process);
      if (iter == process_labels.end())
        LOG_ERROR_P(PARAM_INSERT_PROCESSES) << " process " << target_process << " does not exist in time step " << time_step->label();
      process_labels.insert(iter, new_process);
    }

    time_step->SetInitialisationProcessLabels(label_, process_labels);
  }

  // handle the excludes we've specified
  for (string exclude : exclude_processes_) {
    unsigned count = 0;
    for (auto time_step : time_steps_) {
      vector<string> process_labels = time_step->initialisation_process_labels(label_);
      unsigned size_before = process_labels.size();
      process_labels.erase(std::remove_if(process_labels.begin(), process_labels.end(), [exclude](string& ex) { return exclude == ex; }), process_labels.end());
      unsigned diff = size_before - process_labels.size();

      time_step->SetInitialisationProcessLabels(label_, process_labels);
      count += diff;
    }

    if (count == 0)
      LOG_ERROR_P(PARAM_EXCLUDE_PROCESSES) << " process " << exclude << " does not exist in any time steps to be excluded. Please check your spelling";
  }

  if (convergence_years_.size() != 0) {
    std::sort(convergence_years_.begin(), convergence_years_.end());
    if ((*convergence_years_.rbegin()) != years_)
      convergence_years_.push_back(years_);
  }


  // Build our partition
  vector<string> categories = model_->categories()->category_names();
  partition_.Init(categories);
  cached_partition_.Init(categories);

  // Find any BH_recruitment process in the annual cycle
  unsigned i = 0;
  for (auto time_step : model_->managers().time_step()->ordered_time_steps()) {
    for (auto process : time_step->processes()) {
      if (process->process_type() == ProcessType::kRecruitment && process->type() == PARAM_RECRUITMENT_BEVERTON_HOLT) {
        LOG_FINEST() << "Found a BH process!!!!";
        recruitment_process_.push_back(dynamic_cast<RecruitmentBevertonHolt*>(process));
        if (!recruitment_process_[i])
          LOG_CODE_ERROR() << "BevertonHolt Recruitment exists but dynamic cast pointer cannot be made, if (!recruitment) ";
        i++;
      } else if (process->process_type() == ProcessType::kRecruitment && process->type() == PARAM_RECRUITMENT_BEVERTON_HOLT_WITH_DEVIATIONS) {
        LOG_FINEST() << "Found a BH process!!!!";
        recruitment_process_with_devs_.push_back(dynamic_cast<RecruitmentBevertonHoltWithDeviations*>(process));
        if (!recruitment_process_with_devs_[i])
          LOG_CODE_ERROR() << "BevertonHolt Recruitment with deviations exists but dynamic cast pointer cannot be made, if (!recruitment) ";
        i++;
      }
    }
  }
  */
}

/**
 * Execute our iterative initialisation phases.
 */
void Iterative::Execute() {
/*  if (convergence_years_.size() == 0) {
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

      cached_partition_.BuildCache();
      time_step_manager.ExecuteInitialisation(label_, 1);
      ++total_years;

      if (CheckConvergence()) {
        LOG_FINEST() << " year Convergence was reached = " << years;
        break;
      }
      LOG_FINEST() << "Initial year = " << years;
    }
  }

  LOG_FINEST() << "Number of Beverton-Holt recruitment processes in annual cycle = " << recruitment_process_.size();
  LOG_FINEST() << "Number of Beverton-Holt recruitment processes with deviations in annual cycle = " << recruitment_process_with_devs_.size();
  // We are at Equilibrium state here
  // Check if we have B0 initialised or R0 initialised recruitment
  bool B0_intial_recruitment = false;
  for (auto recruitment_process : recruitment_process_) {
    if (recruitment_process->bo_initialised()) {
      LOG_FINEST() << PARAM_B0 << " has been defined for process labelled " << recruitment_process->label();
      recruitment_process->ScalePartition();
      B0_intial_recruitment = true;
    }
  }
  for (auto recruitment_process_with_devs : recruitment_process_with_devs_) {
    if (recruitment_process_with_devs->bo_initialised()) {
      LOG_FINEST() << PARAM_B0 << " has been defined for process labelled " << recruitment_process_with_devs->label();
      recruitment_process_with_devs->ScalePartition();
      B0_intial_recruitment = true;
    }
  }
  if (B0_intial_recruitment) {
    // Calculate derived quanitities in the right space if we have a B0 initialised model
    timesteps::Manager& time_step_manager = *model_->managers().time_step();
    time_step_manager.ExecuteInitialisation(label_, 1);
  }*/
}

/**
 * Check for convergence on our partition and return true if it exceeds the
 * lambda threshold so we can quit early and save time.
 *
 * @return True if convergence, False otherwise
 */
bool Iterative::CheckConvergence() {
  LOG_TRACE();
/*

  Double variance = 0.0;

  auto cached_category = cached_partition_.begin();
  auto category = partition_.begin();

  for (; category != partition_.end(); ++cached_category, ++category) {
    Double sum = 0.0;
    for (Double value : (*category)->data_)
      sum += value; // Can't use std::accum because of AutoDiff
    if (sum == 0.0)
      return false;

    for (unsigned i = 0; i < (*category)->data_.size(); ++i)
      variance += fabs(cached_category->data_[i] - (*category)->data_[i]) / sum;
  }

  if (variance < lambda_)
    return true;
*/

  return false;
}


} /* namespace initialisationphases */
} /* namespace niwa */
