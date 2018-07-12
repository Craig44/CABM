/**
 * @file Model.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 6/12/2012
 * @section LICENSE
 *
 * Copyright NIWA Science ©2012 - www.niwa.co.nz
 *
 * @modified 12/7/2018 by C.Marsh for IBM usage
 */

// Headers
#include "Model.h"

#include <iostream>
#include <thread>
#include <chrono>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>

#include "Factory.h"
#include "Managers.h"
#include "Objects.h"
#include "Agents/Agent.h"
#include "GlobalConfiguration/GlobalConfiguration.h"
#include "InitialisationPhases/Manager.h"
#include "Logging/Logging.h"
#include "Observations/Manager.h"
#include "World/WorldView.h"
#include "Reports/Manager.h"
#include "TimeSteps/Manager.h"
#include "TimeVarying/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/To.h"

// Namespaces
namespace niwa {

using std::cout;
using std::endl;

/**
 * Default Constructor
 */
Model::Model() {
  LOG_TRACE();
  parameters_.Bind<unsigned>(PARAM_START_YEAR, &start_year_, "Define the first year of the model, immediately following initialisation", R"(Defines the first year of the model, $\ge 1$, e.g. 1990)");
  parameters_.Bind<unsigned>(PARAM_FINAL_YEAR, &final_year_, "Define the final year of the model, excluding years in the projection period", "Defines the last year of the model, i.e., the model is run from start_year to final_year");
  parameters_.Bind<unsigned>(PARAM_MIN_AGE, &min_age_, "Minimum age of individuals in the population", R"($0 \le$ age\textlow{min} $\le$ age\textlow{max})", 0);
  parameters_.Bind<unsigned>(PARAM_MAX_AGE, &max_age_, "Maximum age of individuals in the population", R"($0 \le$ age\textlow{min} $\le$ age\textlow{max})", 0);
  parameters_.Bind<bool>(PARAM_AGE_PLUS, &age_plus_, "Define the oldest age or extra length midpoint (plus group size) as a plus group", "true, false", false);
  parameters_.Bind<string>(PARAM_INITIALISATION_PHASE_LABELS, &initialisation_phases_, "Define the labels of the phases of the initialisation", R"(A list of valid labels defined by \texttt{@initialisation_phase})", true);
  parameters_.Bind<string>(PARAM_TIME_STEPS, &time_steps_, "Define the labels of the time steps, in the order that they are applied, to form the annual cycle", R"(A list of valid labels defined by \texttt{@time_step})");
  parameters_.Bind<unsigned>(PARAM_LENGTH_BINS, &length_bins_, "", "", true);
  parameters_.Bind<bool>(PARAM_LENGTH_PLUS, &length_plus_, "Is the last bin a plus group", "", false);
  parameters_.Bind<string>(PARAM_BASE_LAYER_LABEL, &base_layer_, "Label for the base layer", "");
  parameters_.Bind<unsigned>(PARAM_NROWS, &world_height_, "number of rows in spatial domain", "");
  parameters_.Bind<unsigned>(PARAM_NCOLS, &world_width_, "number of columns in spatial domain", "");



  global_configuration_ = new GlobalConfiguration();
  managers_ = new Managers(this);
  objects_ = new Objects(this);
  factory_ = new Factory(this);
  world_view_ = new WorldView(this);
}

/**
 * Destructor
 */
Model::~Model() {
  delete global_configuration_;
  delete managers_;
  delete objects_;
  delete factory_;
  delete world_view_;
}

/**
 *
 */
vector<unsigned> Model::years() const {
  vector<unsigned> years;
  unsigned year;
  for (year = start_year_; year <= final_year_; ++year)
    years.push_back(year);

  return years;
}

/**
 *
 */
unsigned Model::year_spread() const {
  return (final_year_ - start_year_) + 1;
}

/**
 *
 */
Managers& Model::managers() {
  return *managers_;
}

Objects& Model::objects() {
  return *objects_;
}

Factory& Model::factory() {
  return *factory_;
}

WorldView* Model::world_view() {
  return world_view_;
}

/**
 * Start our model. This is the entry point method for the model
 * after being called from the main() method.
 *
 * This method will start stepping through the states and verifying
 * each step.
 */
bool Model::Start(RunMode::Type run_mode) {
  LOG_TRACE();
  Logging& logging = Logging::Instance();
  if (logging.errors().size() > 0) {
    logging.FlushErrors();
    return false;
  }

  run_mode_ = run_mode;

  if (state_ != State::kStartUp)
    LOG_CODE_ERROR() << "Model state should always be startup when entering the start method";

  managers_->report()->Execute(State::kStartUp);

  LOG_FINE() << "Model: State Change to Validate";
  state_ = State::kValidate;
  Validate();
  if (logging.errors().size() > 0) {
    logging.FlushErrors();
    return false;
  }

  managers_->report()->Execute(state_);

  LOG_FINE() << "Model: State Change to Build";
  state_ = State::kBuild;
  Build();
  if (logging.errors().size() > 0) {
    logging.FlushErrors();
    return false;
  }
  managers_->report()->Execute(state_);

  LOG_FINE() << "Model: State Change to Verify";
  state_ = State::kVerify;
  Verify();
  if (logging.errors().size() > 0) {
    logging.FlushErrors();
    return false;
  }
  managers_->report()->Execute(state_);

  // prepare all reports
  LOG_FINE() << "Preparing Reports";
  managers_->report()->Prepare();

  switch(run_mode_) {
  case RunMode::kBasic:
    RunBasic();
    break;

  default:
    LOG_ERROR() << "Invalid run mode has been specified. This run mode is not supported: " << run_mode_;
    break;
  }

  // finalise all reports
  LOG_FINE() << "Finalising Reports";
  state_ = State::kFinalise;
  for (auto executor : executors_[state_])
    executor->Execute();
  managers_->report()->Execute(state_);
  managers_->report()->Finalise();
  return true;
}

/**
 * Populate the loaded parameters
 */
void Model::PopulateParameters() {
  LOG_TRACE();

  // Check that we've actually defined a @model block
  if (block_type_ == "")
    LOG_ERROR() << "The @model block is missing from configuration file. This block is required.";
  /**
   * Validate the parameters
   */
  parameters_.Populate(this);
}

/**
 * First we will do the local validations. Then we will call validation on the other objects
 */
void Model::Validate() {
  // Check that we've actually defined a @model block
  if (block_type_ == "")
    LOG_ERROR() << "The @model block is missing from configuration file. This block is required.";

  if (!parameters_.has_been_populated())
    parameters_.Populate(this);

  // Call validation for the other objects required by the model
  world_view_->Validate();

  managers_->Validate();

  /**
   * Do some more sanity checks
   */
  initialisationphases::Manager& init_phase_mngr = *managers_->initialisation_phase();
  for (const string& phase : initialisation_phases_) {
    if (!init_phase_mngr.IsPhaseDefined(phase))
      LOG_ERROR_P(PARAM_INITIALISATION_PHASE_LABELS) << "(" << phase << ") has not been defined. Please ensure you have defined it";
  }

  timesteps::Manager& time_step_mngr = *managers_->time_step();
  for (const string time_step : time_steps_) {
    if (!time_step_mngr.GetTimeStep(time_step))
      LOG_ERROR_P(PARAM_TIME_STEPS) << "(" << time_step << ") has not been defined. Please ensure you have defined it";
  }

  /**
   * Do some simple checks
   * e.g Validate that the length_bins are strictly increasing order
   */
  for(unsigned length = 0; length < (length_bins_.size() - 1); ++length) {
    if(length_bins_[length] < 0.0)
    if(length_bins_[length] > length_bins_[length + 1])
      LOG_ERROR_P(PARAM_LENGTH_BINS) << ": Length bins must be strictly increasing " << length_bins_[length] << " is greater than " << length_bins_[length +1];
  }
}

/**
 *
 */
void Model::Build() {
  LOG_TRACE();
  world_view_->Build();
  managers_->Build();

}

/**
 *
 */
void Model::Verify() {
  LOG_TRACE();
  for (auto executor : executors_[state_])
    executor->Execute();
}

/**
 *
 */
void Model::Reset() {
  LOG_TRACE();

  world_view_->Reset();
  managers_->Reset();
}

/**
 *
 */
void Model::RunBasic() {
  LOG_TRACE();
  vector<string> single_step_addressables;
  vector<Double*> estimable_targets;
  // Create an instance of all categories

  // Model is about to run
  for (unsigned i = 0; i < adressable_values_count_; ++i) {
    //if (addressable_values_file_) {
    // estimables.LoadValues(i);
    //  Reset();
    //}

    /**
     * Running the model now
     */
    LOG_FINE() << "Model: State change to Execute";
    state_ = State::kInitialise;
    current_year_ = start_year_;
    // Iterate over all partition members and UpDate Mean Weight for the inital weight calculations

    initialisationphases::Manager& init_phase_manager = *managers_->initialisation_phase();
    init_phase_manager.Execute();
    managers_->report()->Execute(State::kInitialise);

    state_ = State::kExecute;


    timesteps::Manager& time_step_manager = *managers_->time_step();
    timevarying::Manager& time_varying_manager = *managers_->time_varying();
    for (current_year_ = start_year_; current_year_ <= final_year_; ++current_year_) {
      LOG_FINE() << "Iteration year: " << current_year_;
      time_varying_manager.Update(current_year_);
      LOG_FINEST() << "finishing update time varying now Update Category mean length and weight before beginning annual cycle";


      time_step_manager.Execute(current_year_);
    }

    managers_->observation()->CalculateScores();

    for (auto executor : executors_[State::kExecute])
      executor->Execute();

    // Model has finished so we can run finalise.

    LOG_FINE() << "Model: State change to Iteration Complete";
    managers_->report()->Execute(State::kIterationComplete);
  }
}


/**
 * This method will do a single iteration of the model. During
 * a basic run it'll only run once, but during the other run modes i.e. estiamtion and MCMC
 * it'll run multiple times.
 */
void Model::Iterate() {
  LOG_TRACE();
  // Create an instance of all categories
  state_ = State::kInitialise;
  current_year_ = start_year_;
  // Iterate over all partition members and UpDate Mean Weight for the inital weight calculations

  initialisationphases::Manager& init_phase_manager = *managers_->initialisation_phase();
  init_phase_manager.Execute();
  managers_->report()->Execute(State::kInitialise);

  state_ = State::kExecute;
  timesteps::Manager& time_step_manager = *managers_->time_step();
  timevarying::Manager& time_varying_manager = *managers_->time_varying();
  for (current_year_ = start_year_; current_year_ <= final_year_; ++current_year_) {
    LOG_FINE() << "Iteration year: " << current_year_;
    time_varying_manager.Update(current_year_);
    // Iterate over all partition members and UpDate Mean Weight for the inital weight calculations

    time_step_manager.Execute(current_year_);
  }

  managers_->observation()->CalculateScores();

  for (auto executor : executors_[State::kExecute])
    executor->Execute();

  current_year_ = final_year_;
}

/**
 *
 */
void Model::FullIteration() {
  Reset();
  Iterate();
}
} /* namespace niwa */
