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
#include <omp.h>

#include "Factory.h"
#include "Managers.h"
#include "Objects.h"
#include "Agents/Agent.h"
#include "Addressables/Addressables.h"
#include "ConfigurationLoader/AddressableValuesLoader.h"
#include "GlobalConfiguration/GlobalConfiguration.h"
#include "InitialisationPhases/Manager.h"
#include "Logging/Logging.h"
#include "Observations/Manager.h"
#include "Selectivities/Manager.h"
#include "World/WorldView.h"
#include "Reports/Manager.h"
#include "Processes/Manager.h"
#include "TimeSteps/Manager.h"
#include "TimeVarying/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/To.h"
#include <RInside.h>                    // for the embedded R via RInside

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
  parameters_.Bind<string>(PARAM_INITIALISATION_PHASE_LABELS, &initialisation_phases_, "Define the labels of the phases of the initialisation", R"(A list of valid labels defined by \texttt{@initialisation_phase})", false);
  parameters_.Bind<string>(PARAM_TIME_STEPS, &time_steps_, "Define the labels of the time steps, in the order that they are applied, to form the annual cycle", R"(A list of valid labels defined by \texttt{@time_step})");
  parameters_.Bind<unsigned>(PARAM_LENGTH_BINS, &length_bins_, "", "", false);
  parameters_.Bind<bool>(PARAM_LENGTH_PLUS, &length_plus_, "Is the last bin a plus group", "", false);
  parameters_.Bind<string>(PARAM_BASE_LAYER_LABEL, &base_layer_, "Label for the base layer", "");
  parameters_.Bind<float>(PARAM_LATITUDE_BOUNDS, &lat_bounds_, "Latitude bounds for the spatial domain, should include lower and upper bound, so there should be rows + 1 values", "", true);
  parameters_.Bind<float>(PARAM_LONGITUDE_BOUNDS, &lon_bounds_, "Longitude bounds for the spatial domain, should include lower and upper bound, so there should be columns + 1 values", "", true);
  parameters_.Bind<unsigned>(PARAM_NROWS, &world_height_, "number of rows in spatial domain", "")->set_lower_bound(1,true);
  parameters_.Bind<unsigned>(PARAM_NCOLS, &world_width_, "number of columns in spatial domain", "")->set_lower_bound(1,true);
  parameters_.Bind<bool>(PARAM_SEXED, &sex_, "Is sex an attribute of then agents?", "", false);
  //parameters_.Bind<float>(PARAM_PROPORTION_MALE, &proportion_male_, "what proportion of the generated agents should be male?", "", 1.0)->set_range(0.0, 1.0, true, true); //TODO move to recruitment so it can be time-varying.
  parameters_.Bind<string>(PARAM_MATRUITY_OGIVE_LABEL, &maturity_ogives_, "Maturity ogive label for each sex (male female) order is important", "", false);
  parameters_.Bind<string>(PARAM_GROWTH_PROCESS_LABEL, &growth_process_label_, "Label for the growth process in the annual cycle", "");
  parameters_.Bind<string>(PARAM_NATURAL_MORTALITY_PROCESS_LABEL, &natural_mortality_label_, "Label for the natural mortality process in the annual cycle", "");
  //parameters_.Bind<unsigned>(PARAM_MAX_THREADS_TO_USE, &max_threads_, "The maxiumum threads you want to give access to this program", "",1);


  //RegisterAsAddressable(PARAM_PROPORTION_MALE, &proportion_male_); // can make this time-varying

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
  LOG_MEDIUM();
  Logging& logging = Logging::Instance();
  if (logging.errors().size() > 0) {
    logging.FlushErrors();
    return false;
  }

  run_mode_ = run_mode;

  if (state_ != State::kStartUp)
    LOG_CODE_ERROR() << "Model state should always be startup when entering the start method";
  if (global_configuration_->addressable_value_file() != "") {
    configuration::AddressableValuesLoader loader(this);
    loader.LoadValues(global_configuration_->addressable_value_file());
  }

  managers_->report()->Execute(State::kStartUp);

  LOG_MEDIUM() << "Model: State Change to Validate";
  state_ = State::kValidate;
  Validate();
  if (logging.errors().size() > 0) {
    logging.FlushErrors();
    return false;
  }

  managers_->report()->Execute(state_);

  LOG_MEDIUM() << "Model: State Change to Build";
  state_ = State::kBuild;
  Build();
  if (logging.errors().size() > 0) {
    logging.FlushErrors();
    return false;
  }
  managers_->report()->Execute(state_);

  LOG_MEDIUM() << "Model: State Change to Verify";
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
  case RunMode::kSimulation:
    RunBasic();
    break;
  default:
    LOG_ERROR() << "Invalid run mode has been specified. This run mode is not supported: " << run_mode_;
    break;
  }

  // finalise all reports
  LOG_MEDIUM() << "Finalising Reports";
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

  if (sex_ & (maturity_ogives_.size() != 2))
    LOG_FATAL_P(PARAM_MATRUITY_OGIVE_LABEL) << "if you have specified a model with sex, you need to supply two maturity ogives one for each sex (they can be the same)";

  /**
   * Do some simple checks
   * e.g Validate that the length_bins are strictly increasing order
   */
  for(unsigned length = 0; length < (length_bins_.size() - 1); ++length) {
    if(length_bins_[length] < 0.0)
      LOG_ERROR_P(PARAM_LENGTH_BINS) << "the length bin " <<  length_bins_[length] << "is less than zero this is not allowed";
    if(length_bins_[length] > length_bins_[length + 1])
      LOG_ERROR_P(PARAM_LENGTH_BINS) << ": Length bins must be strictly increasing " << length_bins_[length] << " is greater than " << length_bins_[length +1];
  }

  // Calculate length bin midpoints
  for (unsigned length_ndx = 1; length_ndx < length_bins_.size(); ++length_ndx) {
    length_bin_mid_points_.push_back((length_bins_[length_ndx - 1] + ((float)length_bins_[length_ndx] - (float)length_bins_[length_ndx - 1]) / 2));
    LOG_FINE() << "upper = " << length_bins_[length_ndx] << " lower = " << length_bins_[length_ndx - 1] << " mid = " << length_bin_mid_points_[length_ndx - 1];
  }
  length_bin_number_ = length_bin_mid_points_.size();
  LOG_FINE() << "length bins " << length_bins_.size() << " length bin midpoints = " << length_bin_mid_points_.size();



  // Call validation for the other objects required by the model
  LOG_FINE() << "About to validate world";
  world_view_->Validate();
  LOG_FINE() << "About to validate the rest of the managers";
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
  LOG_FINE() << "Exit validation";
}

/**
 *
 */
void Model::Build() {
  LOG_FINE();
  // Build lat long stuff
  if (parameters_.Get(PARAM_LATITUDE_BOUNDS)->has_been_defined() && parameters_.Get(PARAM_LONGITUDE_BOUNDS)->has_been_defined()) {
    if (lon_bounds_.size() != (world_width_ + 1)) {
      LOG_FATAL_P(PARAM_LONGITUDE_BOUNDS) << "longitude bounds must have a minimum and maximum value for each cell, e.g ncol 2 lon_bounds 150 151 152. You supplied '" << lon_bounds_.size() << " but we wanted " << (world_width_ + 1);
    }

    if (lat_bounds_.size() != (world_height_ + 1)) {
      LOG_FATAL_P(PARAM_LATITUDE_BOUNDS) << "latitude bounds must have a minimum and maximum value for each cell, e.g ncol 2 lat_bounds -45 -44 -43. You supplied '" << lat_bounds_.size() << " but we wanted " << (world_height_ + 1);
    }
    for (unsigned lat_ndx = 1; lat_ndx < lat_bounds_.size(); ++lat_ndx) {
      lat_mid_points_.push_back((float)(lat_bounds_[lat_ndx - 1] + ((lat_bounds_[lat_ndx] - lat_bounds_[lat_ndx - 1]) / 2)));
      if (lat_bounds_[lat_ndx] < lat_bounds_[lat_ndx - 1])
        LOG_ERROR_P(PARAM_LATITUDE_BOUNDS) << "must be in ascending order";
    }

    for (unsigned lon_ndx = 1; lon_ndx < lon_bounds_.size(); ++lon_ndx) {
      lon_mid_points_.push_back((float)(lon_bounds_[lon_ndx - 1] + ((lon_bounds_[lon_ndx] - lon_bounds_[lon_ndx - 1]) / 2)));
      if (lon_bounds_[lon_ndx] < lon_bounds_[lon_ndx - 1])
         LOG_ERROR_P(PARAM_LONGITUDE_BOUNDS) << "must be in ascending order";
    }
    min_lon_ = lon_bounds_[0];
    max_lon_ = lon_bounds_[lon_bounds_.size() - 1];
    min_lat_ = lat_bounds_[0];
    max_lat_ = lat_bounds_[lat_bounds_.size() - 1];
    LOG_FINE() << "min lat = " << min_lat_ << " max lat = " << max_lat_ << " min long = " << min_lon_ << " max lon " << max_lon_;
  }



  /*
   * An important sequence in the code, if you cannot obtain pointers at build the order of managers will be important
  */
  managers_->BuildPreWorldView();

  // Do a quick check that we can obtain pointers to mortality and growth process
  processes::Manager& process_manager = *managers_->process();
  auto growth_ptr = process_manager.GetGrowthProcess(growth_process_label_);
  if (growth_ptr == nullptr) {
    LOG_FATAL_P(PARAM_GROWTH_PROCESS_LABEL) << "Does the growth process " << growth_process_label_ << " exist? please check that it does, and is a growth process";
  }
  if (!process_manager.GetMortalityProcess(natural_mortality_label_)) {
    LOG_FATAL_P(PARAM_NATURAL_MORTALITY_PROCESS_LABEL) << "Does the mortality process " << natural_mortality_label_ << " exist? please check that it does, and is a mortality process";
  }
  if (maturity_ogives_.size() > 2) {
    LOG_ERROR_P(PARAM_NATURAL_MORTALITY_PROCESS_LABEL) << "You have specified '" << maturity_ogives_.size() << "' maturity ogives, we only use one if it is unsexed or two if it is sexed";
  }
  // Check maturity Ogive exists and flag to user
  bool first = true;
  bool selectivity_length_based;
  for (unsigned i = 0; i < maturity_ogives_.size(); ++i) {
    Selectivity* selec = managers_->selectivity()->GetSelectivity(maturity_ogives_[i]);
    if (!selec)
      LOG_ERROR_P(PARAM_MATRUITY_OGIVE_LABEL) << "Could not find maturity ogive (selectivity) '" << maturity_ogives_[i] << "', make sure it is defined";
    if (first) {
      first = false;
      selectivity_length_based = selec->is_length_based();
    } else {
      if (selectivity_length_based != selec->is_length_based()) {
        LOG_ERROR_P(PARAM_MATRUITY_OGIVE_LABEL) << "The selectivity  " << maturity_ogives_[i] << " was not the same type (age or length based) as the previous selectivity label";
      }
    }
  }

  world_view_->Build(); // This needs processes to be built, but others want world to be built by DoBuild to do checks, hmmm
  managers_->Build();




  // Check thread logic
  unsigned procs = omp_get_max_threads() - 2; // Default to number availble threads less two one for other stuff and another for reports
  if ( max_threads_ > 0 && max_threads_ < procs )
    procs = max_threads_;
  if ( procs < 1)
    procs = 1;

  LOG_MEDIUM() << "setting number of threads = " << procs;
  omp_set_dynamic(0);
  omp_set_num_threads(procs);


  Addressables& addressables = *managers_->addressables();
  if (addressables.GetValueCount() > 0) {
    addressable_values_file_ = true;
    adressable_values_count_ = addressables.GetValueCount();
  }
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
void Model::RunMSE() {

}
/**
 *
 */
void Model::Reset() {
  LOG_MEDIUM();
  world_view_->Reset();
  managers_->Reset();
}

/**
 *
 */
void Model::RunBasic() {
  LOG_MEDIUM();
  // Model is about to run
  /*
   * - iterate over -i file
   * - iterate over -s values
   */
  // set up R instance and packages
  int argc = 1;
  char* argv[0];
  RInside R(argc, argv);          // create an embedded R instance
  try {
     std::string setup_R = "source(file.path('R','SetUpR.R'))";
     //std::string setup_R = "suppressMessages(source(file.path('R','SetUpR.R')))";
     R.parseEvalQ(setup_R);              // load library, no return value
  } catch(std::exception& ex) {
      std::cerr << "Exception caught: " << ex.what() << std::endl;
  } catch(...) {
      std::cerr << "Unknown exception caught" << std::endl;
  }
  // Run C++ algorithm up to current time-step
  R.parseEvalQ("sim = 1");


  Addressables& addressables = *managers_->addressables();
  LOG_MEDIUM() << "Multi line value = " << adressable_values_count_;

  int simulation_candidates = global_configuration_->simulation_candidates();
  if (simulation_candidates < 1) {
    LOG_FATAL() << "The number of simulations specified at the command line parser must be at least one";
  }
  unsigned suffix_width = (unsigned)(floor(log10((double) (simulation_candidates) * (adressable_values_count_))) + 1);
  LOG_MEDIUM() << "suffix width = " << suffix_width << " value = " << (simulation_candidates) * (adressable_values_count_);
  unsigned suffix_counter = 0;
  for (unsigned i = 0; i < adressable_values_count_; ++i) {
    if (addressable_values_file_) {
      addressables.LoadValues(i);
      LOG_MEDIUM() << "about to reset";
      Reset();
     }
    /**
     * Running the model now
     */
    LOG_MEDIUM() << "Model: State change to Initialisation";
    state_ = State::kInitialise;
    current_year_ = start_year_ - 1;
    // Iterate over all partition members and UpDate Mean Weight for the inital weight calculations
    initialisationphases::Manager& init_phase_manager = *managers_->initialisation_phase();
    timevarying::Manager& time_varying_manager = *managers_->time_varying();
    if (i == 0) {
      LOG_MEDIUM() << "first initialisation phase.";
      init_phase_manager.Execute();
      if (not re_run_initialisation_) {
        LOG_MEDIUM() << "Cache initialsiation";
        world_view_->CachedWorldForInit();
      }
    } else if (re_run_initialisation_) {
      LOG_MEDIUM() << "Re-run initialisation";
      init_phase_manager.Execute();
    } else {
      LOG_MEDIUM() << "resetting initial world view";
      world_view_->MergeWorldForInit();
    }

    managers_->report()->Execute(State::kInitialise);

    LOG_MEDIUM() << "Model: State change to Execute";

    state_ = State::kExecute;

    timesteps::Manager& time_step_manager = *managers_->time_step();
    for (current_year_ = start_year_; current_year_ <= final_year_; ++current_year_) {
      LOG_FINE() << "Iteration year: " << current_year_;
      LOG_FINE() << "Update time varying params";
      time_varying_manager.Update(current_year_);
      LOG_FINE() << "finishing update time varying now Update Category mean length and weight before beginning annual cycle";
      world_view_->rebuild_agent_time_varying_params();
      LOG_FINE() << "rebuild agent values continue to execute the year";
      time_step_manager.Execute(current_year_);
      LOG_FINE() << "finished year exectution";
    }

    LOG_MEDIUM() << " Heading into simulation mode";
    for (int s = 0; s < simulation_candidates; ++s) {
      suffix_counter++;
      string report_suffix = ".";

      unsigned iteration_width = (unsigned)(floor(log10((i + 1) + (s + 1))) + 1);
      LOG_FINE() << "iteration_width = " << iteration_width << " suffix_width = " << suffix_width;
      unsigned diff = suffix_width - iteration_width;
      LOG_FINE() << "diff = " << diff;
      LOG_FINE() << "(i + 1) + (s * simulation_candidates)) = " << (i + 1) + (s * simulation_candidates);
      report_suffix.append(diff,'0');
      report_suffix.append(utilities::ToInline<unsigned, string>(suffix_counter));
      LOG_MEDIUM() << "i = " << i + 1 << " s = " << s + 1 <<  " suffix = " << report_suffix << " what i think it should be doing " << (i + 1) + (s * simulation_candidates) << " diff = " << diff << " iteration width = " << iteration_width;
      managers_->report()->set_report_suffix(report_suffix);

      managers_->observation()->SimulateData();

      // Not convinced this is doing anything
      //for (auto executor : executors_[State::kExecute])
      //  executor->Execute();

      if (s != (simulation_candidates - 1)) { // Only need to execute this s - 1 times as the last run will be done at line 485
        managers_->report()->Execute(State::kIterationComplete);
        managers_->report()->WaitForReportsToFinish();
      }
    }
    // Model has finished so we can run finalise.
    LOG_MEDIUM() << "Model: State change to Iteration Complete";
    managers_->report()->Execute(State::kIterationComplete);
  }
}


/**
 * This method will do a single iteration of the model. During
 * a basic run it'll only run once, but during the other run modes i.e. estiamtion and MCMC
 * it'll run multiple times.
 */
void Model::Iterate() {
  LOG_MEDIUM();
  // Create an instance of all categories
  state_ = State::kInitialise;
  current_year_ = start_year_;
  // Iterate over all partition members and UpDate Mean Weight for the inital weight calculations

  initialisationphases::Manager& init_phase_manager = *managers_->initialisation_phase();
  init_phase_manager.Execute();
  managers_->report()->Execute(State::kInitialise);

  state_ = State::kExecute;
  timesteps::Manager& time_step_manager = *managers_->time_step();
  for (current_year_ = start_year_; current_year_ <= final_year_; ++current_year_) {
    LOG_FINE() << "Iteration year: " << current_year_;
    time_step_manager.Execute(current_year_);
  }

  managers_->observation()->SimulateData();

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

bool Model::lat_and_long_supplied() {
  if (parameters_.Get(PARAM_LATITUDE_BOUNDS)->has_been_defined() && parameters_.Get(PARAM_LONGITUDE_BOUNDS)->has_been_defined())
    return true;
  else
    return false;
}




} /* namespace niwa */
