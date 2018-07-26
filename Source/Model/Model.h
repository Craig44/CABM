/**
 * @file Model.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 16/11/2012
 * @section LICENSE
 *
 * Copyright NIWA Science ©2012 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This class is the primary representation of our model and it's states
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef MODEL_H_
#define MODEL_H_

// Headers
#include "BaseClasses/Executor.h"
#include "BaseClasses/Object.h"
#include "GlobalConfiguration/GlobalConfiguration.h"
#include "Utilities/Math.h"
#include "Utilities/RunMode.h"

// Namespaces
namespace niwa {
using base::Executor;
class Managers;
class Objects;
class Agents;
class Factory;
class Partition;
class ObjectiveFunction;
class WorldView;

namespace State {
enum Type {
  kStartUp, // system is loading from configuration file
  kValidate, // validating user supplied values for variables
  kBuild, // building and checking relationships between objects
  kVerify, // verifying business rules (not yet implemented)
  kInitialise, // running through the initialisation phases
  kExecute, // execute the object
  kIterationComplete, // a single iteration of the model is complete
  kReset, // called between iterations to ensure objects caches are reset
  kInputIterationComplete, // a single run of the mode is complete using an input file to set estimables
  kFinalise // the model is finished
};
}

namespace Units {
enum Type {
  kGrams,
  kKilograms,
  kTonnes
};
} /* namespace Units */

/**
 * Class definition
 */
class Model : public base::Object {
public:
  // Methods
  Model();
  virtual                     ~Model();
  bool                        Start(RunMode::Type run_mode);
  void                        FullIteration();
  void                        Subscribe(State::Type state, Executor* executor) { executors_[state].push_back(executor); }
  void                        PopulateParameters();

  // Accessors
  RunMode::Type               run_mode() const { return run_mode_; }
  State::Type                 state() const { return state_; }
  virtual unsigned            start_year() const { return start_year_; }
  virtual unsigned            final_year() const { return final_year_; }
  unsigned                    projection_final_year() const { return projection_final_year_;}
  bool                        projection_final_phase() {return projection_final_phase_;}
  void                        set_projection_final_phase(bool phase) {projection_final_phase_ = phase;}
  virtual vector<unsigned>    years() const;
  unsigned                    year_spread() const;
  void                        increment_current_year_for_initialisation() {++current_year_ ;}
  void                        set_current_year_in_initialisation(unsigned value) {current_year_ = value;}
  virtual unsigned            current_year() const { return current_year_; } // pass it by reference so I can change it during initialisation
  virtual unsigned            min_age() const { return min_age_; }
  virtual unsigned            max_age() const { return max_age_; }
  virtual unsigned            age_spread() const { return (max_age_ - min_age_) + 1; }
  virtual bool                age_plus() const { return age_plus_; }

  // The following are recruitment specific quantities/accessors these are stock characteristics that are shared around the model
  void                        set_b0(string recruitment_label, float b0) {b0_[recruitment_label] = b0; }
  virtual const map<string, float>&          get_b0s() const {return b0_; }
  void                        set_r0(string recruitment_label, unsigned r0) {r0_[recruitment_label] = r0; }
  unsigned                    get_r0(string recruitment_label)  {return r0_[recruitment_label]; }
  void                        set_ssb(string recruitment_label, unsigned ssb) {ssb_[recruitment_label] = ssb; }
  float                       get_ssb(string recruitment_label)  {return ssb_[recruitment_label]; }
  void                        set_scalar(float value) {global_scalar_ = value; }
  float                       get_scalar()  {return global_scalar_; }
  void                        set_m(float value) {m_ = value; }
  float                       get_m()  {return m_; }
  void                        set_n_agents(unsigned value) {n_agents_ = value; }
  virtual const unsigned      get_n_agents() const {return n_agents_; }

  virtual const vector<string>& time_steps() const { return time_steps_; }
  const vector<string>&       initialisation_phases() const { return initialisation_phases_; }
  virtual const vector<unsigned>&  length_bins() const { return length_bins_; }
  virtual const vector<float>&     length_bin_mid_points() const { return length_bin_mid_points_; }

  virtual bool                length_plus() const { return length_plus_; }
  string&                     get_base_layer_label() {return base_layer_;}
  string&                     get_lat_layer_label() {return lat_layer_label_;}
  string&                     get_long_layer_label() {return lon_layer_label_;}
  unsigned                    get_height() {return world_height_;}
  unsigned                    get_width() {return world_width_;}
  string                      get_mortality_process() {return natural_mortality_label_;}
  string                      get_growth_process() {return growth_process_label_;}
  vector<string>              get_maturity_ogive() {return maturity_ogives_;};
  bool                        get_sexed() const {return sex_;};
  float                       get_male_proportions() const {return proportion_male_;};
  unsigned                    get_max_threads() const {return max_threads_;}


  // manager accessors
  virtual Managers&           managers();
  virtual Objects&            objects();
  GlobalConfiguration&        global_configuration() { return *global_configuration_; }
  virtual Factory&            factory();
  WorldView*                  world_view();

protected:
  // Methods
  void                        Validate();
  void                        Build();
  void                        Verify();
  void                        Iterate();
  void                        Reset();
  void                        RunBasic();

  // Members
  RunMode::Type               run_mode_ = RunMode::kInvalid;
  State::Type                 state_    = State::kStartUp;
  unsigned                    start_year_ = 0;
  unsigned                    final_year_ = 0;
  unsigned                    projection_final_year_ = 0;
  unsigned                    current_year_ = 0;
  unsigned                    min_age_ = 0;
  unsigned                    max_age_ = 0;
  string                      base_weight_units_;
  bool                        age_plus_ = true;
  vector<string>              initialisation_phases_;
  vector<string>              time_steps_;
  vector<unsigned>            length_bins_;
  bool                        length_plus_ = true;
  vector<float>               length_bin_mid_points_;
  bool                        addressable_values_file_ = false;
  unsigned                    adressable_values_count_ = 1;
  bool                        sexed_ = false;
  string                      base_layer_;
  string                      lat_layer_label_; // make it optional
  string                      lon_layer_label_; // make it optional
  unsigned                    world_height_;
  unsigned                    world_width_;
  map<string, float>          b0_;
  map<string, unsigned>       r0_;
  unsigned                    n_agents_;
  string                      natural_mortality_label_;
  string                      growth_process_label_;
  float                       m_; // natural mortality
  map<string, float>          ssb_;
  float                       global_scalar_ = 1.0;
  bool                        sex_;
  float                       proportion_male_;
  vector<string>              maturity_ogives_;
  unsigned                    max_threads_;


  Managers*                   managers_ = nullptr;
  Objects*                    objects_ = nullptr;
  GlobalConfiguration*        global_configuration_ = nullptr;
  Factory*                    factory_ = nullptr;
  WorldView*                  world_view_ = nullptr;
  bool                        projection_final_phase_ = false; // this parameter is for the projection classes. most of the methods are in the reset but they don't need to be applied
  // if the model is in the first iteration and storeing values.
  map<State::Type, vector<Executor*>> executors_;
};

} /* namespace niwa */
#endif /* MODEL_H_ */
