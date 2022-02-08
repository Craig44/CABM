/**
 * @file MortalityEffortBasedWithCovar.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 10/09/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This is a child mortality that applies a fishing event
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_EFFORT_BASED_WITH_COVAR_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_EFFORT_BASED_WITH_COVAR_H_

// headers
#include "Processes/Children/Mortality.h"

#include "Layers/Children/NumericLayer.h"
#include "PreferenceFunctions/PreferenceFunction.h"
#include <omp.h>

#include <ctime>

// Includes for the GammaDiff minimiser
#include "Minimisers/Minimiser.h"

// namespaces
namespace niwa {
class Selectivity;
class Minimiser;
namespace processes {
using std::string;
/**
 * Class definition
 */
class MortalityEffortBasedWithCovar : public Mortality {
public:
  // methods
  explicit MortalityEffortBasedWithCovar(Model* model);
  virtual                     ~MortalityEffortBasedWithCovar() = default;
  virtual void                        DoValidate() override final;
  virtual void                        DoBuild() override final;
  virtual void                        DoReset() override final;
  virtual void                        DoExecute() override final;
  virtual double                      SolveBaranov() override final;

  void                                draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) override final { };
  void                                FillReportCache(ostringstream& cache) override final;

protected:
  vector<double>                       catches_;
  vector<string>                      selectivity_label_;
  vector<Selectivity*>                selectivity_;
  bool                                selectivity_length_based_ = false;
  map<unsigned, double>                catches_by_year_;
  double                               catchability_;
  double                              actual_catch_;
  map<unsigned,double>                actual_catch_by_year_;
  map<unsigned,double>                vulnerable_biomass_by_year_;

  vector<string>                      preference_function_labels_;
  vector<string>                      preference_layer_labels_;
  string                              closure_layer_label_ = "";
  vector<double>                      preference_weights_;
  bool                                calculate_on_the_fly_ = false;
  vector<unsigned>                    non_static_layer_ndx_;
  vector<PreferenceFunction*>         preference_functions_;
  vector<layers::NumericLayer*>       preference_layers_;
  layers::NumericLayer*               closure_layer_ = nullptr;

  vector<vector<double>>              base_preference_;
  vector<vector<double>>              vary_preference_;
  vector<vector<double>>              temp_preference_;

  map<unsigned, vector<vector<double>>> preference_by_year_;

  double                              catch_based_on_baranov_;
  map<unsigned,double>                catch_based_on_baranov_by_year_;

  map<unsigned,double>                lambda_by_year_;
  // objects for thread safety of rng
  vector<double>                       random_numbers_;
  vector<double>                       discard_random_numbers_;
  vector<double>                       selectivity_random_numbers_;

  unsigned                            n_agents_;
  vector<vector<unsigned>>            cell_offset_;

  vector<vector<double>>              effort_by_cell_;
  vector<vector<double>>              vulnerable_by_cell_;
  vector<vector<double>>              F_by_cell_;
  vector<vector<double>>              removals_by_cell_;
  map<unsigned,vector<vector<double>>>  actual_removals_by_year_and_cell_;
  map<unsigned,vector<vector<double>>>  F_by_year_and_cell_;
  map<unsigned,vector<vector<double>>>  vulnerable_by_year_and_cell_;
  map<unsigned,vector<vector<double>>>   effort_by_year_and_cell_;
  vector<size_t>                      effort_index_;
  vector<double>                       vulnerable_biomass_vector_format_;
  vector<double>                       effort_organised_vector_format_;
  vector<double>                       pref_organised_vector_format_;

  vector<double>                       age_freq_census_;
  map<unsigned,vector<double>>         age_freq_census_by_year_;

  vector<double>                       effort_input_;
  string                              effort_layer_label_;
  layers::NumericLayer*               effort_layer_ = nullptr;

  // For reporting
  map<unsigned, double>                actual_removals_by_year_;
  map<unsigned, double>                removals_by_year_;
  Minimiser*                          minimiser_ = nullptr;
  string                              minimiser_label_;
  vector<double>                      start_value_for_lambda_;

  time_t                              start_time_;
  map<unsigned, time_t>               time_by_year_;
  map<unsigned, string>               message_by_year_;
  bool                                first_execute_ = true;
  double                               max_vulnerable_;

};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_EFFORT_BASED_WITH_COVAR_H_ */
