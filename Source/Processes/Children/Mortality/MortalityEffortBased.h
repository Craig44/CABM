/**
 * @file MortalityEffortBased.h
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
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_EFFORT_BASED_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_EFFORT_BASED_H_

// headers
#include "Processes/Children/Mortality.h"

#include "Layers/Children/NumericLayer.h"
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
class MortalityEffortBased : public Mortality {
public:
  // methods
  explicit MortalityEffortBased(Model* model);
  virtual                     ~MortalityEffortBased() = default;
  virtual void                        DoValidate() override final;
  virtual void                        DoBuild() override final;
  virtual void                        DoReset() override final;
  virtual void                        DoExecute() override final;
  virtual double                      SolveBaranov() override final;

  void                                draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) override final { };
  void                                FillReportCache(ostringstream& cache) override final;

protected:
  vector<float>                       catches_;
  vector<string>                      selectivity_label_;
  vector<Selectivity*>                selectivity_;
  bool                                selectivity_length_based_ = false;
  map<unsigned, float>                catches_by_year_;
  float                               catchability_;
  double                              actual_catch_;
  map<unsigned,double>                actual_catch_by_year_;
  map<unsigned,double>                vulnerable_biomass_by_year_;

  double                              catch_based_on_baranov_;
  map<unsigned,double>                catch_based_on_baranov_by_year_;

  map<unsigned,double>                lambda_by_year_;
  // objects for thread safety of rng
  vector<float>                       random_numbers_;
  vector<float>                       discard_random_numbers_;
  vector<float>                       selectivity_random_numbers_;

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
  vector<size_t>                      biomass_index_;
  vector<float>                       vulnerable_biomass_vector_format_;
  vector<float>                       effort_organised_vector_format_;

  vector<float>                       age_freq_census_;
  map<unsigned,vector<float>>         age_freq_census_by_year_;

  vector<float>                       effort_input_;
  string                              effort_layer_label_;
  layers::NumericLayer*               effort_layer_ = nullptr;

  // For reporting
  map<unsigned, float>                actual_removals_by_year_;
  map<unsigned, float>                removals_by_year_;
  Minimiser*                          minimiser_ = nullptr;
  string                              minimiser_label_;
  vector<double>                      start_value_for_lambda_;

  time_t                              start_time_;
  map<unsigned, time_t>               time_by_year_;
  map<unsigned, string>               message_by_year_;
  bool                                effort_layer_based_;
  bool                                first_execute_ = true;


};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_EFFORT_BASED_H_ */
