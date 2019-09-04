/**
 * @file ProcessRemovalsByAge.h
 * @author  C Marsh
 * @version 1.0
 * @date 25/08/15
 * @section LICENSE
 *
 * Copyright NIWA Science 2016 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This observation is a specific process observation class. It is associated with the process type mortality_instantaneous.
 * It calls a catch at object that is created from the process which represents the numbers at age halfway trough the mortality process
 * This class then applies ageing error and converts to a proportion which then gets sent to a likelihood for evaluation.
 *
 */
#ifndef OBSERVATIONS_REMOVEALS_BY_AGE_H_
#define OBSERVATIONS_REMOVEALS_BY_AGE_H_

// Headers
#include "Observations/Observation.h"

#include "Processes/Children/Mortality.h"
#include "AgeingErrors/AgeingError.h"
#include "Layers/Children/CategoricalLayer.h"

// Namespace
namespace niwa {
namespace observations {

using processes::Mortality;

/**
 * Class Definition
 */
class ProcessRemovalsByAge : public niwa::Observation {
public:
  // Methods
  explicit ProcessRemovalsByAge(Model* model);
  virtual                     ~ProcessRemovalsByAge();
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final {first_simualtion_run_ = true;};
  void                        PreExecute() override final;
  void                        Execute() override final;
  void                        Simulate() override final;
  bool                        HasYear(unsigned year) const override final { return std::find(years_.begin(), years_.end(), year) != years_.end(); }

protected:
  // Members
  vector<unsigned>              years_;
  unsigned                      min_age_ = 0;
  unsigned                      max_age_ = 0;
  bool                          plus_group_ = false;
  unsigned                      age_spread_ = 0;
  parameters::Table*            error_values_table_ = nullptr;
  AgeingError*                  ageing_error_ = nullptr;
  string                        ageing_error_label_;
  Mortality*                    mortality_process_ = nullptr;
  vector<float>                 age_results_;
  string                        process_label_;

  unsigned                      time_step_to_execute_;
  map<unsigned, vector<float>>  error_values_by_year_;
  map<unsigned,vector<float>>   proportions_;
  map<unsigned, vector<float>>  error_values_;
  vector<string>                cells_;
  layers::CategoricalLayer*     layer_ = nullptr;
  string                        layer_label_;
  bool                          first_simualtion_run_ = true;
  string                        sexed_;
  unsigned                      sex_match_;
  bool                          sexed_flag_;
  vector<float>                 accumulated_age_frequency_;


};

} /* namespace observations */
} /* namespace niwa */

#endif /* OBSERVATIONS_REMOVEALS_BY_AGE_H_ */
