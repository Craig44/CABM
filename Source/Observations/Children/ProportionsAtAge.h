/**
 * @file ProportionsAtAge.h
 * @author C.Marsh
 * @date 28/09/2013
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * calcualte proportions at age over a mortality block
 */
#ifndef OBSERVATIONS_PROPORTIONS_AT_AGE_H_
#define OBSERVATIONS_PROPORTIONS_AT_AGE_H_

// headers
#include "Observations/Observation.h"
#include "Layers/Children/CategoricalLayer.h"
#include "AgeingErrors/AgeingError.h"

//#include "Catchabilities/Catchability.h"

// namespaces
namespace niwa {
class Selectivity;
namespace observations {

/**
 * class definition
 */
class ProportionsAtAge : public niwa::Observation {
public:
  // methods
  ProportionsAtAge(Model* model);
  virtual                       ~ProportionsAtAge();
  void                          DoValidate() override final;
  virtual void                  DoBuild() override;
  void                          DoReset() override final { };
  void                          PreExecute() override final;
  void                          Execute() override final;
  void                          Simulate() override final;
  bool                          HasYear(unsigned year) const override final { return std::find(years_.begin(), years_.end(), year) != years_.end(); }

protected:
  vector<unsigned>              years_;
  unsigned                      min_age_ = 0;
  unsigned                      max_age_ = 0;
  bool                          plus_group_ = false;
  unsigned                      age_spread_ = 0;
  parameters::Table*            error_values_table_ = nullptr;
  AgeingError*                  ageing_error_ = nullptr;
  string                        ageing_error_label_;
  vector<float>                 age_results_;
  string                        process_label_;
  vector<string>                selectivity_labels_;
  vector<Selectivity*>          selectivities_;
  bool                          selectivity_length_based_ = false;
  string                        time_step_label_ = "";
  float                         time_step_proportion_;
  WorldView*                    world_ = nullptr;
  bool                          sexed_;

  unsigned                      time_step_to_execute_;
  map<unsigned, vector<float>>  error_values_by_year_;
  map<unsigned,vector<float>>   proportions_;
  map<unsigned, vector<float>>  error_values_;
  vector<string>                cells_;
  layers::CategoricalLayer*     layer_ = nullptr;
  string                        layer_label_;
  bool                          are_obs_props_ = true;
  map<string, vector<float>>     pre_age_freq_;
  map<string, vector<float>>     age_freq_;
  map<string, vector<float>>     final_age_freq_;


};

} /* namespace observations */
} /* namespace niwa */

#endif /* OBSERVATIONS_PROPORTIONS_AT_AGE_H_ */
