/**
 * @file MortalityEventBiomassAgeClusters.h
 * @author  C.Marsh github.com/Craig44
 * @date 4/11/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This class takes census information from a mortality process, generates clusters (tow/trip level tags) then randomly selects a sub sample of ages and lengths
 *  to then generate an age frequency either via age length key or direct ageing.
 */
#ifndef OBSERVATIONS_MORTALITY_EVENT_COMPOSITION_H_
#define OBSERVATIONS_MORTALITY_EVENT_COMPOSITION_H_

// headers
#include "Observations/Observation.h"
#include "Layers/Children/CategoricalLayer.h"

#include "AgeingErrors/AgeingError.h"
//#include "Processes/Children/Mortality.h"
#include "Processes/Children/Mortality.h"
#include <boost/math/distributions/normal.hpp>


// namespaces
namespace niwa {
namespace observations {

using processes::Mortality;

/**
 * class definition
 */
class MortalityEventComposition : public niwa::Observation {
public:
  // methods
  MortalityEventComposition(Model* model);
  virtual                     ~MortalityEventComposition();
  void                        DoValidate() override final;
  virtual void                DoBuild() override;
  void                        DoReset() override final { };
  void                        PreExecute() override final;
  void                        Execute() override final;
  void                        Simulate() override final;
  bool                        HasYear(unsigned year) const override final { return std::find(years_.begin(), years_.end(), year) != years_.end(); }
  virtual void                FillReportCache(ostringstream& cache);

protected:
  // Methods
  void                            ResetPreSimulation();

  // members
  vector<unsigned>                years_;
  vector<string>                  cells_;
  layers::CategoricalLayer*       layer_ = nullptr;
  string                          layer_label_;
  AgeingError*                    ageing_error_ = nullptr;
  string                          ageing_error_label_;
  string                          ageing_allocation_ = PARAM_RANDOM;
  AllocationType                  allocation_type_ = AllocationType::kRandom;

  bool                            sexed_flag_;
  unsigned                        n_bins_;
  unsigned                        n_unsexed_bins_;
  // Error values
  map<unsigned, map<string, unsigned>>   error_value_by_year_and_stratum_;  // year x stratum
  //
  string                          obs_type_;
  string                          comp_type_;
  bool                            are_obs_props_;
  bool                            is_age_;

  unsigned                        min_age_ = 0;
  unsigned                        max_age_ = 0;
  bool                            plus_group_ = false;
  unsigned                        age_spread_ = 0;

  string                          fishery_label_;
  Mortality*                      mortality_process_ = nullptr;
  vector<vector<processes::composition_data>>* fishery_comp_data_ = nullptr;
  string                          process_label_;
  string                          stratum_weight_method_ = PARAM_NONE;
  // TODO change from string -> unsigned int for a little speed up
  map<string,vector<unsigned>>    stratum_rows_;
  map<string,vector<unsigned>>    stratum_cols_;
  map<string,float>               stratum_area_;
  vector<float>                   stratum_comp_;
  map<string,float>               stratum_biomass_;
  WorldView*                      world_ = nullptr;
  vector<unsigned>                fishery_years_;
  map<string,vector<float>>       stratum_age_frequency_;
  vector<vector<float>>           ageing_mis_matrix_;
  parameters::Table*              error_values_table_ = nullptr;



};

} /* namespace observations */
} /* namespace niwa */

#endif /* OBSERVATIONS_MORTALITY_EVENT_COMPOSITION_H_ */
