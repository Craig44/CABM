/**
 * @file MortalityEventBiomassScaledAgeFrequency.h
 * @author  C.Marsh github.com/Craig44
 * @date 4/11/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This class takes census information from a mortality process, generates an age an age length key
 * with ageing error and scales the length frequency through the age length key to generate a final age
 * frequency
 */
#ifndef OBSERVATIONS_MORTALITY_EVENT_BIOMASS_SCALED_AGE_FREQUENCY_H_
#define OBSERVATIONS_MORTALITY_EVENT_BIOMASS_SCALED_AGE_FREQUENCY_H_

// headers
#include "Observations/Observation.h"
#include "Layers/Children/CategoricalLayer.h"

#include "AgeingErrors/AgeingError.h"
//#include "Processes/Children/Mortality.h"
#include "Processes/Children/Mortality/MortalityEventBiomass.h"


// namespaces
namespace niwa {
namespace observations {

using processes::MortalityEventBiomass;

/**
 * class definition
 */
class MortalityEventBiomassScaledAgeFrequency : public niwa::Observation {
public:
  // methods
  MortalityEventBiomassScaledAgeFrequency(Model* model);
  virtual                     ~MortalityEventBiomassScaledAgeFrequency();
  void                        DoValidate() override final;
  virtual void                DoBuild() override;
  void                        DoReset() override final { };
  void                        PreExecute() override final;
  void                        Execute() override final;
  void                        Simulate() override final;
  bool                        HasYear(unsigned year) const override final { return std::find(years_.begin(), years_.end(), year) != years_.end(); }

protected:
  // members
  vector<unsigned>                years_;
  vector<string>                  cells_;
  layers::CategoricalLayer*       layer_ = nullptr;
  string                          layer_label_;
  AgeingError*                    ageing_error_ = nullptr;
  string                          ageing_error_label_;
  string                          ageing_allocation_ = PARAM_RANDOM;
  AllocationType                  allocation_type_ = AllocationType::kRandom;
  string                          sexed_;
  unsigned                        sex_match_;
  bool                            sexed_flag_;


  unsigned                        number_of_bootstraps_ = 0;

  string                          fishery_label_;
  MortalityEventBiomass*          mortality_process_ = nullptr;
  string                          process_label_;
  string                          stratum_weight_method_ = PARAM_NONE;
  // TODO change from string -> unsigned int for a little speed up
  map<string,vector<unsigned>>    stratum_rows_;
  map<string,vector<unsigned>>    stratum_cols_;
  map<string,float>               stratum_area_;
  map<string,float>               stratum_biomass_;
  WorldView*                      world_ = nullptr;
  vector<unsigned>                fishery_years_;
  map<string,vector<float>>       stratum_age_frequency_;
  vector<vector<float>>           age_length_key_;


  parameters::Table*              sample_table_ = nullptr;
  parameters::Table*              lf_sample_table_ = nullptr;

  map<unsigned, map<string, unsigned>>  samples_by_year_and_stratum_;
  map<unsigned, map<string, float>>  prop_lf_by_year_and_stratum_;


  map<unsigned,map<string,vector<vector<float>>>>     age_length_key_by_year_stratum_;
  map<unsigned,map<string,vector<float>>>              lf_by_year_stratum_;
};

} /* namespace observations */
} /* namespace niwa */

#endif /* OBSERVATIONS_MORTALITY_EVENT_BIOMASS_SCALED_AGE_FREQUENCY_H_ */
