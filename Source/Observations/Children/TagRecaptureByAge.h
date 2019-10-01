/**
 * @file TagRecaptureByLength.h
 * @author  C.Marsh github.com/Craig44
 * @date 4/11/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This class takes tag-racapture object from mortality which, is responsible for undertaking recapture, and bundles it up so that we can simulated
 * many permutations to identigy the effect of trap shyness and other factors that generate overdispersion, or disrupt conventional assumptions
 */
#ifndef OBSERVATIONS_TAG_RECAPTURE_BY_LENGTH_H_
#define OBSERVATIONS_TAG_RECAPTURE_BY_LENGTH_H_

// headers
#include "Observations/Observation.h"
#include "Layers/Children/CategoricalLayer.h"

#include "AgeingErrors/AgeingError.h"
#include "Processes/Children/Mortality.h"


// namespaces
namespace niwa {
namespace observations {

using processes::Mortality;



/**
 * class definition
 */
class TagRecaptureByLength : public niwa::Observation {
public:
  // methods
  TagRecaptureByLength(Model* model);
  virtual                     ~TagRecaptureByLength();
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

  unsigned                        number_of_bootstraps_ = 0;

  Mortality*                      mortality_process_ = nullptr;
  string                          process_label_;
  string                          stratum_weight_method_ = PARAM_NONE;
  // TODO change from string -> unsigned int for a little speed up
  map<string,vector<unsigned>>    stratum_rows_;
  map<string,vector<unsigned>>    stratum_cols_;
  map<string,float>               stratum_area_;
  map<string,float>               stratum_biomass_;
  WorldView*                      world_ = nullptr;

  map<string,vector<float>>    stratum_age_frequency_;
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

#endif /* OBSERVATIONS_TAG_RECAPTURE_BY_LENGTH_H_ */
