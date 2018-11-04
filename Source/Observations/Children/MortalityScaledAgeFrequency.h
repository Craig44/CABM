/**
 * @file MortalityScaledAgeFrequency.h
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
#ifndef OBSERVATIONS_MORTALITY_SCALED_AGE_FREQUENCY_H_
#define OBSERVATIONS_MORTALITY_SCALED_AGE_FREQUENCY_H_

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
class MortalityScaledAgeFrequency : public niwa::Observation {
public:
  // methods
	MortalityScaledAgeFrequency(Model* model);
  virtual                     ~MortalityScaledAgeFrequency();
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
  Mortality*                      mortality_process_ = nullptr;
  string                          process_label_;
  string                          stratum_weight_method_;
  map<string,vector<unsigned>>    stratum_rows_;
  map<string,vector<unsigned>>    stratum_cols_;

  parameters::Table*              sample_table_ = nullptr;
  map<unsigned, map<string,vector<unsigned>>>  samples_by_year_and_stratum_;



};

} /* namespace observations */
} /* namespace niwa */

#endif /* OBSERVATIONS_MORTALITY_SCALED_AGE_FREQUENCY_H_ */
