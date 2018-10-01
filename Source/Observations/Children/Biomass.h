/**
 * @file Biomass.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 7/01/2014
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * << Add Description >>
 */
#ifndef OBSERVATIONS_BIOMASS_H_
#define OBSERVATIONS_BIOMASS_H_

// headers
#include "Observations/Observation.h"
#include "Layers/Children/CategoricalLayer.h"

//#include "Catchabilities/Catchability.h"

// namespaces
namespace niwa {
class Selectivity;
namespace observations {

/**
 * class definition
 */
class Biomass : public niwa::Observation {
public:
  // methods
  Biomass(Model* model);
  virtual                     ~Biomass() = default;
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
  float                           catchability_value_;
  map<unsigned, float>           error_values_by_year_;
  map<unsigned,map<string, float>> pre_obs_values_by_year_;
  map<unsigned,map<string, float>> obs_values_by_year_;
  vector<float>                  error_values_;
  string                          catchability_label_;
  //Catchability*                   catchability_ = nullptr;
  vector<string>                  selectivity_labels_;
  vector<Selectivity*>            selectivities_;
  string                          time_step_label_ = "";
  vector<string>                  cells_;
  layers::CategoricalLayer*       layer_ = nullptr;
  string                          layer_label_;
  WorldView*                       world_ = nullptr;
  bool                             selectivity_length_based_ = false;
  float                           time_step_proportion_;

};

} /* namespace observations */
} /* namespace niwa */

#endif /* AGE_OBSERVATIONS_BIOMASS_H_ */
