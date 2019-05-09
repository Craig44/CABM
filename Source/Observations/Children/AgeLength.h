/**
 * @file AgeLength.h
 * @author  C.Marsh
 * @date 1/08/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This observation take a random sample in a set of cells and reports the age and length of a fish, could be useful when trying to investigate age-length protocols and sources of bias
 */
#ifndef OBSERVATIONS_AGE_LENGTH_H_
#define OBSERVATIONS_AGE_LENGTH_H_

// headers
#include "Observations/Observation.h"
#include "Layers/Children/CategoricalLayer.h"


// namespaces
namespace niwa {
class Selectivity;
namespace observations {

/**
 * class definition
 */
class AgeLength : public niwa::Observation {
public:
  // methods
  AgeLength(Model* model);
  virtual                     ~AgeLength() = default;
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
  vector<string>                  selectivity_labels_;
  vector<Selectivity*>            selectivities_;
  string                          time_step_label_ = "";
  vector<string>                  cells_;
  map<string,vector<unsigned>>    col_index_for_each_cell_;
  vector<unsigned>                n_samples_;
  map<string,vector<unsigned>>    row_index_for_each_cell_;
  layers::CategoricalLayer*       layer_ = nullptr;
  string                          layer_label_;
  WorldView*                      world_ = nullptr;
  bool                            selectivity_length_based_ = false;


};

} /* namespace observations */
} /* namespace niwa */

#endif /* OBSERVATIONS_AGE_LENGTH_H_ */
