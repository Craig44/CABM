/**
 * @file MovementPreference.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 26/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 * User defines a matrix of probabilities with a selectivity, that describes the probability of moving to another cell.
 *
 *
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MOVEMENT_PREFERENCE_H_
#define SOURCE_PROCESSES_CHILDREN_MOVEMENT_PREFERENCE_H_

// headers
#include "Processes/Process.h"

#include "PreferenceFunctions/PreferenceFunction.h"
#include "Layers/Children/NumericLayer.h"

// namespaces
namespace niwa {
//class Selectivity;
namespace processes {

/**
 * A movement struct that stores movement information
 */


/**
 * Class definition
 */
class MovementPreference : public Process {
public:
  // methods
  explicit MovementPreference(Model* model);
  virtual                     ~MovementPreference() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final;
  void                        FillReportCache(ostringstream& cache);
protected:
  // Methods
  void                        calculate_gradients();
  void                        calculate_diffusion_parameter(float& preference_value, float& diffusion_parameter);  // zeta and d_max are assumed to be available to this function

  //Members
  vector<unsigned>            years_;
//  string                      selectivity_label_;
//  Selectivity*                selectivity_ = nullptr;
  vector<string>              origin_cell_;
  vector<unsigned>            origin_rows_;
  vector<unsigned>            origin_cols_;

  vector<unsigned>            possible_rows_;
  vector<unsigned>            possible_cols_;
  bool                        brownian_motion_;

  // Report containers
  float                       diffusion_parameter_;
  float                       d_max_;
  float                       zeta_;
  vector<string>              preference_function_labels_;
  vector<string>              preference_layer_labels_;

  vector<PreferenceFunction*> preference_functions_;
  vector<layers::NumericLayer*> preference_layers_;

  map<unsigned, vector<vector<float>>> meridonal_gradient_;
  map<unsigned, vector<vector<float>>> zonal_gradient_;
  map<unsigned, vector<vector<float>>> preference_by_year_;
  vector<vector<float>>                initialisation_meridonal_gradient_;
  vector<vector<float>>                initialisation_zonal_gradient_;
  vector<vector<float>>                initialisation_preference_value_;

  // objects for thread safety of rng
  vector<float>                       lat_random_numbers_;
  vector<float>                       lon_random_numbers_;
  unsigned                            n_agents_;
  vector<vector<unsigned>>            cell_offset_;

};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MOVEMENT_PREFERENCE_H_ */
