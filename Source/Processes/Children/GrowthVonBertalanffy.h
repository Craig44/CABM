/**
 * @file GrowthVonBertalanffy.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process is responsible for growing individuals.
 */
#ifndef SOURCE_PROCESSES_CHILDREN_GROWTH_VON_BERTALANFFY_H_
#define SOURCE_PROCESSES_CHILDREN_GROWTH_VON_BERTALANFFY_H_

// headers
#include "Processes/Process.h"

#include "Layers/Children/Numeric/Base/NumericLayer.h"
// namespaces
namespace niwa {
namespace processes {

class NumericLayer;
/**
 * Class definition
 */
class GrowthVonBertalanffy : public Process {
public:
  // methods
  explicit GrowthVonBertalanffy(Model* model);
  virtual                     ~GrowthVonBertalanffy() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final { };
protected:
  std::string				  l_inf_layer_label_;
  std::string				  k_layer_label_;
  niwa::layers::NumericLayer* L_inf_layer_;
  niwa::layers::NumericLayer* k_layer_;  
  vector<double>			  time_step_proportions_;
  double					  cv_;
  string					  distribution_;

  vector<double>			  get_growth_params_by_cell(unsigned row, unsigned col);
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_GROWTH_VON_BERTALANFFY_H_ */
