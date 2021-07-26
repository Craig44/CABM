/**
 * @file GrowthSchnuteWithBasic.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process is responsible for growing individuals length via the Schnute growth relationship and weight by the basic w = a * L^b method, so variaility is in the von-bert function not, a and b
 */
#ifndef SOURCE_PROCESSES_CHILDREN_GROWTH_SCHNUTE_WITH_BASIC_H_
#define SOURCE_PROCESSES_CHILDREN_GROWTH_SCHUNTE_WITH_BASIC_H_

// headers
#include "Processes/Children/Growth.h"

#include "Layers/Children/NumericLayer.h"
#include "Agents/Agent.h"
// namespaces
namespace niwa {
namespace processes {

class NumericLayer;
/**
 * Class definition
 */
class GrowthSchnuteWithBasic : public Growth {
public:
  // methods
  explicit GrowthSchnuteWithBasic(Model* model);
  virtual                     ~GrowthSchnuteWithBasic() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final;
  void                        RebuildCache() override final;
  void                        FillReportCache(ostringstream& cache) override final;
  void                        ApplyStochasticGrowth(vector<Agent>& agents);

  void                        draw_growth_param(unsigned row, unsigned col, unsigned number_of_draws, vector<vector<float>>& vec, unsigned sex) override final;

protected:

  vector<string>				              alpha_layer_label_;
  vector<string>				              beta_layer_label_;
  vector<string>                      t0_layer_label_;
  vector<niwa::layers::NumericLayer*> alpha_layer_;
  vector<niwa::layers::NumericLayer*> beta_layer_;
  vector<niwa::layers::NumericLayer*> t0_layer_;
  vector<string>                      a_layer_label_;
  vector<string>                      b_layer_label_;
  vector<niwa::layers::NumericLayer*> a_layer_;
  vector<niwa::layers::NumericLayer*> b_layer_;
  vector<double>                       tau1_;
  vector<double>                       tau2_;
  vector<double>                       y1_;
  vector<double>                       y2_;
  vector<double>                       beta_;
  vector<double>                       alpha_;
  vector<double>                       a_;
  vector<double>                       b_;
  vector<double>                       t0_;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_GROWTH_SCHNUTE_H_ */
