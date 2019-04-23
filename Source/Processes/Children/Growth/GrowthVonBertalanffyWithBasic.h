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
 * This process is responsible for growing individuals length via the von-bert relationship and weight by the basic w = a * L^b method, so variaility is in the von-bert function not, a and b
 */
#ifndef SOURCE_PROCESSES_CHILDREN_GROWTH_VON_BERTALANFFY_WITH_BASIC_H_
#define SOURCE_PROCESSES_CHILDREN_GROWTH_VON_BERTALANFFY_WITH_BASIC_H_

// headers
#include "Processes/Children/Growth.h"

#include "Layers/Children/NumericLayer.h"

// namespaces
namespace niwa {
namespace processes {

class NumericLayer;
/**
 * Class definition
 */
class GrowthVonBertalanffyWithBasic : public Growth {
public:
  // methods
  explicit GrowthVonBertalanffyWithBasic(Model* model);
  virtual                     ~GrowthVonBertalanffyWithBasic() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final;
  void                        RebuildCache() override final;
  void                        FillReportCache(ostringstream& cache) override final;

  void                        draw_growth_param(unsigned row, unsigned col, unsigned number_of_draws, vector<vector<float>>& vec, unsigned sex) override final;

protected:
  vector<string>				              l_inf_layer_label_;
  vector<string>				              k_layer_label_;
  vector<string>                      t0_layer_label_;
  vector<niwa::layers::NumericLayer*> L_inf_layer_;
  vector<niwa::layers::NumericLayer*> k_layer_;
  vector<niwa::layers::NumericLayer*> t0_layer_;
  vector<string>                      a_layer_label_;
  vector<string>                      b_layer_label_;
  vector<niwa::layers::NumericLayer*> a_layer_;
  vector<niwa::layers::NumericLayer*> b_layer_;
  vector<float>                       l_inf_;
  vector<float>                       k_;
  vector<float>                       a_;
  vector<float>                       b_;
  vector<float>                       t0_;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_GROWTH_VON_BERTALANFFY_H_ */
