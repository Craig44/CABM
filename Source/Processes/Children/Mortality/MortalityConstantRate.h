/**
 * @file MortalityConstantRate.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This is a child mortality that applies an M in each time step it is associated with.
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_CONSTANT_RATE_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_CONSTANT_RATE_H_

// headers
#include "Processes/Children/Mortality.h"

#include "Layers/Children/NumericLayer.h"
#include <omp.h>

// namespaces
namespace niwa {
class Selectivity;
namespace processes {
using std::string;
/**
 * Class definition
 */
class MortalityConstantRate : public Mortality {
public:
  // methods
  explicit MortalityConstantRate(Model* model);
  virtual                     ~MortalityConstantRate() = default;
  virtual void                        DoValidate() override final { };
  virtual void                        DoBuild() override final;
  virtual void                        DoReset() override final { };
  virtual void                        DoExecute() override final;
  void                                draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) override final;
  void                                FillReportCache(ostringstream& cache) override final;

protected:
  string                              m_layer_label_;
  layers::NumericLayer*               m_layer_ = nullptr;
  vector<string>                      selectivity_label_;
  vector<Selectivity*>                selectivity_;
  float                               m_;
  vector<float>                       time_step_proportions_;
  float                               cv_;
  string                              distribution_;
  bool                                selectivity_length_based_;
  // objects for thread safety of rng
  vector<float>                       random_numbers_;
  unsigned                            n_agents_;
  vector<vector<unsigned>>            cell_offset_;
  vector<vector<vector<float>>>       cell_offset_for_selectivity_;



  // For reporting
  map<unsigned, unsigned>             removals_by_year_;

};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_CONSTANT_RATE_H_ */
