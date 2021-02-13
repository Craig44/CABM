/**
 * @file MortalityCull.h
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
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_CULL_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_CULL_H_

// headers
#include "Processes/Children/Mortality.h"

#include "Layers/Children/IntLayer.h"
#include "Agents/Agent.h"
#include <omp.h>

// namespaces
namespace niwa {
class Selectivity;
namespace processes {
using std::string;
/**
 * Class definition
 */
class MortalityCull : public Mortality {
public:
  // methods
  explicit MortalityCull(Model* model);
  virtual                     ~MortalityCull() = default;
  virtual void                        DoValidate() override final { };
  virtual void                        DoBuild() override final;
  virtual void                        DoReset() override final;
  virtual void                        DoExecute() override final;
  void                                draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) override final;
  void                                FillReportCache(ostringstream& cache) override final;
  void                                RebuildCache() override final;
  void                                ApplyStochasticMortality(vector<Agent>& agents);


protected:
  string                              layer_label_;
  layers::IntLayer*                   layer_ = nullptr;
  vector<unsigned>                    years_;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_CULL_H_ */
