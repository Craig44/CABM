/**
 * @file TagShedding.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * Tagged fish will loose tag's as the model continues, so this transfers tagged fish to untagged partition based on some rate
 */
#ifndef SOURCE_PROCESSES_CHILDREN_TAG_LOSS_H_
#define SOURCE_PROCESSES_CHILDREN_TAG_LOSS_H_

// headers
#include "Processes/Process.h"

#include "Layers/Children/NumericLayer.h"
#include "Agents/Agent.h"
#include <omp.h>

// namespaces
namespace niwa {
class Selectivity;
namespace processes {
/**
 * Class definition
 */
class TagShedding : public Process {
public:
  // methods
  explicit TagShedding(Model* model);
  virtual                     ~TagShedding() = default;
  virtual void                        DoValidate() override final { };
  virtual void                        DoBuild() override final;
  virtual void                        DoReset() override final;
  virtual void                        DoExecute() override final;
  void                                FillReportCache(ostringstream& cache) override final;
  void                                RebuildCache() override final;

protected:
  vector<string>                      selectivity_label_;
  vector<Selectivity*>                selectivity_;
  vector<float>                       shedding_rate_;
  vector<unsigned>                    tag_release_year_;
  vector<float>                       time_step_proportions_;
  vector<float>                       ratios_;

  map<unsigned, float>                time_step_ratios_;
  unsigned                            agents_removed_ = 0;
  bool 								  selectivity_length_based_;
  // For reporting
  map<unsigned, unsigned>             tags_loss_by_year_;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_TAG_LOSS_H_ */
