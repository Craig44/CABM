/**
 * @file Growth.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process is the parent growth class, so the manager can differentiate processes via dynamic casts
 * Growth in this sense is both length and weight of an agent
 */
#ifndef SOURCE_PROCESSES_CHILDREN_GROWTH_H_
#define SOURCE_PROCESSES_CHILDREN_GROWTH_H_

// headers
#include "Processes/Process.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * Class definition
 */
class Growth : public Process {
public:
  // methods
  explicit Growth(Model* model);
  virtual                     ~Growth() = default;
  virtual void                        DoValidate(){ };
  virtual void                        DoBuild() { };
  virtual void                        DoReset() { };
  virtual void                        DoExecute() { };

  virtual void                        draw_growth_param(unsigned row, unsigned col, unsigned number_of_draws, vector<vector<float>>& vec) = 0;

protected:
  vector<float>			  time_step_proportions_;
  float					      cv_;
  string					    distribution_;

};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_GROWTH_VON_BERTALANFFY_H_ */
