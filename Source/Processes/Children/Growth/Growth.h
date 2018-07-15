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
 * This process is the parent growth class, so the manager can differentiate processes
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

  virtual void                       draw_growth_param(unsigned row, unsigned col, unsigned number_of_draws, vector<vector<double>>& vec) = 0;

protected:
  vector<double>			  time_step_proportions_;
  double					  cv_;
  string					  distribution_;

};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_GROWTH_VON_BERTALANFFY_H_ */
