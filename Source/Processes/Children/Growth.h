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
 * Growth in this sense is both length and weight of an agent, any child growth class should inherit this as a parent
 */
#ifndef SOURCE_PROCESSES_CHILDREN_GROWTH_H_
#define SOURCE_PROCESSES_CHILDREN_GROWTH_H_

// headers
#include "Processes/Process.h"
#include "Utilities/Distribution.h"
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
  virtual void                        DoBuild();
  virtual void                        DoReset() { };
  virtual void                        DoExecute() { };

  virtual void                        draw_growth_param(unsigned row, unsigned col, unsigned number_of_draws, vector<vector<float>>& vec, unsigned sex) = 0;
  bool                                update_growth() {return update_growth_parameters_;};

protected:
  vector<double>			  time_step_proportions_;
  double					      cv_;
  string					    distribution_label_;
  bool                update_growth_parameters_;
  map<unsigned, double>   time_step_ratios_;
  Distribution        distribution_ = Distribution::kNone;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_GROWTH_VON_BERTALANFFY_H_ */
