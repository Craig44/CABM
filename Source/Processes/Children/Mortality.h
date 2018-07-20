/**
 * @file Mortality.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process is the parent mortlity class, so the manager can generate dynamicasts<> processes
 * any child mortality class should inherit this as a parent
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_H_

// headers
#include "Processes/Process.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * Class definition
 */
class Mortality : public Process {
public:
  // methods
  explicit Mortality(Model* model);
  virtual                     ~Mortality() = default;
  virtual void                        DoValidate(){ };
  virtual void                        DoBuild() { };
  virtual void                        DoReset() { };
  virtual void                        DoExecute() { };

  virtual void                        draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) = 0;
  map<unsigned, vector<unsigned>>&    get_removals_by_age() {return removals_by_age_;};
  map<unsigned, vector<unsigned>>&    get_removals_by_length() {return removals_by_length_;};

protected:
  map<unsigned, vector<unsigned>>     removals_by_age_;
  map<unsigned, vector<unsigned>>     removals_by_length_;


};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_H_ */
