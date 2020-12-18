/**
 * @file Ageing.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 15/10/2020
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process does nothing. It's used for debugging time steps
 */
#ifndef SOURCE_PROCESSES_CHILDREN_AGEING_H_
#define SOURCE_PROCESSES_CHILDREN_AGEING_H_

// headers
#include "Processes/Process.h"

// namespaces
namespace niwa {
namespace processes {

/**
 * Class definition
 */
class Ageing : public Process {
public:
  // methods
  explicit Ageing(Model* model);
  virtual                     ~Ageing() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final { };
  void                        DoReset() override final { };
  void                        DoExecute() override final;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_AGEING_H_ */
