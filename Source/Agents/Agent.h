/**
 * @file Agent.h
 * @author  C.Marsh
 * @version 1.0
 * @date 11/06/2018
 * @section LICENSE
 *
 * Copyright
 *
 * @section The Agent class has describes all the agents characteristics, and has the methods that control a single agents fate. For example natural mortality and movement.
 * Other processes such as recruitment are defined at the partition level.
 *
 */
#ifndef AGENT_H_
#define AGENT_H_

// Headers
#include "BaseClasses/Object.h"

// Namespaces
namespace niwa {
class Model;

/**
 * Class Definition
 */
class Agent : public niwa::base::Object {
public:
  // Methods
  virtual                       ~Agent() = default;
  Agent()                       = delete;
  explicit                      Agent(Model* model);
  void                          Validate();
  void                          Build();
  void                          Reset() { };


  // Accessors

protected:
  // Methods

  // Members
  Model*                      model_ = nullptr;

};
} /* namespace niwa */

#endif /* AGENT_H_ */
