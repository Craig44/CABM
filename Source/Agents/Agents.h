/**
 * @file Agents.h
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
#ifndef AGENTS_H_
#define AGENTS_H_

// Headers
#include "BaseClasses/Object.h"

// Namespaces
namespace niwa {
class Model;

/**
 * Class Definition
 */
class Agents : public niwa::base::Object {
  friend class Model;
public:
  // Methods
  virtual                       ~Agents() = default;
  void                          Validate();
  void                          Build();
  void                          Reset() { };


  // Accessors

protected:
  // Methods
  Agents() = delete;
  explicit Agents(Model* model);


  // Members
  Model*                      model_ = nullptr;

};
} /* namespace niwa */

#endif /* AGENTS_H_ */
