/**
 * @file WorldCell.h
 * @author  C.Marsh
 * @version 1.0
 * @date 12/07/2018
 * @section LICENSE
 * Description :  This class represents 1 cell in our World-Grid. Each
 *                WorldSquare maintains a list of categories and ages; As well as
 *                The number of fish in each location
 *
 *
 */
#ifndef WORLD_CELL_H_
#define WORLD_CELL_H_

// Headers
#include <map>
#include <list>
#include <vector>
#include <string>
#include <memory>

#include "Agents/Agent.h"
#include "Utilities/Types.h"

// Namespaces
namespace niwa {

using std::list;

class Model;
class Agent;
/**
 * Class Definition
 */
class WorldCell {
public:
  // Methods
  WorldCell() = delete;
  explicit WorldCell(Model* model) : model_(model) { };
  virtual                     ~WorldCell() {};
  void                        Validate();
  void                        Build();
  void                        Reset();


protected:
  // Methods

  // Members
  Model*                       model_ = nullptr;
  list<Agent>                  agents_;
  list<Agent>                  cached_agents_;

};

} /* namespace niwa */
#endif /* PARTITION_H_ */
