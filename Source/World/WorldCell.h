/**
 * @file WorldCell.h
 * @author  C.Marsh
 * @version 1.0
 * @date 12/07/2018
 * @section LICENSE
 * Description :  This class represents 1 cell in our World-Grid. Each
 *                WorldSquare maintains a list of Agents
 *                The number of 'entity' in each cell
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
  WorldCell() = default;
  virtual                     ~WorldCell() = default;
  void                        Validate();
  void                        Build();
  void                        Reset();
  void                        set_enabled(bool enabled) {enabled_ = enabled; };
  bool                        is_enabled() {return enabled_; };
  void                        set_area(float area) {area_ = area;}
  void                        seed_agents(unsigned number_agents_to_seed);
  list<Agent>&                get_agents() {return agents_;};
protected:
  // Methods

  // Members
  Model*                       model_ = nullptr;
  list<Agent>                  agents_;
  bool                         enabled_;
  float                        area_;

};

} /* namespace niwa */
#endif /* PARTITION_H_ */
