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
#include "Processes/Children/Growth.h"
#include "Processes/Children/Mortality.h"

// Namespaces
namespace niwa {

using std::list;

class Model;
class Agent;
class Selectivity;
/**
 * Class Definition
 */
class WorldCell {
public:
  // Methods
  WorldCell() = default;
  virtual                     ~WorldCell() = default;
  void                        Validate();
  void                        Build(unsigned row, unsigned col, float lat, float lon, Model* model);
  void                        Reset();
  void                        set_enabled(bool enabled) {enabled_ = enabled; };
  bool                        is_enabled() {return enabled_; };
  void                        set_area(float area) {area_ = area;}
  void                        seed_agents(unsigned number_agents_to_seed, const float& seed_z);
  void                        birth_agents(unsigned number_agents_to_birth);
  //list<Agent>&                get_agents() {return agents_;};
  void                        get_age_frequency(vector<unsigned>& age_freq);
  float                       get_abundance();
  float                       get_biomass();
  float                       get_mature_biomass();
  list<Agent>                 agents_;


protected:
  // Methods

  // Members
  bool                         enabled_ = true;
  float                        area_;
  unsigned                     row_;
  unsigned                     col_;
  // These are model attributes but because I don't know how to give the constructor parameters
  // when its being build as an array by the WorldView I can't give this class a model pointer, pretty annoying.
  float                        lon_ = 0.0;
  float                        lat_ = 0.0;
  Model*                       model_ = nullptr;
  processes::Mortality*        mortality_ = nullptr;
  processes::Growth*           growth_ = nullptr;
  vector<Selectivity*>         selectivity_;
};

} /* namespace niwa */
#endif /* PARTITION_H_ */
