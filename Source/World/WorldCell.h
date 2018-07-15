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
  void                        Build(unsigned row, unsigned col, float lat, float lon, unsigned min_age, unsigned max_age);
  void                        Reset();
  void                        set_enabled(bool enabled) {enabled_ = enabled; };
  bool                        is_enabled() {return enabled_; };
  void                        set_area(float area) {area_ = area;}
  void                        seed_agents(unsigned number_agents_to_seed, const vector<double>&  mort_par, const vector<vector<double>>&  growth_par, const double& seed_z);
  list<Agent>&                get_agents() {return agents_;};
  void                        get_age_frequency(vector<unsigned>& age_freq);
  float                       get_abundance();
  float                       get_biomass();


protected:
  // Methods

  // Members
  list<Agent>                  agents_;
  bool                         enabled_;
  float                        area_;
  unsigned                     row_;
  unsigned                     col_;
  // These are model attributes but because I don't know how to give the constructor parameters
  // when its being build as an array by the WorldView I can't give this class a model pointer, pretty annoying.
  float                        lon_ = 0.0;
  float                        lat_ = 0.0;
  unsigned                     min_age_;
  unsigned                     max_age_;
  unsigned                     age_spread_;

};

} /* namespace niwa */
#endif /* PARTITION_H_ */
