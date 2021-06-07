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
 * - TODO
 *  - ADD Slot agent function, the code is repeasted quite often.
 *  - Think about how, best to integtate Tagged partition, and perhaps expand it, larvae etc.
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
  float                       get_area() {return area_;}

  // Functions for keeping track of individuals to agents when calculating probability of encounter etc.
  void                        calculate_individuals_alive();
  void                        remove_agent_alive(float& scalar) {total_individuals_alive_-= scalar;};
  void                        add_agent_alive(float& scalar) {total_individuals_alive_ += scalar;};
  const double&               get_total_individuals_alive() {return total_individuals_alive_;}
  void                        set_total_individuals_alive(double& val) {total_individuals_alive_ = val;}

  void                        seed_agents(unsigned number_agents_to_seed, const float& seed_z, float scalar);
  void                        birth_agents(unsigned number_agents_to_birth, float scalar);
  void                        set_scalar();
  void                        get_age_frequency(vector<float>& age_freq, bool& is_age);
  void                        get_female_frequency(vector<float>& age_freq, bool& is_age);
  void                        get_male_frequency(vector<float>& age_freq, bool& is_age);
  void                        get_age_agent_frequency(vector<float>& age_freq, bool& is_age);
  void                        get_female_agent_frequency(vector<float>& age_freq, bool& is_age);
  void                        get_male_agent_frequency(vector<float>& age_freq, bool& is_age);


  float                       get_abundance();
  float                       get_biomass();
  float                       get_mature_biomass();
  void                        update_agent_parameters();
  void                        update_mortality_params();
  void                        update_growth_params();
  void                        apply_growth_time_varying();
  void                        apply_mortality_time_varying();


  vector<Agent>               agents_;
  vector<Agent>               tagged_agents_;

  string                      get_cell_label() {return cell_label_;};
  float                       get_lat() {return lat_;};
  float                       get_lon() {return lon_;};


protected:
  // Methods

  // Members
  bool                         enabled_ = true;
  float                        area_;
  unsigned                     row_;
  unsigned                     col_;
  double                       total_individuals_alive_ = 0.0;

  // These are model attributes but because I don't know how to give the constructor parameters
  // when its being build as an array by the WorldView I can't give this class a model pointer, pretty annoying.
  float                        lon_ = 0.0;
  float                        lat_ = 0.0;
  Model*                       model_ = nullptr;
  processes::Mortality*        mortality_ = nullptr;
  processes::Growth*           growth_ = nullptr;
  vector<Selectivity*>         selectivity_;
  string                       cell_label_;
};

} /* namespace niwa */
#endif /* PARTITION_H_ */
