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
class Agent { // Don't make this inherit from BaseClasses/Object.h
public:
  // Methods
  virtual                       ~Agent() = default;
  Agent(float lat, float lon, float first_age_length_par, float second_age_length_par, float M, unsigned age, float first_length_weight_par, float second_length_weigth_par);
  virtual void                  Reset() {};
  // Accessors
  unsigned                     age() const {return age_;};
  bool                         is_alive() const {return alive_ ;};
  float                        get_scalar() const {return scalar_;};
  float                        get_weight() const {return weight_;};
  float                        get_length() const {return length_;};

  //Dynamices
  void                         survival(float& selectivity);  // TODO consider moving these to the process and give processes access
  void                         maturity(float& selectivity);
protected:
  // Methods
  void                        growth_init();
  // Members
  float                       lat_;  // Current location
  float                       lon_;
  float                       first_age_length_par_;  // L_inf for von bert
  float                       second_age_length_par_; // k for von bert
  float                       survival_; // natural mortality
  unsigned                    age_;
  bool                        alive_ = true;
  bool                        mature_ = false;
  bool                        sex_; // 1 = male, 0 = female
  float                       scalar_ = 1.0;
  float                       length_ = 0.0;
  float                       weight_ = 1.0;
  float                       first_length_weight_par_;   // a
  float                       second_length_weight_par_;   // b
  // TODO link an agent to its home for natal homing dynamics


  //


};
} /* namespace niwa */

#endif /* AGENT_H_ */
