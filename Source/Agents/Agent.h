/**
 * @file Agent.h
 * @author  C.Marsh
 * @version 1.0
 * @date 11/06/2018
 * @section LICENSE
 *
 * Copyright
 *
 * @section The Agent class has describes all the agents characteristics, these are manipulated by processes but controlled by their own parameters.
 *
 */
#ifndef AGENT_H_
#define AGENT_H_

// Headers
#include "BaseClasses/Object.h"

// Namespaces
namespace niwa {

/**
 * Class Definition
 */
class Agent { // Don't make this inherit from BaseClasses/Object.h
public:
  // Methods
  virtual                       ~Agent() = default;
  Agent(float lat, float lon, float first_age_length_par, float second_age_length_par, float M, unsigned birth_year, float first_length_weight_par, float second_length_weigth_par, Model* model);
  virtual void                  Reset() {};
  // Accessors
  unsigned                     age();
  virtual const bool&          is_mature() const {return mature_ ;};
  virtual const float&         get_scalar() const {return scalar_;};
  virtual const float&         get_weight() const {return weight_;};
  virtual const float&         get_length() const {return length_;};
  virtual const float&         get_m() const {return M_;};
  virtual const float&         get_first_age_length_par() const {return first_age_length_par_;};
  virtual const float&         get_second_age_length_par() const {return second_age_length_par_;};
  virtual const float&         get_first_length_weight_par() const {return first_length_weight_par_;};
  virtual const float&         get_second_length_weight_par() const {return second_length_weight_par_;};

  void                         set_length(float new_length) {length_ = new_length;}
  void                         set_weight(float new_weight) {weight_ = new_weight;}
  void                         set_scalar(float scalar) {scalar_ = scalar;}
  void                         set_maturity(bool mature) {mature_ = mature;}



  //Dynamices


protected:
  // Methods
  void                        growth_init();
  // Members
  float                       lat_;  // Current location
  float                       lon_;
  float                       first_age_length_par_;  // L_inf for von bert
  float                       second_age_length_par_; // k for von bert
  float                       M_; // natural mortality
  unsigned                    birth_year_;
  bool                        mature_ = false;
  bool                        sex_; // 1 = male, 0 = female TODO
  float                       scalar_ = 1.0;
  float                       length_ = 0.0;
  float                       weight_ = 1.0;
  float                       first_length_weight_par_;   // a
  float                       second_length_weight_par_;   // b
  // TODO link an agent to its home for natal homing dynamics
  Model*                      model_ = nullptr;

private:

  //


};
} /* namespace niwa */

#endif /* AGENT_H_ */
