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
 * This is the core agent class and holds all the attributes for an individual\agent over its life.
 *
 */
class Agent {
public:
  // Methods
  virtual                       ~Agent() = default;
  Agent(float lat, float lon, float first_age_length_par, float second_age_length_par, float third_age_length_par, float M, unsigned birth_year, float first_length_weight_par,
      float second_length_weigth_par, Model* model, bool mature, unsigned sex, float scalar, unsigned home_row, unsigned home_col, unsigned tag);
  // Overload the constructor with another agent, used in tagging to mirror an existing agent.
  //Agent(Agent& agent_to_copy);


  virtual void                  Reset() {};
  // Accessors
  bool                         is_alive() const {return alive_;}
  unsigned                     get_age();
  unsigned                     get_age_index();
  unsigned                     get_length_bin_index();
  virtual const bool&          get_maturity() const {return mature_ ;};
  virtual const unsigned&      get_sex() const {return sex_ ;};
  virtual float&               get_scalar() {return scalar_;};
  virtual const float&         get_weight() const {return weight_;};
  virtual const float&         get_length() const {return length_;};
  virtual const float&         get_m() const {return M_;};
  virtual const float&         get_first_age_length_par() const {return first_age_length_par_;}; // L_inf // alpha
  virtual const float&         get_second_age_length_par() const {return second_age_length_par_;}; // K // beta
  virtual const float&         get_third_age_length_par() const {return third_age_length_par_;}; // t0 // t0
  virtual const float&         get_first_length_weight_par() const {return first_length_weight_par_;}; // a // a
  virtual const float&         get_second_length_weight_par() const {return second_length_weight_par_;}; // b //n
  virtual const unsigned&      get_home_row() const {return home_row_ ;};
  virtual const unsigned&      get_home_col() const {return home_col_ ;};
  virtual const float&         get_lat() const {return lat_;};
  virtual const float&         get_lon() const {return lon_;};
  virtual const unsigned&      get_number_tags() const {return tag_;};
  virtual const unsigned&      get_tag_row() const {return tag_row_;};
  virtual const unsigned&      get_tag_col() const {return tag_col_;};
  virtual const unsigned&      get_tag_release_year() const {return tag_year_;};
  void                         shed_a_tag() {tag_--;} // anything that calls this will need to check get_number_tags() > 1, else this will cause a negative unsigned == issues

  float                        get_length_increment_since_tag()  {return length_ - length_at_tag_;};
  unsigned                     get_time_at_liberty(unsigned current_time_step)  {return current_time_step - tag_time_step_;};

  void                         set_length(float new_length) {length_ = new_length;}
  void                         set_weight(float new_weight) {weight_ = new_weight;}
  void                         set_scalar(float scalar) {scalar_ = scalar;}
  void                         set_maturity(bool mature) {mature_ = mature;}
  void                         set_growth(float scalar) {scalar_ = scalar;}
  void                         set_first_age_length_par(float value) {first_age_length_par_ = value;}
  void                         set_second_age_length_par(float value) {second_age_length_par_ = value;}
  void                         set_third_age_length_par(float value) {third_age_length_par_ = value;}
  void                         set_first_length_weight_par(float value) {first_length_weight_par_ = value;}
  void                         set_second_length_weight_par(float value) {second_length_weight_par_ = value;}
  void                         set_m(float value) {M_ = value;}
  void                         set_lat(float new_lat) {lat_ = new_lat;}
  void                         set_lon(float new_lon) {lon_ = new_lon;}
  void                         apply_tagging_event(unsigned tags, unsigned row, unsigned col);
  void                         dies() {alive_ = false;}
protected:
  // Methods
  void                        growth_init();
  // Members
  float                       lat_;  // Current location
  float                       lon_;
  float                       first_age_length_par_;  // L_inf for von bert
  float                       second_age_length_par_; // k for von bert
  float                       third_age_length_par_;  // t0 for von bert
  float                       M_; // natural mortality
  unsigned                    birth_year_;
  float                       length_ = 0.0;
  float                       weight_ = 1.0;
  float                       first_length_weight_par_;   // a
  float                       second_length_weight_par_;   // b
  // TODO link an agent to its home for natal homing dynamics
  Model*                      model_ = nullptr; // might not need
  bool                        mature_ = false;  // if maturity is in the partition, are all these agents mature?
  unsigned                    sex_ = 0; // 1 = male, 0 = female TODO
  float                       scalar_ = 1.0; // number of individuals this agent represents
  unsigned                    home_row_; // Row agent born in 
  unsigned                    home_col_; // col agent born in 
  unsigned                    tag_ = 0; // number of tags, fish can be double tagged
  unsigned					  tag_year_;
  float                       length_at_tag_ = 0.0;
  unsigned                    tag_time_step_ = 0; // number of tags, fish can be double tagged
  bool                        alive_ = true;
  unsigned                    tag_row_;			   // row agent tagged
  unsigned                    tag_col_;				// col agent was tagged
  unsigned                    life_stage_;  // 1 = larvae, 2 = immature, 3 = mature. // TODO perhaps this could be dealt with in the WorldCell level partition the agents for each cell, then processes can call
                                            // The life history of interest.

private:

  //


};
} /* namespace niwa */

#endif /* AGENT_H_ */
