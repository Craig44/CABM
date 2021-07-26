/**
 * @file DerivedQuantity.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 6/06/2013
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * A derived quantity is a value that is derived from the model
 * at some point during execution. This can be anything
 * from the amount of objects in the partition, to the weight
 * of objects in the partition etc.
 */
#ifndef DERIVEDQUANTITY_H_
#define DERIVEDQUANTITY_H_

// headers
#include "BaseClasses/Executor.h"

// namespaces
namespace niwa {

class Model;
class Selectivity;
class WorldView;

// classes
class DerivedQuantity : public niwa::base::Executor {
public:
  // methods
  DerivedQuantity() = delete;
  explicit DerivedQuantity(Model* model);
  virtual                     ~DerivedQuantity() = default;
  void                        Validate();
  void                        Build();
  void                        Reset();
  double                      GetValue(unsigned year);
  double                      GetValue(unsigned year, unsigned row, unsigned col);
  double                      GetInitialisationValue(unsigned phase = 0, unsigned index = 0);
  double                      GetInitialisationValue(unsigned row, unsigned col, unsigned phase, unsigned index);

  double                      GetLastValueFromInitialisation(unsigned phase);
  double                      GetLastValueFromInitialisation(unsigned phase, unsigned row, unsigned col);

  // pure methods
  virtual void                DoValidate() = 0;
  virtual void                DoBuild() = 0;

  // accessors
  const string&               time_step() { return time_step_label_; }
  vector<vector<double> >&     initialisation_values() { return initialisation_values_; }
  const map<unsigned, double>& values() { return values_; }
  bool                        is_spatial() {return spatial_;};
protected:
  // Members
  Model*                      model_ = nullptr;
  WorldView*                  world_ = nullptr;
  string                      time_step_label_ = "";
  unsigned                    current_initialisation_phase_ = 0;
  vector<vector<double>>       initialisation_values_;
  vector<vector<vector<vector<double>>>>  initialisation_values_by_space_;  //[phase][row][col][value]

  map<unsigned, double>        values_;
  map<unsigned, vector<vector<double>>> values_by_space_;  // year x row x col
  double                       cache_value_;
  double                       value_;

  string                      proportion_method_;
  double                       time_step_proportion_;
  bool                        spatial_ = false;

  // objects for thread safety of rng
  vector<double>                       random_numbers_;
  unsigned                            n_agents_;
  vector<vector<unsigned>>            cell_offset_;


};
} /* namespace niwa */
#endif /* DERIVEDQUANTITY_H_ */
