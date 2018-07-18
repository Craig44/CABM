/**
 * @file Observation.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 6/03/2013
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * The time class represents a moment of time.
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef OBSERVATION_H_
#define OBSERVATION_H_

// Headers
#include "BaseClasses/Executor.h"
#include "Likelihoods/Likelihood.h"
#include "Observations/Comparison.h"
#include "Utilities/Types.h"

// Namespaces
namespace niwa {
namespace obs = niwa::observations;

class Model;


/**
 * Class Definition
 */
class Observation : public niwa::base::Executor {
public:
  // methods
  Observation() = delete;
  explicit Observation(Model* model);
  virtual                     ~Observation() = default;
  void                        Validate();
  void                        Build();
  void                        Reset();

  // pure methods
  virtual void                DoValidate() = 0;
  virtual void                DoBuild() = 0;
  virtual void                DoReset() = 0;
  virtual void                Simulate() = 0;
  virtual bool                HasYear(unsigned year) const = 0;

  // accessors
  string&                       likelihood() { return simulation_likelihood_label_; }
  vector<obs::Comparison>&      comparisons(unsigned year) { return comparisons_[year]; }
  map<unsigned, vector<obs::Comparison> >& comparisons() { return comparisons_; }

protected:
  // methods
  void                        SaveComparison(unsigned age, float length, float expected, float simulated, float error_value);
  // members
  Model*                      model_ = nullptr;
  float                       proportion_of_time_ = 0;
  string                      simulation_likelihood_label_ = "";
  Likelihood*                 likelihood_ = nullptr;
  vector<string>              allowed_likelihood_types_;
  map<unsigned, vector<obs::Comparison> > comparisons_;

};
} /* namespace niwa */
#endif /* OBSERVATION_H_ */
