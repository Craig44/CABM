/**
 * @file Managers.h
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 11/08/2015
 * @section LICENSE
 *
 * Copyright NIWA Science ©2015 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This class holds accessors for the managers in the application. I've gone down
 * this path because I want to be able to easily mock them for use in the unit
 * test suite.
 */
#ifndef SOURCE_MODEL_MANAGERS_H_
#define SOURCE_MODEL_MANAGERS_H_

// namespaces
namespace niwa {

// forward decs
namespace ageingerrors { class Manager; }
namespace asserts { class Manager; }
namespace derivedquantities { class Manager; }
namespace initialisationphases { class Manager; }
namespace layers { class Manager; }
namespace likelihoods { class Manager; }
namespace observations { class Manager; }
namespace preference_functions { class Manager; }
namespace processes { class Manager; }
namespace reports { class Manager; }
namespace selectivities { class Manager; }
namespace timesteps { class Manager; }
namespace minimisers { class Manager; }

class Estimables;
class Model;

// classes
class Managers {
  friend class Model;
  friend class MockManagers;
public:
  // accessors
  virtual ageingerrors::Manager*          ageing_error() { return ageing_error_; }
  virtual asserts::Manager*               assertx() { return assert_; }
  virtual derivedquantities::Manager*     derived_quantity() { return derived_quantity_; }
  virtual initialisationphases::Manager*  initialisation_phase() { return initialisation_phase_; }
  virtual layers::Manager*                layer() { return layer_; }
  virtual likelihoods::Manager*           likelihood() { return likelihood_; }
  virtual observations::Manager*          observation() { return observation_; }
  virtual preference_functions::Manager*  preference_function() { return preference_function_; }
  virtual processes::Manager*             process() { return process_; }
  virtual reports::Manager*               report() { return report_; }
  virtual selectivities::Manager*         selectivity() { return selectivity_; }
  virtual timesteps::Manager*             time_step() { return time_step_; }
  virtual minimisers::Manager*            minimiser() { return minimiser_; }

protected:
  // methods
  Managers(Model* model);
  virtual                     ~Managers();
  void                        Validate();
  void                        Build();
  void                        Reset();
  void                        BuildPreWorldView();


  // members
  Model*                              model_;
  ageingerrors::Manager*              ageing_error_;
  asserts::Manager*                   assert_;
  derivedquantities::Manager*         derived_quantity_;
  initialisationphases::Manager*      initialisation_phase_;
  layers::Manager*                    layer_;
  likelihoods::Manager*               likelihood_;
  observations::Manager*              observation_;
  preference_functions::Manager*      preference_function_;
  processes::Manager*                 process_;
  reports::Manager*                   report_;
  selectivities::Manager*             selectivity_;
  timesteps::Manager*                 time_step_;
  minimisers::Manager*                minimiser_;

};

} /* namespace niwa */

#endif /* SOURCE_MODEL_MANAGERS_H_ */
