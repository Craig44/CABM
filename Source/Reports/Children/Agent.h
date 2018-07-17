/**
 * @file Agent.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 17/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This report will take a sample of agents in the world and summarise their attributes, lengths, weights, age and maturity.
 */
#ifndef SOURCE_REPORTS_CHILDREN_AGENT_H_
#define SOURCE_REPORTS_CHILDREN_AGENT_H_

// headers
#include "Reports/Report.h"

// namespaces
namespace niwa {
class Agent;
namespace reports {

// classes
class Agent : public Report {
public:
  Agent(Model* model);
  virtual                     ~Agent() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final;
  void                        DoExecute() override final;


private:
  unsigned                    n_agents;
  bool                        first_run_ = true;

};

} /* namespace reports */
} /* namespace niwa */

#endif /* SOURCE_REPORTS_CHILDREN_AGENT_H_ */
