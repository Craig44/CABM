/**
 * @file Agent.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 17/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "Agent.h"

#include <boost/algorithm/string/join.hpp>

#include "Model/Managers.h"
#include "Model/Model.h"
#include "Agents/Agent.h"

// namespaces
namespace niwa {
namespace reports {

/**
 * Default constructor
 *
 * @param model Pointer to the current model context
 */
Agent::Agent(Model* model) : Report(model) {
  model_state_ = State::kIterationComplete;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);

  parameters_.Bind<unsigned>(PARAM_NUMBER_OF_AGENTS, &n_agents, "Number of agents to summarise", "");
}

/**
 * Build our relationships between this object and other objects
 */
void Agent::DoBuild() {

}

/**
 * Execute this report
 */
void Agent::DoExecute() {
  LOG_FINE() <<" printing report " << label_;
  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";

  ready_for_writing_ = true;
}


} /* namespace reports */
} /* namespace niwa */
