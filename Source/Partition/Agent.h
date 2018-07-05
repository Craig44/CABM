/**
 * @file Agents.h
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 18/02/2015
 * @section LICENSE
 *
 * Copyright NIWA Science ©2015 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This object represents one category in the model.
 * Because stuff in this is built by the partition there is no need
 * to hide the members behind accessors.
 */
#ifndef PARTITION_AGENT_H_
#define PARTITION_AGENT_H_

// headers
#include <map>
#include <vector>
#include <string>

#include "Utilities/Types.h"
#include "Selectivities/Selectivity.h"

// namespaces
namespace niwa {
class AgeLength;
class LengthWeight;
class Model;

namespace partition {
using std::string;
using std::map;
using std::vector;
using std::pair;
using niwa::utilities::Double;

/**
 * Class definition
 */
class Agent {
public:
  // methods
  Agent(Model* model) : model_(model) { };
  virtual                     ~Agent() = default;

  // members
  unsigned                    age_ = 0;
  float                       length_ = 0;
  unsigned                    area_ = 0;
  bool                        alive_ = true;
  bool                        mature_ = false;

private:
  // members
  Model*                      model_ = nullptr;
};

} /* namespace partitions */
} /* namespace niwa */

#endif /* PARTITION_AGENT_H_ */
