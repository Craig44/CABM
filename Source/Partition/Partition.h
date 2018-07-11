/**
 * @file World.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 18/09/2012
 * @section LICENSE
 *
 * Copyright NIWA Science �2012 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This class holds the partition for our model. It's responsible for ensuring
 * the World::Accessors can access the partition properly.
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef PARTITION_H_
#define PARTITION_H_

// Headers
#include <map>
#include <vector>
#include <string>
#include <memory>

#include "Partition/Agent.h"
#include "Utilities/Types.h"

// Namespaces
namespace niwa {
class Model;

using std::string;
using std::map;
using std::vector;

/**
 * Class Definition
 */
class Partition {
  friend class Model;
public:
  // Methods
  Partition() = delete;
  explicit Partition(Model* model) : model_(model) { };
  virtual                     ~Partition() {};
  void                        Validate();
  void                        Build();
  void                        Reset();
  // TODO add functions like Casal2 that checks to see if attributes are valid such as maturity, sex, area

  // Accessors
  vector<agents::Agent>&   get_partition() {return partition_;}; // accessor for the vector of agents by reference, this is a const so only used for summarising partition e.g SSB, observations

protected:
  // Methods


  // Members
  Model*                       model_ = nullptr;
  vector<agents::Agent>		 partition_;  // A vector of Agent instances....

};

} /* namespace niwa */
#endif /* PARTITION_H_ */
