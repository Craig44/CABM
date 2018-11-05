/*
 * Selectivity.h
 *
 *  Created on: 21/12/2012
 *      Author: Admin
 */

#ifndef SELECTIVITY_H_
#define SELECTIVITY_H_

// Headers
#include <map>

#include "BaseClasses/Object.h"
#include "Model/Model.h"
#include "Utilities/Types.h"

// Namespaces
namespace niwa {

// Using
using std::map;

/**
 * Class Definition
 */
class Selectivity : public niwa::base::Object {
public:
  // Methods
  Selectivity() = delete;
  explicit Selectivity(Model* model);
  virtual                     ~Selectivity() = default;
  void                        Validate();
  virtual void                Build() { RebuildCache(); };
  void                        Reset();
  virtual float               GetResult(unsigned age_or_length); // Age and or Length INDEX not the actual value
  bool                        is_length_based() {return length_based_;}
  bool                        include_zero_age_values_;

protected:
  // pure methods
  virtual void                DoValidate() = 0;
  // Members
  Model*                      model_ = nullptr;
  bool                        length_based_;
//  map<unsigned, float>       values_;
  vector<float>               values_;
  vector<float>               length_values_;
  unsigned                    min_index_;
};
} /* namespace niwa */
#endif /* SELECTIVITY_H_ */
