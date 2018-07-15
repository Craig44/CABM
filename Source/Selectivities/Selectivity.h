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
using niwa::utilities::Double;

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
  virtual Double              GetResult(unsigned age_or_length);

protected:
  // pure methods
  virtual void                DoValidate() = 0;
  // Members
  Model*                      model_ = nullptr;
  bool                        length_based_ = false;
//  map<unsigned, Double>       values_;
  vector<Double>              values_;
  vector<Double>              length_values_;
  unsigned                    min_index_;
};
} /* namespace niwa */
#endif /* SELECTIVITY_H_ */
