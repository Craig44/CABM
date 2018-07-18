/**
 * @file Binomial.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 25/03/2013
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers
#include "Binomial.h"

#include "Utilities/Math.h"
#include "Utilities/RandomNumberGenerator.h"

// Namespaces
namespace niwa {
namespace likelihoods {

namespace math = niwa::utilities::math;


/**
 * Simulate observed values
 *
 * @param comparisons A collection of comparisons passed by the observation
 */
void Binomial::SimulateObserved(map<unsigned, vector<observations::Comparison> >& comparisons) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  float error_value = 0.0;
  auto iterator = comparisons.begin();
  for (; iterator != comparisons.end(); ++iterator) {
    LOG_FINE() << "Simulating values for year: " << iterator->first;
    for (observations::Comparison& comparison : iterator->second) {
      error_value = ceil(comparison.error_value_);

      if (comparison.expected_ <= 0.0 || error_value <= 0.0)
        comparison.simulated_ = 0.0;
      else
        comparison.simulated_ = rng.binomial(comparison.expected_, error_value) / error_value;
    }
  }
}

} /* namespace likelihoods */
} /* namespace niwa */
