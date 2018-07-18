//============================================================================
// Name        : Dirichlet.cpp
// Author      : C.Marsh
// Date        : 21/07/2015
// Copyright   : Copyright NIWA Science ©2009 - www.niwa.co.nz
// Description :
//============================================================================

// Global headers
#include <cmath>

// Local headers
#include "Dirichlet.h"
#include <cmath>
#include <set>

#include "Utilities/Math.h"
#include "Utilities/RandomNumberGenerator.h"


// Namespaces
namespace niwa {
namespace likelihoods {

using std::set;
namespace math = niwa::utilities::math;


/**
 * Simulate observed values
 *
 * @param comparisons A collection of comparisons passed by the observation
 */

void Dirichlet::SimulateObserved(map<unsigned, vector<observations::Comparison> >& comparisons) {
  // instance the random number generator
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  float totals = 0.0;

  auto iterator = comparisons.begin();
  for (; iterator != comparisons.end(); ++iterator) {
    LOG_FINE() << "Simulating values for year: " << iterator->first;
    for (observations::Comparison& comparison : iterator->second) {
      float error_value = comparison.error_value_;
      if (comparison.expected_ <= 0.0 || error_value <= 0.0)
        comparison.simulated_ = 0.0;
      else
        comparison.simulated_ = rng.gamma(comparison.expected_ * error_value);

      totals += comparison.simulated_;
    }

    for (observations::Comparison& comparison : iterator->second)
      comparison.simulated_ /= totals;
  }
}

} /* namespace likelihoods */
} /* namespace niwa */
