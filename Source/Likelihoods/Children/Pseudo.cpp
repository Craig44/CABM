/**
 * @file Pseudo.cpp
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
#include "Pseudo.h"

// Namespaces
namespace niwa {
namespace likelihoods {


/**
 * Simulate observed values
 *
 * @param comparisons A collection of comparisons passed by the observation
 */
void Pseudo::SimulateObserved(map<unsigned, map<string, vector<observations::Comparison> > >& comparisons) {
  auto iterator = comparisons.begin();
  for (; iterator != comparisons.end(); ++iterator) {
    LOG_FINE() << "Simulating values for year: " << iterator->first;
    for (auto second_iter = iterator->second.begin(); second_iter != iterator->second.end(); ++second_iter) {
      LOG_FINE() << "Simulating values for cell: " << second_iter->first;
      for (observations::Comparison& comparison : second_iter->second) {
        comparison.simulated_ = 0.0;
      }
    }
  }
}


} /* namespace likelihoods */
} /* namespace niwa */
