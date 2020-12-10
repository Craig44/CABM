/**
 * @file LogNormalWithQ.cpp
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
#include "LogNormalWithQ.h"

#include <cmath>

#include "Utilities/RandomNumberGenerator.h"
#include "Model/Model.h"

// Namespaces
namespace niwa {
namespace likelihoods {


/**
 * Simulate observed values
 *
 * @param comparisons A collection of comparisons passed by the observation
 */
void LogNormalWithQ::SimulateObserved(map<unsigned, map<string, vector<observations::Comparison> > >& comparisons) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();

  float error_value = 0.0;
  auto iterator = comparisons.begin();
  if (model_->run_mode() == (RunMode::Type)(RunMode::kMSE)) {
    vector<unsigned> sim_years = model_->simulation_years();
    for (; iterator != comparisons.end(); ++iterator) {
      LOG_FINE() << "Simulating values for year: " << iterator->first;
      if((iterator->first >= sim_years[0]) & (iterator->first <= sim_years[sim_years.size() - 1])) {
        for (auto second_iter = iterator->second.begin(); second_iter != iterator->second.end(); ++second_iter) {
          LOG_FINE() << "Simulating values for cell: " << second_iter->first;
          for (observations::Comparison& comparison : second_iter->second) {
            error_value = comparison.error_value_;

            if (comparison.expected_ <= 0.0 || error_value <= 0.0)
              comparison.simulated_ = 0.0;
            else
              comparison.simulated_ = rng.lognormal(comparison.expected_, error_value);
          }
        }
      }
    }
  } else {
    for (; iterator != comparisons.end(); ++iterator) {
      LOG_FINE() << "Simulating values for year: " << iterator->first;
      for (auto second_iter = iterator->second.begin(); second_iter != iterator->second.end(); ++second_iter) {
        LOG_FINE() << "Simulating values for cell: " << second_iter->first;
        for (observations::Comparison& comparison : second_iter->second) {
          error_value = comparison.error_value_;

          if (comparison.expected_ <= 0.0 || error_value <= 0.0)
            comparison.simulated_ = 0.0;
          else
            comparison.simulated_ = rng.lognormal(comparison.expected_, error_value);
        }
      }
    }
  }
}

} /* namespace likelihoods */
} /* namespace niwa */
