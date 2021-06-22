/**
 * @file Multinomial.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 25/03/2013
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers
#include "Multinomial.h"

#include <cmath>
#include <set>

#include "Utilities/Math.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Model/Model.h"

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
void Multinomial::SimulateObserved(map<unsigned, map<string, vector<observations::Comparison> > >& comparisons) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  auto iterator = comparisons.begin();
  double error_val = 1;
  double resid_prop_sum = 1.0;
  if (model_->run_mode() == (RunMode::Type)(RunMode::kMSE)) {
    vector<unsigned> sim_years = model_->simulation_years();
    for (; iterator != comparisons.end(); ++iterator) {
      LOG_FINE() << "Simulating values for year: " << iterator->first;
      if((iterator->first >= sim_years[0]) & (iterator->first <= sim_years[sim_years.size() - 1])) {
        for (auto second_iter = iterator->second.begin(); second_iter != iterator->second.end(); ++second_iter) {
          enter_cell_ = true;
          LOG_FINE() << "Simulating values for cell: " << second_iter->first;
          resid_prop_sum = 1.0;
          // Now simulate need to adapt error_vals_props
          for (observations::Comparison& comparison : second_iter->second) {
            // For the multinomial we need a single N per composition
            if(enter_cell_) {
              error_val = comparison.error_value_;
              enter_cell_ = false;
            }
            if (comparison.expected_ <= 0.0 || error_val <= 0.0) {
              comparison.simulated_ = 0.0;
            } else {
              if((resid_prop_sum <= 0.0) | (error_val <= 0.0)) {
                comparison.simulated_ = 0.0;
              } else {
                comparison.simulated_ = rng.binomial((comparison.expected_ / resid_prop_sum), error_val);
                error_val -= comparison.simulated_;
                resid_prop_sum -= comparison.expected_;
              }
              LOG_FINEST() << "expected = " << comparison.expected_ <<  " Simulated = " << comparison.simulated_;
            }
          }
        }
      }
    }
  } else {
    for (; iterator != comparisons.end(); ++iterator) {
      LOG_FINE() << "Simulating values for year: " << iterator->first;
      for (auto second_iter = iterator->second.begin(); second_iter != iterator->second.end(); ++second_iter) {
           enter_cell_ = true;
          LOG_FINE() << "Simulating values for cell: " << second_iter->first;
          resid_prop_sum = 1.0;
          // Now simulate need to adapt error_vals_props
          for (observations::Comparison& comparison : second_iter->second) {
            // For the multinomial we need a single N per composition
            if(enter_cell_) {
              error_val = comparison.error_value_;
              enter_cell_ = false;
            }
            if (comparison.expected_ <= 0.0 || error_val <= 0.0) {
              comparison.simulated_ = 0.0;
            } else {
              if((resid_prop_sum <= 0.0) | (error_val <= 0.0)) {
                comparison.simulated_ = 0.0;
              } else {
                comparison.simulated_ = rng.binomial((comparison.expected_ / resid_prop_sum), error_val);
                error_val -= comparison.simulated_;
                resid_prop_sum -= comparison.expected_;
              }
              LOG_FINEST() << "expected = " << comparison.expected_ <<  " Simulated = " << comparison.simulated_;
            }
          }    
       }
    }
  }
}
} /* namespace likelihoods */
} /* namespace niwa */
