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
  LOG_MEDIUM();
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  auto iterator = comparisons.begin();
  double error_val = 1;
  double previous_error_val = 1;

  //double resid_prop_sum = 1.0;
  if (model_->run_mode() == (RunMode::Type)(RunMode::kMSE)) {
    vector<unsigned> sim_years = model_->simulation_years();
    double rng_uniform = 0.0;
    double cumulative_expect = 0.0;
    for (; iterator != comparisons.end(); ++iterator) {
      LOG_MEDIUM() << "Simulating values for year: " << iterator->first;
      if((iterator->first >= sim_years[0]) & (iterator->first <= sim_years[sim_years.size() - 1])) {
        for (auto second_iter = iterator->second.begin(); second_iter != iterator->second.end(); ++second_iter) {
          enter_cell_ = true;
          LOG_MEDIUM() << "Simulating values for cell: " << second_iter->first;
          //resid_prop_sum = 1.0;
          // Now simulate need to adapt error_vals_props
          vector<double> expected_vals(second_iter->second.size(), 0.0);
          unsigned counter = 0;
          for (observations::Comparison& comparison : second_iter->second) {
            expected_vals[counter] = comparison.expected_;
            if(counter == 0) {
              error_val = comparison.error_value_;
              previous_error_val = error_val;
            }
            if(previous_error_val != comparison.error_value_)
              LOG_WARNING() << "Found N for multinomial in year " << iterator->first  << " and cell " << second_iter->first << " which differend with theprevious N. When simualting from multinomial need a single N for all bins in the composition. We are just using the first N and ignoreing the others";
            
            counter++;
          }
          // Now simulate using the uniform random variable
          vector<double> sim_obs(expected_vals.size(),0.0); 
          while(error_val > 0) {
            rng_uniform = rng.chance();
            cumulative_expect = 0.0;
            for (unsigned i = 0; i < expected_vals.size(); ++i) {
              cumulative_expect += expected_vals[i];
              if (cumulative_expect >= rng_uniform) {
                sim_obs[i]++;
                break;
              }
            }
            error_val--;
          }
          // Save back into the comparisons
          counter = 0;
          for (observations::Comparison& comparison : second_iter->second) {
            comparison.simulated_ = sim_obs[counter];
            counter++;
          }

          /*
          for (observations::Comparison& comparison : second_iter->second) {
            // For the multinomial we need a single N per composition
            if(enter_cell_) {
              error_val = comparison.error_value_;
              previous_error_val = comparison.error_value_;
              enter_cell_ = false;
            }
            if(previous_error_val != comparison.error_value_)
              LOG_WARNING() << "Found N for multinomial in year " << iterator->first  << " and cell " << second_iter->first << " which differend with theprevious N. When simualting from multinomial need a single N for all bins in the composition. We are just using the first N and ignoreing the others";
            
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
          */
        }
      }
    }
  } else {
    double rng_uniform = 0.0;
    double cumulative_expect = 0.0;
    for (; iterator != comparisons.end(); ++iterator) {
      LOG_MEDIUM() << "Simulating values for year: " << iterator->first;
      for (auto second_iter = iterator->second.begin(); second_iter != iterator->second.end(); ++second_iter) {
          enter_cell_ = true;
          LOG_MEDIUM() << "Simulating values for cell: " << second_iter->first;
          //resid_prop_sum = 1.0;
          // Now simulate need to adapt error_vals_props
          vector<double> expected_vals(second_iter->second.size(), 0.0);
          unsigned counter = 0;
          for (observations::Comparison& comparison : second_iter->second) {
            expected_vals[counter] = comparison.expected_;
            if(counter == 0) {
              error_val = comparison.error_value_;
              previous_error_val = error_val;
            }
            if(previous_error_val != comparison.error_value_)
              LOG_WARNING() << "Found N for multinomial in year " << iterator->first  << " and cell " << second_iter->first << " which differend with theprevious N. When simualting from multinomial need a single N for all bins in the composition. We are just using the first N and ignoreing the others";
            
            counter++;
          }
          // Now simulate using the uniform random variable
          vector<double> sim_obs(expected_vals.size(),0.0); 
          LOG_MEDIUM() << "error value: " << error_val;

          while(error_val > 0) {
            rng_uniform = rng.chance();
            cumulative_expect = 0.0;
            for (unsigned i = 0; i < expected_vals.size(); ++i) {
              cumulative_expect += expected_vals[i];
              if (cumulative_expect >= rng_uniform) {
                sim_obs[i]++;
                break;
              }
            }
            error_val--;
          }
          // Save back into the comparisons
          counter = 0;
          for (observations::Comparison& comparison : second_iter->second) {
            comparison.simulated_ = sim_obs[counter];
            counter++;
          }
          /*
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
          */
       }
    }
  }
}
} /* namespace likelihoods */
} /* namespace niwa */
