//============================================================================
// Name        : MultinomialDirichletMultinomial.cpp
// Author      : C.Marsh
// Date        : 11/05/22
// Copyright   : Copyright NIWA Science 2020 - www.niwa.co.nz
// Description :
//============================================================================

// Global headers
#include "DirichletMultinomial.h"

#include <cmath>

// Local headers
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

DirichletMultinomial::DirichletMultinomial(Model* model) : Likelihood(model) {
  parameters_.Bind<string>(PARAM_LABEL, &label_, "Label for the Dirichlet-multinomial distribution", "");
  parameters_.Bind<string>(PARAM_TYPE, &type_, "Type of likelihood", "");
  parameters_.Bind<double>(PARAM_THETA, &theta_, "Theta parameter (account for overdispersion)", "")->set_lower_bound(0.0);
  //parameters_.Bind<string>(PARAM_OVERDISPERSION_TYPE, &theta_model_, "Is theta linear or saturated", "", PARAM_LINEAR)->set_allowed_values({PARAM_LINEAR,PARAM_SATURATED});
  RegisterAsAddressable(PARAM_THETA, &theta_);
}

/*
* Validate user hasn't supplied a label that is a type from another likelihood.
*/
void DirichletMultinomial::DoValidate() {
}


/**
 * Simulate observed values
 *
 * @param comparisons A collection of comparisons passed by the observation
 */

void DirichletMultinomial::SimulateObserved(map<unsigned, map<string, vector<observations::Comparison> > >& comparisons) {
  // instance the random number generator
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  double error_val = 1;
  double previous_error_val = 1;
  if (model_->run_mode() == (RunMode::Type)(RunMode::kMSE)) {

  } else {
    auto iterator = comparisons.begin();
    for (; iterator != comparisons.end(); ++iterator) {
      LOG_FINEST() << "Simulating values for year: " << iterator->first;
      double rng_uniform = 0.0;
      double cumulative_expect = 0.0;
      for (auto second_iter = iterator->second.begin(); second_iter != iterator->second.end(); ++second_iter) {
        LOG_FINEST() << "Simulating values for cell: " << second_iter->first << " next loop = " << second_iter->second.size();
        double total = 0.0;
        double max_neff = 0;
        for (observations::Comparison& comparison : second_iter->second) {
          comparison.simulated_ = rng.gamma(comparison.expected_ * theta_ * comparison.error_value_ );
          if(max_neff < comparison.error_value_)
            max_neff = comparison.error_value_;
          total += comparison.simulated_;
        }
        for (observations::Comparison& comparison : second_iter->second)
          comparison.simulated_ /= total;
        // Now draw from a multinomial
        enter_cell_ = true;
        //resid_prop_sum = 1.0;
        // Now simulate need to adapt error_vals_props
        vector<double> expected_vals(second_iter->second.size(), 0.0);
        unsigned counter = 0;
        for (observations::Comparison& comparison : second_iter->second) {
          expected_vals[counter] = comparison.simulated_;
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
      }
    } // for (; iterator != comparisons.end(); ++iterator) {
  } // if (model_->run_mode() == (RunMode::Type)(RunMode::kMSE)) {
} // SimulateObserved()

} /* namespace likelihoods */
} /* namespace niwa */
