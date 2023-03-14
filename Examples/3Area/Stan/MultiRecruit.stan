/*
 * @file AgeStructuredModel_1.0.stan
 * @author C. Marsh
 *
 * A stan module for a simple age structured dynamics model
 * vector<> is a column vector!! remember this when structuring matrices used in multiplication
 */
 
 
functions {
  // logistic ogive function
  vector logistic_ogive(vector ages, real sel_50, real sel_95) {
    int n_ages = rows(ages);
    vector[n_ages] logis;
    for (age in 1:n_ages) {
      logis[age] = 1 / (1 + pow(19, (sel_50 - ages[age]) / sel_95));
    }
    return logis;
  }
  
  // Casal's lognormal prior calculation
  real lnorm_prior(real X, real Expectation, real sigma) {
    real result = log(X) + log(sigma) + 0.5 * (log(X/ Expectation) / sigma + sigma * 0.5)^2;
  return result;
  }
  
  // Multinomial objective calculation, from CASAL
  // returns negative log-likelihood.
  real obj_multinomal_casal_lpdf(matrix obs, matrix Exp, vector error, int[] indicator) {
	  real obj = 0.0;
	  for (y in 1:cols(Exp)) {
	    if (indicator[y] == 1) {
	      obj -= lgamma(error[y] + 1.0); // lgamma(n + 1) = lfactorial(n)
        obj += sum(lgamma((error[y] * obs[,y]) + 1.0) - error[y] * obs[,y] .* log(Exp[,y]));
	    }
	  }
    return obj;
  }
  
  // negative log likelihood for the lognormal distribution
  real obj_lognomal_lpdf(vector observed, vector expected, vector error, int[] indicator) {
    int years = rows(observed);
    vector[years] sigma;
    vector[years] score;
    real neg_ll = 0.0;
    for (y in 1:years) {
	    if (indicator[y] == 1) {
        sigma[y] = sqrt(log(1.0 + error[y] * error[y]));
        score[y] = log(observed[y] / expected[y])/ sigma[y] + 0.5 * sigma[y];
        neg_ll += log(sigma[y]) + 0.5 * score[y] * score[y];
	    }
    }
    return neg_ll;
  }  
  // Weight at age
  vector mean_weight_at_age(vector L_a, real a, real b) {
    int n_ages = rows(L_a);
    vector[n_ages] W_a;
    for (age in 1:n_ages) {
      W_a[age] = a*pow(L_a[age],b);
    }
    return W_a;
  }
  
  // Length at age
  vector VonBertalanffy(vector ages, real L_inf, real k, real t0) {
    int n_ages = rows(ages);
    vector[n_ages] L_a;
    for (age in 1:n_ages) {
      L_a[age] = L_inf * (1-exp(-k*(age -t0)));
    }
    return L_a;
  }  
  
  // Beverton-Holt SR relationship function
  real BevertonHolt(real SSB, real B0, real h) {
    real ssb_ratio = SSB / B0;
    real part_2 = (1 - ((5*h - 1) / (4*h)) * ( 1 - ssb_ratio));
    return ssb_ratio / part_2;
  }
}


data {
  // Model dimensions
  int<lower=1> Y; // number of years
  int<lower=1> A; // number of ages
  int<lower=1> R; // number of regions
  int<lower=1> T; // number of Tag release events - represent tagged partitions.

  int years[Y] ; // model years
  int tag_years[T] ;  // tag release years
  
  // Tag release is expected to occur post ageing pre-movement, so first age shouldn't have tags.
  matrix[A,T] tag_release_by_age[R]; 

  // Observational inputs
  int<lower = 0, upper = 1> fishery_obs_indicator[R,Y]; // 0 no observation in this year and fishery
  int<lower = 0, upper = 1> biomass_indicator[R,Y];     // 0 no observation in this year and fishery
  int<lower = 0, upper = 1> tag_recapture_indicator[R,T,Y];     // 0 no observation in this year and fishery


  // Observational inputs
  matrix[A,Y] fishery_at_age_obs[R]; // Observations
  vector[Y] biomass_obs[R]; // Observations
  matrix[R,Y] tag_recapture_proportions[T]; // Don't have to sum = 1, we calculate the NR group during observation calculation.
  int tag_recapture_obs[T,R,Y]; // Int's for the negative binomial

  vector[Y] fishery_at_age_error[R]; // effective sample size
  vector[Y] biomass_error[R]; // CV's for a lognormal
  vector[T] tag_recapture_eff; // effective sample size for Tag-recapture observations.

  
  real<lower = 0, upper = 1> proportion_mortality_spawning;
  real<lower = 0, upper = 1> proportion_mortality_survey;
  real<lower = 0, upper = 1> u_max;
  real<lower = 0> catch_penalty;
  int<lower = 0, upper = 1> penalty_in_log_space; // 0 = no, 1 = yes

  // Catch history
  vector[Y] catches[R]; // catches biomass
  
  // Biological parameters
  real<lower = 0> M; // M
  real<lower = 0> a; // a in the mean weight calculation
  real<lower = 0> b; // b in the mean weight calculation
  real<lower = 0> L_inf; // L-inf in the Von Bertallanffy equation
  real<lower = 0> k; // k in the Von Bertallanffy equation
  real<lower = -5> t0; // t0 in the Von Bertallanffy equation
  real<lower = 0> h; // steepness
  real<lower = 1, upper = 20> mat_a50; // Maturity Parameters
  real<lower = 1, upper = 20> mat_ato95; // Maturity Parameters
  int<lower = 0, upper = 1> apply_prior;// 0 = don't apply prior, 1 = apply priors
  int<lower = 0, upper = 1> ycs_prior_applies_to_standardised;// 0 = estimated Y's ~lognormal, 1 = YCS (standardised)~lognromal()
  vector<lower = 0, upper = 1>[R] proportion_recruitment_by_region;
  //matrix[R,R] movement_matrix;

  simplex[R] prob_move[R];  
  int<lower = 0, upper = 1> debug;
  //vector[L] length_bins; // The last bin is a plus group, need to add empirical length at age.
  /*
  simplex[Y] unity_YCS; // YCS
  real<lower = 10, upper = 20> ln_R0; // R0
  real<lower = 0> Q; // survey_catchability

  real<lower = 1, upper = 20> s_a50;
  real<lower = 1, upper = 20> s_ato95;
  real<lower = 1, upper = 20> f_a50;
  real<lower = 1, upper = 20> f_ato95;
  */
}

// This section is equivalent to dobuild and do validate, for operations that don't need to be
// re executed they should go here. For example, truncating observations etc
transformed data {
  int total_groups = T + 1; // need to increment it so that we have untagged group accounted for when interating over the partition
}

/*
 * The estimated parameters
*/
parameters {
  // Can expand this later for area specific vulnerability
  /*
  real<lower = 1, upper = 20> s_a50;
  real<lower = 1, upper = 20> s_ato95;
  real<lower = 1, upper = 20> f_a50;
  real<lower = 1, upper = 20> f_ato95;
  real<lower = 0.01, upper = 1> Q;
    vector <lower = 0.01, upper = 20>[Y] YCS; // YCS
  real<lower = 0> R0; // R0
  
  real<lower = 0, upper = 1> p;
*/
  //simplex[R] prob_move[R];
  
  simplex[Y] unity_YCS; // YCS
  real<lower = 13, upper = 20> ln_R0; // R0
  real<lower = 0.000001, upper = 1> Q; // survey_catchability

  
  real<lower = 1, upper = 20> f_a50;
  real<lower = 1, upper = 20> f_ato95;
  //real<lower = 0.000001> overdispersion_tag_ll;
  
  //real<lower = 0.001, upper = 2> sigma_R;
/*
  real<lower = 1, upper = 20> s_a50;
  real<lower = 1, upper = 20> s_ato95;
  real<lower = 2.5, upper = 3.5> s_a50;
  real<lower = 2.5, upper = 3.5> s_ato95;
  real<lower = 2.5, upper = 3.5> f_a50;
  real<lower = 2.5, upper = 3.5> f_ato95;
  real<lower = 0.2, upper = 0.25> Q;
  vector <lower = 0.999, upper = 1.0001>[Y] YCS; // YCS
  real<lower = 50000000, upper = 51000000> R0; // R0*/

}

/*
 * Transformed Data
 * The purpose of this section of code is for ...
*/

transformed parameters{
  matrix[A,Y+1] N[R,total_groups]; // an array representing each region with a matrix of numbers at age for each year
  vector[A] cache_N[R,total_groups];

  vector[A] ages; // vector of ages
  // model expectations
  matrix[A,Y] pre_age_expectations[R]; // Expected survey comp values
  vector[Y] standardised_ycs;
  matrix[A,Y] age_expectations[R]; // Expected survey comp values
  matrix[A,Y] fishery_age_expectations[R];// Expected survey comp values
  vector[Y] biomass_expectations[R]; // expected biomass
  vector[A] numbers_at_age_with_error;
  matrix[A,Y] recapture_expectations[R,T]; // Expected survey comp values
  vector[Y] recapture_expected_props[R,T]; // Expected survey comp values

  // model quantities
  vector[A] maturity_at_age; // maturity schedule
  vector[A] fish_select_at_age; // Fishery selectivity
  vector[A] length_at_age; // length at age
  vector[A] weight_at_age; // weight at age
  vector[Y] SSB;
  vector[Y] pre_SSB;
  vector[Y] recruits;
  vector[Y] exploitation[R];
  vector[Y] vulnerable[R];
  vector[Y] actual_catches[R];
  
  vector[A] temp_partition;

  vector[A] zero_vec; // used to broadcast arrays to 0
  matrix[A,Y] zero_matrix;
  matrix[A,Y+1] alt_zero_matrix;

  vector[A - 1] temp_ageing_partition;

  vector[R] u_obs;
  real plus_group;
  real B0 = 0.0;
  real pre_B0 = 0.0;
  real post_SSB = 0.0;
  real flag_catch_penalty = 0.0;
  real ll_catch_penalty = 0.0;
  real ll_tag_penalty = 0.0;
  real neg_ll_for_reporting = 0;
  vector[R] neg_ll_fishery_age;
  vector[R] neg_ll_bio;
  vector[T] neg_ll_tag_recapture;

  real ssb_ratio = 0.0;
  real part_2 = 0.0;
  real R0 = exp(ln_R0);
  real tag_exp_NR_group = 0.0;
  real tag_obs_NR_group = 0.0;

  for (age in 1:A) {
    ages[age] = age;
  }
  // Do preliminary calculations
  length_at_age = VonBertalanffy(ages, L_inf, k, t0);
  weight_at_age = mean_weight_at_age(length_at_age, a, b);
  fish_select_at_age = logistic_ogive(ages, f_a50, f_ato95);
  maturity_at_age = logistic_ogive(ages, mat_a50, mat_ato95);
  // Initialise arrays and containers to be 0.0, have to do this because main function does a lot of incremental += stuff so need to set to 0
  zero_vec = rep_vector(0.0, A);
  neg_ll_tag_recapture = rep_vector(0.0, T);
  zero_matrix = rep_matrix(0.0,A,Y);
  alt_zero_matrix = rep_matrix(0.0,A,Y + 1);
  u_obs = rep_vector(0.0,R);
  SSB = rep_vector(0.0, Y);
  pre_SSB = rep_vector(0.0, Y);
  vulnerable = rep_array(SSB,R);
  exploitation = rep_array(SSB,R);
  actual_catches = rep_array(SSB,R);
  biomass_expectations = rep_array(SSB, R);
  N = rep_array(alt_zero_matrix, R,total_groups);
  pre_age_expectations =  rep_array(zero_matrix, R);
  age_expectations =  rep_array(zero_matrix, R);
  fishery_age_expectations =  rep_array(zero_matrix, R);
  recapture_expectations =  rep_array(zero_matrix, R,T);
  recapture_expected_props =  rep_array(SSB, R,T);
  for (t in 1:Y) 
    standardised_ycs[t] = unity_YCS[t] * Y;
  /* 
  *  ======================
  *  Initialise Partition
  *  ======================
  */   
  //
  // Approximate equilibrium age-structure without movement.
  for (reg in 1:R) {
    N[reg, 1,1,1] = R0 * proportion_recruitment_by_region[reg];
    for(age in 2:(A - 1)) {
      N[reg, 1, age, 1] = N[reg, 1, age - 1, 1] * exp(-M);
    }
    N[reg, 1,A,1] =  N[reg, 1,A - 1, 1] * exp(-M)/(1-exp(-M));
    temp_partition = N[reg, 1,,1];
    N[reg, 1,,1] = temp_partition * exp(-M);
  }
  
  // Apply movement, Rec, M, A times
  for(init_years in 1:(A * 2)) {
    // Movement
    cache_N = rep_array(zero_vec,R,total_groups);
    for (from_reg in 1:R) {
      // Ageing
      plus_group = N[from_reg, 1, A, 1];
      temp_ageing_partition = N[from_reg, 1, 1:(A-1), 1]; // To stop deep copying alisaing
      N[from_reg, 1,2:A, 1] = temp_ageing_partition;
      N[from_reg, 1,A,1] += plus_group;
      N[from_reg, 1,1,1] = 0.0; 
      // Now do movement
      for (to_reg in 1:R) {
        for (age in 1:A) {
          cache_N[to_reg, 1, age] += N[from_reg, 1, age, 1] * prob_move[from_reg, to_reg];
        }
      }
    }
    // Re-Merge, Cached N
    pre_B0 = 0.0;
    for (reg in 1:R) {
      N[reg, 1,,1] = cache_N[reg, 1,];
      N[reg, 1,1,1] = R0 * proportion_recruitment_by_region[reg];
      pre_B0 += sum(N[reg, 1,,1] .* maturity_at_age .* weight_at_age);
      for (age in 1:A) 
        N[reg, 1,age,1] *= exp(-M);
    }
  }
  // Calculate B0
  post_SSB = 0.0;
  for (reg in 1:R) {
    post_SSB += sum(N[reg, 1,,1] .* maturity_at_age .* weight_at_age);
  }
  B0 = pre_B0 + (post_SSB - pre_B0) * proportion_mortality_spawning;
  
  /* 
  *  ======================
  *     Main function
  *  ======================
  * - Annual cycle
  * -- Ageing
  * -- Move
  * -- Recruit
  * -- Z 
  */ 
  for(y in 2:(Y + 1)) {
    post_SSB = 0.0;
    // Ageing
    for (reg in 1:R) {
      for (tag in 1:total_groups) {
        plus_group = N[reg, tag, A, y - 1];
        temp_ageing_partition = N[reg, tag, 1:(A-1), y - 1]; // To stop deep copying alisaing
        N[reg, tag, 2:A, y] = temp_ageing_partition;
        N[reg, tag, A, y] += plus_group;
        N[reg, tag, 1, y] = 0.0; 
      }
    }
    // Check if we are releasing tags this year
    // Apply it using an exploitation so that we don't end up with negative partition during estimation
    // of crappy parameter space.
    for(tag_year in 1:T) {
      if (tag_years[tag_year] == years[y - 1]) {
        for (reg in 1:R) {
          for (age in 1:A) {
            if ((tag_release_by_age[reg, age, tag_year] / N[reg,1,age,y]) > u_max) {
              N[reg, tag_year + 1, age, y] =  N[reg,1,age,y] * u_max;
              N[reg,1,age,y] -= N[reg,1,age,y] * u_max;
              ll_tag_penalty += ((N[reg,1,age,y] * u_max) - tag_release_by_age[reg, age, tag_year]) * ((N[reg,1,age,y] * u_max) - tag_release_by_age[reg, age, tag_year]);
            } else {
              N[reg, tag_year + 1, age, y] = tag_release_by_age[reg, age, tag_year];
              N[reg,1,age,y] -= tag_release_by_age[reg, age, tag_year];
            }

          }
        }

      }
    }
    
    // Movement
    cache_N = rep_array(zero_vec,R,total_groups);
    
    
    for (from_reg in 1:R) {
      for (to_reg in 1:R) {
        for (tag in 1:total_groups) {
          for (age in 1:A) {
            cache_N[to_reg, tag, age] += N[from_reg, tag, age, y] * prob_move[from_reg, to_reg];
          }
        }
      }
    }
    // Re-Merge, Cached N
    for (reg in 1:R) {
      for (tag in 1:total_groups) {
        N[reg, tag, ,y] = cache_N[reg, tag ,];
      }
    }
    
    // Apply recruitment
    // B-H calcualtion
    if (y == 2) {
      recruits[y - 1] = R0 * standardised_ycs[y - 1];
    } else {
      ssb_ratio = SSB[y - 2] / B0; // Because of how I indexed the yearly loop we need to go back 2 years. 
      part_2 = (1 - ((5*h - 1) / (4*h)) * ( 1 - ssb_ratio));
      recruits[y - 1] = R0 * ssb_ratio / part_2 * standardised_ycs[y - 1];
    }
    // recruitment only occurs with untagged fish
    for (reg in 1:R) {
      N[reg, 1, 1,y] = recruits[y - 1] * proportion_recruitment_by_region[reg];
    }
    
    ///----------------
    //  Pre-Execute
    //-----------------
    for (reg in 1:R) {
      for (tag in 1:total_groups) {
        pre_SSB[y-1] += sum(N[reg, tag, ,y] .* maturity_at_age .* weight_at_age);
        //for (age in 1:A) 
        for(age in 1:A)
          pre_age_expectations[reg, age ,y-1] += N[reg, tag, age ,y] * fish_select_at_age[age];
      }
    }

    // Calculate F for each area
    for (reg in 1:R) {
      for (tag in 1:total_groups)
        vulnerable[reg,y-1] += sum(N[reg, tag, , y] * exp(-0.5 * M) .* fish_select_at_age .* weight_at_age);
      
      u_obs[reg] = max(catches[reg,y-1] / vulnerable[reg,y-1] * fish_select_at_age);
      // Uobs is just for reporting and comparing with u_max
      if (u_obs[reg] > u_max) {
        exploitation[reg, y - 1] = (catches[reg, y - 1] / vulnerable[reg, y - 1]) * (u_max / u_obs[reg] );
        flag_catch_penalty = 1.0;
        actual_catches[reg, y - 1] = exploitation[reg, y - 1] * vulnerable[reg, y - 1];
        u_obs[reg] = u_max;
        if (penalty_in_log_space == 1) {
          ll_catch_penalty += pow(log(catches[reg, y - 1]) - log(actual_catches[reg, y - 1]),2); 
        } else {
          ll_catch_penalty += pow(catches[reg, y - 1] - actual_catches[reg, y - 1],2); 
        }
      } else {
        exploitation[reg, y - 1] = u_obs[reg] ;
        actual_catches[reg, y - 1] = catches[reg, y - 1];
        u_obs[reg] = catches[reg, y - 1] / vulnerable[reg, y - 1];
      }   
    }
    
    // Calculate Fishery observations
    for (reg in 1:R) {
      if (fishery_obs_indicator[reg, y-1] == 1) {
        temp_partition = rep_vector(0.0, A);
        for(tag in 1:total_groups) {
          for(age in 1:A)
            temp_partition[age] += N[reg, tag, age ,y] * u_obs[reg] .* fish_select_at_age[age] * exp(-M * 0.5);
        }
        /*
        // Ageing error TODO add at later date
        temp_partition = fishery_age_expectations[reg,,y-1];
        // A bit annoying just reset the container to be 0's
        numbers_at_age_with_error = rep_vector(A,0.0);
          
        for (j_ndx in 1:A) {
          for (k_ndx in 1:A) {
            numbers_at_age_with_error[k_ndx] += temp_partition[j_ndx] * ageing_error[j_ndx,k_ndx];
          }
        }
        temp_partition = numbers_at_age_with_error;
        */
        
        // convert Numbers at age to proportions at age
        for(age in 1:A)
          fishery_age_expectations[reg,age,y-1] += (temp_partition[age] / sum(temp_partition));        
        
      }
    }
    
    // Tag recapture observations
    for (reg in 1:R) {
      for(tag in 1:T) {
        if (tag_recapture_indicator[reg,tag,y-1] > 0) {
          // We need to predict recaptures this year
          recapture_expectations[reg,tag,,y-1] = N[reg, tag + 1, ,y] * u_obs[reg] .* fish_select_at_age * exp(-M * 0.5);
          
        }
      } 
    }
    
    
    // Apply mortality
    for (reg in 1:R) {
      for(tag in 1:total_groups) {
        temp_partition = N[reg, tag, , y];
        N[reg, tag, , y] = temp_partition * exp(-M) .* (1.0 - u_obs[reg] * fish_select_at_age);
      }
    } 
    //----------------
    //  Post-Execute
    //----------------
    // Calculate Survey observations
    for (reg in 1:R) {
      temp_partition = rep_vector(0.0, A);
      for(tag in 1:total_groups) {
        post_SSB += sum(N[reg, tag, , y] .* maturity_at_age .* weight_at_age);
        for(age in 1:A)
          temp_partition[age]  += N[reg, tag, age ,y] .* fish_select_at_age[age];
      }
      age_expectations[reg, ,y-1] = pre_age_expectations[reg, ,y-1] + (temp_partition - pre_age_expectations[reg, ,y-1]) * proportion_mortality_survey;
      temp_partition = age_expectations[reg, ,y-1];
      age_expectations[reg, ,y-1] = temp_partition / sum(temp_partition);
      biomass_expectations[reg, y - 1] += sum(temp_partition .* weight_at_age) * Q;
      
    }
    SSB[y - 1] = pre_SSB[y-1] + (post_SSB - pre_SSB[y-1]) * proportion_mortality_spawning;

  }
  
  // Evaluate likelihood
  for (reg in 1:R) {
    neg_ll_fishery_age[reg] = obj_multinomal_casal_lpdf(fishery_at_age_obs[reg] | fishery_age_expectations[reg],  fishery_at_age_error[reg],fishery_obs_indicator[reg]);
  }
  /*
  for (reg in 1:R) {
    neg_ll_survey_age[reg] = obj_multinomal_casal_lpdf(survey_at_age_obs[reg] | survey_age_expectations[reg],  survey_at_age_error[reg],survey_age_indicator[reg]);
  }
  */
  for (reg in 1:R) {
    neg_ll_bio[reg] = obj_lognomal_lpdf(biomass_obs[reg] | biomass_expectations[reg], biomass_error[reg], biomass_indicator[reg]);
  }
  // Custom Tag-recapture observations 
  // multinomial first
  for(tag in 1:T) {
    tag_exp_NR_group = 0.0;
    tag_obs_NR_group = 0.0;
    for(y in 1:Y) {
      for (reg in 1:R) {
        if (tag_recapture_indicator[reg,tag,y] > 0) {
          // Calculate proportions
          recapture_expected_props[reg,tag,y] = sum(recapture_expectations[reg,tag,,y]) / sum(tag_release_by_age[reg, , tag]);
          neg_ll_tag_recapture[tag] += tag_recapture_proportions[tag,reg,y] * log(recapture_expected_props[reg,tag,y]);
          
          tag_exp_NR_group += recapture_expected_props[reg,tag,y];
          tag_obs_NR_group += tag_recapture_proportions[tag,reg,y];
        }
      }
    }
    // NR group.
    neg_ll_tag_recapture[tag] += (1 - tag_obs_NR_group) * log(1 - tag_exp_NR_group);
    neg_ll_tag_recapture[tag] *= tag_recapture_eff[tag];
  }
  
  // Negative binomial
}

model {
  // turn it into a positive because stan does maximum likelihood not neg ll
  target += -sum(neg_ll_fishery_age) - sum(neg_ll_bio) - ll_tag_penalty - ll_catch_penalty;// - sum(neg_ll_tag_recapture);
    
    // Tag recapture likelihood contribution negative binomial.
    for (reg in 1:R) {
      for(tag in 1:T) {
        for(y in 1:Y) {
          if (tag_recapture_indicator[reg,tag,y] > 0) {
            // We need to predict recaptures this year
            tag_recapture_obs[tag, reg, y] ~ poisson(recapture_expectations[reg,tag,,y]);
            //tag_recapture_obs[tag, reg, y] ~ neg_binomial_2(recapture_expectations[reg,tag,,y], overdispersion_tag_ll);
          }
        }
      } 
    }
  //target += (-1.0);
}


generated quantities {
   // ... declarations ... statements ...
   
  int sim_tag_recapture[T,R,Y];
  vector[Y + 1] expectations;
  int temp[Y + 1];
  real NR_group = 0.0;
  for(tag in 1:T) {
    for (reg in 1:R) {
      expectations[1:Y] = recapture_expected_props[reg,tag,];
      expectations[Y + 1] = 1.0 - sum(recapture_expected_props[reg,tag,]);
      temp = multinomial_rng(expectations, 2000);
      sim_tag_recapture[tag,reg,1:Y] = temp[1:Y];
    }
  }
  
  
  
}
