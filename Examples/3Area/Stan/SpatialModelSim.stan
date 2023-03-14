/*
 * @file AgeStructuredModel_1.0.stan
 * @author C. Marsh
 *
 * A stan module for a simple age structured dynamics model
 * vector<> is a column vector!! remember this when structuring matrices used in multiplication
 */
 
 
functions {
  // logistic ogive function
  real gmean(row_vector x) {
    real tot = 1.0;
    for (i in 1:num_elements(x))
      tot *= x[i];
    return tot^(1.0/num_elements(x));
  }
  /* compute the cholesky factor of an AR1 correlation matrix
   * @author paul.buerkner@gmail.com
   * Args: 
   *   ar: AR1 autocorrelation 
   *   nrows: number of rows of the covariance matrix 
   * Returns: 
   *   A nrows x nrows matrix 
   */ 
   matrix cholesky_cor_ar1(real ar, int nrows) { 
     matrix[nrows, nrows] mat; 
     vector[nrows - 1] gamma; 
     mat = diag_matrix(rep_vector(1, nrows)); 
     for (i in 2:nrows) { 
       gamma[i - 1] = pow(ar, i - 1); 
       for (j in 1:(i - 1)) { 
         mat[i, j] = gamma[i - j]; 
         mat[j, i] = gamma[i - j]; 
       } 
     } 
     return cholesky_decompose(mat); 
   }
  // Logisitic likelihood negative loglikelihood
  real calculate_logistic_lpdf(matrix Observations, matrix expectations, matrix covariance, vector sample_size_by_year) {
    int A = cols(expectations);
    int Y = rows(expectations);
    matrix[A - 1, A - 1] v_mat;
    matrix[A - 1, A - 1] inv_v_mat;
    matrix[A - 1, A] k_mat;
    matrix[A, A - 1] t_k_mat;
  
    matrix[Y, A-1]  ww;
    vector[Y] weights_by_year;  // Sample sizes to 
    real log_det;
    real temp_val = 0.0;
    real log_total_obs = 0.0;
    
    real neg_ll;
    // initialisa
    k_mat = rep_matrix(0.0, A - 1, A);
    v_mat = rep_matrix(0.0, A - 1, A - 1);
    ww = rep_matrix(0.0, Y, A - 1);
    k_mat[,A] = rep_vector(-1.0, A - 1);
    for (i in 1:(A-1)) 
      k_mat[i,i] = 1.0;
    t_k_mat = covariance * (k_mat');
    v_mat = k_mat * t_k_mat;
    inv_v_mat = inverse(v_mat);
  
    log_det = log_determinant(v_mat);
    log_total_obs = sum(log(Observations));
    neg_ll = 0.5 * Y * (A - 1)* log(2.0*pi()) + log_total_obs + 0.5 * Y * log_det;
    weights_by_year = sqrt(mean(sample_size_by_year) ./ sample_size_by_year);
    neg_ll += (A - 1) * sum(log(weights_by_year));
    for (y in 1:Y) 
      ww[y,] = log(Observations[y, 1:(A-1)] ./ Observations[y, A]) - log(expectations[y, 1:(A-1)] ./ expectations[y, A]);
    
    for (y in 1:Y) 
      neg_ll += sum((0.5/(weights_by_year[y]*weights_by_year[y])) * (ww[y,] * inv_v_mat) .* ww[y,]);
    return neg_ll;
  }
  
  matrix standardised_residuals(matrix Observations, matrix expectations, matrix covariance, vector sample_size_by_year, int centered) {
    int A = cols(expectations);
    int Y = rows(expectations);
    matrix[A - 1, A - 1] v_mat;
    matrix[A - 1, A - 1] inv_v_mat;
    matrix[A - 1, A] k_mat;
    matrix[A - 1, A] FF;
    matrix[A, A - 1] t_k_mat;
    matrix[A-1, A-1] HH;
    matrix[Y, A-1]  ww;
    matrix[A, A-1] FHinv;
    matrix[A, A] Gam;
    matrix[Y, A] sdmat;
    matrix[Y, A - 1] sdmat_non_centered;
    matrix[Y, A] sres;
    matrix[Y, A - 1] sres_non_centered;

    vector[Y] gmean_obs;  // Sample sizes to 
    vector[Y] gmean_exp;
    vector[Y] weights_by_year;  // Sample sizes to 
    real log_det;
    real temp_val = 0.0;
    real log_total_obs = 0.0;
    
    // initialisa
    k_mat = rep_matrix(0.0, A - 1, A);
    v_mat = rep_matrix(0.0, A - 1, A - 1);
    ww = rep_matrix(0.0, Y, A - 1);
    weights_by_year = sqrt(mean(sample_size_by_year) ./ sample_size_by_year);

    k_mat[,A] = rep_vector(-1.0, A - 1);
    for (i in 1:(A-1)) 
      k_mat[i,i] = 1.0;

    
    t_k_mat = covariance * (k_mat');
    v_mat = k_mat * t_k_mat;
    inv_v_mat = inverse(v_mat);
    
    if(centered) {
        HH = rep_matrix(1.0, A - 1, A - 1);
        for (i in 1:cols(HH))
          HH[i,i] = 2.0;
          
        FF  = k_mat;
        FF[,A] = rep_vector(1.0, A - 1);
        FHinv = FF' * inverse(HH);
        Gam = FHinv * (v_mat * FHinv');
        for(y in 1:Y) 
          sdmat[y,] = (sqrt(diagonal(Gam)) * weights_by_year[y])';
        for(y in 1:Y)  {
          gmean_obs[y] = gmean(Observations[y,]);
          gmean_exp[y] = gmean(expectations[y,]);
          sres[y, ] = (log(Observations[y,] / gmean_obs[y]) - log(expectations[y,] / gmean_exp[y])) ./ sdmat[y,];
        }
        return sres;
    } else {
      for(y in 1:Y) 
        sdmat_non_centered[y,] = (sqrt(diagonal(v_mat)) * weights_by_year[y])';
      for(y in 1:Y)  {
        sres_non_centered[y, ] = (log(Observations[y,1:(A - 1)] / Observations[y,A]) - log(expectations[y,1:(A - 1)] / expectations[y,A])) ./ sdmat_non_centered[y,];
      }
      return sres_non_centered;
    }
    
    // Should have exited before here, but compiler will complain.
    return sres;
  }

  matrix ar1_covar(real sigma, real rho, int bins) {
    matrix[bins,bins] covariance_matrix;
    matrix[bins,bins] cholesky_cor = cholesky_cor_ar1(rho, bins);
    covariance_matrix = diag_pre_multiply(rep_vector(sigma,bins), cholesky_cor) * diag_pre_multiply(rep_vector(sigma,bins), cholesky_cor)';
    return covariance_matrix;
  }
  
  
  // logistic ogive function
  vector logistic_ogive(vector ages, real sel_50, real sel_95) {
    int n_ages = rows(ages);
    vector[n_ages] logis;
    for (age in 1:n_ages) {
      logis[age] = 1 / (1 + pow(19, (sel_50 - ages[age]) / sel_95));
    }
    return logis;
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
  
  
  // Logisitc likelihood from "replacing the multinomial"
  // EQ A.9
  /*
  real logisitic_normal_lpdf(matrix obs, matrix Exp, matrix covar, int[] indicator) {
	  // calculate v_inverse
  }
  */
}


data {
  // Model dimensions
  int<lower=1> Y; // number of years
  int<lower=1> A; // number of ages
  int<lower=1> R; // number of regions
  int<lower=1> T; // number of Tag release events - represent tagged partitions.

  int years[Y] ; // model years
  int tag_years[T] ;  // tag release years
  int N_sims; // number of simulations

  // Tag release is expected to occur post ageing pre-movement, so first age shouldn't have tags.
  matrix[A,T] tag_release_by_age[R]; 

  // Observational inputs
  // Indicators
  int<lower = 0, upper = 1> fishery_obs_indicator[R,Y];     // 0 = no, 1 = yes observation in this year and fishery
  int<lower = 0, upper = 1> biomass_indicator[R,Y];         // 0 = no, 1 = yes observation in this year and fishery
  int<lower = 0, upper = 1> tag_recapture_indicator[R,T,Y]; // 0 = no, 1 = yesobservation in this year and fishery
  // observed values
  int fishery_at_age_obs[R,A,Y];      // multinomial counts, observed proportions and error values
  //matrix[A,Y] fishery_at_age_obs[R]; // Proportions at age for logistiic normal
  

  vector[Y] biomass_obs[R];
  // Error values
  vector[Y] biomass_error[R];         // CV's for a lognormal
  
  real<lower = 0, upper = 1> proportion_mortality_spawning;
  real<lower = 0, upper = 1> proportion_mortality_biomass;
  real<lower = 0, upper = 1> u_max; // Maxium exploitaton rate, flag penalty if greater than this
  real<lower = 0> catch_penalty;  // Penalty multiplier to add to likelihood
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
  vector<lower = 0, upper = 1>[R] proportion_recruitment_by_region;
  simplex[R] prob_move[R];  // rows = from, cols = to
  
  simplex[Y] unity_YCS; // YCS
  real<lower = 13, upper = 20> ln_R0; // R0
  real<lower = 0.000001, upper = 1> Q; // survey_catchability

  real<lower = 1, upper = 20> f_a50;
  real<lower = 1, upper = 20> f_ato95;
  //vector[L] length_bins; // The last bin is a plus group, need to add empirical length at age.

  // Logistic likelihood stuff
  real rho;
  real sigma;
  
  real<lower = 0> tag_neg_bin_overdispersion;
  
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

  real<lower = 0, upper = 1> p;
 
  //simplex[R] prob_move[R];
  

  //real<lower = 0.000001> overdispersion_tag_ll;
  
  //real<lower = 0.001, upper = 2> sigma_R;

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
  matrix[A,Y] pre_age_expectations[R];      // Expected survey comp values
  vector[Y] standardised_ycs;               // Once converted from Simplex
  matrix[A,Y] age_expectations[R];          // Expected survey comp values
  matrix[A,Y] fishery_age_expectations[R];  // Expected survey comp values
  vector[Y] biomass_expectations[R];        // expected biomass
  vector[A] numbers_at_age_with_error;      // If ageing error is applied
  matrix[A,Y] recapture_expectations[R,T];  // Expected survey comp values
  vector[Y] recapture_expected_props[R,T];  // Expected survey comp values

  // model quantities
  vector[A] maturity_at_age;                // maturity schedule
  vector[A] fish_select_at_age;             // Fishery selectivity
  vector[A] length_at_age;                  // length at age
  vector[A] weight_at_age;                  // weight at age
  vector[Y] SSB;
  vector[Y] pre_SSB;
  vector[Y] recruits;
  vector[Y] exploitation[R];
  vector[Y] vulnerable[R];
  matrix[A,Y] tagged_vulnerable[R];
  matrix[A,Y] total_vulnerable[R];

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
  
  // Evaluate likelihood logistic likelihood
  matrix[A,A] covariance_matrix;
  covariance_matrix = ar1_covar(sigma, rho, A);
	
	
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
  
  tagged_vulnerable =  rep_array(zero_matrix, R);
  total_vulnerable =  rep_array(zero_matrix, R);

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
  for(init_years in 1:(A * 4)) {
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
  * -- Tagging
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
      for (tag in 1:total_groups) {
        vulnerable[reg,y-1] += sum(N[reg, tag, , y] * exp(-0.5 * M) .* fish_select_at_age .* weight_at_age);
        for (age in 1:A) {
          if (tag > 1) 
            tagged_vulnerable[reg, age, y-1] += N[reg, tag, age, y] * exp(-0.5 * M) * fish_select_at_age[age];
          total_vulnerable[reg, age, y-1] += N[reg, tag, age, y] * exp(-0.5 * M) * fish_select_at_age[age];
          
        }
      }
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
      biomass_expectations[reg, y - 1] = sum(pre_age_expectations[reg, ,y-1] .* weight_at_age) + (sum(temp_partition .* weight_at_age) - sum(pre_age_expectations[reg, ,y-1] .* weight_at_age)) * proportion_mortality_biomass;
      biomass_expectations[reg, y - 1] *= Q;
      
      age_expectations[reg, ,y-1] = pre_age_expectations[reg, ,y-1] + (temp_partition - pre_age_expectations[reg, ,y-1]) * proportion_mortality_biomass;
      temp_partition = age_expectations[reg, ,y-1];
      age_expectations[reg, ,y-1] = temp_partition / sum(temp_partition);
      //biomass_expectations[reg, y - 1] += sum(temp_partition .* weight_at_age) * Q;
      
    }
    SSB[y - 1] = pre_SSB[y-1] + (post_SSB - pre_SSB[y-1]) * proportion_mortality_spawning;

  }

	  

}

model {
  // Evaluate likelihood
  /*
  real lognormal_sigmas;
  for (y in 1:Y) {
    for (reg in 1:R) {
      // Multinomial age dist first
      if (fishery_obs_indicator[reg, y] == 1)
        fishery_at_age_obs[reg, , y] ~ multinomial(fishery_age_expectations[reg,,y]);
      if (biomass_indicator[reg, y] == 1) {
        lognormal_sigmas = sqrt(log(1.0 + biomass_error[reg, y] * biomass_error[reg,y]));
        
        biomass_obs[reg,y] ~ lognormal(log(biomass_expectations[reg,y]) - 0.5*lognormal_sigmas*lognormal_sigmas,lognormal_sigmas);
      }
    }
  }
  // turn it into a positive because stan does maximum likelihood not neg ll
  target += -ll_tag_penalty -ll_catch_penalty;// - sum(neg_ll_tag_recapture);
  */
    /*
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
    }*/
  target += (-1.0);
}


generated quantities {
   // ... declarations ... statements ...
  int sim_tag_recapture_poisson[R,T,Y,N_sims];
  real lambda_tag_recapture_poisson[R,T,Y];
  int sim_tag_recapture_negative_binomial[R,T,Y,N_sims];
  matrix[Y,A] sim_fishery_at_age_obs_logis[R,N_sims];

  vector[Y] sim_biomass_obs[R,N_sims];
  real lognormal_sigmas;

  lambda_tag_recapture_poisson = rep_array(0.0,R,T,Y);
  //sim_tag_recapture_poisson = rep_array(0,R,T,Y,N_sims);

  for (reg in 1:R) {
    for (y in 1:Y) {
      for(tag in 1:T) {
        sim_tag_recapture_poisson[reg, tag, y,] = rep_array(0,N_sims);
        sim_tag_recapture_negative_binomial[reg, tag, y,] = rep_array(0,N_sims);
        if (sum(recapture_expectations[reg,tag,,y]) > 0) {
          lambda_tag_recapture_poisson[reg, tag, y] = sum(recapture_expectations[reg,tag,,y]);
          for(sim in 1:N_sims) {
            sim_tag_recapture_poisson[reg, tag, y,sim] = poisson_rng(lambda_tag_recapture_poisson[reg, tag, y] );
            sim_tag_recapture_negative_binomial[reg, tag, y,sim] = neg_binomial_2_rng(lambda_tag_recapture_poisson[reg, tag, y], tag_neg_bin_overdispersion);
          }
        }
        
      }
    
      // Multinomial age dist first
      if (fishery_obs_indicator[reg, y] == 1) {
        for(sim in 1:N_sims) {
          sim_fishery_at_age_obs_logis[reg, sim, y, ] = softmax(multi_normal_rng(log(fishery_age_expectations[reg,,y]), covariance_matrix))';
        }
      }
      if (biomass_indicator[reg, y] == 1) {
        lognormal_sigmas = sqrt(log(1.0 + biomass_error[reg, y] * biomass_error[reg,y]));
        for(sim in 1:N_sims) 
          sim_biomass_obs[reg,sim, y] = lognormal_rng(log(biomass_expectations[reg,y]) - 0.5*lognormal_sigmas*lognormal_sigmas,lognormal_sigmas);
      }
    }
  }
}
