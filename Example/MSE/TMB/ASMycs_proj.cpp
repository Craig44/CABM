#include <TMB.hpp>
/*
 * isNA
 */
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}
/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}
/*
 * 
 */
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}
/*
 * 
 */
template <class Type> 
Type square(Type x){return x*x;}

template <class Type> 
vector<Type> square(vector<Type>& x) {
  return x*x;
}
// transform Y -Inf-Inf -> X bound lb - ub
template <class Type> 
Type invlogit_general(Type& Y, Type& lb, Type& ub) {
  return(lb + (ub - lb) * (1 / (1 + exp(-Y))));
}


// logistic ogive function
template <class Type> 
vector<Type> logistic_ogive(vector<Type> ages, Type sel_50, Type sel_95) {
  std::cout << "logistic_ogive\n";
  int n_ages = ages.size();
  vector<Type> logis(n_ages);
  for (int age = 0;  age < n_ages; ++age) {
    logis[age] = Type(1.0) / (Type(1.0) + pow(Type(19.0), (sel_50 - ages[age]) / sel_95));
  }
  return logis;
}
/*
 * Geometric mean
 */
template <class Type> 
Type geo_mean(vector<Type>& x){
  return exp((log(x).sum())/x.size());
}

template <class Type>
void get_covar(vector<Type> theta, matrix<Type>& covar, int n, int type) {
  // https://github.com/glmmTMB/glmmTMB/blob/9cb23d77afc042be5bced3713c40529023d9ef7e/glmmTMB/src/glmmTMB.cpp#L361
  if (type == 1) {
    // case: diag_covstruct
    //std::cerr << "diag_covstruct" << std::endl;
    Type sd = exp(theta(0));
    for(int i = 0; i < n; ++i)
      covar(i,i) = sd * sd;
  } else if (type == 2) {
    // case: us_covstruct
    //std::cerr << "us_covstruct" << std::endl;
    vector<Type> sd = exp(theta.head(n));
    vector<Type> corr_transf = theta.tail(theta.size() - n);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    // for covariance
    //density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
        covar(i,j) = nldens.cov()(i,j) * sd(i) * sd(j);
      }
    }
  } else if (type == 3){
    // case: cs_covstruct Compound Symmetry: Heterogenous.
    //std::cerr << "cs_covstruct" << std::endl;
    
    vector<Type> sd = exp(theta.head(n));
    Type corr_transf = theta(n);
    Type a = Type(1) / (n - Type(1)); // Correcting for violations of sphericity 
    Type rho = invlogit(corr_transf) * (Type(1) + a) - a;
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        covar(i,j) = (i==j ? Type(1) * sd(i) * sd(i) : rho * sd(i) * sd(j)) ;
    
  } else if (type == 4){
    // case: toep_covstruct Toeplitz: Heterogenous.
    //std::cerr << "toep_covstruct" << std::endl;
    vector<Type> sd = exp(theta.head(n));
    vector<Type> parms = theta.tail(n-1);              // Corr parms
    parms = parms / sqrt(Type(1.0) + parms * parms );  // Now in (-1,1)
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        covar(i,j) = (i==j ? Type(1) * sd(i) * sd(i) : parms( (i > j ? i-j : j-i) - 1 ) * sd(i) * sd(j));
  } else if (type == 5){
    //std::cerr << "ar1_covstruct" << std::endl;
    // case: ar1_covstruct
    //  * NOTE: Valid parameter space is phi in [-1, 1]
    //  * NOTE: 'times' not used as we assume unit distance between consecutive time points.
    vector<Type> sd(n);
    sd.fill(exp(theta(0)));
    Type corr_transf = theta(1);
    Type phi = corr_transf / sqrt(1.0 + pow(corr_transf, 2));
    for(int i=0; i < n; i++){
      for(int j=0; j < n; j++){
        covar(i,j) = pow(phi, abs(i-j)) * sd(i) * sd(i);
      }
    }
  }
  
  return;
}
/*
 * centred log transform
 */
template <class Type> 
vector<Type> crl(vector<Type>& x) { 
  return log(x / geo_mean(x));
}
/*
 * rmultinomm - for simulate call
 */
template <class Type> 
vector<Type> rmultinom(vector<Type> prob, Type N) { 
  vector<Type> sim_X(prob.size());
  sim_X.setZero();
  for(int i = 0; i < prob.size(); ++i)
    sim_X(i) = rbinom(N, prob(i));
  return (sim_X);
}

/*
 *  Simulate a single draw from a multinomial-dirichlet distribution
 */
template <class Type> 
vector<Type> rdirichletmulti(vector<Type> fitted_props, Type& n_eff, Type& theta) {
  vector<Type> dirichlet_draw(fitted_props.size()); 
  for(int ndx = 0; ndx < fitted_props.size(); ndx++) 
    dirichlet_draw(ndx) = rgamma(fitted_props(ndx) * theta * n_eff, (Type)1.0);// shape, rate = 1.0
  
  Type dirich_total = dirichlet_draw.sum();
  dirichlet_draw /= dirich_total;
  return(rmultinom(dirichlet_draw, n_eff));
}
/*
 * inverse centred log transform
 * Up to a constant so standardise
 */
template <class Type> 
vector<Type> inv_crl(vector<Type>& y){ 
  return exp(y) / (exp(y)).sum();
}

// Beverton-Holt SR relationship function
template<class Type>
Type BevertonHolt(Type SSB, Type B0, Type h) {
  Type ssb_ratio = SSB / B0;
  Type part_2 = (1 - ((5*h - 1) / (4*h)) * ( 1 - ssb_ratio));
  return (ssb_ratio / part_2);
}

// Beverton-Holt SR relationship function without equilibrium assumptions
template<class Type>
Type BevertonHoltNoEquil(Type a, Type b, Type SSB) {
  Type Rt = (a + SSB) / (SSB * b);
  return Rt;
}



/*
 * A simple age structure stock assessment in TMB, that has two fisheries
 * 
 */

template<class Type>
Type objective_function<Type>::operator() () {
  /*
   * Declare Namespace
   */
  using namespace density;
  // model dimensions
  DATA_VECTOR(ages);    // assumes min(ages) >= 1
  DATA_VECTOR(years);   // annual years
  int nyears = years.size();
  int nages = ages.size();
  DATA_INTEGER(maxAgePlusGroup);  // 1 = yes, 0 = no
  
  // Observation info
  DATA_MATRIX(ageing_error_matrix);     // nages * nages
  DATA_IVECTOR(survey_year_indicator);  // 1 = calculate, 0 = ignore nyears
  DATA_VECTOR(survey_obs);              // Relative index (this will be ignored if survey_comp_type = 1)
  DATA_VECTOR(survey_cv);               // CV for relative index (this will be ignored if survey_comp_type = 1)
  DATA_VECTOR_INDICATOR(survey_index_keep, survey_obs);  // For one-step predictions
  
  DATA_ARRAY(survey_trans_comp);        // ALN transformed observations dim = length()survey_min_age:survey_max_age) - 1  NB x sum(survey_year_indicator) 
  // fishery_comp_type = 1, log(N), fishery_comp_type = 0 & fishery_comp_likelihood = 0 alr(N), fishery_comp_likelihood = 1 N_eff * prop(N)
  DATA_ARRAY_INDICATOR(keep_survey_comp, survey_trans_comp);
  DATA_SCALAR(survey_min_age);         // survey_age_indicator.size() = nages, 
  DATA_SCALAR(survey_max_age);          // survey_age_indicator.size() = nages, 
  DATA_INTEGER(survey_age_plus_group);  // 0 no plus group, 1 = an accumlation group
  DATA_INTEGER(survey_comp_type);       // 1 = numbers at age, 0 = proportions at age;
  DATA_INTEGER(survey_covar_structure); // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(survey_comp_likelihood); // 0 = MVN (can be applied to both comp_type), 1 = Multinomial, 2 = dirichlet-multinomial
  
  DATA_IVECTOR(fishery_year_indicator);  // 1 = calculate, 0 = ignore
  DATA_SCALAR(fishery_min_age);          // survey_age_indicator.size() = nages, 
  DATA_SCALAR(fishery_max_age);          // survey_age_indicator.size() = nages, 
  DATA_INTEGER(fishery_age_plus_group);  // 0 no plus group, 1 = an accumlation group
  DATA_ARRAY(fishery_trans_comp);        // sum(fishery_year_indicator) x fishery_nages NB!!!! different dims to most other matrices..
  // fishery_comp_type = 1, log(N), fishery_comp_type = 0 & fishery_comp_likelihood = 0 alr(N), fishery_comp_likelihood = 1 N_eff * prop(N)
  DATA_INTEGER(fishery_comp_type);       
  // numbers at age = 1, proportions at age = 0
  DATA_INTEGER(fishery_comp_likelihood); // 0 = MVN (can be applied to both comp_type), 1 = Multinomial, 2 = dirichlet-multinomial
  DATA_INTEGER(fishery_covar_structure); // 1 = iid (1), 5 = AR(1) (2), 2 = Unstructured (A_f(A_f - 1)/2).
  DATA_ARRAY_INDICATOR(keep_fishery_comp, fishery_trans_comp);
  
  DATA_IVECTOR(ycs_estimated);    // 1 = estimated, 0 = ignore
  DATA_INTEGER(standardise_ycs);  // 1 = yes, 0 = No
  
  DATA_INTEGER(n_proj_years);     // 1 = yes, 0 = No
  DATA_INTEGER(proj_method);      // 0 = F, 1 = Catch
  DATA_SCALAR(future_vals);       // either constant catch or F depends on proj_method
  
  // Catch info
  DATA_VECTOR(catches);               // length = nyears
  DATA_VECTOR(propZ_ssb);             // proportion of Z for SSB, length nyears
  DATA_VECTOR(propZ_survey);          // proportion of Z for SSB, length nyears
  // Biological info
  DATA_ARRAY(propMat);                // Proportion Mature (used in SSB)    dim = nages x nyears
  DATA_ARRAY(stockMeanLength);        // Stock Mean weight used in SSB + survey calculation  dim = nages x nyears
  DATA_ARRAY(catchMeanLength);        // Stock Mean weight used in Catch equation calculation dim = nages x nyears
  DATA_SCALAR(natMor);                // Instantaneous natural mortality        
  DATA_SCALAR(steepness);             // Instantaneous natural mortality     
  DATA_SCALAR(mean_weight_a);         // a parameter mean weight from length 
  DATA_SCALAR(mean_weight_b);         // b parameter mean weight from length    
  
  DATA_INTEGER(stockRecruitmentModelCode); // SR relationship 0 = RW, 1 = Ricker, 2 = BH

  // Bounds for selectivity parameters for model stability
  DATA_VECTOR(sel_ato95_bounds);  // length 2
  DATA_VECTOR(sel_a50_bounds);    // length 2

  
/*
 * Parameters to estiamte.
 */
  PARAMETER(ln_R0); 
  PARAMETER_VECTOR(ln_ycs_est);           // length(recruit_devs) = sum(ycs_estimated)
  PARAMETER(ln_sigma_r);                    // logistic fishery ogive
  PARAMETER( ln_extra_survey_cv );          // Additional survey cv.
  PARAMETER_VECTOR(trans_survey_error);     // 
  PARAMETER_VECTOR(trans_fishery_error);    // 

  PARAMETER(logit_f_a50);                   // logistic fishery ogive paramss
  PARAMETER(logit_f_ato95);                 // logistic fishery ogive paramss
  
  PARAMETER(logit_survey_a50);              // logistic survey ogive paramss
  PARAMETER(logit_survey_ato95);            // logistic survey ogive paramss
  // have trouble estiming these parameters with the Q so, need to bound them.
  PARAMETER(logit_surveyQ);                 // logit transformes bound between 0-1
  PARAMETER_VECTOR(ln_F);                   // length n_years
  PARAMETER(ln_catch_sd);
  /*
   * Parameter transformations
  */
  int year_ndx, age_ndx, iter, age_iter;
  Type extra_survey_cv = exp(ln_extra_survey_cv);
  Type R0 = exp(ln_R0);
  Type sigma_r = exp(ln_sigma_r);
  Type B0 = 0.0;
  Type catch_sd = exp(ln_catch_sd);
  array<Type> N(nages, nyears + 1);
  N.fill(0.0);
  // Convert mean length to mean weight
  array<Type> stockMeanWeight(stockMeanLength.dim[0], stockMeanLength.dim[1]);
  stockMeanWeight.fill(0.0);
  array<Type> catchMeanWeight(catchMeanLength.dim[0], catchMeanLength.dim[1]);
  stockMeanWeight.fill(0.0);
  
  for(iter = 0; iter < catchMeanLength.dim(0); ++iter) {
    for(age_iter = 0; age_iter < catchMeanLength.dim(1); ++age_iter) {
      stockMeanWeight(iter, age_iter) = mean_weight_a * pow(stockMeanLength(iter, age_iter), mean_weight_b);
      catchMeanWeight(iter, age_iter) = mean_weight_a * pow(catchMeanLength(iter, age_iter), mean_weight_b);
    }
  }
   
  // deal with YCS
  vector<Type> ycs(nyears);
  iter = 0;
  Type recruit_nll = 0.0;
  
  
  for(year_ndx = 0; year_ndx < nyears; ++year_ndx) {
    if (ycs_estimated[year_ndx] == 1) {
      ycs(year_ndx) = exp(ln_ycs_est(iter));
      ++iter;
    } else {
      ycs(year_ndx) = exp(rnorm(-0.5 * sigma_r * sigma_r, sigma_r));
    }
  }
  if (standardise_ycs == 1) {
    ycs /= ycs.mean();
  } 
  // Note this contains constants (non estimated ycs values), and probably needs a jacombian for the transformation.
  // mean of random variables 
  vector<Type> future_ycs(n_proj_years);
  for(year_ndx = 0; year_ndx < n_proj_years; ++year_ndx) 
    future_ycs(year_ndx) = exp(rnorm(-0.5 * sigma_r * sigma_r, sigma_r));
    
    
  Type survey_a50 = invlogit_general(logit_survey_a50, sel_a50_bounds(0), sel_a50_bounds(1));
  Type survey_ato95 = invlogit_general(logit_survey_ato95, sel_ato95_bounds(0), sel_ato95_bounds(1));
  Type f_a50 = invlogit_general(logit_f_a50, sel_a50_bounds(0), sel_a50_bounds(1));
  Type f_ato95 = invlogit_general(logit_f_ato95, sel_ato95_bounds(0), sel_ato95_bounds(1));
  Type survey_Q = invlogit(logit_surveyQ);
  /*
   * Set up container storage
   */
  vector<Type> ssb(nyears + 1);
  ssb.setZero();
  vector<Type> annual_F = exp(ln_F);
  array<Type> F_ay(nages, nyears);
  F_ay.fill(0.0);
  array<Type> Z_ay(nages, nyears);
  Z_ay.fill(0.0);
  // Fitted value containers
  
  vector<Type> predlogN(nages); 
  vector<Type> temp_partition(nages); 

  vector<Type> pred_catches(nyears);
  pred_catches.setZero();
   
  vector<Type> fishery_selectivity(nages);
  vector<Type> survey_selectivity(nages);
  vector<Type> survey_sd(survey_cv.size());
  
  // Calculate vulnerable biomass and U
  Type delta = 0.0001;
  /*
   * Build Covariance's for Observations and states currently just iid with different covariances
   */

  /*
   * Deal with F stuff
   */
  for(age_ndx = 0; age_ndx < nages; ++age_ndx) {
    fishery_selectivity(age_ndx) = Type(1)/(Type(1) + pow(Type(19),((f_a50 - ages(age_ndx))/f_ato95)));
    survey_selectivity(age_ndx) = Type(1)/(Type(1) + pow(Type(19),((survey_a50 - ages(age_ndx))/survey_ato95)));
  }
  for(year_ndx = 0; year_ndx < nyears; year_ndx++) {
    F_ay.col(year_ndx) = fishery_selectivity * annual_F(year_ndx);
    Z_ay.col(year_ndx) = F_ay.col(year_ndx) + natMor;
  }
  /*
   * Initialise first year
   */
  // Initialise Partition
  for(age_ndx = 0; age_ndx < nages; ++age_ndx) 
    N(age_ndx, 0) = R0 * exp(-(ages(age_ndx)) * natMor);
  if(maxAgePlusGroup == 1)
    N(nages - 1, 0) = R0 * exp(- ages(nages - 1) * natMor) / (1.0 - exp(-natMor));
  // Applying ageing
  temp_partition = N.col(0);
  N((nages - 1), 0) += N(nages - 2, 0);
  for(age_ndx = 1; age_ndx < (nages - 1); ++age_ndx) 
    N(age_ndx, 0) =  temp_partition(age_ndx - 1, 0);
  N(0, 0) = R0;
  ssb(0) =  sum(N.col(0) * exp(-natMor * propZ_ssb(0))  * stockMeanWeight.col(0) * propMat.col(0));
  B0 = ssb(0);
  // Now the rest M
  //for(age_ndx = 0; age_ndx < nages; ++age_ndx) 
  //  N(age_ndx, 0) *= exp(-natMor);
  
  /*
   * Start Model
   */
  for(year_ndx = 1; year_ndx <= nyears; year_ndx++) {
    //----------------
    //Recuritment
    //----------------
    if(stockRecruitmentModelCode == 0) { // straight RW 
      predlogN(0) = log(N(0,year_ndx - 1) * ycs(year_ndx - 1));
    } else if(stockRecruitmentModelCode == 2) { //BH
      predlogN(0) = log(R0 * BevertonHolt(ssb(year_ndx - 1), B0, steepness) * ycs(year_ndx - 1));
    } else{
      error("SR model code not recognized");
    }
    
    //----------------
    // F + M + ageing
    //----------------
    for(age_ndx = 1; age_ndx < nages; ++age_ndx) 
      predlogN(age_ndx) = log(N(age_ndx - 1, year_ndx - 1)) - Z_ay(age_ndx - 1, year_ndx - 1);
    
    if(maxAgePlusGroup == 1) {
      predlogN(nages - 1) = log(N(nages - 2, year_ndx - 1) * exp(- Z_ay(nages - 2, year_ndx - 1)) +
        N(nages - 1, year_ndx - 1) * exp(-Z_ay(nages - 1, year_ndx - 1)));
    }  

    // Transform from log space
    N.col(year_ndx) = exp(predlogN);
    
    // Calculate SSBs an interpolation bewtween the year, starting with previous years Paritition
    for(age_ndx = 0; age_ndx < nages; ++age_ndx) {
      ssb(year_ndx) += N(age_ndx, year_ndx - 1) * exp(-Z_ay(age_ndx, year_ndx - 1) * propZ_ssb(year_ndx - 1)) * propMat(age_ndx, year_ndx - 1) * stockMeanWeight(age_ndx, year_ndx - 1);
      pred_catches(year_ndx - 1) += exp(log(N(age_ndx, year_ndx - 1)) + log(F_ay(age_ndx, year_ndx - 1)) - log(Z_ay(age_ndx, year_ndx - 1)) + log(1 - exp(-Z_ay(age_ndx, year_ndx - 1)))) * catchMeanWeight(age_ndx, year_ndx - 1);
    }
  }
  // Projections
  vector<Type> partition = N.col(nyears);
  
  vector<Type> ssb_proj(n_proj_years + 1);
  vector<Type> proj_pred_catches(n_proj_years);
  vector<Type> proj_F(n_proj_years);
  
  proj_pred_catches.setZero();
  ssb_proj.setZero();
  ssb_proj(0) = ssb(nyears - 1);
  vector<Type> predN(nages);
  predN.setZero();
  
  if(proj_method == 0) {
    vector<Type> proj_F_by_age = future_vals * fishery_selectivity;
    vector<Type> proj_Z_by_age = future_vals * fishery_selectivity + natMor;
    for(year_ndx = 1; year_ndx <= n_proj_years; year_ndx++) {
      proj_F(year_ndx - 1) = future_vals;
      //----------------
      //Recuritment
      //----------------
      if(stockRecruitmentModelCode == 0) { // straight RW 
        predN(0) = (partition(0) * future_ycs(year_ndx - 1));
      } else if(stockRecruitmentModelCode == 2) { //BH
        predN(0) = (R0 * BevertonHolt(ssb_proj(year_ndx - 1), B0, steepness) * future_ycs(year_ndx - 1));
      } else{
        error("SR model code not recognized");
      }
      
      //----------------
      // F + M + ageing
      //----------------
      for(age_ndx = 1; age_ndx < nages; ++age_ndx) 
        predN(age_ndx) = partition(age_ndx - 1) * exp( - proj_Z_by_age(age_ndx - 1));
      
      if(maxAgePlusGroup == 1) {
        predN(nages - 1) = partition(nages - 2) * exp(-proj_Z_by_age(nages - 2)) + partition(nages - 1) * exp(-proj_Z_by_age(nages - 1));
      }  
      
      // Calculate SSBs an interpolation bewtween the year, starting with previous years Paritition
      for(age_ndx = 0; age_ndx < nages; ++age_ndx) {
        ssb_proj(year_ndx) += partition(age_ndx) * exp(-proj_Z_by_age(age_ndx) * propZ_ssb(nyears - 1)) * propMat(age_ndx, nyears - 1) * stockMeanWeight(age_ndx, nyears - 1);
        proj_pred_catches(year_ndx - 1) += proj_F_by_age(age_ndx) / proj_Z_by_age(age_ndx) * partition(age_ndx) * (1 - exp(-proj_Z_by_age(age_ndx))) * catchMeanWeight(age_ndx, nyears - 1);
      }
      partition = predN;
    }
  } else if (proj_method == 1) {
    vector<Type> vulnerable_proj(n_proj_years);
    vulnerable_proj.setZero();
    vector<Type> proj_F_by_age;
    vector<Type> proj_Z_by_age;
    Type u_proj;
    Type plus_group = 0.0;
    // constant catch
    for(year_ndx = 1; year_ndx <= n_proj_years; year_ndx++) {
      //----------------
      // F + M + ageing
      //----------------
      for(age_ndx = 0; age_ndx < nages; ++age_ndx)
        vulnerable_proj(year_ndx - 1) += partition(age_ndx) * exp( -0.5 * natMor) * catchMeanWeight(age_ndx, nyears - 1) * fishery_selectivity(age_ndx);
      u_proj = future_vals / vulnerable_proj(year_ndx - 1);
      if (u_proj > 0.95)
        u_proj = 0.95;
      proj_F(year_ndx - 1) = -log(1 - u_proj);
      proj_pred_catches(year_ndx - 1) = future_vals;
      proj_F_by_age = proj_F(year_ndx - 1) * fishery_selectivity;
      proj_Z_by_age = proj_F(year_ndx - 1) * fishery_selectivity + natMor;
      //----------------
      //Recuritment
      //----------------
      if(stockRecruitmentModelCode == 0) { // straight RW 
        predN(0) = (partition(0) * future_ycs(year_ndx - 1));
      } else if(stockRecruitmentModelCode == 2) { //BH
        predN(0) = (R0 * BevertonHolt(ssb_proj(year_ndx - 1), B0, steepness) * future_ycs(year_ndx - 1));
      } else{
        error("SR model code not recognized");
      }
      //----------------
      // F + M + ageing
      //----------------
      for(age_ndx = 1; age_ndx < nages; ++age_ndx) 
        predN(age_ndx) = partition(age_ndx - 1) * exp( - proj_Z_by_age(age_ndx - 1));
      
      if(maxAgePlusGroup == 1) {
        predN(nages - 1) = partition(nages - 2) * exp(-proj_Z_by_age(nages - 2)) + partition(nages - 1) * exp(-proj_Z_by_age(nages - 1));
      }  
      // Calculate SSBs an interpolation bewtween the year, starting with previous years Paritition
      for(age_ndx = 0; age_ndx < nages; ++age_ndx) {
        ssb_proj(year_ndx) += partition(age_ndx) * exp(-proj_Z_by_age(age_ndx) * propZ_ssb(nyears - 1)) * propMat(age_ndx, nyears - 1) * stockMeanWeight(age_ndx, nyears - 1);
        proj_pred_catches(year_ndx - 1) += proj_F_by_age(age_ndx) / proj_Z_by_age(age_ndx) * partition(age_ndx) * (1 - exp(-proj_Z_by_age(age_ndx))) * catchMeanWeight(age_ndx, nyears - 1);
      }
      partition = predN;
    }
    REPORT( vulnerable_proj );
  }
  
 
 
 
 
 
  //std::cout << "Joint ll = " << joint_ll << " catch pen1 = " << catch_pen << " catch pen2 = " << catch_pen2 <<"\n";
  REPORT( ssb );
  REPORT( R0 );
  REPORT( sigma_r );
  REPORT( B0 );
  //REPORT( annual_Fs );
  REPORT( pred_catches );

  REPORT( fishery_selectivity );
  REPORT( survey_selectivity );

  REPORT( N );
  REPORT( ycs );
  REPORT( future_ycs );
  
  REPORT( survey_a50 );
  REPORT( survey_ato95 ); 
  REPORT( survey_Q );
  REPORT( f_a50 );
  REPORT( f_ato95 );
  REPORT( annual_F );
  
  REPORT( F_ay );  
  REPORT( Z_ay );
  ADREPORT( ssb );

  REPORT(ssb_proj);
  REPORT(proj_pred_catches);
  REPORT( proj_F );
  
  REPORT(stockMeanWeight);
  REPORT(catchMeanWeight);
  
  // ADREPORT(logN);
  // ADREPORT(logF);
  return 0.0; // we are minimising with nlminb
  //return 0.0;
}
