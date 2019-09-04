/**
 * @file Math.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 25/03/2013
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * The time class represents a moment of time.
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef UTILITIES_MATH_H_
#define UTILITIES_MATH_H_


// Headers
#include <cmath>

#include "Utilities/DoubleCompare.h"
#include "Utilities/Types.h"
#include <boost/math/distributions/normal.hpp>
#include <numeric>
#define PI 3.14159265358979

// Namespaces
namespace niwa {
namespace utilities {
namespace math {

namespace dc = niwa::utilities::doublecompare;
using niwa::utilities::Double;

/**
 * LnGamma
 */
inline double LnGamma(double t) {
  double x, y, tmp, ser;
  double cof[6] = {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  y = x = t;
  tmp = x + 5.5 - (x + 0.5) * log(x + 5.5);
  ser = 1.000000000190015;
  for (unsigned j = 0; j <= 5; j++) {
    y = y + 1.0;
    ser += (cof[j] / y);
  }
  return(log(2.5066282746310005 * ser / x) - tmp);
}

/**
 * LnFactorial
 */
inline double LnFactorial(double t) {
  return niwa::utilities::math::LnGamma(t + 1.0);
}
/*
 *
*/

/**
 * Normal Distribution CDF Method
 *
 * @param x X value
 * @param mu Mu value
 * @param sigma Sigma value
 * @return Normal CDF
 */

inline float NormalCDF(float x, float mu, float sigma) {
  if (sigma <= 0.0 && x < mu)
    return 0;
  else if (sigma <= 0.0 && x >= mu)
    return 1;

  boost::math::normal s(mu, sigma);
  return cdf(s, (x));
}


/**
 * dnorm: return the pdf for the normal
 */
inline double dnorm(const double& x, const double& mu, const double& sigma) {
  double z = 1 / (sigma * sqrt(2 * PI)) * exp(-((x - mu) * (x - mu))/(2 * sigma * sigma));
  return(z);
}

//************************************
// These are distributional functions taken from CASAL, some will wont to be updated
//************************************
/*
 * dlognorm: return the pdf for the log-normal
 */
inline double dlognorm(const double& x, const double& mu = 0.0, const double& sigma = 1.0) {
  // Parameterised by the mean and standard deviation of the (normal) distribution
  //  of log(x), NOT by those of the (lognormal) distribution of x.
  if (x <= 0)
    return 0.0;
  else
    return dnorm(log(x),mu,sigma);
}

/**
 * pnorm: return the cdf for the normal
 */
inline double pnorm(const double& x, const double& mu = 0.0, const double& sigma = 1.0) {
  // Abramowitz & Stegun eqn 26.2.18
  // Equations: z = fabs((x-mu)/sigma);
  //            p = 1-0.5*pow((1+0.196854*z+0.115194*z*z+0.000344*z*z*z+0.019527*z*z*z*z),-4);
  //            if (x<mu) p=1-p;
  double z = fabs((x - mu)/sigma);
  double p = 1 - 0.5*pow((1+0.196854*z+0.115194*z*z+0.000344*z*z*z+0.019527*z*z*z*z),-4);
  if (x < mu)
    p = 1 - p;
  return(p);
}

/**
 * given class bins (lower and upper bin definitions) iterate over all bins and apportion the CDF into each bin
 */
inline vector<float> block_cdf(const vector<float>& class_bins, const double& mu = 0.0, const double& sigma = 1.0) {
  vector<float> prob(class_bins.size(), 0.0);
  vector<float> prob_bins(class_bins.size() - 1, 0.0);
  unsigned index = 0;
  bool check_max_of_one = false;
  for(auto& val : class_bins) {
    prob[index] = pnorm(val, mu, sigma);
    //LOG_FINE() << "prob = " << prob[index] << " val = " << val;
    if (prob[index] >= 1.0)
      check_max_of_one = true;
    ++index;
  }
  if (not check_max_of_one)
    prob[class_bins.size() - 1] += (1.0 - prob[class_bins.size() - 1]);
  // calculate difference between bins
  for(unsigned i = 0; i < (class_bins.size() - 1); ++i)
    prob_bins[i] = prob[i + 1] - prob[i];
  if (prob[0] > 0.0)
    prob_bins[0] += prob[0];
  return(prob_bins);
}
/**
 * pnorm2: return the cdf for the normal that may be more expensive (computationally) but is a better approximation.
 */
inline double pnorm2(const double& x, const double& mu = 0.0, const double& sigma = 1.0) {
  // Ian Doonan's code, A better approximation of the normal CDF as there is no closed form
  double norm, ttt, p;
  double z = fabs((x - mu)/sigma);
  double tt = 1.0 / (1.0 + 0.2316419 * z);
  norm = 1.0 / sqrt(2.0 * PI) * exp(-0.5 * z * z);
  ttt = tt;
  p = 0.319381530 * ttt;
  ttt = ttt * tt;
  p = p - 0.356563782 * ttt;
  ttt = ttt * tt;
  p = p + 1.781477937 * ttt;
  ttt = ttt * tt;
  p = p - 1.821255978 * ttt;
  ttt = ttt * tt;
  p = p + 1.330274429 * ttt;
  p *= norm;
  if (x < mu)
    p = 1-p;
  return(p);
}

/**
 * plnorm: return the cdf for the normal
 */
inline double plognorm(const double& x, const double& mu, const double& sigma) {
  // Parameterised by the mean and standard deviation of the (normal) distribution
  //  of log(x), NOT by those of the (lognormal) distribution of x.
  if (x <= 0)
    return 0.0;
  else
    return pnorm(log(x),mu,sigma);
}

// Distribution: this is taken from CASAL, will be useful for Samik and his length structured stuff.
// Probability distribution of a random variable over classes.
// distribution(class_mins,plus_group,dist,mean,stdev)[i]
// is the probability that a random variable with distribution 'dist', 'mean' and 'stdev'
// falls between class_mins[i] and class_mins[i+1]
// (unless i = class_mins.indexmax() and plus_group!=0, in which case
// distribution(...)[i] is the probability that the random variable exceeds class_mins[i]).
//
// We use an approximation: P(X is more than 5 std.devs away from its mean) = 0.
//  Almost true for the normal distribution, but may be problematic if you use something more skewed.
inline vector<double> distribution(const vector<double>& class_mins, int plus_group = 0, string dist = "normal", const double& mean = 0.0, const double& stdev = 1.0) {
  int n_bins = class_mins.size() - (plus_group ? 0 : 1);
  vector<double> result(n_bins, 0.0);
  double so_far;
  double mu, sigma;
  if (dist == PARAM_LOGNORMAL){
    sigma = sqrt(log(1+pow(stdev/mean,2)));
    mu = log(mean) - sigma*sigma/2;
  }
  LOG_TRACE();
  if (class_mins[0] < mean - 5 * stdev){
    so_far = 0;
  } else {
    if (dist == PARAM_NORMAL){
      so_far = pnorm(class_mins[0],mean,stdev);
    } else if (dist == PARAM_LOGNORMAL){
      so_far = plognorm(class_mins[0],mu,sigma);
    } else
      LOG_CODE_ERROR() << "unknown distribution supplies '" << dist;
  }
  LOG_TRACE();

  int c;
  for (c = 0; c < (n_bins - 1); c++){
    if (class_mins[c + 1] > mean + 5 * stdev){
      result[c] = 1-so_far;
      so_far = 1;
    } else if (class_mins[c + 1] < mean - 5 * stdev){
      result[c] = 0;
      so_far = 0;
    } else {
      if (dist == PARAM_NORMAL){
        result[c] = pnorm(class_mins[c + 1], mean ,stdev) - so_far;
      } else if (dist == PARAM_LOGNORMAL){
        result[c] = plognorm(class_mins[c + 1], mu, sigma) - so_far;
      }
      so_far += result[c];
    }
    if (result[c] < 0 || result[c]!=result[c]) {
      LOG_CODE_ERROR() << "bad result in distribution, got " << result[c];
    }
  }
  LOG_TRACE();

  c = n_bins - 1;
  if (plus_group){
    result[c] = 1 - so_far;
  } else {
    if (class_mins[c + 1] > mean + 5 * stdev){
      result[c] = 1 - so_far;
    } else {
      if (dist == PARAM_NORMAL){
        result[c] = pnorm(class_mins[c + 1], mean, stdev) - so_far;
      } else if (dist == PARAM_LOGNORMAL){
        result[c] = plognorm(class_mins[c + 1], mu, sigma) - so_far;
      }
    }
    if (result[c] < 0 || result[c] != result[c]){
      LOG_CODE_ERROR() << "bad result in distribution, got " << result[c];
    }
  }
  LOG_TRACE();
  return result;
}
//**********************************************************************
// void Engine::condassign( double &res, const double &cond, const double &arg1, const double &arg2 ) {
// Conditional Assignment
//**********************************************************************
inline void cond_assign(double &res, const double &cond, const double &arg1, const double &arg2) {
  res = (cond) > 0 ? arg1 : arg2;
}

//**********************************************************************
// void Engine::condassign( double &res, const double &cond, const double &arg)
// Conditional Assignment
//**********************************************************************
inline void cond_assign(double &res, const double &cond, const double &arg) {
  res = (cond) > 0 ? arg : res;
}

/**
 * double Engine::boundpin(double y, double fmin, double fmax)
 * Boundary Pin
 */

inline double scale_value(double value, double min, double max) {
  if (dc::IsEqual(value, min))
    return -1;
  else if (dc::IsEqual(value, max))
    return 1;

  return asin(2 * (value - min) / (max - min) - 1) / 1.57079633;
}

/**
 *
 */
inline double unscale_value(const double& value, double& penalty, double min, double max) {
  // courtesy of AUTODIF - modified to correct error -
  // penalty on values outside [-1,1] multiplied by 100 as of 14/1/02.
  double t = 0.0;
  double y = 0.0;

  t = min + (max - min) * (sin(value * 1.57079633) + 1) / 2;
  cond_assign(y, -.9999 - value, (value + .9999) * (value + .9999), 0);
  penalty += y;
  cond_assign(y, value - .9999, (value - .9999) * (value - .9999), 0);
  penalty += y;
  cond_assign(y, -1 - value, 1e5 * (value + 1) * (value + 1), 0);
  penalty += y;
  cond_assign(y, value - 1, 1e5 * (value - 1) * (value - 1), 0);
  penalty += y;

  return (t);
}


//**********************************************************************
//    General math utilities
//**********************************************************************

// Check if a vecotr contains all ones
inline bool all_ones(const vector<double>& x) {
  for(auto num : x) {
    if (num != 1.0)
      return false;
  }
  return true;
}

// Return the mean for a vector
inline double mean(const vector<double>& Values){
  double mu = 0.0;
  double total = 0.0;
  for (const auto& value : Values)
    total += value;
  double n = Values.size();
  mu = total / n;
  return mu;
}

// Return the mean for a vector
inline double mean(const vector<unsigned>& Values){
  unsigned total = 0.0;
  for (const auto& value : Values)
    total += value;
  double n = Values.size();
  return (double)total / n;
}

// Return the mean for an unsigned map
inline double mean(const map<unsigned, double>& Values){
  double mu = 0.0;
  double total = 0.0;
  for (const auto& value : Values)
    total += value.second;
  double n = Values.size();
  mu = total / n;
  return mu;
}

// Return the Variance for a vector
inline double Var(const vector<double>& Values){
  double mean_ = math::mean(Values);
  double variance = 0;
  for (const auto& value : Values)
    variance += (value - mean_) * (value - mean_);
  double n = Values.size();
  double var = variance / (n - 1.0);
  return var;
}


// Return the Variance for a vector
inline double Var(const vector<unsigned>& Values){
  double mean_ = math::mean(Values);
  double variance = 0;
  for (const auto& value : Values)
    variance += (double)((value - mean_) * (value - mean_));
  double n = Values.size();
  double var = variance / (n - 1.0);
  return var;
}

// Return the Variance for an unsigned map
inline double Var(const map<unsigned, double>& Values){
  double mean_ = math::mean(Values);
  double variance = 0;
  for (const auto& value : Values)
    variance += (value.second - mean_) * (value.second - mean_);
  double n = Values.size();
  double var = variance / (n - 1.0);
  return var;
}

// Return the Standard Deviation for a vector
inline double std_dev(const vector<double>& Values){
  double sd;
  sd = sqrt(math::Var(Values));
  return sd;
}

// Return the Standard Deviation for an unsigned map
inline double std_dev(const map<unsigned, double>& Values){
  double sd;
  sd = sqrt(math::Var(Values));
  return sd;
}

// Return the Sum of a vector
inline double Sum(const vector<double>& Values){
  double total = 0;
  for (auto val : Values)
    total += val;
  return total;
}

// Return the Sum of a vector
inline float Sum(const vector<float>& Values){
  float total = 0;
  for (auto val : Values)
    total += val;
  return total;
}

// Return the maximum value for a vector
inline double Max(const vector<double>& Values){
  double max = 0.0;
  unsigned iter = 1;
  for (auto value : Values) {
    if (iter == 1)
      max = value;
    else if (max < value)
      max = value;
   ++iter;
  }
  return max;
}


// Return the row index and column index of a matrix given an index. This follows R functionality.
// Assumes the index follows byrow = T.
inline vector<int> get_mat_index(int rows, int cols, double index) {
  //unsigned  row_index, col_index;
  vector<int> result;
  int row_ind= floor(index / rows);
  int col_index = index - (row_ind * cols);
  result.push_back(row_ind);
  result.push_back(col_index);
  return result;
}

// Calculate the Sum of a vector
inline vector<double> elem_prod(vector<double>& x,vector<double>& y) {
  if(x.size() != y.size())
    LOG_FATAL() << "method  elem_prod() can only work for equal length vectors";
  vector<double> result(y.size(),0.0);
    for (unsigned i = 0; i< x.size(); ++i) {
      result[i] = x[i] * y[i];
    }
    return result;
}


template <typename T>
inline vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

// Cumulatively sum a vector
template <typename T>
inline vector<T> cumsum(const vector<T> &v) {
  vector<T> result(v.size());
  partial_sum(v.begin(), v.end(), result.begin(), std::plus<T>());
  return result;
}


} /* namespace math */
} /* namespace utilities */
} /* namespace niwa */
#endif /* UTILITIES_MATH_H_ */
