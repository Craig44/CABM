/*
 * LogisticNormal.h
 *
 *  Created on: Oct 26, 2016
 *      Author: Zaita
 */

#ifndef SOURCE_LIKELIHOODS_CHILDREN_LOGISTICNORMAL_H_
#define SOURCE_LIKELIHOODS_CHILDREN_LOGISTICNORMAL_H_

#include "Likelihoods/Likelihood.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

namespace niwa {
namespace likelihoods {
namespace ublas = boost::numeric::ublas;


class LogisticNormal : public niwa::Likelihood {
public:
  // Methods
  LogisticNormal(Model* model);
  virtual                     ~LogisticNormal();
  void                        DoValidate() override final;
  void                        SimulateObserved(map<unsigned, map<string, vector<observations::Comparison> > >& comparisons) override final;
  virtual void                FillReportCache(ostringstream& cache);  // If we want to store and report more information within a process use this method

protected:
  // Estimable parameters
  float                       sigma_;
  vector<unsigned>            bins_;
  vector<float>               rho_;
  bool                        arma_;
  bool                        robust_;
  unsigned                    n_bins_;
  bool                        sep_by_sex_;
  bool                        sex_lag_;
  bool                        sexed_;
  unsigned                    unique_bins_;

  // Covariance containers
  ublas::matrix<float>       covariance_matrix_;
  ublas::matrix<float>       covariance_matrix_lt;
  parameters::Table*         covariance_table_ = nullptr;
  // Methods
  void                       calculate_covariance();
  vector<float>              GetRho(vector<float>& Phi, unsigned nBin, bool ARMA);
  vector<float>              RecursiveFilter(vector<float>& ar_coef, unsigned nBins, vector<float>& initial_vals);
  bool                       DoCholeskyDecmposition();
  bool                       InvertMatrix(const ublas::matrix<float>& input, ublas::matrix<float>& inverse);
  float                      det_fast(const ublas::matrix<float>& matrix);
};

} /* namespace likelihoods */
} /* namespace niwa */
#endif /* SOURCE_LIKELIHOODS_CHILDREN_LOGISTICNORMAL_H_ */
