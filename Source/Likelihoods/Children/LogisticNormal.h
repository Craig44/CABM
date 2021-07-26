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
  double                       sigma_;
  vector<unsigned>            bins_;
  vector<double>               rho_;
  bool                        arma_;
  bool                        robust_;
  unsigned                    n_bins_;
  bool                        sep_by_sex_;
  bool                        sex_lag_;
  bool                        sexed_;
  unsigned                    unique_bins_;

  // Covariance containers
  ublas::matrix<double>       covariance_matrix_;
  ublas::matrix<double>       covariance_matrix_lt;
  parameters::Table*         covariance_table_ = nullptr;
  // Methods
  void                       calculate_covariance();
  vector<double>              GetRho(vector<double>& Phi, unsigned nBin, bool ARMA);
  vector<double>              RecursiveFilter(vector<double>& ar_coef, unsigned nBins, vector<double>& initial_vals);
  bool                       DoCholeskyDecmposition();
  bool                       InvertMatrix(const ublas::matrix<double>& input, ublas::matrix<double>& inverse);
  double                      det_fast(const ublas::matrix<double>& matrix);
};

} /* namespace likelihoods */
} /* namespace niwa */
#endif /* SOURCE_LIKELIHOODS_CHILDREN_LOGISTICNORMAL_H_ */
