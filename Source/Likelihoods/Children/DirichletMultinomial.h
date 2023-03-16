//============================================================================
// Name        : DirichletMultinomial.h
// Author      : C.Marsh
// Date        : 11/05/22
// Copyright   : Copyright NIWA Science 2020 - www.niwa.co.nz
// Description :
//============================================================================

// DirichletMultinomial  error distribution

#ifndef LIKELIHOODS_DIRICHLET_MULTINOMIAL_H_
#define LIKELIHOODS_DIRICHLET_MULTINOMIAL_H_

// Headers
#include "../Likelihood.h"
// Namespaces
namespace niwa {
namespace likelihoods {

/**
 * Class definition
 */
class DirichletMultinomial : public niwa::Likelihood {
public:
  // Methods
  DirichletMultinomial(Model* model);
  virtual ~DirichletMultinomial() = default;
  void   DoValidate() override final;
  void   SimulateObserved(map<unsigned, map<string, vector<observations::Comparison> > >& comparisons) override final;

protected:
  double                      theta_;
  bool                        enter_cell_ = true;

};

} /* namespace likelihoods */
} /* namespace niwa */
#endif /* LIKELIHOODS_DIRICHLET_MULTINOMIAL_H_ */
