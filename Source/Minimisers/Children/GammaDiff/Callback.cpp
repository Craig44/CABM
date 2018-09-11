/**
 * @file BaranovCallBack.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 17/04/2013
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
// headers
#include "Callback.h"


// namespaces
namespace niwa {
namespace minimisers {
namespace gammadiff {

/**
 * Default Constructor
 */
CallBack::CallBack(Model* model) : model_(model) {
}

//**********************************************************************
// double CGammaDiffCallback::operator()(const vector<double>& Parameters)
// Operatior() for Minimiser CallBack
//**********************************************************************
double CallBack::operator()(const vector<double>& Parameters) {
  // I have done a simple linear regression problem to help me understand how the minimiser works
  // solves a and b in the equation y = a + bX;
  params_[0] = Parameters[0];
  params_[1] = Parameters[1];

  data_hat_.resize(data_.size(),0);
  double SSE = 0;
  for(unsigned i = 0; i < data_.size(); ++i) {
    data_hat_[i]  = Parameters[0] +  Parameters[1] * covariate_[i];
    SSE += pow(data_[i] - data_hat_[i], 2);
  }

  LOG_FINE() << "a = " <<  Parameters[0] << " B = " <<  Parameters[1] << " SSE = "  << SSE;
  return SSE;
}

} /* namespace gammadiff */
} /* namespace minimiser */
} /* namespace niwa */
