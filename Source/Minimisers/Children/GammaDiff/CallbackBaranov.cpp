/**
 * @file BaranovCallBack.cpp
 * @author  C.Marsh
 * @version 1.0
 * @date 17/04/2013
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
// headers
#include "CallbackBaranov.h"


// namespaces
namespace niwa {
namespace minimisers {
namespace gammadiff {

/**
 * Default Constructor
 */
CallbackBaranov::CallbackBaranov(Mortality* mortality) : mortality_(mortality) {

}

//**********************************************************************
// double CGammaDiffCallback::operator()(const vector<double>& Parameters)
// Operatior() for Minimiser CallBack
//**********************************************************************
double CallbackBaranov::operator()(const vector<double>& Parameters) {
  LOG_FINE() << "Setup Baranov equation";
  mortality_->set_lambda(Parameters[0]);
  double catch_SSE = mortality_->SolveBaranov();
  return catch_SSE;
}

} /* namespace gammadiff */
} /* namespace minimiser */
} /* namespace niwa */
