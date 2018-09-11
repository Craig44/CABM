/**
 * @file CGammaDiffCallback.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 17/04/2013
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
#ifndef MINIMISERS_GAMMADIFF_CALLBACK_H_
#define MINIMISERS_GAMMADIFF_CALLBACK_H_

// Headers
#include <vector>

#include "Model/Model.h"

// namespaces
namespace niwa {
namespace minimisers {
namespace gammadiff {

using std::vector;

/**
 * Class definition
 */
class CallBack {
public:
  CallBack(Model* model);
  virtual                     ~CallBack() = default;
  double                      operator()(const vector<double>& Parameters);

private:
  Model*                    model_;
  vector<double>            data_ = {32.61641, 33.90679, 44.67287, 35.51160, 35.33105, 46.71992, 38.16378, 25.83634, 31.43997, 31.94780};
  vector<double>            covariate_ = {38.79049, 45.39645, 81.17417, 51.41017, 52.58575, 84.30130, 59.21832, 24.69878, 36.26294, 41.08676};
  vector<double>            data_hat_;
  vector<double>            params_ = {0.5,0.5};

};

} /* namespace gammadiff */
} /* namespace minimiser */
} /* namespace niwa */

#endif /* MINIMISERS_GAMMADIFF_CALLBACK_H_ */
