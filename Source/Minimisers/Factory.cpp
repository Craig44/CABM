/**
 * @file Factory.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 28/02/2013
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers
#include "Factory.h"

#include "Model/Model.h"
#include "Model/Managers.h"
#include "Minimisers/Manager.h"
#include "Minimisers/Children/Dummy/Dummy.h"

#include "Minimisers/Children/GammaDiff.h"
//#include "Minimisers/Children/STANBFGS.h"

// Namespaces
namespace niwa {
namespace minimisers {

/**
 * Create the instance of our object as defined by the two parameters
 * object_type and sub_type.
 *
 * @param object_type The type of object to create (e.g age_size, process)
 * @param sub_type The child type of the object to create (e.g ageing, schnute)
 * @return shared_ptr to the object we've created
 */
Minimiser* Factory::Create(Model* model, const string& object_type, const string& sub_type) {
  Minimiser* result = nullptr;

  if (object_type == PARAM_MINIMIZER) {
    if (sub_type == PARAM_GAMMADIFF)
      result = new GammaDiff(model);

    if (result)
      model->managers().minimiser()->AddObject(result);
  }
  return result;
}

} /* namespace minimisers */
} /* namespace niwa */
