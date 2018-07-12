/**
 * @file Factory.cpp
 * @author  C.Marsh
 * @date 12/7/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "Factory.h"

#include "Model/Model.h"
#include "Model/Managers.h"
#include "Layers/Manager.h"

// namespaces
namespace niwa {
namespace layers {

/**
 * Create the instance of our object as defined by the two parameters
 * object_type and sub_type.
 *
 * @param object_type The type of object to create (e.g age_size, process)
 * @param sub_type The child type of the object to create (e.g ageing, schnute)
 * @return shared_ptr to the object we've created
 */
Layer* Factory::Create(Model* model, const string& object_type, const string& sub_type) {
  Layer* result = nullptr;


  if (result)
    model->managers().layer()->AddObject(result);

  return result;
}


} /* namespace layers */
} /* namespace niwa */
