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
#include "Layers/Children/IntLayer.h"
#include "Layers/Children/NumericLayer.h"
#include "Layers/Children/CategoricalLayer.h"
#include "Layers/Children/Numeric/NumericMetaLayer.h"


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
  if (object_type == PARAM_LAYER) {
    if (sub_type == PARAM_INTEGER)
      result = new IntLayer(model);
    else if (sub_type == PARAM_NUMERIC)
      result = new NumericLayer(model);
    else if (sub_type == PARAM_CATEGORICAL)
      result = new CategoricalLayer(model);
    else if (sub_type == PARAM_NUMERIC_META)
      result = new NumericMetaLayer(model);

    if (result)
      model->managers().layer()->AddObject(result);
  }




  return result;
}


} /* namespace layers */
} /* namespace niwa */
