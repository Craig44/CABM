/**
 * @file Manager.cpp
 * @author  C.Marsh
 * @date 12/7/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "Manager.h"

// namespaces
namespace niwa {
namespace layers {

/**
 * Return a named layer from our collection
 *
 * @param label The label of the layer
 * @return pointer to quantity, or empty pointer if not found
 */
Layer* Manager::GetLayer(const string& label) {
  for (auto layer : objects_) {
    if (layer->label() == label)
      return layer;
  }
  return nullptr;
}

} /* namespace layers */
} /* namespace niwa */
