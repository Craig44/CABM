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

IntLayer* Manager::GetIntLayer(const string& label) {
  IntLayer* pPtr = nullptr;
  for (auto layer : objects_) {
    if (layer->label() == label && layer->type() == PARAM_INTEGER)
      pPtr = dynamic_cast<IntLayer*>(layer);
    return pPtr;
  }
  return nullptr;
}

NumericLayer* Manager::GetNumericLayer(const string& label) {
  NumericLayer* pPtr = nullptr;
  for (auto layer : objects_) {
    if (layer->label() == label && layer->type() == PARAM_NUMERIC)
      pPtr = dynamic_cast<NumericLayer*>(layer);
    return pPtr;
  }
  return nullptr;
}

} /* namespace layers */
} /* namespace niwa */
