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

/**
 * Return a int layer from our collection
 *
 * @param label The label of the layer
 * @return pointer to quantity, or empty pointer if not found
 */
IntLayer* Manager::GetIntLayer(const string& label) {
  IntLayer* pPtr = nullptr;
  for (auto layer : objects_) {
    if (layer->label() == label && (layer->layer_type() == LayerType::kInteger)) {
      pPtr = dynamic_cast<IntLayer*>(layer);
      return pPtr;
    }
  }
  return pPtr;
}

/**
 * Return a numeric layer from our collection
 *
 * @param label The label of the layer
 * @return pointer to quantity, or empty pointer if not found
 */
NumericLayer* Manager::GetNumericLayer(const string& label) {
  LOG_TRACE();
  NumericLayer* pPtr = nullptr;
  for (auto layer : objects_) {
    if ((layer->label() == label) && (layer->layer_type() == LayerType::kNumeric)) {
      pPtr = dynamic_cast<NumericLayer*>(layer);
      return pPtr;
    }
  }
  return pPtr;
}

CategoricalLayer* Manager::GetCategoricalLayer(const string& label) {
  LOG_TRACE();
  CategoricalLayer* pPtr = nullptr;
  for (auto layer : objects_) {
    if ((layer->label() == label) && (layer->layer_type() == LayerType::kCategorical)) {
      pPtr = dynamic_cast<CategoricalLayer*>(layer);
      return pPtr;
    }
  }
  return pPtr;
}


/**
 * override Base classes and Build all non integer and numeric layers
 */
void Manager::BuildPostWorldLayers() {
  LOG_FINEST() << "Starting Build... with " << objects_.size() << " objects";
  for(auto stored_object : objects_) {
    if ((stored_object->layer_type() != LayerType::kInteger) || ((stored_object->layer_type() == LayerType::kNumeric) && (stored_object->type() == PARAM_BIOMASS))) {
      LOG_FINEST() << "Building process = " << stored_object->label();
      stored_object->Build();
    }
  }
  LOG_FINEST() << "Build Finished";
}

/**
 * override Base classes and Build mortality and growth processes
 */
void Manager::BuildPreWorldLayers() {
  LOG_FINEST() << "Starting Build... with " << objects_.size() << " objects";
  for(auto stored_object : objects_) {
    if ((stored_object->layer_type() == LayerType::kInteger) || ((stored_object->layer_type() == LayerType::kNumeric) && (stored_object->type() != PARAM_BIOMASS))) {
      stored_object->Build();
      LOG_FINEST() << "Building process = " << stored_object->label();
    }
  }
  LOG_FINEST() << "Build Finished";
}


} /* namespace layers */
} /* namespace niwa */
