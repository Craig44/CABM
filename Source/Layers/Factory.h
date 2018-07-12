/**
 * @file Factory.h
 * @author  C.Marsh
 * @date 12/7/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * Standard factory class for the layers
 */
#ifndef LAYERS_FACTORY_H_
#define LAYERS_FACTORY_H_

// namespaces
#include "Layers/Layer.h"

// namespaces
namespace niwa {
class Model;

namespace layers {

// classes
class Factory {
public:
  // methods
  static Layer*     Create(Model* model, const string& object_type, const string& sub_type);

private:
  // methods
  Factory() = delete;
  virtual ~Factory() = delete;
};

} /* namespace layers */
} /* namespace niwa */
#endif /* LAYERS_FACTORY_H_ */
