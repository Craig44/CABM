/**
 * @file Manager.h
 * @author  C.Marsh
 * @date 12/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This class is the manager for the Layers
 */
#ifndef LAYERS_MANAGER_H_
#define LAYERS_MANAGER_H_

// headers
#include "BaseClasses/Manager.h"

#include "Layers/Layer.h"
#include "Children/IntLayer.h"
#include "Children/NumericLayer.h"
#include "Children/CategoricalLayer.h"


#include "Model/Managers.h"

// namespaces
namespace niwa {
namespace layers {

// classes
class Manager : public niwa::base::Manager<layers::Manager, niwa::Layer> {
  friend class niwa::base::Manager<layers::Manager, niwa::Layer>;
  friend class niwa::Managers;
public:
  // methods
  virtual                     ~Manager() noexcept(true) { };
  Layer*                      GetLayer(const string& label);
  IntLayer*                   GetIntLayer(const string& label);
  NumericLayer*               GetNumericLayer(const string& label);
  CategoricalLayer*           GetCategoricalLayer(const string& label);

  void                        Build() override final { };
  void                        BuildPreWorldLayers();
  void                        BuildPostWorldLayers();

};

} /* namespace layers */
} /* namespace niwa */
#endif /* LAYERS_MANAGER_H_ */
