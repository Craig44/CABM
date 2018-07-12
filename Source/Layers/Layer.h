/**
 * @file Layer.h
 * @author  C.Marsh
 * @date 12/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 *       This class is responsible for maintaining a list of our
 *       Layers. Like all classes it is responsible for the validate and
 *       build calls to the children.
 *
 *       Layers are a very dynamic concept. During the execution of a process
 *       a layer can be defined. This layer is overlayed against the World
 *       grid and a check against it's parameters will determine if this
 *       World-Cell is a valid one for this process or not.
 *
 *       i.e. We have a layer called Depth. This tells us the depth of water
 *       at each World-Square. We know that the fish we are looking for
 *       only hang out in Water Depth 250-300M. So we can put a layer_min
 *       and layer_max of 250 and 300. Then we can apply this layer to a
 *       process. We know this process (e.g Recruitment) will only affect the
 *       World-Squares with a depth of 250-300 (Inclusive).
 *
 *       Different types of layers have different functionalities. You can
 *       check individual functionalites in
 */
#ifndef LAYER_H_
#define LAYER_H_

// headers
#include "BaseClasses/Object.h"
#include "Model/Model.h"
// namespaces
namespace niwa {

class Model;

// classes
class Layer : public niwa::base::Object {
public:
  // methods
  explicit Layer(Model* model);
  virtual                     ~Layer() = default;
  void                        Validate();
  void                        Build();
  void                        Reset();

  // pure methods
  virtual void                DoValidate() = 0;
  virtual void                DoBuild() = 0;

  // accessors


protected:
  // Members
  Model*                      model_ = nullptr;
  unsigned                    height_;
  unsigned                    width_;
};
} /* namespace niwa */
#endif /* LAYER_H_ */
