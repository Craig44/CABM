/**
 * @file Categorical.h
 * @author  C.Marsh
 * @date 26/07/2018
 * @section LICENSE
 * @description
 * All numeric categorical will inherit from this class
 *
 */
#ifndef CATEGORICAL_LAYER_H_
#define CATEGORICAL_LAYER_H_

// Local Headers
#include "Layers/Layer.h"

namespace niwa {
namespace layers {

class CategoricalLayer : public niwa::Layer {
public:
  // Functions
  explicit CategoricalLayer(Model* model);
  virtual                     ~CategoricalLayer();
  virtual void                 set_value(unsigned RowIndex, unsigned ColIndex, string Value);  // TODO once we have children make this pure virtual
  virtual string               get_value(unsigned RowIndex, unsigned ColIndex);


protected:
  //member
  virtual void                DoValidate();
  virtual void                DoBuild();
  // Variables
  string                      **grid_;
  parameters::Table*          data_table_ = nullptr;
  bool                        proportion_;

};

} /* namespace layers */
} /* namespace niwa */

#endif /*CATEGORICAL_LAYER_H_*/
