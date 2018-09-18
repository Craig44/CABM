/**
 * @file NumericLayer.h
 * @author  C.Marsh
 * @date 12/07/2018
 * @section LICENSE
 * @description
 * All numeric layers will inherit from this class
 *
 */
#ifndef NUMERIC_LAYER_H_
#define NUMERIC_LAYER_H_

// Local Headers
#include "Layers/Layer.h"

namespace niwa {
namespace layers {

class NumericLayer : public niwa::Layer {
public:
  // Functions
  explicit NumericLayer(Model* model);
  virtual                     ~NumericLayer();
  virtual void                set_value(unsigned RowIndex, unsigned ColIndex, float Value);  // TODO once we have children make this pure virtual
  virtual float               get_value(unsigned RowIndex, unsigned ColIndex);
  virtual float               get_value(unsigned RowIndex, unsigned ColIndex, unsigned year);
  bool                        is_static() {return static_layer_;};
  bool                        is_proportion() {return proportion_;};


protected:
  //member
  virtual void                DoValidate();
  virtual void                DoBuild();
  // Variables
  float                       **grid_;
  parameters::Table*          data_table_ = nullptr;
  bool                        proportion_;
  bool                        static_layer_; // is this layer pre-defined (user input) or calculated on the fly

};

} /* namespace layers */
} /* namespace niwa */

#endif /*INTLAYER_H_*/
