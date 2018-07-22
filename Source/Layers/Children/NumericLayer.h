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
  int                         countValidSpaces();
  void                        set_value(unsigned RowIndex, unsigned ColIndex, float Value);  // TODO once we have children make this pure virtual
  float                       get_value(unsigned RowIndex, unsigned ColIndex);


protected:
  //member
  virtual void                DoValidate() override final;
  virtual void                DoBuild() override final;
  // Variables
  float                       **grid_;
  parameters::Table*          data_table_ = nullptr;
  bool                        proportion_;

};

} /* namespace layers */
} /* namespace niwa */

#endif /*INTLAYER_H_*/
