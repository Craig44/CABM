/**
 * @file IntLayer.h
 * @author  C.Marsh
 * @date 12/07/2018
 * @section LICENSE
 * @description
 * Other Layer types will have a Base class, but as they will have more childrent, I can only think of one base class (Base) perhaps Abundance could be TODO consider later
 *
 */
#ifndef INTLAYER_H_
#define INTLAYER_H_

// Local Headers
#include "Layers/Layer.h"

namespace niwa {
namespace layers {

class IntLayer : public niwa::Layer {
public:
  // Functions
  explicit IntLayer(Model* model);
  virtual                     ~IntLayer();
  int                         countValidSpaces();
  void                        set_value(int RowIndex, int ColIndex, unsigned Value);
  unsigned                    get_value(int RowIndex, int ColIndex);


protected:
  //member
  virtual void                DoValidate() override final;
  virtual void                DoBuild() override final;
  // Variables
  unsigned                    **grid_;
  parameters::Table*          int_table_ = nullptr;

};

} /* namespace layers */
} /* namespace niwa */

#endif /*INTLAYER_H_*/
