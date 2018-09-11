/**
 * @file Abundance.h
 * @author  C.Marsh
 * @date 10/09/2018
 * @section LICENSE
 * @description
 *
 */

#ifndef NUMERIC_ABUNDANCE_H_
#define NUMERIC_ABUNDANCE_H_

// headers
#include "../NumericLayer.h"

namespace niwa {
class WorldView;
namespace layers {

//**********************************************************************
//
//
//**********************************************************************
class Abundance : public NumericLayer {
public:
  // Functions
  Abundance(Model* model);
  virtual                     ~Abundance() {};
  virtual void                set_value(unsigned RowIndex, unsigned ColIndex, float Value) {};  // TODO once we have children make this pure virtual
  virtual float               get_value(unsigned RowIndex, unsigned ColIndex);
  virtual float               get_value(unsigned RowIndex, unsigned ColIndex, unsigned year);

protected:
  virtual void                DoValidate() override final;
  virtual void                DoBuild() override final;
  // Variables
  bool                        mature_;
  WorldView*                  world_ = nullptr;
};
} /* namespace layers */
} /* namespace niwa */
#endif /* NUMERIC_ABUNDANCE_H_ */
