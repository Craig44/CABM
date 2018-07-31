/**
 * @file Biomass.h
 * @author  C.Marsh
 * @date 31/07/2018
 * @section LICENSE
 * @description
 *
 */

#ifndef NUMERIC_BIOMASS_H_
#define NUMERIC_BIOMASS_H_

// headers
#include "../NumericLayer.h"

namespace niwa {
class WorldView;
namespace layers {

//**********************************************************************
//
//
//**********************************************************************
class Biomass : public NumericLayer {
public:
  // Functions
  Biomass(Model* model);
  virtual                     ~Biomass() {};
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
#endif /* NUMERIC_BIOMASS_H_ */
