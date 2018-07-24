//============================================================================
// Name        : NumericMetaLayer.h
// Author      : S.Rasmussen
// Date        : 16/01/2009
// Copyright   : Copyright NIWA Science ©2009 - www.niwa.co.nz
// Description :
// $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
//============================================================================
#ifndef NUMERIC_META_LAYER_H_
#define NUMERIC_META_LAYER_H_

// headers
#include "../NumericLayer.h"

namespace niwa {
namespace layers {

//**********************************************************************
//
//
//**********************************************************************
class NumericMetaLayer : public NumericLayer {
public:
  // Functions
  NumericMetaLayer(Model* model);
  virtual                     ~NumericMetaLayer() {};
  virtual void                set_value(unsigned RowIndex, unsigned ColIndex, float Value) {};  // TODO once we have children make this pure virtual
  virtual float               get_value(unsigned RowIndex, unsigned ColIndex);
protected:
  virtual void                DoValidate() override final;
  virtual void                DoBuild() override final;
  // Variables
  string                       default_Layer_label_;
  vector<unsigned>             years_;
  vector<string>               layer_names_;
  map<unsigned, NumericLayer*> years_layer_;
  NumericLayer*                default_layer_;
  bool                         has_years_;
};
} /* namespace layers */
} /* namespace niwa */
#endif /* NUMERIC_META_LAYER_H_ */
