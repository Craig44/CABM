/**
 * @file CategoricalMetaLayer.h
 * @author  C Marsh
 * @version 1.0
 * @date 4/11/2018
 * @section LICENSE
 *
 * @description a class for nicely dealing with time varying categorical layers.
 */

#ifndef CATEGORICAL_META_LAYER_H_
#define CATEGORICAL_META_LAYER_H_

// headers
#include "../CategoricalLayer.h"

namespace niwa {
namespace layers {

//**********************************************************************
//
//
//**********************************************************************
class CategoricalMetaLayer : public CategoricalLayer {
public:
  // Functions
  CategoricalMetaLayer(Model* model);
  virtual                     ~CategoricalMetaLayer() {};
  virtual void                set_value(unsigned RowIndex, unsigned ColIndex, string Value) {};  // TODO once we have children make this pure virtual
  virtual string              get_value(unsigned RowIndex, unsigned ColIndex);
  virtual string              get_value(unsigned RowIndex, unsigned ColIndex, unsigned year);

protected:
  virtual void                DoValidate() override final;
  virtual void                DoBuild() override final;
  // Variables
  string                       default_Layer_label_;
  vector<unsigned>             years_;
  vector<string>               layer_names_;
  map<unsigned, CategoricalLayer*> years_layer_;
  CategoricalLayer*                default_layer_;
  bool                         has_years_;
};
} /* namespace layers */
} /* namespace niwa */
#endif /* CATEGORICAL_META_LAYER_H_ */
