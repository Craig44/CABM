/**
 * @file NumericLayer.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 31/07/208
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This report will take a target numeric layer and print all of
 * the values
 */
#ifndef SOURCE_REPORTS_CHILDREN_NUMERIC_LAYER_H_
#define SOURCE_REPORTS_CHILDREN_NUMERIC_LAYER_H_

// headers
#include "Reports/Report.h"

#include "Layers/Children/NumericLayer.h"

// namespaces
namespace niwa {
namespace reports {

// classes
class NumericLayer : public Report {
public:
  NumericLayer(Model* model);
  virtual                     ~NumericLayer() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final;
  void                        DoExecute() override final;
  void                        DoReset() override final  { };


private:
  string                      layer_label_ = "";
  niwa::layers::NumericLayer* layer_ = nullptr;
  bool                        first_run_ = true;

};

} /* namespace reports */
} /* namespace niwa */

#endif /* SOURCE_REPORTS_CHILDREN_NUMERIC_LAYER_H_ */
