/**
 * @file ModelAttributes.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 18/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This report will take print model attributes
 */
#ifndef SOURCE_REPORTS_CHILDREN_MODEL_ATTRIBUTES_H_
#define SOURCE_REPORTS_CHILDREN_MODEL_ATTRIBUTES_H_

// headers
#include "Reports/Report.h"

// namespaces
namespace niwa {
class WorldView;

namespace reports {

// classes
class ModelAttributes : public Report {
public:
  ModelAttributes(Model* model);
  virtual                     ~ModelAttributes() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final;
  void                        DoExecute() override final;


private:
  bool                        first_run_ = true;
  WorldView*                  world_;
};

} /* namespace reports */
} /* namespace niwa */

#endif /* SOURCE_REPORTS_CHILDREN_MODEL_ATTRIBUTES_H_ */
