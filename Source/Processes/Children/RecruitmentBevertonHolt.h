/**
 * @file RecruitmentBevertonHolt.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 12/07/2018
 * @section LICENSE
 *
 * Copyright NIWA Science �2015 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This process is responsible for seeding new recruits into the world
 */
#ifndef SOURCE_PROCESSES_CHILDREN_RECRUITMENT_BEVERTON_HOLT_H_
#define SOURCE_PROCESSES_CHILDREN_RECRUITMENT_BEVERTON_HOLT_H_

// headers
#include "Processes/Process.h"
#include "Layers/Children/Numeric/Base/NumericLayer.h"

// namespaces
namespace niwa {
namespace processes {
using std::string;
/**
 * Class definition
 */
class RecruitmentBevertonHolt : public Process {
public:
  // methods
  explicit RecruitmentBevertonHolt(Model* model);
  virtual                     ~RecruitmentBevertonHolt() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final { };
protected:
  string                      ssb_label_;
  string                      recruitment_layer_label_;
  vector<double>              ycs_values_;
  double                      b0_;
  double                      steepness_;
  layers::NumericLayer*       recruitment_layer_ = nullptr;
  double                      initial_scalar_;


};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_RECRUITMENT_BEVERTON_HOLT_H_ */
