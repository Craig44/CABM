/**
 * @file RecruitmentBevertonHolt.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 12/07/2018
 * @section LICENSE
 *
 * Copyright NIWA Science ©2015 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This process is responsible for seeding new recruits into the world
 * TODO you could make this a random process, but perhaps it doesn't matter, that is randomly seed an individual in one of the areas
 *
 */
#ifndef SOURCE_PROCESSES_CHILDREN_RECRUITMENT_BEVERTON_HOLT_H_
#define SOURCE_PROCESSES_CHILDREN_RECRUITMENT_BEVERTON_HOLT_H_

// headers
#include "Processes/Children/Recruitment.h"
#include "Layers/Children/NumericLayer.h"
#include "DerivedQuantities/DerivedQuantity.h"
// namespaces
namespace niwa {
namespace processes {
using std::string;
/**
 * Class definition
 */
class RecruitmentBevertonHolt : public Recruitment {
public:
  // methods
  explicit RecruitmentBevertonHolt(Model* model);
  virtual                     ~RecruitmentBevertonHolt() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final;
  void                        DoExecute() override final;
  void                        FillReportCache(ostringstream& cache) override final;

protected:
  vector<float>               ycs_values_;
  float                       steepness_;

  // Reporting containers that will be printed in FillReportCache() method
  map<unsigned, float>        ycs_values_by_year_;
  map<unsigned, float>        ssb_ratio_;
  map<unsigned, float>        SR_;
  map<unsigned, float>        true_ycs_;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_RECRUITMENT_BEVERTON_HOLT_H_ */
