/**
 * @file RecruitmentConstant.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 12/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process is responsible for seeding new recruits into the world
 * TODO you could make this a random process, but perhaps it doesn't matter, that is randomly seed an individual in one of the areas
 */
#ifndef SOURCE_PROCESSES_CHILDREN_RECRUITMENT_CONSTANT_H_
#define SOURCE_PROCESSES_CHILDREN_RECRUITMENT_CONSTANT_H_

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
class RecruitmentConstant : public Recruitment {
public:
  // methods
  explicit RecruitmentConstant(Model* model);
  virtual                     ~RecruitmentConstant() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final;
  void                        FillReportCache(ostringstream& cache) override final;

protected:


};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_RECRUITMENT_CONSTANT_H_ */
