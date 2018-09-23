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
#ifndef SOURCE_PROCESSES_CHILDREN_RECRUITMENT_H_
#define SOURCE_PROCESSES_CHILDREN_RECRUITMENT_H_

// headers
#include "Processes/Process.h"
#include "Layers/Children/NumericLayer.h"
#include "DerivedQuantities/DerivedQuantity.h"
// namespaces
namespace niwa {
namespace processes {
using std::string;
/**
 * Class definition
 */
class Recruitment : public Process {
public:
  // methods
  explicit Recruitment(Model* model);
  virtual                     ~Recruitment() = default;
  virtual void                DoValidate() { };
  virtual void                DoBuild();
  virtual void                DoReset()  { };
  virtual void                DoExecute() { };
  virtual void                FillReportCache(ostringstream& cache) { };

protected:
  string                      ssb_label_;
  string                      recruitment_layer_label_;
  float                       b0_;
  layers::NumericLayer*       recruitment_layer_ = nullptr;
  DerivedQuantity*            derived_quantity_ = nullptr;
  float                       initial_scalar_;

  // Reporting containers that will be printed in FillReportCache() method
  map<unsigned, float>        recruits_by_year_;
  unsigned                    initial_recruits_;
  bool                        first_enter_execute_ = true;
  map<unsigned, float>        ssb_by_year_;
  float                       scalar_;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_RECRUITMENT_H_ */
