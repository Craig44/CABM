/**
 * @file AgeFrequencyByCell.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 23/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This report will take print the state of the partition at two points in the initialisation phase.
 * After the approximation, and after the entire initialisation phase.
 */
#ifndef SOURCE_REPORTS_CHILDREN_AGE_FEQUENCY_H_
#define SOURCE_REPORTS_CHILDREN_AGE_FEQUENCY_H_

// headers
#include "Reports/Report.h"

// namespaces
namespace niwa {
class WorldView;

namespace reports {

// classes
class AgeFrequencyByCell : public Report {
public:
  AgeFrequencyByCell(Model* model);
  virtual                     ~AgeFrequencyByCell() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final;
  void                        DoExecute() override final;


private:
  bool                        first_run_ = true;
  WorldView*                  world_;
  bool                        call_number_ = true; // initialisation phase reports are called twice once in the middle of initialisation and one once it is complete.
};

} /* namespace reports */
} /* namespace niwa */

#endif /* SOURCE_REPORTS_CHILDREN_AGE_FEQUENCY_H_ */
