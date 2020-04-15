/**
 * @file InitialisationPartition.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 15/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This report will take print the state of the partition at two points in the initialisation phase.
 * After the approximation, and after the entire initialisation phase.
 */
#ifndef SOURCE_REPORTS_CHILDREN_INITIALISATION_PARTITION_H_
#define SOURCE_REPORTS_CHILDREN_INITIALISATION_PARTITION_H_

// headers
#include "Reports/Report.h"

// namespaces
namespace niwa {
class WorldView;

namespace reports {

// classes
class InitialisationPartition : public Report {
public:
  InitialisationPartition(Model* model);
  virtual                     ~InitialisationPartition() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final;
  void                        DoReset() override final {call_number_ = true;};
  void                        DoExecute() override final;


private:
  bool                        do_length_frequency_;
  bool                        do_age_ = true;
  bool                        first_run_ = true;
  WorldView*                  world_;
  bool                        call_number_ = true; // initialisation phase reports are called twice once in the middle of initialisation DoExecute and one once it is complete in RunBasic Model.cpp
};

} /* namespace reports */
} /* namespace niwa */

#endif /* SOURCE_REPORTS_CHILDREN_INITIALISATION_PARTITION_H_ */
