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
 * This report will take print the a bunch of tagging partition information
 */
#ifndef SOURCE_REPORTS_CHILDREN_TAGGING_INFO_H_
#define SOURCE_REPORTS_CHILDREN_TAGGING_INFO_H_

// headers
#include "Reports/Report.h"

// namespaces
namespace niwa {
class WorldView;

namespace reports {

// classes
class TaggingInfo : public Report {
public:
  TaggingInfo(Model* model);
  virtual                     ~TaggingInfo() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final;
  void                        DoExecute() override final;
  void                        DoReset() override final  { };


private:
  bool                        first_run_ = true;
  WorldView*                  world_;
  bool                        call_number_ = true; // initialisation phase reports are called twice once in the middle of initialisation and one once it is complete.
  vector<unsigned>            release_year_;
  vector<string>              release_region_;
  vector<unsigned>            origin_rows_;
  vector<unsigned>            origin_cols_;
  vector<float>               age_freq_;
};

} /* namespace reports */
} /* namespace niwa */

#endif /* SOURCE_REPORTS_CHILDREN_TAGGING_INFO_H_ */
