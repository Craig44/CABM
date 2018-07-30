/**
 * @file Maturity.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 15/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process does nothing. It's used for debugging time steps
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MATURITY_H_
#define SOURCE_PROCESSES_CHILDREN_MATURITY_H_

// headers
#include "Processes/Process.h"
#include <omp.h>

// namespaces
namespace niwa {
class Selectivity;
namespace processes {

/**
 * Class definition
 */
class Maturity : public Process {
public:
  // methods
  explicit Maturity(Model* model);
  virtual                     ~Maturity() = default;
  void                        DoValidate() override final { };
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final;
  void                        FillReportCache(ostringstream& cache) override final;

protected:
  vector<string>             selectivity_label_;
  vector<Selectivity*>       selectivity_;
  map<unsigned, unsigned>    mature_individuals_by_year_;

  // objects for thread safety of rng
  vector<float>                       random_numbers_;
  unsigned                            n_agents_;
  vector<vector<unsigned>>            cell_offset_;

  vector<vector<vector<float>>>       cell_offset_for_selectivity_;
  bool                                length_based_selectivity_ = false;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MATURITY_H_ */
