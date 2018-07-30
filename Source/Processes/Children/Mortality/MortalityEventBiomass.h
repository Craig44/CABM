/**
 * @file MortalityEventBiomass.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 18/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This is a child mortality that applies a fishing event
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_EVENT_BIOMASS_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_EVENT_BIOMASS_H_

// headers
#include "Processes/Children/Mortality.h"

#include "Layers/Children/NumericLayer.h"
#include <omp.h>


// namespaces
namespace niwa {
class Selectivity;
namespace processes {
using std::string;
/**
 * Class definition
 */
class MortalityEventBiomass : public Mortality {
public:
  // methods
  explicit MortalityEventBiomass(Model* model);
  virtual                     ~MortalityEventBiomass() = default;
  virtual void                        DoValidate() override final;
  virtual void                        DoBuild() override final;
  virtual void                        DoReset() override final { };
  virtual void                        DoExecute() override final;
  void                                draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) override final { };
  void                                FillReportCache(ostringstream& cache) override final;

protected:
  vector<string>                      catch_layer_label_;
  vector<layers::NumericLayer*>       catch_layer_;
  vector<string>                      selectivity_label_;
  vector<Selectivity*>                selectivity_;
  vector<unsigned>                    years_;
  bool                                selectivity_length_based_ = false;
  float                               discard_mortality_;
  float                               mls_;
  // objects for thread safety of rng
  vector<float>                       random_numbers_;
  vector<float>                       discard_random_numbers_;
  vector<float>                       selectivity_random_numbers_;
  vector<vector<vector<float>>>       cell_offset_for_selectivity_;

  unsigned                            n_agents_;
  vector<vector<unsigned>>            cell_offset_;

  // For reporting
  map<unsigned, float>                actual_removals_by_year_;
  map<unsigned, float>                removals_by_year_;

};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_EVENT_BIOMASS_H_ */
