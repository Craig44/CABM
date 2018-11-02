/**
 * @file MortalityBaranov.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 21/09/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This is a child mortality that applies a fishing event
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_BARANOV_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_BARANOV_H_

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
class MortalityBaranov : public Mortality {
public:
  // methods
  explicit MortalityBaranov(Model* model);
  virtual                     ~MortalityBaranov() = default;
  virtual void                        DoValidate() override final;
  virtual void                        DoBuild() override final;
  virtual void                        DoReset() override final { };
  virtual void                        DoExecute() override final;
  void                                draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) override final;
  void                                FillReportCache(ostringstream& cache) override final;
  void                                RebuildCache() override final;
  bool                                check_years(vector<unsigned> years_to_check_);
protected:
  // F variables
  vector<string>                      f_layer_label_;
  vector<layers::NumericLayer*>       f_layer_;
  vector<string>                      fishery_selectivity_label_;
  vector<Selectivity*>                fishery_selectivity_;
  vector<unsigned>                    years_;
  float                               discard_mortality_;
  float                               mls_;
  bool                                f_selectivity_length_based_ = false;
  map<unsigned,float>                 catch_by_year_;
  bool                                print_extra_info_ = false;
  // For reporting
   map<unsigned, float>                actual_removals_by_year_;
   map<unsigned, float>                removals_by_year_;
   vector<vector<vector<float>>>       F_by_cell_year_;
  // For Tag-recaptures
  vector<unsigned>                    scanning_years_;
  vector<float>                       scanning_proportion_;
  vector<vector<float>>               scanning_prop_year_by_space_;
  vector<float>                       scanning_random_numbers_;

  // M variables
  vector<string>                      natural_mortality_selectivity_label_;
  string                              m_layer_label_;
  layers::NumericLayer*               m_layer_ = nullptr;
  vector<Selectivity*>                natural_mortality_selectivity_;
  float                               cv_;
  string                              distribution_;
  bool                                m_selectivity_length_based_ = false;
  float                               m_;

  // Threading variables
  vector<vector<unsigned>>            cell_offset_;
  vector<vector<unsigned>>            model_length_bins_;
  vector<vector<unsigned>>            model_age_bins_;
  vector<vector<float>>               mls_by_space_;
  vector<vector<float>>               discard_by_space_;
  vector<vector<unsigned>>            current_year_by_space_;
  vector<vector<unsigned>>            current_time_step_by_space_;
  unsigned                            n_agents_;



};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_BARANOV_H_ */
