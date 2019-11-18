/**
 * @file Tagging.h
 * @author C>Marsh
 * @github https://github.com/Craig44
 * @date 26/7/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process does nothing. It's used for debugging time steps
 */
#ifndef SOURCE_PROCESSES_CHILDREN_TAGGING_H_
#define SOURCE_PROCESSES_CHILDREN_TAGGING_H_

// headers
#include "Processes/Process.h"

#include "Layers/Children/IntLayer.h"

// namespaces
namespace niwa {
class Selectivity;
namespace processes {

/**
 * Class definition
 */
class Tagging : public Process {
public:
  // methods
  explicit                            Tagging(Model* model);
  virtual                             ~Tagging();
  void                                DoValidate() override final;
  void                                DoBuild() override final;
  void                                DoReset() override final;
  void                                DoExecute() override final;
  void                                FillReportCache(ostringstream& cache) override final;

protected:
  vector<unsigned>                    years_;
  vector<string>                      selectivity_labels_;
  vector<Selectivity*>                selectivities_;
  bool                                selectivity_length_based_;
  bool                                apply_using_proportions_ = false;
  // objects for thread safety of rng
  vector<float>                       random_numbers_;
  vector<float>                       selectivity_random_numbers_;
  vector<float>                       handling_mortality_random_numbers_;
  vector<vector<vector<float>>>       cell_offset_for_selectivity_;

  unsigned                            n_agents_;
  vector<vector<unsigned>>            cell_offset_;
  vector<vector<unsigned>>            model_length_bins_;
  vector<vector<unsigned>>            model_age_bins_;
  vector<vector<unsigned>>            current_year_by_space_;
  vector<vector<unsigned>>            handling_mort_by_space_;

  vector<string>                      tag_layer_label_;
  vector<layers::IntLayer*>           tag_layer_;
  vector<unsigned>                    scanning_years_;
  vector<float>                       scanning_proportion_;
  float                               handling_mortality_;
  // Reporting
  map<unsigned,vector<unsigned>>      age_distribution_of_tagged_fish_by_year_;
  map<unsigned,vector<unsigned>>      length_distribution_of_tagged_fish_by_year_;

  vector<vector<vector<vector<unsigned>>>>    length_distribution_of_tagged_fish_by_year_cell_; // year * row * col * length_bins
  vector<vector<vector<vector<unsigned>>>>    age_distribution_of_tagged_fish_by_year_cell_; // year * row * col * age_bins

  parameters::Table*                  proportions_table_ = nullptr;
  vector<vector<float>>               proportions_data_; // n_rows x n_length_bins
  vector<unsigned>                    table_rows_;
  vector<unsigned>                    table_cols_;
  vector<unsigned>                    table_years_;


};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_TAGGING_H_ */
