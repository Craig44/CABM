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
  bool                                selectivity_length_based_ = false;
  bool                                apply_using_proportions_ = false;
  // objects for thread safety of rng
  vector<double>                       random_numbers_;
  vector<double>                       selectivity_random_numbers_;
  vector<double>                       handling_mortality_random_numbers_;
  unsigned                            n_agents_;

  vector<unsigned>                    age_freq_;
  vector<unsigned>                    length_freq_;

  vector<string>                      tag_layer_label_;
  vector<layers::IntLayer*>           tag_layer_;
  vector<unsigned>                    scanning_years_;
  vector<double>                       scanning_proportion_;
  double                               handling_mortality_;
  // Reporting
  map<unsigned,vector<unsigned>>      age_distribution_of_tagged_fish_by_year_;
  map<unsigned,vector<unsigned>>      length_distribution_of_tagged_fish_by_year_;

  vector<vector<vector<vector<double>>>>       age_length_param1_of_tagged_fish_by_year_cell_; // year * row * col * agents
  vector<vector<vector<vector<double>>>>       age_length_param2_of_tagged_fish_by_year_cell_; // year * row * col * agents
  vector<vector<vector<vector<vector<double>>>>>  age_length_key_by_release_event; // year * row * col * age * length

  vector<vector<vector<vector<unsigned>>>>    length_observed_tag_of_tagged_fish_by_year_cell_; // year * row * col * length_bins
  vector<vector<vector<vector<unsigned>>>>    length_distribution_of_tagged_fish_by_year_cell_; // year * row * col * length_bins
  vector<vector<vector<vector<unsigned>>>>    age_distribution_of_tagged_fish_by_year_cell_; // year * row * col * age_bins

  parameters::Table*                  proportions_table_ = nullptr;
  vector<vector<double>>               proportions_data_; // n_rows x n_length_bins
  vector<unsigned>                    table_rows_;
  vector<unsigned>                    table_cols_;
  vector<unsigned>                    table_years_;


};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_TAGGING_H_ */
