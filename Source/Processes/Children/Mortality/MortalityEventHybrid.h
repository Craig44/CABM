/**
 * @file MortalityEventHybrid.h
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
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_EVENT_HYBRID_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_EVENT_HYBRID_H_

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
class MortalityEventHybrid : public Mortality {
public:
  // methods
  explicit MortalityEventHybrid(Model* model);
  virtual                        ~MortalityEventHybrid();
  virtual void                        DoValidate() override final;
  virtual void                        DoBuild() override final;
  virtual void                        DoReset() override final;
  virtual void                        DoExecute() override final;
  void                                draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) override final { };
  void                                FillReportCache(ostringstream& cache) override final;
  virtual void                        set_HCR(map<unsigned, map<string, float>> future_catches)  override final; // years * fishery label * catch


protected:
  parameters::Table*                  f_table_ = nullptr;
  parameters::Table*                  method_table_ = nullptr;
  vector<vector<string>>              fishery_f_layer_labels_; // n_fishery * n_years
  vector<vector<float>>               fishery_actual_catch_taken_; // n_fishery * n_years
  vector<vector<float>>               fishery_f_to_take_; // n_fishery * n_years
  vector<vector<vector<vector<vector<float>>>>>  F_by_year_bin_; // n_years * n_col *  * n_sex * n_bins (either age or length)
  vector<vector<vector<vector<vector<float>>>>>  prop_F_fishery_and_bin_; //n_fishery * n_sex * n_bins (either age or length)
  vector<vector<layers::NumericLayer*>>  fishery_f_layer_; // n_fishery * n_years
  vector<vector<string>>              fishery_selectivity_label_; // n_fishery * n_sex
  vector<vector<Selectivity*>>        fishery_selectivity_; // n_fishery * n_sex
  vector<float>                       fishery_mls_; // fishery specific
  vector<float>                       fishery_hand_mort_; // fishery specific
  bool                                selectivity_length_based_ = false;
  vector<float>                       f_to_take_by_fishery_;
  // For reporting
  map<unsigned, float>                actual_removals_by_year_;
  map<unsigned, float>                removals_by_year_;
  bool                                print_census_info_ = false;
  vector<unsigned>                    harvest_control_years_;
  map<unsigned, map<string, float>>   harvest_control_Catches_;
  bool                                use_HCR_vals_ = false;
  // need to be reset for any multi run format in DoReset
  vector<vector<vector<vector<float>>>> actual_catch_by_area_; // n_years * n_fishery * rows * cols

  vector<unsigned>                    cell_ndx_;
  unsigned                            n_bins_;

  // event biomass stuff
  vector<float>                       catch_to_take_by_fishery_;

};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_EVENT_BIOMASS_H_ */
