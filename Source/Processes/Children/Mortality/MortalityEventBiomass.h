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
  virtual                        ~MortalityEventBiomass();
  virtual void                        DoValidate() override final;
  virtual void                        DoBuild() override final;
  virtual void                        DoReset() override final;
  virtual void                        DoExecute() override final;
  void                                draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) override final { };
  void                                FillReportCache(ostringstream& cache) override final;
  virtual void                        set_HCR(map<unsigned, map<string, double>> future_catches)  override final; // years * fishery label * catch




protected:
  parameters::Table*                  catch_table_ = nullptr;
  parameters::Table*                  method_table_ = nullptr;
  vector<vector<string>>              fishery_catch_layer_labels_; // n_fishery * n_years
  vector<vector<double>>               fishery_actual_catch_taken_; // n_fishery * n_years
  vector<vector<double>>               fishery_catch_to_take_; // n_fishery * n_years
  vector<vector<layers::NumericLayer*>>  fishery_catch_layer_; // n_fishery * n_years
  vector<vector<string>>               fishery_selectivity_label_; // n_fishery * n_sex
  vector<vector<Selectivity*>>         fishery_selectivity_; // n_fishery * n_sex
  vector<double>                       fishery_mls_; // fishery specific
  vector<double>                       fishery_hand_mort_; // fishery specific
  bool                                 selectivity_length_based_ = false;
  vector<double>                       catch_to_take_by_fishery_;
  vector<double>                       vulnerable_biomass_;

  // For scanning and Tag-recaptures optional
  parameters::Table*                  scanning_table_ = nullptr;
  vector<unsigned>                    scanning_years_;
  vector<vector<double>>               scanning_proportion_by_fishery_; // n_fishery * n_years;
  bool                                scanning_ = false;
  bool                                account_for_tag_population_ = false;
  vector<bool>                        scanning_this_year_; // n_fishery
  vector<bool>                        only_mature_partition_;
  // For reporting
  map<unsigned, double>                actual_removals_by_year_;
  map<unsigned, double>                removals_by_year_;
  bool                                print_census_info_ = false;
  bool                                print_tag_recap_info_ = false;

  vector<unsigned>                    harvest_control_years_;
  map<unsigned, map<string, double>>   harvest_control_Fs_;
  bool                                use_HCR_vals_ = false;

  // need to be reset for any multi run format in DoReset
  vector<vector<vector<vector<double>>>> actual_catch_by_area_; // n_years * n_fishery * rows * cols


  vector<unsigned>                    cell_ndx_;


};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_EVENT_BIOMASS_H_ */
