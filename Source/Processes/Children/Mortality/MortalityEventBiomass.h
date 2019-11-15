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
  virtual bool                        check_years(vector<unsigned> years_to_check_) override final;
  bool                                check_fishery_exists(string fishery_label);
  vector<unsigned>&                   get_fishery_years() {return catch_year_;};
  vector<vector<census_data>>         get_fishery_census_data(string fishery_label);


protected:
  parameters::Table*                  catch_table_ = nullptr;
  parameters::Table*                  method_table_ = nullptr;
  vector<unsigned>                    fishery_index_; // used for look up on all the vectors specific fishery objects, better than maps for random access
  vector<string>                      fishery_label_;
  vector<unsigned>                    catch_year_;
  vector<vector<string>>              fishery_catch_layer_labels_; // n_fishery * n_years
  vector<vector<float>>               fishery_actual_catch_taken_; // n_fishery * n_years
  vector<vector<float>>               fishery_catch_to_take_; // n_fishery * n_years
  vector<vector<layers::NumericLayer*>>  fishery_catch_layer_; // n_fishery * n_years
  vector<vector<string>>              fishery_selectivity_label_; // n_fishery * n_sex
  vector<vector<Selectivity*>>        fishery_selectivity_; // n_fishery * n_sex
  vector<float>                       fishery_mls_; // fishery specific
  vector<float>                       fishery_hand_mort_; // fishery specific
  bool                                selectivity_length_based_ = false;
  vector<float>                       catch_to_take_by_fishery_;
  // For scanning and Tag-recaptures optional
  parameters::Table*                  scanning_table_ = nullptr;
  vector<unsigned>                    scanning_years_;
  vector<vector<float>>               scanning_proportion_by_fishery_; // n_fishery * n_years;
  bool                                scanning_ = false;
  bool                                account_for_tag_population_ = false;
  vector<bool>                        scanning_this_year_; // n_fishery

  // For reporting
  map<unsigned, float>                actual_removals_by_year_;
  map<unsigned, float>                removals_by_year_;
  bool                                print_extra_info_ = false;
  vector<vector<vector<composition_data>>>    age_comp_by_fishery_; // n_fishery * n_years * n_cells
  vector<vector<vector<composition_data>>>   length_comp_by_fishery_; // n_fishery * n_years * n_cells
  vector<vector<vector<census_data>>>   fishery_census_data_; // n_fishery * n_years * n_cells
  vector<vector<vector<vector<float>>>> actual_catch_by_area_; // n_years * n_fishery * rows * cols

  vector<unsigned>                    cell_ndx_;


};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_EVENT_BIOMASS_H_ */
