/**
 * @file Mortality.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process is the parent mortlity class, so the manager can generate dynamicasts<> processes
 * any child mortality class should inherit this as a parent
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_H_

// headers
#include "Processes/Process.h"

// namespaces
namespace niwa {
namespace processes {


// Aggregated compositional data
struct composition_data {
  string type_;
  unsigned year_;
  unsigned row_;
  unsigned col_;
  vector<float> frequency_;
  vector<float> female_frequency_;
  composition_data(string type, unsigned year, unsigned row, unsigned col, unsigned size) : type_(type), year_(year),
      row_(row), col_(col)
  {
    frequency_.resize(size);
    female_frequency_.resize(size);
  }
};

// Census data that links every age and length to one another of every agent caught, may only want this for a few scenerios could be memory\computational hungry
struct census_data {
  unsigned year_;
  unsigned row_;
  unsigned col_;
  float    biomass_; // Total biomass if unsexed otherwise male biomass in sexed model
  float    female_biomass_; // Total female biomass
  vector<unsigned> fishery_ndx_;
  vector<unsigned> sex_;
  vector<unsigned> age_ndx_;
  vector<unsigned> length_ndx_;
  vector<float> scalar_;
  vector<float> weight_;
  census_data(unsigned year, unsigned row, unsigned col) : year_(year),
      row_(row), col_(col) { biomass_ = 0.0, female_biomass_ = 0.0; }
};

// Tag Recapture class
struct tag_recapture {
  unsigned year_;
  unsigned row_;
  unsigned col_;
  unsigned time_step_;
  unsigned scanned_fish_ = 0;  // individuals
  unsigned scanned_agents_ = 0;// agents
  double   expected_scanned_ = 0;
  unsigned tag_draws_ = 0;
  vector<float> tagged_fish_available_;
  vector<float> all_fish_available_;
  unsigned agents_sampled_ = 0;
  unsigned agents_caught_ = 0;
  double individuals_caught_ = 0.0;
  double proportion_inidividuals_tagged_ = 0;
  double prob_sample_tagged_agents_ = 0;
  vector<unsigned> age_;
  vector<unsigned> sex_;
  vector<unsigned> fishery_ndx_;
  vector<unsigned> length_ndx_;
  vector<unsigned> time_at_liberty_;
  vector<unsigned> tag_row_;
  vector<unsigned> tag_col_;
  vector<unsigned> tag_release_year_;
  vector<float> length_increment_;
  tag_recapture(unsigned year, unsigned row, unsigned col, unsigned time_step) : year_(year),
      row_(row), col_(col), time_step_(time_step) { }
};
/**
 * Class definition
 */
class Mortality : public Process {
public:
  // methods
  explicit Mortality(Model* model);
  virtual                     ~Mortality() = default;
  virtual void                        DoValidate(){ };
  virtual void                        DoBuild() { };
  virtual void                        DoReset();
  virtual void                        DoExecute() { };

  virtual void                        draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) = 0;
  vector<composition_data>&           get_removals_by_age() {return removals_by_age_and_area_;};
  vector<composition_data>&           get_removals_by_length() {return removals_by_length_and_area_;};
  vector<census_data>&                get_census_data() {return removals_census_;};
  vector<tag_recapture>&              get_tag_recapture_info() {return removals_tag_recapture_;};


  virtual bool                        update_mortality() {return update_natural_mortality_parameters_;};
  virtual double                      SolveBaranov() { return 1.0;};
  void                                set_lambda(double lambda) {lambda_ = lambda;};
  virtual bool                        check_years(vector<unsigned> years_to_check_) {return false;};

protected:
  vector<composition_data>            removals_by_age_and_area_;
  vector<composition_data>            removals_by_length_and_area_;
  vector<census_data>                 removals_census_;
  vector<tag_recapture>               removals_tag_recapture_;
  map<unsigned, vector<unsigned>>     removals_by_age_;
  map<unsigned, vector<unsigned>>     removals_by_length_;
  bool                                update_natural_mortality_parameters_;
  double                              lambda_;



};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_H_ */
