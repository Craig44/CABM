/**
 * @file Mortality.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 */

// headers
#include "Mortality.h"

#include "Layers/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
Mortality::Mortality(Model* model) : Process(model) {
  process_type_ = ProcessType::kMortality;

}

void Mortality::DoReset() {
  LOG_FINE() << "clearing containers";
  removals_by_age_and_area_.clear();
  removals_by_length_and_area_.clear();
  removals_census_.clear();
  removals_tag_recapture_.clear();
}


// true means all years are found, false means there is a mismatch in years
bool Mortality::check_years(vector<unsigned>& years_to_check_) {
  LOG_FINE();
  for (unsigned year_ndx = 0; year_ndx < years_to_check_.size(); ++year_ndx) {
    if (find(years_.begin(), years_.end(), years_to_check_[year_ndx]) == years_.end())
      return false;
  }

  return true;

}
// Does this fishery exist?
// yes = true
// no = false
bool Mortality::check_fishery_exists(string fishery_label) {
  LOG_FINE();
  if (find(fishery_label_.begin(), fishery_label_.end(), fishery_label) != fishery_label_.end())
    return true;

  return false;

}

// find and return age comp data for this fishery.
// TODO: change this, make void and pass a pointer of class vector<vector<census_data>> and fishery label,
// and assign the address of the pointer to this fisher label, that way no copies are being made.
vector<vector<census_data>> Mortality::get_fishery_census_data(string fishery_label) {
  LOG_FINE();
  vector<vector<census_data>> fish_comp;
  for (unsigned i = 0; i < fishery_label_.size(); ++i) {
    if (fishery_label_[i] == fishery_label)
      return fishery_census_data_[fishery_index_[i]];
  }
  return fish_comp;
}

// Return Fishery specific Length composition
vector<vector<composition_data>>* Mortality::get_fishery_length_comp(string &fishery_label) {
  LOG_FINE();
  vector<vector<composition_data>> *fish_comp = nullptr;
  for (unsigned i = 0; i < fishery_label_.size(); ++i) {
    if (fishery_label_[i] == fishery_label)
      return &length_comp_by_fishery_[fishery_index_[i]];
  }
  LOG_FATAL()
  << "This should not reach this point";
  return fish_comp;
}

// Return Fishery specific Age composition
vector<vector<composition_data>>* Mortality::get_fishery_age_comp(string &fishery_label) {
  LOG_FINE();
  vector<vector<composition_data>> *fish_comp = nullptr;
  for (unsigned i = 0; i < fishery_label_.size(); ++i) {
    if (fishery_label_[i] == fishery_label)
      return &age_comp_by_fishery_[fishery_index_[i]];
  }
  LOG_FATAL()
  << "This should not reach this point";
  return fish_comp;
}

} /* namespace processes */
} /* namespace niwa */
