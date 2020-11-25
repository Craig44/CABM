/**
 * @file ModelAttributes.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 18/07/2018
 * @section LICENSE
 *
 *
 */

// headers
#include "ModelAttributes.h"

#include <boost/algorithm/string/join.hpp>

#include "Model/Managers.h"
#include "Model/Model.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"

// namespaces
namespace niwa {
namespace reports {

/**
 * Default constructor
 *
 * @param model Pointer to the current model context
 */
ModelAttributes::ModelAttributes(Model* model) : Report(model) {
  model_state_ = State::kIterationComplete;;
  run_mode_    = (RunMode::Type)(RunMode::kBasic);
}

/**
 * Build our relationships between this object and other objects
 */
void ModelAttributes::DoBuild() {
  LOG_TRACE();
}

/**
 * Execute this report
 */
void ModelAttributes::DoExecute() {
  LOG_MEDIUM() <<" printing report " << label_;

  auto scalars = model_->get_scalars();
  cache_ << "*"<< type_ << "[" << label_ << "]" << "\n";
  for (auto scalar : scalars)
    cache_ << scalar.first << ": " << scalar.second << "\n";

  cache_ << "model_years: ";
  for (auto year : model_->years())
    cache_ << year << " ";
  cache_ << "\n";

  cache_ << "ages: ";
  for (auto age = model_->min_age(); age <= model_->max_age(); ++age)
    cache_ << age << " ";
  cache_ << "\n";

  cache_ << "length_mid_points: ";
  vector<float> length_mid_point = model_->length_bin_mid_points();
  for (unsigned length_ndx = 0; length_ndx < length_mid_point.size(); ++length_ndx)
    cache_ << length_mid_point[length_ndx] << " ";
  cache_ << "\n";
  cache_ << "length_bins: ";
  vector<unsigned> length_bins = model_->length_bins();
  for (unsigned length_ndx = 0; length_ndx < length_bins.size(); ++length_ndx)
    cache_ << length_bins[length_ndx] << " ";
  cache_ << "\n";

  if(model_->get_lat_mid_points().size() > 0) {
    cache_ << "latitude_mid_points: ";
    for (unsigned ndx = 0; ndx < model_->get_lat_mid_points().size(); ++ndx)
      cache_ << model_->get_lat_mid_points()[ndx] << " ";
    cache_ << "\n";

    cache_ << "latitude_matrix "  << REPORT_R_MATRIX << "\n";
    for (unsigned lat_ndx = 0; lat_ndx < model_->get_lat_mid_points().size(); ++lat_ndx) {
      for (unsigned lon_ndx = 0; lon_ndx < model_->get_lon_mid_points().size(); ++lon_ndx) {
        cache_ <<  model_->get_lat_mid_points()[lat_ndx] << " ";
      }
      cache_ << "\n";
    }
  }

  if(model_->get_lon_mid_points().size() > 0) {
    cache_ << "longitude_mid_points: ";
    for (unsigned ndx = 0; ndx < model_->get_lon_mid_points().size(); ++ndx)
      cache_ << model_->get_lon_mid_points()[ndx] << " ";
    cache_ << "\n";
    cache_ << "longitude_matrix "  << REPORT_R_MATRIX << "\n";
    for (unsigned lon_ndx = 0; lon_ndx < model_->get_lon_mid_points().size(); ++lon_ndx) {
      for (unsigned lat_ndx = 0; lat_ndx < model_->get_lat_mid_points().size(); ++lat_ndx) {
        cache_ <<  model_->get_lon_mid_points()[lon_ndx] << " ";
      }
      cache_ << "\n";
    }
  }


  cache_ << "growth_model: ";
  switch(model_->get_growth_model()) {

  case Growth::Type::kVonbert:
    cache_ << "VonBertalanffy";
    break;
  case Growth::Type::kSchnute:
    cache_ << "Schnute";
    break;
  case Growth::Type::kInvalid:
    cache_ << "Invalid";
    break;
  }

  cache_ << "\n";
  ready_for_writing_ = true;

}

} /* namespace reports */
} /* namespace niwa */
