/**
 * @file CategoricalMetaLayer.cpp
 * @author  C Marsh
 * @version 1.0
 * @date 4/11/2018
 * @section LICENSE
 *
 */

// Local headers
#include "CategoricalMetaLayer.h"

#include "Layers/Manager.h"

namespace niwa {
namespace layers {

// Constructor
CategoricalMetaLayer::CategoricalMetaLayer(Model* model) : CategoricalLayer(model) {
  // Default Variables
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "the years to apply the following layer", "");
  parameters_.Bind<string>(PARAM_DEFAULT_LAYER, &default_Layer_label_, "The default layer to apply in initialisation phase", "");
  parameters_.Bind<string>(PARAM_LAYER_LABELS, &layer_names_, "The layer label associated with each year", "");
  layer_type_ = LayerType::kCategorical;

}

//**********************************************************************
// void CNumericMetaLayer::validate()
// Validate the layer
//**********************************************************************
void CategoricalMetaLayer::DoValidate() {
  LOG_TRACE();
  if (years_.size() != layer_names_.size())
    LOG_ERROR_P(PARAM_LAYER_LABELS) << "there needs to be a year for each layer label, you supplied '" << years_.size() << "' but '" << layer_names_.size() << "' layer labels. Please sort this out, chairs =)";
}

//**********************************************************************
// void CNumericMetaLayer::build()
// Build the layer
//**********************************************************************
void CategoricalMetaLayer::DoBuild() {
  LOG_TRACE();
  for (auto year : years_) {
    if (find(model_->years().begin(), model_->years().end(), year) == model_->years().end())
      LOG_ERROR_P(PARAM_YEARS) << "the year " << year << " could not be found in model years, for obvious reasons the year must exist between start and final year in the @model block";
  }
  CategoricalLayer* layer = nullptr;
  for (unsigned i = 0; i < years_.size(); ++i) {
    layer = model_->managers().layer()->GetCategoricalLayer(layer_names_[i]);
    if (!layer)
      LOG_FATAL_P(PARAM_LAYER_LABELS) << "the layer = " << layer_names_[i] << " can not be found, does it exist? if it does make sure it is type categorical";

    if (layer->layer_type() != LayerType::kCategorical) {
      LOG_FATAL_P(PARAM_LAYER_LABELS) << "the layer = " << layer_names_[i] << " is not categorical, can you please sort that out";

    }
    years_layer_[years_[i]] = layer;
  }

  default_layer_ = model_->managers().layer()->GetCategoricalLayer(default_Layer_label_);
  if (!default_layer_)
    LOG_FATAL_P(PARAM_DEFAULT_LAYER) << "the default layer = " << default_Layer_label_ << " can not be found, does it exist? if it does make sure it is type categorical";

}


//**********************************************************************
// double getValue(int RowIndex, int ColIndex, int TargetRow, int TargetCol)
// get value
//**********************************************************************
string CategoricalMetaLayer::get_value(unsigned RowIndex, unsigned ColIndex) {
  LOG_FINEST();
  string value = "";
  if ((model_->current_year() == 0) || (model_->state() == State::kInitialise)) {
    LOG_FINEST() << "return default value";
    return default_layer_->get_value(RowIndex, ColIndex);
  }
  LOG_FINEST() << "return value from year " << model_->current_year();
  value = years_layer_[model_->current_year()]->get_value(RowIndex, ColIndex);
  return value;
}

string CategoricalMetaLayer::get_value(unsigned RowIndex, unsigned ColIndex, unsigned year) {
  LOG_FINEST();
  string value = "";
  if ((year == 0) || (find(years_.begin(),years_.end(),year) == years_.end()))
    return default_layer_->get_value(RowIndex, ColIndex);
  value = years_layer_[year]->get_value(RowIndex, ColIndex);
  return value;
}

} /* namespace layers */
} /* namespace niwa */
