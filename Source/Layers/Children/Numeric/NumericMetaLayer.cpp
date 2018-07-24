//============================================================================
// Name        : CNumericMetaLayer.cpp
// Author      : S.Rasmussen
// Date        : 16/01/2009
// Copyright   : Copyright NIWA Science ©2009 - www.niwa.co.nz
// Description :
// $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
//============================================================================

// Local headers
#include "NumericMetaLayer.h"

#include "Layers/Manager.h"

namespace niwa {
namespace layers {

// Constructor
NumericMetaLayer::NumericMetaLayer(Model* model) : NumericLayer(model) {
  // Default Variables
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "the years to apply the following layer", "");
  parameters_.Bind<string>(PARAM_DEFAULT_LAYER, &default_Layer_label_, "The default layer to apply in initialisation phase", "");
  parameters_.Bind<string>(PARAM_LAYER_LABELS, &layer_names_, "The layer label associated with each year", "");
  layer_type_ = LayerType::kNumeric;

}

//**********************************************************************
// void CNumericMetaLayer::validate()
// Validate the layer
//**********************************************************************
void NumericMetaLayer::DoValidate() {
  LOG_TRACE();
  if (years_.size() != layer_names_.size())
    LOG_ERROR_P(PARAM_LAYER_LABELS) << "there needs to be a year for each layer label, you supplied '" << years_.size() << "' but '" << layer_names_.size() << "' layer labels. Please sort this out, chairs =)";
}

//**********************************************************************
// void CNumericMetaLayer::build()
// Build the layer
//**********************************************************************
void NumericMetaLayer::DoBuild() {
  LOG_TRACE();
  for (auto year : model_->years()) {
    if (find(years_.begin(), years_.end(), year) == years_.end())
      LOG_ERROR_P(PARAM_YEARS) << "the model year " << year << " could not be found in your years parameter, there must be a year for every model year unfortunatley, feel free to change this if you want.";
  }
  NumericLayer* layer = nullptr;
  for (unsigned i = 0; i < years_.size(); ++i) {
    layer = model_->managers().layer()->GetNumericLayer(layer_names_[i]);
    if (!layer)
      LOG_FATAL_P(PARAM_LAYER_LABELS) << "the layer = " << layer_names_[i] << " can not be found, does it exist? if it does make sure it is type numeric";

    if (layer->layer_type() != LayerType::kNumeric) {
      LOG_FATAL_P(PARAM_LAYER_LABELS) << "the layer = " << layer_names_[i] << " is not numeric, can you please sort that out";

    }
    years_layer_[years_[i]] = layer;
  }

  default_layer_ = model_->managers().layer()->GetNumericLayer(default_Layer_label_);
  if (!default_layer_)
    LOG_FATAL_P(PARAM_DEFAULT_LAYER) << "the default layer = " << default_Layer_label_ << " can not be found, does it exist? if it does make sure it is type numeric";


}


//**********************************************************************
// double getValue(int RowIndex, int ColIndex, int TargetRow, int TargetCol)
// get value
//**********************************************************************
float NumericMetaLayer::get_value(unsigned RowIndex, unsigned ColIndex) {
  float value = 0.0;
  if (model_->state() == State::kInitialise)
    return default_layer_->get_value(RowIndex, ColIndex);
  value = years_layer_[model_->current_year()]->get_value(RowIndex, ColIndex);
  return value;
}


} /* namespace layers */
} /* namespace niwa */
