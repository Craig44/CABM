/**
 * @file GrowthVonBertalanffy.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This file exists to keep documentation generator happy
 */

// headers
#include "GrowthVonBertalanffy.h"

#include "Layers/Manager.h"
// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
GrowthVonBertalanffy::GrowthVonBertalanffy(Model* model) : Process(model) {
  process_type_ = ProcessType::kGrowth;
  parameters_.Bind<string>(PARAM_L_INF_LAYER_LABEL, &l_inf_layer_label_, "Label for the numeric layer that describes mean L_inf through space", "");
  parameters_.Bind<string>(PARAM_K_LAYER_LABEL, &k_layer_label_, "Label for the numeric layer that describes mean k through space", "");
  parameters_.Bind<double>(PARAM_TIME_STEP_PROPORTION, &time_step_proportions_, "A vector of proportions that describe how much of the annual increment to add to the lenght of an agent in each time step this process is applied", "");
  parameters_.Bind<string>(PARAM_DISTRIBUTION, &distribution_, "the distribution to allocate the parameters to the agents", "");
  parameters_.Bind<double>(PARAM_CV, &cv_, "The cv of the distribution", "");
}

// Build relationships between classes
void GrowthVonBertalanffy::DoBuild() {
	// Get the layers
	L_inf_layer_ = model_->managers().layer()->GetNumericLayer(l_inf_layer_label_);
	if (!L_inf_layer_) {
		LOG_ERROR_P(PARAM_L_INF_LAYER_LABEL) << "could not find the layer '" << l_inf_layer_label_ << "', please make sure it exists and that it is type 'numeric'";
	}
	k_layer_ = model_->managers().layer()->GetNumericLayer(k_layer_label_);
	if (!k_layer_) {
		LOG_ERROR_P(PARAM_K_LAYER_LABEL) << "could not find the layer '" << k_layer_label_ << "', please make sure it exists and that it is type 'numeric'";
	}

	// Check time steps are the same as this process is subscribed to

	// Check that the layers are all positive


}

vector<double> GrowthVonBertalanffy::get_growth_params_by_cell(unsigned row, unsigned col) {
	vector<double> result;
	result.push_back(L_inf_layer_->get_value(row, col));
	result.push_back(k_layer_->get_value(row, col));
	return result;
}

} /* namespace processes */
} /* namespace niwa */