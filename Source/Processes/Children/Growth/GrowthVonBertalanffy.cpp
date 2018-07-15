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
 */

// headers
#include "GrowthVonBertalanffy.h"

#include "Layers/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
GrowthVonBertalanffy::GrowthVonBertalanffy(Model* model) : Growth(model) {
  parameters_.Bind<string>(PARAM_L_INF_LAYER_LABEL, &l_inf_layer_label_, "Label for the numeric layer that describes mean L_inf through space", "");
  parameters_.Bind<string>(PARAM_K_LAYER_LABEL, &k_layer_label_, "Label for the numeric layer that describes mean k through space", "");
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


/*
 * This method is called at when ever an agent is created/seeded or moves. Agents will get a new/updated growth parameters
 * based on the spatial cells of the process. This is called in initialisation/Recruitment and movement processes if needed.
*/
void  GrowthVonBertalanffy::draw_growth_param(unsigned row, unsigned col, unsigned number_of_draws, vector<vector<double>>& vec) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
	double mean_linf = L_inf_layer_->get_value(row, col);
	double mean_k = k_layer_->get_value(row, col);
	vec.clear();
	vec.resize(number_of_draws);

	for (unsigned i = 0; i < number_of_draws; ++i) {
	  vec[i].push_back(rng.lognormal(mean_linf, cv_));
	  vec[i].push_back(rng.lognormal(mean_k, cv_));
	}
}

} /* namespace processes */
} /* namespace niwa */
