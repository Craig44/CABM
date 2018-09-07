/**
 * @file GrowthVonBertalanffyWithBasic.cpp
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
#include "GrowthVonBertalanffyWithBasic.h"

#include "Layers/Manager.h"
#include "Utilities/RandomNumberGenerator.h"
#include "World/WorldCell.h"
#include "World/WorldView.h"
#include "TimeSteps/Manager.h"

#include <omp.h>

// namespaces
namespace niwa {
namespace processes {

/**
 * constructor
 */
GrowthVonBertalanffyWithBasic::GrowthVonBertalanffyWithBasic(Model* model) : Growth(model) {
  parameters_.Bind<string>(PARAM_LINF_LAYER_LABEL, &l_inf_layer_label_, "Label for the numeric layer that describes mean L_inf through space", "", "");
  parameters_.Bind<string>(PARAM_K_LAYER_LABEL, &k_layer_label_, "Label for the numeric layer that describes mean k through space", "", "");
  parameters_.Bind<string>(PARAM_A_LAYER_LABEL, &a_layer_label_, "Label for the numeric layer that describes mean a in the weight calcualtion through space", "", "");
  parameters_.Bind<string>(PARAM_B_LAYER_LABEL, &b_layer_label_, "Label for the numeric layer that describes mean b in the weight calcualtion through space", "", "");
  parameters_.Bind<float>(PARAM_LINF, &l_inf_, "Label for the numeric layer that describes mean L_inf through space", "", 0);
  parameters_.Bind<float>(PARAM_K, &k_, "Label for the numeric layer that describes mean k through space", "", 0);
  parameters_.Bind<float>(PARAM_A, &a_, "Label for the numeric layer that describes mean a in the weight calcualtion through space", "", 0);
  parameters_.Bind<float>(PARAM_B, &b_, "Label for the numeric layer that describes mean b in the weight calcualtion through space", "", 0);
}

// check users have defined spatial or non_spatial
void GrowthVonBertalanffyWithBasic::DoValidate() {
  if (parameters_.Get(PARAM_LINF_LAYER_LABEL)->has_been_defined() & parameters_.Get(PARAM_LINF)->has_been_defined())
    LOG_ERROR_P(PARAM_LINF) << "You have specified both a value and a layer label for the 'L_inf' parameter, you must pick one or the other";
  if (parameters_.Get(PARAM_K_LAYER_LABEL)->has_been_defined() & parameters_.Get(PARAM_K)->has_been_defined())
    LOG_ERROR_P(PARAM_K) << "You have specified both a value and a layer label for the 'k' parameter, you must pick one or the other";
  if (parameters_.Get(PARAM_A_LAYER_LABEL)->has_been_defined() & parameters_.Get(PARAM_A)->has_been_defined())
    LOG_ERROR_P(PARAM_A) << "You have specified both a value and a layer label for the 'a' parameter, you must pick one or the other";
  if (parameters_.Get(PARAM_B_LAYER_LABEL)->has_been_defined() & parameters_.Get(PARAM_B)->has_been_defined())
    LOG_ERROR_P(PARAM_B) << "You have specified both a value and a layer label for the 'b' parameter, you must pick one or the other";

  if (!parameters_.Get(PARAM_LINF_LAYER_LABEL)->has_been_defined() & !parameters_.Get(PARAM_LINF)->has_been_defined())
    LOG_ERROR_P(PARAM_LABEL) << "You have not specified a value or a layer label for the 'L_inf' parameter, you must pick one of these options";
  if (!parameters_.Get(PARAM_K_LAYER_LABEL)->has_been_defined() & !parameters_.Get(PARAM_K)->has_been_defined())
    LOG_ERROR_P(PARAM_LABEL) << "You have not specified a value or a layer label for the 'k' parameter, you must pick one of these options";
  if (!parameters_.Get(PARAM_A_LAYER_LABEL)->has_been_defined() & !parameters_.Get(PARAM_A)->has_been_defined())
    LOG_ERROR_P(PARAM_LABEL) << "You have not specified a value or a layer label for the 'a' parameter, you must pick one of these options";
  if (!parameters_.Get(PARAM_B_LAYER_LABEL)->has_been_defined() & !parameters_.Get(PARAM_B)->has_been_defined())
    LOG_ERROR_P(PARAM_LABEL) << "You have not specified a value or a layer label for the 'b' parameter, you must pick one of these options";

}


// Build relationships between classes
void GrowthVonBertalanffyWithBasic::DoBuild() {
	// Get the layers if they have been defined
  if (l_inf_layer_label_ != "") {
    L_inf_layer_ = model_->managers().layer()->GetNumericLayer(l_inf_layer_label_);
    if (!L_inf_layer_) {
      LOG_ERROR_P(PARAM_LINF_LAYER_LABEL) << "could not find the layer '" << l_inf_layer_label_ << "', please make sure it exists and that it is type 'numeric'";
    }
  }
  if (k_layer_label_ != "") {
    k_layer_ = model_->managers().layer()->GetNumericLayer(k_layer_label_);
    if (!k_layer_) {
      LOG_ERROR_P(PARAM_K_LAYER_LABEL) << "could not find the layer '" << k_layer_label_ << "', please make sure it exists and that it is type 'numeric'";
    }
  }
  if (a_layer_label_ != "") {
    a_layer_ = model_->managers().layer()->GetNumericLayer(a_layer_label_);
    if (!a_layer_) {
      LOG_ERROR_P(PARAM_A_LAYER_LABEL) << "could not find the layer '" << a_layer_label_ << "', please make sure it exists and that it is type 'numeric'";
    }
  }
  if (b_layer_label_ != "") {
    b_layer_ = model_->managers().layer()->GetNumericLayer(b_layer_label_);
    if (!b_layer_) {
      LOG_ERROR_P(PARAM_B_LAYER_LABEL) << "could not find the layer '" << b_layer_label_ << "', please make sure it exists and that it is type 'numeric'";
    }
  }

	// Check that the layers are all positive

}

// Execute the process
void GrowthVonBertalanffyWithBasic::DoExecute() {
  LOG_TRACE();
  #pragma omp parallel for collapse(2)
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        float length_prop = time_step_proportions_[ model_->managers().time_step()->current_time_step()];
        for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++iter) {
          //LOG_FINEST() << "length = " << (*iter).get_length() << " weight = " << (*iter).get_weight() << " L-inf " << (*iter).get_first_age_length_par() << " k = " << (*iter).get_second_age_length_par();
          float new_length =  (*iter).get_length() + length_prop * ((*iter).get_first_age_length_par() - (*iter).get_length()) * (1 - exp(-(*iter).get_second_age_length_par()));
          float weight = (*iter).get_first_length_weight_par() * pow(new_length, (*iter).get_second_length_weight_par());
          //LOG_FINEST() << "length = " << new_length << " weight = " << weight;
          (*iter).set_length(new_length);
          (*iter).set_weight(weight);
        }
      }
    }
  }
}

/*
 * This method is called at when ever an agent is created/seeded or moves. Agents will get a new/updated growth parameters
 * based on the spatial cells of the process. This is called in initialisation/Recruitment and movement processes if needed.
*/
void  GrowthVonBertalanffyWithBasic::draw_growth_param(unsigned row, unsigned col, unsigned number_of_draws, vector<vector<float>>& vec) {
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  float mean_linf, mean_k, a, b;
  if (L_inf_layer_)
    mean_linf = L_inf_layer_->get_value(row, col);
  else
    mean_linf = l_inf_;

  if (k_layer_)
	  mean_k = k_layer_->get_value(row, col);
  else
    mean_k = k_;

  if (a_layer_)
    a = a_layer_->get_value(row, col);
  else
    a = a_;

  if (b_layer_)
    b = b_layer_->get_value(row, col);
  else
    b = b_;

  vec.clear();
	vec.resize(number_of_draws);

	for (unsigned i = 0; i < number_of_draws; ++i) {
	  vec[i].push_back(rng.lognormal(mean_linf, cv_));
	  vec[i].push_back(rng.lognormal(mean_k, cv_));
	  //LOG_FINEST() << vec[i][0] << " " << vec[i][1];
	  vec[i].push_back(a);
	  vec[i].push_back(b);
	}
}

} /* namespace processes */
} /* namespace niwa */
