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
  parameters_.Bind<string>(PARAM_LINF_LAYER_LABEL, &l_inf_layer_label_, "Label for the numeric layer that describes mean L inf through space", "", "");
  parameters_.Bind<string>(PARAM_K_LAYER_LABEL, &k_layer_label_, "Label for the numeric layer that describes mean k through space", "", "");
  parameters_.Bind<float>(PARAM_T0, &t0_, "The value for t0 default = 0", "", 0);
  parameters_.Bind<string>(PARAM_A_LAYER_LABEL, &a_layer_label_, "Label for the numeric layer that describes mean a in the weight calcualtion through space", "", "");
  parameters_.Bind<string>(PARAM_B_LAYER_LABEL, &b_layer_label_, "Label for the numeric layer that describes mean b in the weight calcualtion through space", "", "");
  parameters_.Bind<float>(PARAM_LINF, &l_inf_, "Value of mean L inf multiplied by the layer value if supplied", "", 0);
  parameters_.Bind<float>(PARAM_K, &k_, "Value of mean k multiplied by the layer value if supplied", "", 0);
  parameters_.Bind<float>(PARAM_A, &a_, "alpha value for weight at length function", "", 0);
  parameters_.Bind<float>(PARAM_B, &b_, "beta value for weight at length function", "", 0);

  RegisterAsAddressable(PARAM_LINF, &l_inf_);
  RegisterAsAddressable(PARAM_K, &k_);
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
  LOG_FINE();
  vector<TimeStep*> time_steps = model_->managers().time_step()->ordered_time_steps();
  LOG_FINEST() << "time_steps.size(): " << time_steps.size();
  vector<unsigned> active_time_steps;
  for (unsigned i = 0; i < time_steps.size(); ++i) {
    if (time_steps[i]->HasProcess(label_))
      active_time_steps.push_back(i);
  }

  if (time_step_proportions_.size() == 0) {
    for (unsigned i : active_time_steps)
      time_step_proportions_[i] = 1.0;
  } else {
    if (time_step_proportions_.size() != active_time_steps.size())
      LOG_FATAL_P(PARAM_TIME_STEP_PROPORTIONS) << " length (" << time_step_proportions_.size()
          << ") does not match the number of time steps this process has been assigned to (" << active_time_steps.size() << ")";

    for (float value : time_step_proportions_) {
      if (value < 0.0 || value > 1.0)
        LOG_ERROR_P(PARAM_TIME_STEP_PROPORTIONS) << " value (" << value << ") must be between 0.0 (exclusive) and 1.0 (inclusive)";
    }

    for (unsigned i = 0; i < time_step_proportions_.size(); ++i) {
      time_step_ratios_[active_time_steps[i]] = time_step_proportions_[i];
      LOG_FINE() << "setting growth in time step " << active_time_steps[i] << " = " << time_step_proportions_[i];
    }
  }
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
  LOG_MEDIUM();
 // #pragma omp parallel for collapse(2)
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        float length_prop = time_step_ratios_[model_->managers().time_step()->current_time_step()];
        //LOG_FINE() << "length prop = " << length_prop;
        for (auto iter = cell->agents_.begin(); iter != cell->agents_.end(); ++iter) {
          if ((*iter).is_alive()) {
            //LOG_FINEST() << "length = " << (*iter).get_length() << " weight = " << (*iter).get_weight() << " L-inf " << (*iter).get_first_age_length_par() << " k = " << (*iter).get_second_age_length_par() << " prop = " << length_prop;
            float new_length =  (*iter).get_length() + length_prop * ((*iter).get_first_age_length_par() - (*iter).get_length()) * (1 - exp(-(*iter).get_second_age_length_par()));
            float weight = (*iter).get_first_length_weight_par() * pow(new_length, (*iter).get_second_length_weight_par());
            LOG_FINEST() << "old length = " << (*iter).get_length() << " new length = " << new_length << " weight = " << weight << " k = " << (*iter).get_second_age_length_par() << " linf = " << (*iter).get_first_age_length_par();
            (*iter).set_length(new_length);
            (*iter).set_weight(weight);
          }
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
    vec[i].push_back(t0_);
	  vec[i].push_back(a);
	  vec[i].push_back(b);
	}
}
// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void GrowthVonBertalanffyWithBasic::FillReportCache(ostringstream& cache) {
  LOG_TRACE();

  if ((L_inf_layer_ != nullptr) | (k_layer_ != nullptr)) {// | (a_layer_  != nullptr) | (b_layer_ != nullptr)) {
    LOG_FINE() << "Spatially varying growth "<< REPORT_R_DATAFRAME << "\n" << "cell ";
    for (unsigned t = 0; t <= model_->max_age(); ++t)
      cache << t << " ";
    cache << "\n";

    for (unsigned row = 0; row < model_->get_height(); ++row) {
      for (unsigned col = 0; col < model_->get_width(); ++col) {
        vector<float> length;
        cache << row << "-" << col << " " << t0_ << " " ;
        float mean_linf = 0, mean_k = 0;//, a = 0, b = 0;
        if (L_inf_layer_)
          mean_linf = L_inf_layer_->get_value(row, col);
        else
          mean_linf = l_inf_;

        if (k_layer_)
          mean_k = k_layer_->get_value(row, col);
        else
          mean_k = k_;

        length.push_back(mean_linf * (1-exp(-mean_k * (0.001 - t0_))));

/*
        if (a_layer_)
          a = a_layer_->get_value(row, col);
        else
          a = a_;

        if (b_layer_)
          b = b_layer_->get_value(row, col);
        else
          b = b_;*/

        for (unsigned t = 1; t <= model_->max_age(); ++t) {
          length.push_back(length[t - 1]  + (mean_linf - length[t - 1]) * (1 - exp(-mean_k)));
          cache << length[t] << " ";
        }
        cache << "\n";
      }
    }
  } else {
    cache << "time: ";
    vector<float> length;
    length.push_back(l_inf_ * (1-exp(-k_ * (0.001 - t0_))));
    for (unsigned t = 0; t <= model_->max_age(); ++t) {
      cache << t << " ";
      if (t > 0)
        length.push_back(length[t - 1]  + (l_inf_ - length[t - 1]) * (1 - exp(-k_)));
    }
    cache << "\nlength: ";
    for (unsigned i = 0; i < length.size(); ++i) {
      cache << length[i] << " ";
    }
    cache << "\n";
  }
}
// A Method for telling the world we need to redistribute Mortality parmaeters
void GrowthVonBertalanffyWithBasic::RebuildCache() {
  LOG_FINE();
  world_->rebuild_growth_params();
}
} /* namespace processes */
} /* namespace niwa */
