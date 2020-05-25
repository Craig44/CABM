/**
 * @file GrowthSchnuteWithBasic.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 * The schnute Growth increment formula
 */

// headers
#include "GrowthSchnuteWithBasic.h"

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
GrowthSchnuteWithBasic::GrowthSchnuteWithBasic(Model* model) : Growth(model) {
  parameters_.Bind<string>(PARAM_ALPHA_LAYER_LABEL, &alpha_layer_label_, "Label for the numeric layer that describes mean L inf by area", "", true);
  parameters_.Bind<string>(PARAM_BETA_LAYER_LABEL, &beta_layer_label_, "Label for the numeric layer that describes mean k by area", "", true);
  parameters_.Bind<string>(PARAM_T0_LAYER_LABEL, &t0_layer_label_, "Label for the numeric layer that describes mean t0 by area", "", true);
  parameters_.Bind<string>(PARAM_A_LAYER_LABEL, &a_layer_label_, "Label for the numeric layer that describes mean a in the weight calcualtion through space", "", true);
  parameters_.Bind<string>(PARAM_B_LAYER_LABEL, &b_layer_label_, "Label for the numeric layer that describes mean b in the weight calcualtion through space", "", true);

  parameters_.Bind<float>(PARAM_T0, &t0_, "The value for t0 default = 0", "", true);
  parameters_.Bind<float>(PARAM_ALPHA, &alpha_, "alpha value for schnute growth curve", "", true);
  parameters_.Bind<float>(PARAM_BETA, &beta_, "beta value for schnute growth curve", "", true);
  parameters_.Bind<float>(PARAM_A, &a_, "alpha value for weight at length function", "", true);
  parameters_.Bind<float>(PARAM_B, &b_, "beta value for weight at length function", "", true);
  parameters_.Bind<float>(PARAM_TAU1, &tau1_, "reference age for y1", "");
  parameters_.Bind<float>(PARAM_TAU2, &tau2_, "reference age for y2", "");
  parameters_.Bind<float>(PARAM_Y1, &y1_, "mean size at reference ages tau1", "");
  parameters_.Bind<float>(PARAM_Y2, &y2_, "mean size at reference ages tau2", "");

  RegisterAsAddressable(PARAM_ALPHA, &alpha_, addressable::kAll, addressable::kyes);
  RegisterAsAddressable(PARAM_BETA, &beta_, addressable::kAll, addressable::kyes);
  RegisterAsAddressable(PARAM_TAU1, &tau1_, addressable::kAll, addressable::kyes);
  RegisterAsAddressable(PARAM_TAU2, &tau2_, addressable::kAll, addressable::kyes);
  model_->set_growth_model(niwa::Growth::kSchnute);
}

// check users have defined spatial or non_spatial
void GrowthSchnuteWithBasic::DoValidate() {
  // Check user hasn't specified them both
  if (parameters_.Get(PARAM_ALPHA_LAYER_LABEL)->has_been_defined() & parameters_.Get(PARAM_ALPHA)->has_been_defined())
    LOG_ERROR_P(PARAM_ALPHA) << "You have specified both a value and a layer label for the " << PARAM_ALPHA << " parameter, you must pick one or the other";
  if (parameters_.Get(PARAM_BETA_LAYER_LABEL)->has_been_defined() & parameters_.Get(PARAM_BETA)->has_been_defined())
    LOG_ERROR_P(PARAM_BETA) << "You have specified both a value and a layer label for the " << PARAM_BETA << " parameter, you must pick one or the other";
  if (parameters_.Get(PARAM_A_LAYER_LABEL)->has_been_defined() & parameters_.Get(PARAM_A)->has_been_defined())
    LOG_ERROR_P(PARAM_A) << "You have specified both a value and a layer label for the " << PARAM_A << " parameter, you must pick one or the other";
  if (parameters_.Get(PARAM_B_LAYER_LABEL)->has_been_defined() & parameters_.Get(PARAM_B)->has_been_defined())
    LOG_ERROR_P(PARAM_B) << "You have specified both a value and a layer label for the " << PARAM_B << " parameter, you must pick one or the other";
  if(parameters_.Get(PARAM_T0)->has_been_defined() & parameters_.Get(PARAM_T0_LAYER_LABEL)->has_been_defined())
    LOG_ERROR_P(PARAM_LABEL) << "You have specified both a value or a layer label for the " << PARAM_T0 << " parameter, you must pick one or the other";
  // Check user has specified at least one
  if (!parameters_.Get(PARAM_A_LAYER_LABEL)->has_been_defined() & !parameters_.Get(PARAM_A)->has_been_defined())
    LOG_ERROR_P(PARAM_LABEL) << "You have not specified a value or a layer label for the 'a' parameter, you must pick one of these options";
  if (!parameters_.Get(PARAM_B_LAYER_LABEL)->has_been_defined() & !parameters_.Get(PARAM_B)->has_been_defined())
    LOG_ERROR_P(PARAM_LABEL) << "You have not specified a value or a layer label for the 'b' parameter, you must pick one of these options";
  if(!parameters_.Get(PARAM_T0)->has_been_defined() & !parameters_.Get(PARAM_T0_LAYER_LABEL)->has_been_defined())
    t0_.push_back(0.0); // Default value of t0 = 0
  if (!parameters_.Get(PARAM_ALPHA_LAYER_LABEL)->has_been_defined() & !parameters_.Get(PARAM_ALPHA)->has_been_defined())
    LOG_ERROR_P(PARAM_LABEL) << "You have not specified a value or a layer label for the " << PARAM_ALPHA << " parameter, you must pick one or the other";
  if (!parameters_.Get(PARAM_BETA_LAYER_LABEL)->has_been_defined() & !parameters_.Get(PARAM_BETA)->has_been_defined())
    LOG_ERROR_P(PARAM_LABEL) << "You have not specified a value or a layer label for the " << PARAM_BETA << " parameter, you must pick one or the other";

  for(auto& tau1 : tau1_) {
    if(tau1 != model_->min_age())
      LOG_ERROR_P(PARAM_TAU1) << "tau1 should be the same as min_age on the model.";
  }

  if (distribution_label_ == PARAM_NORMAL) {
    distribution_ = Distribution::kNormal;
    LOG_FINE() << "normal length increment";
  } else if (distribution_label_ == PARAM_LOGNORMAL) {
    distribution_ = Distribution::kLogNormal;
    LOG_FINE() << "lognormal length increment";
  } else {
    LOG_ERROR_P(PARAM_DISTRIBUTION) << "unknown distribution";
  }


}


// Build relationships between classes
void GrowthSchnuteWithBasic::DoBuild() {
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
  if (alpha_layer_label_.size() > 0) {
    niwa::layers::NumericLayer* temp_layer = nullptr;
    for (unsigned i = 0; i < alpha_layer_label_.size(); ++i) {
      temp_layer = model_->managers().layer()->GetNumericLayer(alpha_layer_label_[i]);
      if (!temp_layer) {
        LOG_FATAL_P(PARAM_LINF_LAYER_LABEL) << "could not find the layer '" << alpha_layer_label_[i] << "', please make sure it exists and that it is type 'numeric'";
      }
      alpha_layer_.push_back(temp_layer);
    }
  }
  if (t0_layer_label_.size() > 0) {
    niwa::layers::NumericLayer* temp_layer = nullptr;
    for (unsigned i = 0; i < t0_layer_label_.size(); ++i) {
      temp_layer = model_->managers().layer()->GetNumericLayer(t0_layer_label_[i]);
      if (!temp_layer) {
        LOG_FATAL_P(PARAM_T0_LAYER_LABEL) << "could not find the layer '" << t0_layer_label_[i] << "', please make sure it exists and that it is type 'numeric'";
      }
      t0_layer_.push_back(temp_layer);
    }
  }
  if (beta_layer_label_.size() > 0) {
    niwa::layers::NumericLayer* temp_layer = nullptr;
    for (unsigned i = 0; i < beta_layer_label_.size(); ++i) {
      temp_layer = model_->managers().layer()->GetNumericLayer(beta_layer_label_[i]);
      if (!temp_layer) {
        LOG_FATAL_P(PARAM_K_LAYER_LABEL) << "could not find the layer '" << beta_layer_label_[i] << "', please make sure it exists and that it is type 'numeric'";
      }
      alpha_layer_.push_back(temp_layer);
    }
  }


  if (a_layer_label_.size() > 0) {
    niwa::layers::NumericLayer* temp_layer = nullptr;
    for (unsigned i = 0; i < a_layer_label_.size(); ++i) {
      temp_layer = model_->managers().layer()->GetNumericLayer(a_layer_label_[i]);
      if (!temp_layer) {
        LOG_FATAL_P(PARAM_A_LAYER_LABEL) << "could not find the layer '" << a_layer_label_[i] << "', please make sure it exists and that it is type 'numeric'";
      }
      a_layer_.push_back(temp_layer);
    }
  }
  if (b_layer_label_.size() > 0) {
    niwa::layers::NumericLayer* temp_layer = nullptr;
    for (unsigned i = 0; i < b_layer_label_.size(); ++i) {
      temp_layer = model_->managers().layer()->GetNumericLayer(b_layer_label_[i]);
      if (!temp_layer) {
        LOG_FATAL_P(PARAM_B_LAYER_LABEL) << "could not find the layer '" << b_layer_label_[i] << "', please make sure it exists and that it is type 'numeric'";
      }
      b_layer_.push_back(temp_layer);
    }
  }
  // Check that there is consistency in model sex and growth sex
  if (model_->get_sexed()) {
    // LAYERS
    if (b_layer_.size() == 1)
      b_layer_.assign(2, b_layer_[0]);
    if (a_layer_.size() == 1)
      a_layer_.assign(2, a_layer_[0]);
    if (alpha_layer_.size() == 1)
      alpha_layer_.assign(2, alpha_layer_[0]);
    if (beta_layer_.size() == 1)
      beta_layer_.assign(2, beta_layer_[0]);
    if (t0_layer_.size() == 1)
      t0_layer_.assign(2, t0_layer_[0]);
    // Values
    if (tau1_.size() == 1)
      tau1_.assign(2, tau1_[0]);
    if (tau2_.size() == 1)
      tau2_.assign(2, tau2_[0]);
    if (y1_.size() == 1)
      y1_.assign(2, y1_[0]);
    if (y2_.size() == 1)
      y2_.assign(2, y2_[0]);
    if (beta_.size() == 1)
      beta_.assign(2, beta_[0]);
    if (alpha_.size() == 1)
      alpha_.assign(2, alpha_[0]);
    if (a_.size() == 1)
      a_.assign(2, a_[0]);
    if (b_.size() == 1)
      b_.assign(2, b_[0]);
    if (t0_.size() == 1)
      t0_.assign(2, t0_[0]);
    if (t0_.size() == 0)
      t0_.assign(2, 0.0);
  }

	// Check that the layers are all positive
}

// The work horse of the function.
void GrowthSchnuteWithBasic::ApplyStochasticGrowth(vector<Agent>& agents) {
  //utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  float new_length = 0.0;
  float weight = 0.0;
  float length_prop = time_step_ratios_[model_->managers().time_step()->current_time_step()];
  float alpha, beta;
  unsigned sex_ndx;
  for(auto& agent : agents) {
    if (agent.is_alive()) {
      sex_ndx = agent.get_sex();
      alpha = agent.get_first_age_length_par();
      beta = agent.get_second_age_length_par();
      new_length = pow((pow(agent.get_length(), beta) * exp(-alpha * length_prop) + ((pow(y2_[sex_ndx], beta) - pow(y1_[sex_ndx], beta) * exp(-alpha * (tau2_[sex_ndx] - tau1_[sex_ndx]))) * (1 - exp(-alpha * length_prop))) / (1 -  exp(-alpha * (tau2_[sex_ndx] - tau1_[sex_ndx])))),1  /beta);
      weight = agent.get_first_length_weight_par() * pow(new_length, agent.get_second_length_weight_par());
      //LOG_FINEST() << "old length = " << (*iter).get_length() << " new length = " << new_length << " weight = " << weight << " k = " << (*iter).get_second_age_length_par() << " linf = " << (*iter).get_first_age_length_par();
      agent.set_length(new_length);
      agent.set_weight(weight);
    }
  }
}

// Execute the process
void GrowthSchnuteWithBasic::DoExecute() {
  LOG_MEDIUM();
 // #pragma omp parallel for collapse(2)
  for (unsigned row = 0; row < model_->get_height(); ++row) {
    for (unsigned col = 0; col < model_->get_width(); ++col) {
      WorldCell* cell = world_->get_base_square(row, col);
      if (cell->is_enabled()) {
        // Apply growth to elements in a cell.
        ApplyStochasticGrowth(cell->agents_);
        ApplyStochasticGrowth(cell->tagged_agents_);

      }
    }
  }
}

/*
 * This method is called at when ever an agent is created/seeded or moves (WorldCell.cpp) . Agents will get a new/updated growth parameters
 * based on the spatial cells of the process. This is called in initialisation/Recruitment and movement processes if needed.
*/
void  GrowthSchnuteWithBasic::draw_growth_param(unsigned row, unsigned col, unsigned number_of_draws, vector<vector<float>>& vec, unsigned sex) {
  LOG_FINE() << "sex " << sex;
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  float mean_alpha, mean_beta, a, b, t0;
  if (alpha_layer_.size() > 0)
    mean_alpha = alpha_layer_[sex]->get_value(row, col);
  else
    mean_alpha = alpha_[sex];

  if (beta_layer_.size() > 0)
    mean_beta = beta_layer_[sex]->get_value(row, col);
  else
    mean_beta = beta_[sex];

  if (t0_layer_.size() > 0)
    t0 = t0_layer_[sex]->get_value(row, col);
  else
    t0 = t0_[sex];

  if (a_layer_.size() > 0)
    a = a_layer_[sex]->get_value(row, col);
  else
    a = a_[sex];

  if (b_layer_.size() > 0)
    b = b_layer_[sex]->get_value(row, col);
  else
    b = b_[sex];

  vec.clear();
	vec.resize(number_of_draws);
  LOG_FINE() << "mean_alpha " << mean_alpha << " mean_beta " << mean_beta << " a " << a << " b " << b;// << " t0 " << t0_[sex];

  if (distribution_ == Distribution::kNormal) {
    for (unsigned i = 0; i < number_of_draws; ++i) {
      vec[i].push_back(rng.normal(mean_alpha, cv_ * mean_alpha));
      vec[i].push_back(rng.normal(mean_beta, cv_ * mean_beta));
      vec[i].push_back(t0);
      vec[i].push_back(y1_[sex]);
      vec[i].push_back(y2_[sex]);
      vec[i].push_back(tau2_[sex]);
      vec[i].push_back(a);
      vec[i].push_back(b);
    }
  } else if(distribution_ == Distribution::kLogNormal) {
    for (unsigned i = 0; i < number_of_draws; ++i) {
      vec[i].push_back(rng.lognormal(mean_alpha, cv_));
      vec[i].push_back(rng.lognormal(mean_beta, cv_));
      //vec[i].push_back(mean_alpha);
      //vec[i].push_back(mean_beta);
      vec[i].push_back(t0);
      vec[i].push_back(y1_[sex]);
      vec[i].push_back(y2_[sex]);
      vec[i].push_back(tau2_[sex]);
      vec[i].push_back(a);
      vec[i].push_back(b);
    }
  }

}
// FillReportCache, called in the report class, it will print out additional information that is stored in
// containers in this class.
void GrowthSchnuteWithBasic::FillReportCache(ostringstream& cache) {
  LOG_TRACE();

/*  if ((L_inf_layer_.size() > 0) | (k_layer_.size() > 0)) {// | (a_layer_  != nullptr) | (b_layer_ != nullptr)) {
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


        if (a_layer_)
          a = a_layer_->get_value(row, col);
        else
          a = a_;

        if (b_layer_)
          b = b_layer_->get_value(row, col);
        else
          b = b_;

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
  }*/
}
// A Method for telling the world we need to redistribute Mortality parmaeters
void GrowthSchnuteWithBasic::RebuildCache() {
  LOG_FINE();
  world_->rebuild_growth_params();
}
} /* namespace processes */
} /* namespace niwa */
