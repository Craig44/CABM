/**
 * @file RunParameters.h
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date Jan 12, 2016
 * @section LICENSE
 *
 * Copyright NIWA Science ©2016 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This file contains the struct definition for the object we will
 * pass in to the CASAL2 DLL holding all of the optional information
 * we need to give it.
 */
#ifndef RUNPARAMETERS_H_
#define RUNPARAMETERS_H_

// headers
#include <string>

#include "Utilities/RunMode.h"

// namespaces
namespace niwa {
namespace utilities {
using std::string;

// structs
struct RunParameters {
  RunMode::Type run_mode_ = RunMode::kInvalid;
  string        config_file_ = "config.ibm";
  string        addressable_value_input_file_ = "";
  unsigned      random_number_seed_ = 123u;
  bool          debug_mode_ = false;
  bool          no_std_report_ = false;
  string        log_level_ = "warning";
  string        output_ = "";
  unsigned      simulation_candidates_ = 1u;
  bool          override_random_number_seed_ = false;
  unsigned      override_rng_seed_value_ = 123u;
};

} /* namespace niwa */
} /* namespace utilities */

#endif /* RUNPARAMETERS_H_ */
