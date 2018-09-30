/**
 * @file GlobalConfiguration.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 18/09/2012
 * @section LICENSE
 *
 * Copyright NIWA Science ©2012 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This class is a singleton object that holds some global
 * configuration data for the application.
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef GLOBALCONFIGURATION_H_
#define GLOBALCONFIGURATION_H_

// Headers
#include <map>
#include <vector>
#include <string>

#include "BaseClasses/Object.h"
#include "Translations/Translations.h"
#include "Utilities/RunParameters.h"
#include "Utilities/To.h"

// Namespaces
using std::map;
using std::vector;
using std::string;

namespace util = niwa::utilities;

namespace niwa {
class Model;

/**
 * Class Definitiion
 */
class GlobalConfiguration {
public:
  // Methods
  GlobalConfiguration() = default;
  virtual                     ~GlobalConfiguration() = default;
  void                        Clear();
  void                        ParseOptions(Model* model);

  // Accessors and Mutators
  void                  set_command_line_parameters(vector<string> &parameters) { command_line_parameters_ = parameters; }
  vector<string>&       command_line_parameters() { return command_line_parameters_; }
  string                addressable_value_file() const { return options_.addressable_value_input_file_; }
  void                  set_run_parameters(utilities::RunParameters& options) { options_ = options; }
  unsigned              random_seed() { return options_.random_number_seed_; }
  string                config_file() { return options_.config_file_; }
  unsigned              simulation_candidates() const { return options_.simulation_candidates_; }
  bool                  skip_loading_config_file() const { return skip_loading_config_file_; }
  void                  flag_skip_config_file() { skip_loading_config_file_ = true; }
  bool                  disable_standard_report() { return options_.no_std_report_; }

private:
  // Members
  vector<string>              command_line_parameters_;
  utilities::RunParameters    options_;
  bool                        skip_loading_config_file_;
};
} /* namespace niwa */
#endif /* GLOBALCONFIGURATION_H_ */
