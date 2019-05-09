/**
 * @file CommandLineParser.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 16/11/2012
 * @section LICENSE
 *
 * Copyright NIWA Science ©2012 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

 // Headers
#include "CommandLineParser.h"

#include <iostream>
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>

#include "License.h"
#include "Version.h"

// Namespaces
namespace niwa {
namespace utilities {

using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::cout;
using std::endl;
using std::ostringstream;

/**
 * This method will take the raw command line input from the main() method
 * and process them into something more useful.
 *
 * @param argc The number of arguments that have been provided
 * @param argv Pointer to an array containing the arguments
 * @param options The options object to fille with the values
 */
void CommandLineParser::Parse(int argc, char* argv[], RunParameters& options) {
  // Build Options menu
  options_description oDesc("Usage");
  oDesc.add_options()
    ("help,h", "Print help")
    ("license,l", "Display IBM license")
    ("version,v", "Display version information")
    ("simulation,s", value<unsigned>(), "Simulation mode (arg = number of candidates)")
    ("config,c", value<string>(), "Configuration file")
    ("run,r", "Basic model run mode")
    ("input,i", value<string>(), "Load free parameter values from file")
    ("seed,g", value<unsigned>(), "Random number seed")
    ("loglevel", value<string>(), "Set log level: finest, fine, trace, none(default)");


  ostringstream o;
  o << oDesc;
  command_line_usage_ = o.str();

  // Read our Parameters
  variables_map parameters;

  try {
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, oDesc), parameters);
    notify(parameters);

  } catch (boost::program_options::unknown_option &ex) {
    cout << "An error occurred while processing the command line. " << ex.what() << endl;
  }

  /**
   * Load any variables into the global config that need to be available
   * immediately
   */
  if (parameters.count("loglevel")) {
    options.log_level_ = parameters["loglevel"].as<string>();
    Logging::Instance().SetLogLevel(options.log_level_);
  }

  LOG_TRACE();


  if (parameters.count("config"))
    options.config_file_ = parameters["config"].as<string>();
  if (parameters.count("input"))
    options.addressable_value_input_file_ = parameters["input"].as<string>();

  /**
   * Determine what run mode we should be in. If we're
   * in help, version or license then we don't need to continue.
   */
  if ( (parameters.count("help")) || (parameters.size() == 0) ) {
    options.run_mode_ = RunMode::kHelp;
    cout << command_line_usage_ << endl;
    return;

  } else if (parameters.count("version")) {
    options.run_mode_ = RunMode::kVersion;
    cout << SOURCE_CONTROL_VERSION << endl;
    return;

  } else if (parameters.count("license")) {
    options.run_mode_ = RunMode::kLicense;
    cout << license << endl;
    return;

  }

  /**
   * At this point we know we've been asked to do an actual model
   * run. So we need to check to ensure the command line makes
   * sense.
   */
  unsigned run_mode_count = 0;
  run_mode_count += parameters.count("run");
  run_mode_count += parameters.count("simulation");

  if (run_mode_count == 0)
    LOG_ERROR() << "No valid run mode has been specified on the command line. Please specify a valid run mode (e.g -r)";
  if (run_mode_count > 1)
    LOG_ERROR() << "Multiple run modes have been specified on the command line. Only 1 run mode is valid";

  if (parameters.count("run"))
    options.run_mode_ = RunMode::kBasic;
  else if (parameters.count("simulation")) {
     options.simulation_candidates_ = parameters["simulation"].as<unsigned>();
     options.run_mode_ = RunMode::kSimulation;

   } else {
    LOG_ERROR() << "An invalid or unknown run mode has been specified on the command line.";
  }

  /**
   * Now we store any variables we want to use to override global defaults.
   */
  if (parameters.count("seed")) {
    options.override_rng_seed_value_ = parameters["seed"].as<unsigned>();
    options.override_random_number_seed_ = true;
  }

}

} /* namespace utilities */
} /* namespace niwa */








