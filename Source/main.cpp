// defines
#ifndef BOOST_USE_WINDOWS_H
#define BOOST_USE_WINDOWS_H
#endif

// Headers
#include <iostream>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <thread>

#include "ConfigurationLoader/Loader.h"
#include "GlobalConfiguration/GlobalConfiguration.h"
#include "Model/Factory.h"
#include "Model/Managers.h"
#include "Model/Model.h"
#include "Reports/Children/StandardHeader.h"
#include "Reports/Manager.h"
#include "Utilities/CommandLineParser/CommandLineParser.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/RunParameters.h"
#include "Logging/Logging.h"

// Namespaces
using namespace niwa;
using std::cout;
using std::cerr;
using std::endl;

/**
 * Application entry point
 */
int main(int argc, char * argv[]) {
  int return_code = 0;
  bool model_start_return_success = true;

  try {
    Model model;
    reports::StandardHeader standard_report(&model);


    utilities::RunParameters parameters;

    utilities::CommandLineParser parser;
    parser.Parse(argc, argv, parameters);
    model.global_configuration().set_run_parameters(parameters);

    vector<string> cmd_parameters;
    for (int i = 0; i < argc; ++i) cmd_parameters.push_back(argv[i]);
    model.global_configuration().set_command_line_parameters(cmd_parameters);

    RunMode::Type run_mode = parameters.run_mode_;

    // Check if the right exe for this run mode
    #ifndef MSE
      if(run_mode == RunMode::kMSE)
        LOG_FATAL() << "Cannot run MSE runmode with this executable. You will need to recompile with the build command doBuild release_mse";
    #endif
    /**
     * Check the run mode and call the handler.
     */
    switch (run_mode) {
    case RunMode::kInvalid:
      LOG_ERROR() << "Invalid run mode specified.";
      break;
    case RunMode::kVersion:
    case RunMode::kHelp:
    case RunMode::kLicense:
      break;

    case RunMode::kSimulation:
    case RunMode::kMSE:
    case RunMode::kBasic:
      LOG_MEDIUM() << "set up global config";
      if (!model.global_configuration().disable_standard_report()) {
         standard_report.Prepare();
         model.managers().report()->set_std_header(standard_report.header());
       }
      LOG_MEDIUM() << "set up config loader";
       // load our configuration file
       configuration::Loader config_loader(model);
       if (!config_loader.LoadConfigFile()) {
         Logging::Instance().FlushErrors();
         return_code = -1;
         break;
       }
       LOG_MEDIUM() << "set up logging";
       Logging& logging = Logging::Instance();
        config_loader.ParseFileLines();
        if (logging.errors().size() > 0) {
          logging.FlushErrors();
          return_code = -1;
          break;
        }

       // override any config file values from our command line
       model.global_configuration().ParseOptions(&model);
       utilities::RandomNumberGenerator::Instance().Reset(model.global_configuration().random_seed());
       LOG_MEDIUM() << "RNG";

       // Thread off the reports
       reports::Manager* report_manager = model.managers().report();
       std::thread report_thread([&report_manager]() { report_manager->FlushReports(); });
       LOG_MEDIUM() << "Threading";

       // Run the model
       model_start_return_success = model.Start(run_mode);
       LOG_MEDIUM() << "Finished running the model";

       // finish report thread
       report_manager->StopThread();
       report_thread.join();

       if (logging.errors().size() > 0) {
         logging.FlushErrors();
         return_code = -1;
       } else if (not parameters.surpress_warnings_) {
         //cout << "about to flush warnings() " << parameters.surpress_warnings_ << "\n";
         logging.FlushWarnings();
       }

       if (!model.global_configuration().disable_standard_report())
         standard_report.Finalise();
       break;
    } // switch(run_mode)

  } catch (const string &exception) {
    cout << "## ERROR - CABM experienced a problem and has stopped execution" << endl;
    cout << "Error: " << exception << endl;
    return_code = -1;

  } catch (std::exception& e) {
    cout << "## ERROR - CABM experienced a problem and has stopped execution" << endl;
    cout << e.what() << endl;
    return_code = -1;

  } catch(...) {
    cout << "## ERROR - CABM experienced a problem and has stopped execution" << endl;
    cout << "The exception was caught with the catch-all. The type was unknown" << endl;
    cout << "Please contact the application developer" << endl;
    return_code = -1;
  }

  if (!model_start_return_success) {
    LOG_FINEST() << "Done with errors";
    return -1;
  }

  LOG_FINEST() << "Done";
	return return_code;
}


