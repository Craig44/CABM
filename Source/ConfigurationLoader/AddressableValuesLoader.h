/**
 * @file AddressableValuesLoader.h
 * @author Scott Rasmussen (scott.rasmussen@zaita.com) / C Marsh
 * @date 19/12/2014 modified 30/09/2018
 * @section LICENSE
 *
 * Copyright NIWA Science ©2014 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * The estimable values loader is responsible for loading a text file (csv)
 * of estimable values to do multiple iterations.
 *
 * Depending on the run mode the estimable values will be used different.
 * For a standard run the model will do N iterations where N is the amount
 * of values specified for the estimables.
 */
#ifndef CONFIGURATION_ADDRESSABLEVALUESLOADER_H_
#define CONFIGURATION_ADDRESSABLEVALUESLOADER_H_

// headers
#include <string>

#include "Model/Model.h"

// namespaces
namespace niwa {
namespace configuration {

using std::string;

// classes
class AddressableValuesLoader {
public:
  // methods
  AddressableValuesLoader(Model* model) : model_(model) { }
  virtual                     ~AddressableValuesLoader() = default;
  void                        LoadValues(const string& file_name);

private:
  // members
  Model*                    model_ = nullptr;
};

} /* namespace configuration */
} /* namespace niwa */

#endif /* CONFIGURATION_ADDRESSABLEVALUESLOADER_H_ */
