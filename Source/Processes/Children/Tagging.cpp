/**
 * @file Tagging.cpp
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 26/7/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This file exists to keep documentation generator happy
 */

// headers
#include "Tagging.h"

#include "Selectivities/Manager.h"
#include "Layers/Manager.h"


// namespaces
namespace niwa {
namespace processes {

/**
 * Empty constructor
 */
Tagging::Tagging(Model* model) : Process(model) {
  numbers_table_ = new parameters::Table(PARAM_NUMBERS);
  parameters_.BindTable(PARAM_NUMBERS, numbers_table_, "Table of N data by year, cell and length bin", "", true, true);
  parameters_.Bind<string>(PARAM_SELECTIVITIES, &selectivity_labels_, "", "");
  parameters_.Bind<unsigned>(PARAM_YEARS, &years_, "Years to execute the transition in", "");
  process_type_ = ProcessType::kTransition;
}


void Tagging::DoValidate() {

}

void Tagging::DoBuild() {

}

void Tagging::DoExecute() {

}

}
}
