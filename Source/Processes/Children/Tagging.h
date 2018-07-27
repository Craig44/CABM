/**
 * @file Tagging.h
 * @author C>Marsh
 * @github https://github.com/Craig44
 * @date 26/7/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process does nothing. It's used for debugging time steps
 */
#ifndef SOURCE_PROCESSES_CHILDREN_TAGGING_H_
#define SOURCE_PROCESSES_CHILDREN_TAGGING_H_

// headers
#include "Processes/Process.h"

#include "Layers/Children/CategoricalLayer.h"


// namespaces
namespace niwa {
class Selectivity;
namespace processes {

/**
 * Class definition
 */
class Tagging : public Process {
public:
  // methods
  explicit Tagging(Model* model);
  virtual                     ~Tagging() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final;

protected:
  map<unsigned, vector<float>> numbers_;
  parameters::Table*           numbers_table_ = nullptr;
  vector<unsigned>             years_;
  vector<string>               selectivity_labels_;
  vector<Selectivity*>         selectivities_;

};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_TAGGING_H_ */
