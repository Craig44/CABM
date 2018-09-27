/**
 * @file MortalityBaranov.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 21/09/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This is a child mortality that applies a fishing event
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_BARANOV_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_BARANOV_H_

// headers
#include "Processes/Children/Mortality.h"

#include "Layers/Children/NumericLayer.h"
#include <omp.h>


// namespaces
namespace niwa {
class Selectivity;
namespace processes {
using std::string;
/**
 * Class definition
 */
class MortalityBaranov : public Mortality {
public:
  // methods
  explicit MortalityBaranov(Model* model);
  virtual                     ~MortalityBaranov() = default;
  virtual void                        DoValidate() override final;
  virtual void                        DoBuild() override final;
  virtual void                        DoReset() override final { };
  virtual void                        DoExecute() override final;
  void                                draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) override final;
  void                                FillReportCache(ostringstream& cache) override final;
  void                                RebuildCache() override final;

protected:
  // F variables
  vector<string>                      f_layer_label_;
  vector<layers::NumericLayer*>       f_layer_;
  vector<string>                      fishery_selectivity_label_;
  vector<Selectivity*>                fishery_selectivity_;
  vector<unsigned>                    years_;
  float                               discard_mortality_;
  float                               mls_;
  bool                                f_selectivity_length_based_ = false;
  map<unsigned,float>                 catch_by_year_;
  // M variables
  vector<string>                      natural_mortality_selectivity_label_;
  string                              m_layer_label_;
  layers::NumericLayer*               m_layer_ = nullptr;
  vector<Selectivity*>                natural_mortality_selectivity_;
  float                               cv_;
  string                              distribution_;
  bool                                m_selectivity_length_based_ = false;
  float                               m_;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_BARANOV_H_ */
