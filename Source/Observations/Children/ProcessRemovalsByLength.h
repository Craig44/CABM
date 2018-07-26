/**
 * @file ProcessRemovalsByLength.h
 * @author  C Marsh
 * @version 1.0
 * @date 25/08/15
 * @section LICENSE
 *
 * Copyright NIWA Science 2016 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This observation is a specific process observation class. It is associated with the process type mortality_instantaneous.
 *
 */
#ifndef OBSERVATIONS_REMOVEALS_BY_LENGTH_H_
#define OBSERVATIONS_REMOVEALS_BY_LENGTH_H_

// Headers
#include "Observations/Observation.h"

#include "Processes/Children/Mortality.h"
#include "Layers/Children/CategoricalLayer.h"

// Namespace
namespace niwa {
namespace observations {

using processes::Mortality;

/**
 * Class Definition
 */
class ProcessRemovalsByLength : public niwa::Observation {
public:
  // Methods
  explicit ProcessRemovalsByLength(Model* model);
  virtual                     ~ProcessRemovalsByLength();
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        PreExecute() override final;
  void                        Execute() override final;
  void                        Simulate() override final;
  bool                        HasYear(unsigned year) const override final { return std::find(years_.begin(), years_.end(), year) != years_.end(); }

protected:
  // Members
  vector<unsigned>              years_;
  parameters::Table*            error_values_table_ = nullptr;
  Mortality*                    mortality_process_ = nullptr;
  string                        process_label_;

  map<unsigned, vector<float>>  error_values_by_year_;
  map<unsigned, vector<float>>  error_values_;
  vector<string>                cells_;
  vector<unsigned>              cell_rows_;
  vector<unsigned>              cell_cols_;
  layers::CategoricalLayer*     layer_ = nullptr;
  string                        layer_label_;


};

} /* namespace observations */
} /* namespace niwa */

#endif /* OBSERVATIONS_REMOVEALS_BY_LENGTH_H_ */
