/**
 * @file MovementBoxTransfer.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 19/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 * User defines a matrix of probabilities with a selectivity, that describes the probability of moving to another cell.
 *
 *
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MOVEMENT_BOX_TRANSFER_H_
#define SOURCE_PROCESSES_CHILDREN_MOVEMENT_BOX_TRANSFER_H_

// headers
#include "Processes/Children/Movement.h"
#include "Layers/Children/NumericLayer.h"
#include <omp.h>

// namespaces
namespace niwa {
class Selectivity;
namespace processes {


enum class MovementType {
  kUnknown = 0,
  kMarkovian = 1,
  kNatal_homing = 2
};



/**
 * Class definition
 */
class MovementBoxTransfer : public Movement {
public:
  // methods
  explicit MovementBoxTransfer(Model* model);
  virtual                     ~MovementBoxTransfer() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final;
  void                        FillReportCache(ostringstream& cache)  override final;
protected:
  vector<string>              origin_cell_;
  vector<string>              probability_layer_labels_;
  vector<layers::NumericLayer*>  probability_layers_;
  vector<unsigned>            origin_rows_;
  vector<unsigned>            origin_cols_;
  vector<unsigned>            possible_rows_;
  vector<unsigned>            possible_cols_;
  string                      movement_type_string_;
  MovementType                movement_type_ = MovementType::kUnknown;

  // objects for thread safety of rng
  vector<float>                       random_numbers_;
  unsigned                            n_agents_;
  vector<vector<unsigned>>            cell_offset_;

  // Selectivity
  vector<string>                      selectivity_label_;
  vector<Selectivity*>                selectivity_;
  bool                                selectivity_length_based_;



};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MOVEMENT_BOX_TRANSFER_H_ */
