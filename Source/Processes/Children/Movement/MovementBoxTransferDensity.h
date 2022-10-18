/**
 * @file MovementBoxTransferDensity.h
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
#ifndef SOURCE_PROCESSES_CHILDREN_MOVEMENT_BOX_TRANSFER_WITH_DENSITY_H_
#define SOURCE_PROCESSES_CHILDREN_MOVEMENT_BOX_TRANSFER_WITH_DENSITY_H_

// headers
#include "Processes/Children/Movement.h"
#include "Layers/Children/NumericLayer.h"
#include "Agents/Agent.h"
#include "DerivedQuantities/DerivedQuantity.h"

// namespaces
namespace niwa {
class Selectivity;
namespace processes {




/**
 * Class definition
 */
class MovementBoxTransferDensity : public Movement {
public:
  // methods
  explicit MovementBoxTransferDensity(Model* model);
  virtual                     ~MovementBoxTransferDensity() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final;
  void                        ApplyStochasticMovement(vector<Agent>& agents, MovementData& store_infor, bool tagged_partition, unsigned& origin_element, unsigned& row, unsigned& col);

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
  string                      derived_quantity_label_;
  DerivedQuantity*            derived_quantity_ = nullptr;

  // objects for thread safety of rng
  vector<float>                       random_numbers_;
  unsigned                            n_agents_;
  vector<vector<unsigned>>            cell_offset_;

  // Selectivity
  vector<string>                      selectivity_label_;
  vector<Selectivity*>                selectivity_;
  bool                                selectivity_length_based_;

  double                              threshold_;
  string                              direction_;
  bool                                scale_dq_by_initial_value_;
  map<unsigned, double>               dq_by_year_;
  vector<unsigned>                    years_;
  map<unsigned, string>               applied_movement_;

};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MOVEMENT_BOX_TRANSFER_WITH_DENSITY_H_ */
