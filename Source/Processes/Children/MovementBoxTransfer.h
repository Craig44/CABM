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
#include "Processes/Process.h"
#include "Layers/Children/NumericLayer.h"
#include <omp.h>

// namespaces
namespace niwa {
//class Selectivity;
namespace processes {


enum class MovementType {
  kUnknown = 0,
  kMarkovian = 1,
  kNatal_homing = 2
};


/**
 * A movement struct that stores movement information
 */
struct MovementData {
  string origin_cell_;
  unsigned year_;
  unsigned initial_numbers_ = 0;
  vector<vector<unsigned>> destination_of_agents_moved_;
  MovementData(unsigned rows, unsigned cols, string origin_cell, unsigned year) : origin_cell_(origin_cell), year_(year)
  {
    // set up matrix in constructor, save some sloppy run time code
    destination_of_agents_moved_.resize(rows);
    for (unsigned i = 0; i < rows; ++i)
      destination_of_agents_moved_[i].resize(cols);
  }
};

/**
 * Class definition
 */
class MovementBoxTransfer : public Process {
public:
  // methods
  explicit MovementBoxTransfer(Model* model);
  virtual                     ~MovementBoxTransfer() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final;
  void                        FillReportCache(ostringstream& cache);
protected:
  vector<unsigned>            years_;
//  string                      selectivity_label_;
//  Selectivity*                selectivity_ = nullptr;
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

  // Report containers
  vector<MovementData>        moved_agents_by_year_;

};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MOVEMENT_BOX_TRANSFER_H_ */
