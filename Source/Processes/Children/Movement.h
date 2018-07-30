/**
 * @file Movement.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 30/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process does nothing. It's used for debugging time steps
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MOVEMENT_H_
#define SOURCE_PROCESSES_CHILDREN_MOVEMENT_H_

// headers
#include "Processes/Process.h"

// namespaces
namespace niwa {
namespace processes {

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
class Movement : public Process {
public:
  // methods
  explicit Movement(Model* model);
  virtual                     ~Movement() = default;
  virtual void                        DoValidate() {};
  virtual void                        DoBuild() {};
  virtual void                        DoReset() {};
  virtual void                        DoExecute() {};
  virtual void                        FillReportCache(ostringstream& cache) {};

protected:
  // Report containers
  vector<MovementData>        moved_agents_by_year_;
};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MOVEMENT_H_ */
