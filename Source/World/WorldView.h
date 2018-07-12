/**
 * @file WorldView.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 12/07/2018
 * @section LICENSE
 *
 * @section DESCRIPTION
 *
 *    Following SPM - The world maintains a Grid (base_grid_) of WorldCell objects. This
 *    represents the current state of the world as it is now. A second
 *    grid called the Drifference Grid (cached_grid_) is an exact copy
 *    of pGrid in dimensions, but it maintains only the difference
 *    (adjustments) we need to make to the world. After each process is run
 *    the difference grid is merged into the world grid to update the current
 *    state.
 *
 *    The runMode calls allow the application entry (int main()) easy access
 *    to call functions based on each run mode. The World then sets up what
 *    is required to kick-off each type of run.
 *
 *    Note: CWorld is a singleton because we only want one ever. Multiple
 *    worlds are not prohibited.

 */
#ifndef WORLD_H_
#define WORLD_H_

// headers
#include <map>
#include <list>
#include <vector>
#include <string>

// namespaces
namespace niwa {

class Model;
class WorldCell;
class Layer;

using std::string;
/**
 * Class definition
 */
class WorldView {

public:
  // methods
  WorldView(Model* model) : model_(model) { };
  virtual                     ~WorldView() = default;
  void                        Validate() {};
  void                        Build();
  void                        Reset() {};

  // Accessors
  double                      get_abundance(void) {return 1.0; };
  double                      get_biomass(void) {return 1.0; };
  //unsigned                    get_height() { return height_; }  // these can be called off the model
  //unsigned                    get_width() { return width_; }  // these can be called off the model

  //TODO
  //WorldCell*                  get_base_square(int RowIndex, int ColIndex);
  //WorldCell*                  get_cached_square(int RowIndex, int ColIndex);
  void                        MergeCachedAgents() {};

protected:
  // members
  WorldCell                   **base_grid_;
  WorldCell                   **cached_grid_;
  Layer*                      base_layer_ = nullptr;


private:
  // members
  Model*                      model_ = nullptr;

};
} /* namespace niwa */

#endif /* WORLD_H_ */
