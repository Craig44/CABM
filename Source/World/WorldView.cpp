/**
 * @file WorldView.cpp
 * @author C.Marsh
 * @github
 * @date 11/06/2018
 * @section LICENSE
 *
 */

// headers
#include "WorldView.h"

#include "Model/Model.h"
#include "Layers/Manager.h"
#include "Utilities/To.h"
#include "Utilities/Types.h"

// namespaces
namespace niwa {

void WorldView::Build() {
  LOG_TRACE();
  base_layer_ = model_->managers().layer()->GetLayer(model_->get_base_layer());
  if(!base_layer_) {
    LOG_ERROR() << "The base layer '"<< model_->get_base_layer() << "' found in the @model block could not be found, please check there is a @layer defined for this layer";
  }


}


} /* namespace niwa */
