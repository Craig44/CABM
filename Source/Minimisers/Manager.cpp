/**
 * @file Manager.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 14/05/2013
 * @section LICENSE
 *
 * Copyright NIWA Science ©2013 - www.niwa.co.nz
 *
 */

// headers
#include "Manager.h"

// namespaces
namespace niwa {
namespace minimisers {

/**
 * Default constructor
 */
Manager::Manager() {
}

/**
 * Destructor
 */
Manager::~Manager() noexcept(true) {

}

/**
 * Validate the minimisers
 */
void Manager::Validate() {
  LOG_TRACE();
  for (auto minimiser : objects_)
    minimiser->Validate();

  /*
  if (objects_.size() > 1) {
    LOG_FATAL() << "Found more than one @minimiser block currently this is not allowed, please sort out";
  }
  if (objects_.size() == 1) {
    objects_[0]->set_active(true);
    active_minimiser_ = objects_[0];

  }
   */
}

Minimiser* Manager::get_minimiser(string label) {
  for (auto this_minimiser : objects_) {
    if (this_minimiser->label() == label)
      return this_minimiser;
  }
  return nullptr;
}

} /* namespace minimisers */
} /* namespace niwa */
