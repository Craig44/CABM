/**
 * @file Addressables.cpp
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 20/01/2015
 * @section LICENSE
 *
 * Copyright NIWA Science ©2014 - www.niwa.co.nz
 *
 */

// headers
#include "Addressables.h"

#include "GlobalConfiguration/GlobalConfiguration.h"
#include "Model/Model.h"
#include "Model/Objects.h"
#include "Logging/Logging.h"

// namespaces
namespace niwa {

/**
 *
 */
void Addressables::AddValue(const string& addressable_label, float value) {
  addressable_value_[addressable_label].push_back(value);
}

/**
 *
 */
vector<string> Addressables::GetAddressables() const {
  vector<string> result;
  for (auto iter : addressable_value_)
    result.push_back(iter.first);

  return result;
}

/**
 *
 */
unsigned Addressables::GetValueCount() const {
  if (addressable_value_.size() == 0)
    return 0;

  auto iter = addressable_value_.begin();
  return iter->second.size();
}
/**
 *
 */
map<string, float> Addressables::GetValues(unsigned index) const {
  map<string, float> result;
  for (auto iter : addressable_value_)
    result[iter.first] = iter.second[index];
  return result;
}

/**
 *
 */
void Addressables::LoadValues(unsigned index) {
  /**
   * load our estimables if they haven't been loaded already
   */
  if (addressables_.size() == 0) {
    string error = "";
    for (auto iter : addressable_value_) {
      if (!model_->objects().VerfiyAddressableForUse(iter.first, addressable::kInputRun, error)) {
        LOG_FATAL() << "The addressable " << iter.first << " could not be verified for use in -i run. Error was " << error;
      }
      float* ptr = model_->objects().GetAddressable(iter.first);
      addressables_[iter.first] = ptr;
    }
  }

  for (auto iter : addressables_) {
    if (index >= addressable_value_[iter.first].size())
      LOG_CODE_ERROR() << "index >= addressable_value_[iter.first].size()";
    LOG_FINE() << "Changing parameter " << iter.first  << " to = " << addressable_value_[iter.first][index];
    (*addressables_[iter.first]) = addressable_value_[iter.first][index];
  }
}




} /* namespace niwa */
