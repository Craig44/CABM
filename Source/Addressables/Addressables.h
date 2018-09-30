/**
 * @file Addressables.h
 * @author C.Marsh
 * @date 30/09/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * Because we handle the values loaded by an input file differently depending on what
 * run mode we are in this class is responsible for loading the values from the input
 * file and making them available to the sub-systems that require them.
 */
#ifndef ADDRESSABLES_H_
#define ADDRESSABLES_H_

// headers
#include <map>
#include <string>
#include <memory>

#include "Model/Managers.h"
#include "Utilities/Map.h"
#include "Utilities/Types.h"

// namespaces
namespace niwa {

using std::shared_ptr;
using utilities::Double;
using std::string;
using std::vector;
using std::map;
class Model;

// Enumerated Types
enum class AddressableType {
  kInvalid      = 0,
  kSingle       = 1,
  kVector       = 2,
  kStringMap    = 3,
  kUnsignedMap  = 4
};


/**
 * Class definition
 */
class Addressables {
public:
  // methods
  Addressables(Model* model) : model_(model) { };
  virtual                       ~Addressables() = default;
  void                          AddValue(const string& addressable_label, float value);
  vector<string>                GetAddressables() const;
  unsigned                      GetValueCount() const;
  map<string, float>            GetValues(unsigned index) const;
  void                          LoadValues(unsigned index);

private:
  // members
  Model*                        model_ = nullptr;
  map<string, vector<float>>    addressable_value_;
  map<string, float*>           addressables_;

};
} /* namespace niwa */

#endif /* ADDRESSABLES_H_ */
