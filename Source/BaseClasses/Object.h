/**
 * @file BaseObject.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 18/09/2012
 * @section LICENSE
 *
 * Copyright NIWA Science �2012 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This class is the highest level class in our object tree. All objects
 * that have to be created from a configuration block need to inherit from this class.
 *
 * This class creates access to some central global objects within the system
 * including: GlobalConfiguration, TheModel, ParameterList etc
 *
 */
#ifndef BASE_OBJECT_H_
#define BASE_OBJECT_H_

// Headers
#include <string>
#include <memory>

#include "Logging/Logging.h"
#include "ParameterList/ParameterList.h"
#include "Translations/Translations.h"
#include "Utilities/Map.h"
#include "Utilities/NoCopy.h"
#include "Utilities/Types.h"


// Namespaces
namespace niwa {
// Enumerated Types
namespace addressable {
enum Type {
  kInvalid      = 0,
  kSingle       = 1,
  kMultiple     = 2,
  kVector       = 3,
  kStringMap    = 4,
  kUnsignedMap  = 5,
  kVectorStringMap = 6
};

enum Usage {
  kNone         = 0,
  kLookup       = 1, // Assert, Additional Prior, Equation, Reports
  kInputRun     = 4,
  kSimulate     = 32,
  kTimeVarying  = 128,
  kAll          = 255
};

enum rerun_initialisation {
  kno         = 0,
  kyes         = 1
};
};

namespace base {

using std::string;
using utilities::OrderedMap;
using std::vector;
using std::map;

// classes
class Object {
public:
  // Methods
  Object() = default;
  virtual                         ~Object() {};
  bool                            HasAddressable(const string& label) const;
  bool                            HasAddressableUsage(const string& label, const addressable::Usage&) const;
  bool                            IsAddressableAVector(const string& label) const;
  unsigned                        GetAddressableSize(const string& label) const;
  double*                          GetAddressable(const string& label);
  virtual double*                  GetAddressable(const string& label, const string& index);
  vector<double*>*                 GetAddressables(const string& absolute_label, const vector<string> indexes);
  map<unsigned, double>*           GetAddressableUMap(const string& label);
  map<unsigned, double>*           GetAddressableUMap(const string& label, bool& create_missing);
  OrderedMap<string, double>*      GetAddressableSMap(const string& label);
  vector<double>*                  GetAddressableVector(const string& label);
  addressable::Type               GetAddressableType(const string& label) const;
  addressable::Usage              GetAddressableUsage(const string& label) const;
  addressable::rerun_initialisation GetAddressableInit(const string& label) const;

  void                            PrintParameterQueryInfo();
  void                            SubscribeToRebuildCache(Object* subscriber);
  void                            NotifySubscribers();
  virtual void                    RebuildCache();

  // pure virtual methods
  virtual void                    Reset() = 0;


  // Accessors and Mutators
  string                      label() const { return label_; }
  string                      type() const { return type_; }
  ParameterList&              parameters() { return parameters_; }
  string                      location();
  void                        set_block_type(string value) { block_type_ = value; parameters_.set_parent_block_type(value); }
  void                        set_label(string value) { label_ = value;}
  void                        set_defined_file_name(string value) { parameters_.set_defined_file_name(value); }
  void                        set_defined_line_number(unsigned value) { parameters_.set_defined_line_number(value); }
  string                      block_type() const { return block_type_; }

protected:
  // Methods
  void                        RegisterAsAddressable(const string& label, double* variable, addressable::Usage usage = addressable::kAll, addressable::rerun_initialisation re_run_init = addressable::kno);
  void                        RegisterAsAddressable(const string& label, vector<double>* variables, addressable::Usage usage = addressable::kAll, addressable::rerun_initialisation re_run_init = addressable::kno);
  void                        RegisterAsAddressable(const string& label, OrderedMap<string, double>* variables, addressable::Usage usage = addressable::kAll, addressable::rerun_initialisation re_run_init = addressable::kno);
  void                        RegisterAsAddressable(const string& label, map<unsigned, double>* variables, addressable::Usage usage = addressable::kAll, addressable::rerun_initialisation re_run_init = addressable::kno);
  void                        RegisterAsAddressable(map<string, vector<double>>* variables);

  // Members
  string                          block_type_           = "";
  string                          label_                = "";
  string                          type_                 = "";
  ParameterList                   parameters_;
  map<string, double*>             addressables_;
  map<string, bool>               create_missing_addressables_;
  map<string, vector<double>* >    addressable_vectors_;
  map<string, vector<double*> >    addressable_custom_vectors_;
  map<string, addressable::Type>  addressable_types_;
  map<string, addressable::Usage> addressable_usage_;
  map<string, addressable::rerun_initialisation> addressable_initphase_;

  vector<Object*>                 rebuild_cache_subscribers_;

  map<string, map<unsigned, double>* >      addressable_u_maps_;
  map<string, OrderedMap<string, double>* > addressable_s_maps_;
  vector<map<string, vector<double>>* >     unnamed_addressable_s_map_vector_;

  DISALLOW_COPY_AND_ASSIGN(Object);
};

} /* namespace base */
} /* namespace niwa */
#endif /* BASE_OBJECT_H_ */
