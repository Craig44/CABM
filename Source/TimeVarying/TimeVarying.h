/**
 * @file TimeVarying.h
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 27/01/2015
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * << Add Description >>
 */
#ifndef TIMEVARYING_H_
#define TIMEVARYING_H_

// headers
#include "BaseClasses/Object.h"
#include "Model/Model.h"

// namespaces
namespace niwa {

/**
 *
 */
class TimeVarying : public niwa::base::Object {
  typedef void (TimeVarying::*UpdateFunction)(double);
public:
  // methods
  TimeVarying() = delete;
  explicit TimeVarying(Model* model);
  virtual                     ~TimeVarying() = default;
  void                        Validate();
  void                        Build();
  void                        Reset();
  void                        Update(unsigned current_year);

  //accessors
  map<unsigned, double>&      get_parameter_by_year() { return parameter_by_year_; }

protected:
  // methods
  void                        RestoreOriginalValue();

  // settors
  void                        set_single_value(double value);
  void                        set_vector_value(double value);
  void                        set_map_value(double value);

  // pure virtual methods
  virtual void                DoValidate() = 0;
  virtual void                DoBuild() = 0;
  virtual void                DoReset() = 0;
  virtual void                DoUpdate() = 0;

  // function pointers
  UpdateFunction              update_function_ = 0;

  // members
  Model*                      model_;
  base::Object*               target_object_;
  string                      type_ = "";
  vector<unsigned>            years_;
  string                      parameter_;
  double                       original_value_ = 0;
  map<unsigned, double>*       addressable_map_ = 0;
  vector<double>*              addressable_vector_ = 0;
  double*                      addressable_ = 0;
  map<unsigned, double>        parameter_by_year_;
};

typedef std::shared_ptr<TimeVarying> TimeVaryingPtr;

} /* namespace niwa */

#endif /* TIMEVARYING_H_ */
