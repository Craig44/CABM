/**
 * @file Mortality.h
 * @author C.Marsh
 * @github https://github.com/Craig44
 * @date 13/07/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This process is the parent mortlity class, so the manager can generate dynamicasts<> processes
 * any child mortality class should inherit this as a parent
 */
#ifndef SOURCE_PROCESSES_CHILDREN_MORTALITY_H_
#define SOURCE_PROCESSES_CHILDREN_MORTALITY_H_

// headers
#include "Processes/Process.h"

// namespaces
namespace niwa {
namespace processes {

struct composition_data {
  string type_;
  unsigned year_;
  unsigned row_;
  unsigned col_;
  vector<unsigned> frequency_;
  composition_data(string type, unsigned year, unsigned row, unsigned col, unsigned size) : type_(type), year_(year),
      row_(row), col_(col)
  {
    frequency_.resize(size);
  }
};

/**
 * Class definition
 */
class Mortality : public Process {
public:
  // methods
  explicit Mortality(Model* model);
  virtual                     ~Mortality() = default;
  virtual void                        DoValidate(){ };
  virtual void                        DoBuild() { };
  virtual void                        DoReset() { };
  virtual void                        DoExecute() { };

  virtual void                        draw_rate_param(unsigned row, unsigned col, unsigned number_of_draws, vector<float>& vector) = 0;
  vector<composition_data>&           get_removals_by_age() {return removals_by_age_and_area_;};
  vector<composition_data>&           get_removals_by_length() {return removals_by_length_and_area_;};
  virtual bool                        update_mortality() {return update_natural_mortality_parameters_;};
  virtual double                      SolveBaranov() { return 1.0;};
  void                                set_lambda(double lambda) {lambda_ = lambda;};
protected:
  vector<composition_data>            removals_by_age_and_area_;
  vector<composition_data>            removals_by_length_and_area_;
  map<unsigned, vector<unsigned>>     removals_by_age_;
  map<unsigned, vector<unsigned>>     removals_by_length_;
  bool                                update_natural_mortality_parameters_;
  double                              lambda_;



};

} /* namespace processes */
} /* namespace niwa */

#endif /* SOURCE_PROCESSES_CHILDREN_MORTALITY_H_ */
